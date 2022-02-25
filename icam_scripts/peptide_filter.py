#    iCaM2.0 peptide filter script
#    Copyright (C) 2011-2016 TRON gGmbH
#    This software implements procedures protected by the patents:
#    - PCT/EP2013/003559 (WO2014/082729) Individualized Vaccines for Cancer
#    - PCT/EP2014/001232 (WO2014/180569) Immunogenicity Prediction of T Cell Epitope
#    - PCT/EP2015/053021 () Predicting T cell epitopes useful for vaccination
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

import sys
import operator
import math
from nucleotide_to_peptide_classes import mutation_filter

def _get_TROMPA_column_dict():
  d = {}
  d["key"] = 0
  d["mut"] = 1
  d["trid"] = 4
  d["expr"] = 5
  d["codon"] = 9
  d["subst"] = 10
  d["pep"] = 11
  d["mhc"] = 16
  d["mhcII"] = 22
  d["phase"] = 31
  d["rnavaf"] = 41
  d["pep_WT"] = 12
  d["mRNA"] = 13
  d["epi_clI"] = 17
  d["epi_clII"] = 23
  return d

def run(infile,outfile,filter="default",sep="\t",options=[None,48,5],do_squish=True,write_squish_data=True):
  l = []
  with open(infile) as f:
    for line in f:
      l.append(line.strip("\n").split(sep))
  if len(l) <= 1:
    open(outfile, 'a').close()  # Write empty file
    return 0, 0, 0, 0, 0, [], 0, 0, 0
  l1 = len(l) - 1
  if do_squish:
    l = uniq_mutpeptide(l,write_squish_data)
  if filter == "class_II_I":
    l,returned,passed,total,nmut,filterinfo, x, xx, xxx = class_II_I_filter(l,*options)
  else:
    mut = 1
    nmut = len(set(zip(*l[1:])[mut]))
    l2 = len(l)
    returned,passed,total,filterinfo,x,xx,xxx = l2,l2,l2,[],0,0,0
  with open(outfile,"w") as f:
    for line in l:
      f.write(sep.join(line))
      f.write("\n")
  return l1, total, passed, returned, nmut, filterinfo,x,xx,xxx

def main():
  f = "class_II_I"
  if len(sys.argv) > 3:
    f = sys.argv[3]
  n=96
  n2=5
  if len(sys.argv) > 4:
    n = int(sys.argv[4])
  if len(sys.argv) > 5:
    n2 = int(sys.argv[5])
  l1, total, passed, returned, nmut, filterinfo ,x,xx,xxx = run(sys.argv[1],sys.argv[2],f,"\t",[None,n,n2],False,False)
  print l1, total, passed, returned, nmut
  for i in filterinfo:
    print i

def uniq_mutpeptide(l,write_squish_data): 
  """
  collapse equivalency classes, use the class representative with highest transcript expression
  """
  dt = _get_TROMPA_column_dict()
  mut = dt["mut"]
  pep = dt["pep"]
  transcript_exp = dt["expr"]
  trid = dt["trid"]
  KEY = dt["key"]
  d = {}
  c = 0
  for i in l:
    c += 1
    if c == 1:
      header = i
      header.append("other_transcripts_with_this_peptide")
      new1 = len(header)-1
      header.append("peptide_resulting_from_this_mutation")
      new2 = len(header)-1
      header.append("distinct_peptides_resulting_from_this_mutation")
      new3 = len(header)-1
      header.append("keys_of_distinct_peptides_resulting_from_this_mutation")
      new4 = len(header)-1
      continue
    id = i[mut] + "_" + i[pep]
    x = d.get(id,[])
    x.append(i)
    d[id] = x
  r = []
  # get all transcript ids with same peptide, selecte the one with the highest expression (could be also the sum ?)
  for k in d.keys():
    l = d[k]
    l = sorted(l,key=lambda x:to_float(x[transcript_exp]),reverse=True)
    sel = l[0]
    tl = [x[trid] for x in l[1:]]
    sel.append(" ".join(tl))
    r.append(sel)
  ## get info on distinct peptides from one mutation
  d = {}
  for i in r:
    x = d.get(i[mut],[])
    x.append(i)
    d[i[mut]] = x
  for k in d.keys():
    l = d[k]
    keylist = []
    for i,j in enumerate(l):
      j.append(str(i+1))
      j.append(str(len(l)))
      keylist.append(j[KEY])
    for i,j in enumerate(l):
      idx = keylist.index(j[KEY])
      j.append(" ".join(keylist[:idx] + keylist[idx+1:]))
  r.sort(key=lambda x:to_float(x[new2][0]))
  r.sort(key=operator.itemgetter(mut))
  if not write_squish_data: # remove last four columns, if they are not wanted
    header = header[:-4]
    r2 = r
    r = []
    for i in r2:
      r.append(i[:-4])
  return [header] + r

def select_by_mutation_maximum(r,mut,maxmut):
  # select up to maxmut mutations (!= number of entries in list r)
  if (len(r) == 0) or (maxmut == 0):
    return [],0
  nmut = len(set(zip(*r)[mut]))
  newlen = len(r)
  while (nmut > maxmut):
    newlen = newlen - 1
    nmut = len(set(zip(*r[:newlen])[mut]))
  r = r[:newlen]
  return r,nmut

def class_II_I_filter(l_input,mhc_replacement,n,mut_indels_max=5):

  mut_indels = mut_indels_max
  n_wo_indels = n - mut_indels
  n_per_set = (n_wo_indels - (n_wo_indels % 10)) / 2 # we want two sets for the current default filter, one will be priotizrd by MHC class I score, the other one by MHC class II scores
  #mut_indels = n - (n_per_set * 2)
  #if mut_indels > mut_indels_max: # use max. 5 indels -> this could be a parameter in the future
  #  mut_indels = mut_indels_max
  mut_classII = mut_classI = n_per_set

  # add ranks
  l_input = _get_ranks(l_input)

  # this list will hold human readable infromation on the filter steps
  complete_filter_info = []

  # remove obvious non-candidates: phasing problems, stop codons, transcript not expressed
  complete_filter_info.append("## Peptide prefilter:")
  l0,filterinfo_pre,passed1 = _prefilter(l_input)
  header = l0[0]
  complete_filter_info += filterinfo_pre

  # split indels and SNVs
  dt = _get_TROMPA_column_dict()
  subst = dt["subst"]
  l1 = []
  indel_l = []
  for words in l0:
    if words[subst] == "":
      indel_l.append(words)
    else:
      l1.append(words)

  # select indels
  complete_filter_info.append("## Indels:")
  l_indels,n_indel_mutations,filterinfo_indel,passed_indel = _select_indels(indel_l,mut_indels)
  complete_filter_info += filterinfo_indel

  # first run of class_II_I filter, with RNA_VAF > 0
  complete_filter_info.append("## Selecting mutations with RNA VAF > 0")
  snv_slots = n - n_indel_mutations
  l_selected_first_round,remaining_first_round,selected_snvs_first_round,passed_first_round,filterinfo_1st = _class_II_I_filter_SNV(l1,mhc_replacement,snv_slots,mut_classII,mut_classI,0.0,"larger_than")
  complete_filter_info += filterinfo_1st
  # are we done yet ?
  # - any free slots in the available slots for snvs ?
  # - are there any SNVs left ?
  if (selected_snvs_first_round < snv_slots) and (len(remaining_first_round) > 0):
    complete_filter_info.append("## Selected only {0} mutations, up to {1} are possible ".format(selected_snvs_first_round,snv_slots))
    complete_filter_info.append("## Selecting mutations with RNA VAF == 0")
    snv_slots2 = snv_slots - selected_snvs_first_round
    # run second round of selection, selecting only from mhc_class_I binders
    remaining_first_round = [l1[0]] + remaining_first_round 
    l_selected_2nd_round,remaining_2nd_round,selected_snvs_2nd_round,passed_2nd_round,filterinfo_2nd = _class_II_I_filter_SNV(remaining_first_round,mhc_replacement,snv_slots2,0,snv_slots2,0.0,"equals")
    complete_filter_info += filterinfo_2nd
  else:
    l_selected_2nd_round = []
    selected_snvs_2nd_round = 0
    passed_2nd_round = 0
  lr = [header] + l_selected_first_round + l_selected_2nd_round + l_indels
  returned = len(lr) - 1
  passed = passed_first_round + passed_2nd_round + passed_indel
  total = len(indel_l) + len(l1) - 1
  nmut = n_indel_mutations + selected_snvs_first_round + selected_snvs_2nd_round
  return lr,returned,passed,total,nmut,complete_filter_info,len(l_selected_first_round),len(l_selected_2nd_round),len(l_indels)
  
def _select_indels(indel_l,mut_indels,maxlen_AA=50,minlen_AA=14):
  if len(indel_l) == 0:
    return [],0,[],0
  dt = _get_TROMPA_column_dict()
  rnavaf = dt["rnavaf"]
  mut = dt["mut"]
  expr = dt["expr"]
  mhc = dt["mhc"]
  pep = dt["pep"]
  indels_total = len(indel_l)
  rx = []
  for words in indel_l:
    if float(words[rnavaf]) > 0.0:
      words_mod = _cut_indels(words,maxlen_AA)
      if len(words_mod[pep]) >= minlen_AA:
        rx.append(words_mod)
  passed_indel = len(rx)
  r = sorted(sorted(rx,key=lambda x:to_float(x[expr]),reverse=True),key=lambda x:to_float(x[mhc]))
  r,n_indel_mutations = select_by_mutation_maximum(r,mut,mut_indels)
  indel_pep_filtered = len(r)
  filterinfo = []
  filterinfo.append("Selected {0} expressed indel peptides out of {1} (corresponding to {2} mutations), sorted my MHC class I binding".format(indel_pep_filtered,len(indel_l),n_indel_mutations))
  return r,n_indel_mutations,filterinfo,passed_indel

def _cut_indels(line,maxlen_AA):
  dt = _get_TROMPA_column_dict()
  pep = line[dt["pep"]]
  mRNA = line[dt["mRNA"]]
  if len(pep) <= maxlen_AA:
    return line
  epi_I = line[dt["epi_clI"]]
  epi_II = line[dt["epi_clII"]]
  p1 = pep.find(epi_I)
  p2 = pep.find(epi_II)
  positions = [p1,p2,p1 + len(epi_I), p2 + len(epi_II)]
  start = min(positions)
  end = max(positions)
  if (end - start) <= maxlen_AA:
    l = list(line)
    selected_peptide = pep[start:end]
    p1 = start
  else:
    selected_peptide = epi_I  
  diff = int(math.floor(float(maxlen_AA)/2) - math.ceil(float(len(selected_peptide))/2))
  diff = max(0,diff)
  start = max(0,p1 - diff)
  end = min(p1 + len(selected_peptide) + diff,len(pep))
  l = list(line)
  l[dt["pep"]] = pep[start:end]
  l[dt["mRNA"]] = mRNA[start*3:end*3]
  return l

def _get_ranks(l):
  header = l[0]
  r = l[1:]
  dt = _get_TROMPA_column_dict()
  expr = dt["expr"]
  mhc = dt["mhc"]
  mhcII = dt["mhcII"]
  # sort exon expression, then by MHCI
  r = sorted(sorted(r,key=lambda x:to_float(x[expr]),reverse=True),key=lambda x:to_float(x[mhc]))
  # add rank info
  header.append("RANK_MHC")
  for i,j in enumerate(r):
    j.append(str(i+1))
  # sort exon expression, then by MHCI
  r = sorted(sorted(r,key=lambda x:to_float(x[expr]),reverse=True),key=lambda x:to_float(x[mhcII]))
  # add rank info
  header.append("RANK_MHC_II")
  for i,j in enumerate(r):
    j.append(str(i+1))
  # sort MHC, then by transcript expression
  r = sorted(sorted(r,key=lambda x:to_float(x[mhc])),key=lambda x:to_float(x[expr]),reverse=True)
  # add rank info
  header.append("RANK_EXPR")
  for i,j in enumerate(r):
    j.append(str(i+1))
  return [header] + r

def _prefilter(l):
  dt = _get_TROMPA_column_dict()
  phase = dt["phase"]
  expr = dt["expr"]
  subst = dt["subst"]
  codon = dt["codon"]
  mhc = dt["mhc"]
  mhcII = dt["mhcII"]
  header = l[0]
  lx = l[1:]
  filter_list = []
  passed = 0
  filterinfo = []
  r = []
  filter_list.append(mutation_filter(phase,header,str,"does_not_contain","rejected"))
  filter_list.append(mutation_filter(phase,header,str,"does_not_contain","unknown"))
  filter_list.append(mutation_filter(expr,header,float,"larger_than",0.0))
  filter_list.append(mutation_filter(subst,header,str,"does_not_end_with","-"))
  filter_list.append(mutation_filter(codon,header,str,"is_not_equal_to","0"))
  
  mhcI_filter = mutation_filter(mhc,header,float,"is_not_equal_to",999.0)
  filter_list.append(mhcI_filter)
  mhcII_filter = mutation_filter(mhcII,header,float,"is_not_equal_to",999.0)
  # test if all MHC class 2 predictions are invalid and skip mhcII_filter then
  r_test = []
  for words in lx:
    words = mhcII_filter.apply(words)
    if len(words) != 0:
      r_test.append(words)
  if len(r_test) != 0:
    filter_list.append(mhcII_filter)
  else:
    filterinfo.append("No HLA class II prediction, skipping class II prefilter")
  
  for words in lx:
    try:
      for f in filter_list:
        words = f.apply(words)
      if len(words) != 0:
        r.append(words)
        passed += 1
    except Exception as e:
      raise e
  for f in filter_list:
    filterinfo.append(f.get_string())
  return [header] + r,filterinfo,passed


def _class_II_I_filter_SNV(l,mhc_replacement,maxmut,mut_classII,mut_classI,rnaVAF_threshold,rnaVAF_filter):
  dt = _get_TROMPA_column_dict()
  mut = dt["mut"]
  expr = dt["expr"]
  mhc = dt["mhc"]
  mhcII = dt["mhcII"]
  phase = dt["phase"]
  rnavaf = dt["rnavaf"]
  passed = 0
  filter_list = []
  r = []
  header = l[0]
  lx = l[1:]
  remaining = []
  filter_list.append(mutation_filter(rnavaf,header,float,rnaVAF_filter,rnaVAF_threshold))
  for words in lx:
    w_tmp = list(words) # store for later use
    if (mhc_replacement!=None) and words[mhc] == "":
      words[mhc] = str(mhc_replacement)
    try:
      for f in filter_list:
        words = f.apply(words)
      if len(words) != 0:
        r.append(words)
        passed+=1
      else:
        remaining.append(w_tmp) # did not pass filter
    except Exception as e:
      raise e 
 
  # the next steps are just resorting and tracking the differnt set sizes etc
  if mut_classII > 0: 
    nx = mut_classII # maximum number of mutations represented by the first block (mhc_classII binders)
    min_expr = 10.0 # minimum gene expression
    r = [tuple(x) for x in r] # compareable elements
    # sort and select (MHC II)
    r = sorted(sorted(r,key=lambda x:to_float(x[expr]),reverse=True),key=lambda x:to_float(x[mhcII]))
    r_temp = []
    for i in r:
      if float(i[expr]) >=  min_expr:
        r_temp.append(i)
    r1,nmut1 = select_by_mutation_maximum(r_temp,mut,nx)
    maxmut = maxmut - nmut1
    # remove those in r1 from r -> r2
    r2 = []
    for i in r:
      if i not in r1:
        r2.append(i)
  else:
    r2 = r
    r1 = []
    nmut1 = 0
 
  nx = mut_classI # maximum number of mutations represented by the second block (mhc_classI binders)
  min_expr = 1.0 # different for class 1 
  # sort and select (MHC I)
  r2 = sorted(sorted(r2,key=lambda x:to_float(x[expr]),reverse=True),key=lambda x:to_float(x[mhc]))
  r_temp = []
  for i in r2:
    if float(i[expr]) >=  min_expr:
      r_temp.append(i)
  if len(r_temp) > 0:
    r3,nmut3 = select_by_mutation_maximum(r_temp,mut,nx)
    maxmut = maxmut - nmut3
    # remove those in r3 from r2 -> r4
  else:
    r3 = []
    nmut3 = 0

  r4 = []
  for i in r2:
    if i not in r3:
      r4.append(i)

  # sort and select (expression)
  r4 = sorted(sorted(r4,key=lambda x:to_float(x[mhc])),key=lambda x:to_float(x[expr]),reverse=True)
  if len(r4) > 0:
    r5,nmut5 = select_by_mutation_maximum(r4,mut,maxmut)
  else:
    r5 = []
    nmut5 = 0

  # store the ones which are still remaining
  for i in r4:
    if i not in r5:
      remaining.append(i)

  # gather textual information on filters
  filterinfo = []
  for f in filter_list:
    filterinfo.append(f.get_string())
  filterinfo.append(str(len(r1)) + " peptides sorted by transcript expression and then by HLA class II prediction score (transcript expression >= 10)")
  filterinfo.append(str(len(r3)) + " peptides sorted by transcript expression and then by HLA class I prediction score (transcript expression >= 1)")
  filterinfo.append(str(len(r5)) + " peptides sorted by transcript expression")

  selected_snvs = nmut1 + nmut3 + nmut5
  ret = [list(x) for x in r1 + r3 + r5]
  return ret,remaining,selected_snvs,passed,filterinfo

def to_float(x):
  if x == "":
    return -1
  else:
    return float(x)



if __name__ == '__main__':
  main()


