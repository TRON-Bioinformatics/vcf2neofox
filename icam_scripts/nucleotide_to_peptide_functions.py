#    iCaM2.0 NGS data analysis pipeline
#    Copyright (C) 2011-2016 TRON gGmbH
#    See LICENSE for licensing details.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

import sys
import math
import traceback
import operator
import time

# import constant data
from nucleotide_to_peptide_constants import *

from nucleotide_to_peptide_IEDB_mp import run_mhc_prediction3 

import phase

def check_phase(mut_list,bam_files,samtools_bin,phasing_dict,phasing_report,phasing_info,min_reads=1):
  """
  Check for each mutation pair in the input list if the mutations are in phase. 

  Args:
    mut_list: list of objects of calss mutation
    bam_files: list of filenames of BAM files
    samtools_bin: path of samtools binary
    phasing_dict: a dictionary
    phasing_report: a list
    phasing_info: a dictionary
    min_reads: minimum number of read with both mutations to accept phase (default=1)

  Returns:
    A triple (d,report,p_dict) with:
      d: a dictionary, including all data of argument phasing_dict plus new key-value pairs, where the keys are mutation positions and the values are the sets of mutation positions which are in phase with the mutation identified by the key 
      report: a list, including all elements of argument phasing_report, plus appended strings which describe detected mutation in phase and the evidence from the BAM files in human readable form
      p_dict: a dictionary, including all data of argument phasing_info plus new key-value pairs, where the keys are strings consisting of a mutation pair (e.g. "chr1_2314|chr1_24314") and the values are stings which describe the phase of the mutation pair - one of "in_phase" or "phase_unknown"; the phasing check is skipped when a mutation pair is already in this dictionary
"""
  d = phasing_dict
  report = phasing_report
  p_dict = phasing_info
  for bi,b in enumerate(bam_files):
    for m1 in mut_list:
      for m2 in mut_list:
        pos1 = m1.get_pos()
        pos2 = m2.get_pos()
        if pos1 == pos2:
           d = _add_to_set(d,pos1,pos1) # every mutation is in phase with itself
        p12 = "|".join(sorted([pos1,pos2]))
        if p12 in p_dict:
          continue # check was done before, so skip it
        elif not pos1 == pos2:
          p = phase.run_m(m1,m2,b,samtools_bin)
          if (p[0] > 0) and (p[1] > 0):  # is the any coverage on both loci ?
            if (p[2] >= min_reads):
              d = _add_to_set(d,pos1,pos2)
              d = _add_to_set(d,pos2,pos1)
              report.append("{0} is in phase with {1} (coverage is {2} and {3}, {4} reads have both mutations in bam file {5})".format(m1.get_pos(),m2.get_pos(),p[0],p[1],p[2],bi))
              p_dict[p12] = "in_phase"
            else:
              p_dict[p12] = "phase_unknown" # we are conservative here; we cannot be sure that we have sampled sufficent number of reads: if we could be sure, "not_in_phase" would be the correct value
          else:
            p_dict[p12] = "phase_unknown"
  return d,report,p_dict

def _add_to_set(d,k,v):
  x = d.get(k,set())
  x.add(v)
  d[k] = x
  return d

def _set_all_phase_info(all_seqs_mutated,value):
  for i in all_seqs_mutated:
    for j in i:
      j.next_mutationPHASE = value
  return all_seqs_mutated

def modify_mutated_sequence(mutated_seq1,mutated_seq2):
  # modify mutant seq
  seq = list(mutated_seq1.mutated_sequence)
  seq[mutated_seq2.nucl_pos] = mutated_seq2.mutated_sequence[mutated_seq2.nucl_pos]
  mutated_seq1.mutated_sequence = "".join(seq)

  # also modify WT sequence - when there is a SNP it is also copied in the new WT sequence
  seq = list(mutated_seq1.wt_sequence)
  seq[mutated_seq2.nucl_pos] = mutated_seq2.wt_sequence[mutated_seq2.nucl_pos]
  mutated_seq1.wt_sequence = "".join(seq)  
  
  # re-initialize the mutated_seq1 object
  mutated_seq1._init_seqs(mutated_seq1.len_5p,mutated_seq1.len_3p,mutated_seq1.strand)



def merge_mutated_seqs_in_phase(this_transcript_info,all_seqs_mutated,mut,peptide_len_snv,max_dist=None):
  """
  Merge non-synonymous individual SNV and/or SNP effects in differnet codons of the input sequences, or merge a non-synonymous individual SNV and/or SNP effect with a synonymous SNP/SNV in the same codon, when this results in a new mutant sequence. The function will only merge up to two nearby substitution events. Larger substitution groups are rejected, as well as groups of substitutions and frameshifts. The field "next_mutationPHASE" of each mutated_sequence object will be set to a value according to the merging result: "rejected:more_than_two_mutations", "rejected:mutli_allelic", "rejected:snv_indel", "in_phase:merged_two_codons" or "in_phase:merged_same_codon".

  Args:
    this_transcript_info: an object of class transcript_info, representing a transcript
    all_seqs_mutated: list of tuples with objects of class mutated_sequence, representing the individual effects of all mutations associated with the given transcript; includes mutant sequences
    mut: list of objects of class mutation, representing all mutations associated with the given transcript
    peptide_len_snv: max. C- and N-terminal extension length of a peptide around mutant amino acids 
    max_dist: maximum distance of two mutations on a transcript in order to be considered for merging

  Returns:
    r: a list of lists of objects of class mutated_sequence, representing the (possibly) merged sequences
    messages: a list of strings, representing textual information on the merging events
"""
  # all_seqs_mutated contains also syn sequences
  # ... unless they are in an UTR

  if not max_dist:
    max_dist = (peptide_len_snv * 3) + 2
	
  # set id and distance of next mutation, if there is any  
  for i in all_seqs_mutated:
    for j in i:
      j.find_next_mutation()

  if (len(mut) == 1) or (len(all_seqs_mutated) <=1):
    # nothing to do
    return all_seqs_mutated,[]

  # get distances beween mutations and sort mutations whcih are close into clusters
  clusters = []
  tmp_pos = []
  for i in all_seqs_mutated:
    tmp_pos.append(i[0].CDS_nucl_pos)
  all_seqs_mutated = zip(*sorted(zip(tmp_pos,all_seqs_mutated), key=operator.itemgetter(0)))[1]
  cluster = [all_seqs_mutated[0]]
  for i in range(1,len(all_seqs_mutated)):
    mutated_seq1 = all_seqs_mutated[i][0]
    mutated_seq2 = all_seqs_mutated[i-1][0]
    mut1 = mutated_seq1.mutation.get_pos()
    mut2 = mutated_seq2.mutation.get_pos()
    dist = abs(mutated_seq1.CDS_nucl_pos - mutated_seq2.CDS_nucl_pos)
    if dist > max_dist:
      clusters.append(list(cluster))
      cluster = [all_seqs_mutated[i]]
    else:
     cluster.append(all_seqs_mutated[i])
  clusters.append(cluster)

  # set phasing status for next mutation
  for i in all_seqs_mutated:
    mut1 = i[0].mutation.get_pos()
    mut2 = i[0].next_mutationID
    if mut2 != "":
      i[0].next_mutationPHASE = this_transcript_info.phase_info.get("|".join(sorted([mut1,mut2])),"not_tested")

  r = []
  messages = []
  for i,cluster in enumerate(clusters):
    if len(cluster) <= 1:
      # noting to do, just set phase info to an empty value_ the phase of the next mutation is irrelevant - we know that there are more than one mutation in this transcript
      # but the mutation in this respective (single member) cluster is too distant from all other ones
      cluster = _set_all_phase_info(cluster,"")
      r += cluster
      continue
    if len(cluster) > 2:
      messages.append("more than two phased mutations in transcript " + this_transcript_info.id + ", cluster "+ str(i) +", skipping seq merging")
      cluster = _set_all_phase_info(cluster,"rejected:more_than_two_mutations")
      r += cluster
      continue
    if ((len(cluster[0]) > 1) and cluster[0][0].prot_is_different and cluster[0][1].prot_is_different) or ((len(cluster[1]) > 1) and cluster[1][0].prot_is_different and cluster[1][1].prot_is_different):
      messages.append("multi allelic phased mutations in transcript " + this_transcript_info.id + ", cluster "+ str(i) +", skipping seq merging")
      cluster = _set_all_phase_info(cluster,"rejected:mutli_allelic")
      r += cluster
      continue
    if cluster[0][0].prot_is_different or (len(cluster[0]) == 1):
      mutated_seq1 = cluster[0][0]
    else:
       mutated_seq1 = cluster[0][1]
    if cluster[1][0].prot_is_different or (len(cluster[1]) == 1):
      mutated_seq2 = cluster[1][0]
    else:
      mutated_seq2 = cluster[1][1]
    # check for a indels
    if mutated_seq1.mutation.is_indel or mutated_seq2.mutation.is_indel:
      messages.append("indel and snv mutations in transcript (potentially in phase) " + this_transcript_info.id + ", cluster "+ str(i) +", skipping seq merging")
      cluster = _set_all_phase_info(cluster,"rejected:snv_indel")
      r += cluster
      continue
    mut1 = mutated_seq1.mutation.get_pos()
    mut2 = mutated_seq2.mutation.get_pos()
    assert(mut1 != mut2)
    if (not (mut1 in this_transcript_info.phase_dict)) and (not mut2 in this_transcript_info.phase_dict[mut1]):
      # noting to do; not in phase or phasing not observed
      r += cluster
      continue
    # not same codon:
    if mutated_seq1.next_mutationDIST != 0:
      if (not mutated_seq1.prot_is_different) or (not mutated_seq2.prot_is_different):
        # noting to do; zero or only one non-syn. effect
        r += cluster
        continue
      modify_mutated_sequence(mutated_seq1,mutated_seq2)
      modify_mutated_sequence(mutated_seq2,mutated_seq1)
      messages.append("two phased mutations in transcript " + this_transcript_info.id + ", cluster "+ str(i) +", merging effects of mutations "+ mutated_seq1.mutation.get_pos() +"/"+ mutated_seq2.mutation.get_pos())
      cluster = _set_all_phase_info(cluster,"in_phase:merged_two_codons")
      r += cluster
      continue
    # same codon:
    else:
      modify_mutated_sequence(mutated_seq1,mutated_seq2)
      mutated_seq1.next_mutationPHASE = "in_phase:merged_same_codon"
      messages.append("two phased mutations in transcript " + this_transcript_info.id + ", cluster "+ str(i) +", merging effects of mutations (which affect the same codon) "+ mutated_seq1.mutation.get_pos() +"/"+ mutated_seq2.mutation.get_pos())
      r.append([mutated_seq1]) # transcript info object is not valid any more, but this does not matter; the object is not needed anymore, too
  return r,messages

def check_sequence_bounds(seq,start,end):
  if end > len(seq):
    end = len(seq)
  if start < 0:
    start = 0
  return start,end

def reverse_complement(seq,reverse=True):
  ret = []
  for i in seq:
    ret.append(complement_table[i])
  if reverse:
    ret.reverse()
  return "".join(ret)

def translation(seq,include_stop=False):
  if len(seq) != 0:
    ret = "" 
    mod = len(seq) % 3
    re = ""
    if mod != 0:
      re = seq[-mod:]
      seq = seq[:-mod]
    for i in range(0,len(seq),3):
      codon = "".join(seq[i:i+3]).upper()
      ret = ret + codon_table.get(codon,"X")
    if not include_stop:
      if ret[-1] == stop_codon:
        ret = ret[:-1]
    return (ret,re)
  else:
    return ("","")

def get_mutation_transcript_annotation(x,i):
  y = zip(x.exons,x.exon_expression,x.transcripts,x.gene_symbols,x.gene_expression,x.transcript_strand)
  for j,t in enumerate(y):
    if t[2] == i:
      return t  

def find_peptides(prot_seq,unmutated_prot,ws,peptide=""):
  s = prot_seq
  y = ""
  return_l = []
  for i in range(0,(len(s)-ws+1),1):
    window = s[i:i+ws]
    if (window not in peptide) and peptide != "": 
      continue # skip all windows which are not in "peptide" aka the 27mer
    x = unmutated_prot.find(window)
    if x == -1:
      if y == "":
        y = window
        start = i
      else:
        y = y + window[-1]
    else:
      if y != "":
        return_l.append([y,start])
        y = ""
  if y != "":
    return_l.append([y,start])
  return return_l

def parse_mutations(mut_ref,seq,genomic_seq):
  ret = []
  for i in range(len(seq)):
    if seq[i].startswith("("):
      mut = seq[i].replace("(","").replace(")","").split("|")
      assert len(mut) == 4
      if int(mut[3]) == mut_ref:
        prefix = list(genomic_seq[0:i])
        suffix = list(genomic_seq[i+1:])
        middle = list(genomic_seq[i])
        for j in mut[0:2]:
          if len(j) == 0:
            pass 
          elif j[0] == "+": #insertion
            ret.append(( prefix + middle + list(j[1:]) + suffix, i))
          elif j[0] == "-": #deletion
            x = len(j[1:])
            x = x if x <= len(suffix) else len(suffix)
            if j[1:].upper() != "".join(suffix[:x]).upper():
              # may happen when a deletion affects the end of an exon AND the following intron, i.e. crosses an exon-intron boundary
              # not sure what to do here; currenty this is ignored
              continue
            assert j[1:].upper() == "".join(suffix[:x]).upper()
            ret.append(( prefix + list(mut[2]) + suffix[x:], i))
          elif j[0].upper() in 'ACGT':
            # add double mutation here
            ret.append(( prefix + list(j[0].upper()) + suffix, i))
          elif j[0] == "*":
            ret.append(( prefix + list(mut[2]) + suffix, i))
        break
  return ret

def find_stop(pep,ws=None):
  stop_pos = pep.find(stop_codon)
  r = None
  if stop_pos == -1:
    r = pep
  else:
    short_seq = pep[:stop_pos]
    if ws == None:
       r = short_seq
    elif len(short_seq) >= ws:
       r = short_seq
  return r

def get_min_mhc(score_list,id,pos):
  idx = -1
  l = []
  for i,s in enumerate(score_list):
    try:
      p = id[i].index(pos)
      l.append((s[0],s[1],s[2][p],s[3][p],s[4][p]))
    except ValueError:
      pass
  m = 10000000
  for i,s in enumerate(l):
    x = s[2]
    if x <= m:
      idx = i
      m = x
  if idx == -1:
    r = (0,"",999.0,"","")
  else:
    r = l[idx]
  return r,l

def format_list(l,sep,human_readable=False,fl=False,dec="."):
  if human_readable:
    x = []
    for i in l:
      if not (i in x):
        x.append(i)
  else:
    x=l
  if fl:
    x = [format_float(y,dec) for y in x]
  else:
    x = [str(y) for y in x]
  return sep.join(x)
 
def format_float(f,dec):
  return str(f).replace(".",dec)

def format_mhc(l,sep,dec):
  x = l
  if len(x) == 0:  # six empty fields
    x = sep * 5
  else:
    if isinstance(x[5],list) and len(x[5]) != 0:
      x[5] = x[5][0]
    #if isinstance(x[6],list) and len(x[6]) != 0:  
    #  x[6] = x[6][0]
    x = format_list(l,sep,False,True,dec)
  return x
  

def format_line(i,sep="\t",subsep=" ",dec=".",human_readable=False):
  slist = (   i.get_pos(),                   ## fixed position
              i.ref_base,
              i.first_allele,
              i.second_allele,
              i.classification,
              i.is_indel,
              format_list(i.frameshift,subsep,human_readable),
              i.synonymous,
              i.is_in_utr,
              format_list(i.aa_change,subsep,human_readable),
              format_list(i.transcript_strand,subsep,human_readable),
              format_list(i.gene_symbols,subsep,human_readable),
              format_list(i.transcripts,subsep,human_readable),
              format_list(i.mRNA_environment,subsep,human_readable),
              format_list(i.mRNA_environment_u,subsep,human_readable),
              format_list(i.exons,subsep,human_readable),
              format_float(i.avg_expression,dec),
              format_list(i.gene_expression,subsep,human_readable,True,dec),
              format_list(i.exon_expression,subsep,False,True,dec),
              format_list(i.in_repeat,subsep,human_readable),
              i.dbsnp,
              i.dbsnp_val,
              i.chr,
              i.pos + 1,
              i.genomic_environment,                ## fixed position
              i.filename,
              len(i.filename.split(",")) 
              )
  slist = list(slist)
  # fields only in the non human readable output:
  if not human_readable:
    pass
  s = sep.join([str(x) for x in slist])
  return s

def get_header(human_readable):
  h = ["position","ref","alt1","alt2","classification","indel","frameshift","synonymous","UTR","amino_acid_substitution","strand","gene","transcript","mRNA_environment(indel_only)","mRNA_environment_WT(indel_only)","exons","mean_expression","gene_expression","exon_expression","repeat","dbsnp","dbsnp_validation_level","chr","pos_(one_based)","+-300nt","dataset","number_of_datasets"]
  if not human_readable:
    pass
  return h





