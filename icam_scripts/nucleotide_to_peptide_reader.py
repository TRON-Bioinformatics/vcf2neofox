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

import operator

from mutation import *

#################################################### Mutation calls & FDR

def read_mutations(infiles,data_set_names):
  d = {}
  if len(data_set_names) == (len(infiles)/2):
    names = data_set_names
  else:
    names = []
    for i in range(0,len(infiles),2):
      names.append(infiles[i])
  pos = set()
  for idx,file_no in enumerate(range(0,len(infiles),2)):
    fn = infiles[file_no]
    fn_format = infiles[file_no + 1]
    fn_format = fn_format.split(",")
    somatic = True
    if len(fn_format) > 1:
      somatic = fn_format[1].lower() != "germline"
    fn_format = fn_format[0]
    mutation_object_list = get_formats()[fn_format](fn)
    for jdx,i in enumerate(mutation_object_list):
      #check coordinate collision
      if i.get_pos() in pos:
        for j in d[i.chr]:
          if j.pos == i.pos:
            j.filename = j.filename + "," + names[idx]
            break
        continue
      pos.add(i.get_pos())
      # add some class attributes
      i.transcripts = []
      i.exons = []
      i.genomic_environment = ""
      i.gene_symbols = []
      i.gene_expression = []
      i.exon_expression = []
      i.transcript_length = []
      i.mRNA_environment_u = []
      i.mRNA_environment = []
      i.transcript_strand = []
      i.filename = names[idx]
      i.is_in_utr = False
      i.in_repeat = []
      i.dbsnp = ""
      i.dbsnp_val = ""
      i.frameshift = []
      i.aa_change = []
      i.classification = []
      i.synonymous = True
      i.somatic = somatic
      # append muation to mutation dictionary
      l = d.get(i.chr,[])
      l.append(i)
      d[i.chr] = l
 
  x = 0
  for i in d.keys():
    d[i].sort(key=operator.attrgetter('pos'))
    x += len(d[i])
  
  return d,x

#################################################### dbSNP

def add_dbSNP_information(d,dbSNP_file):
  d2 = {}
  for i in d.keys():
    d2[i] = [x.pos for x in d[i]]
  with open(dbSNP_file) as f:
    for line in f:
      if not line[0] == "#":
        words = line[:-1].split("\t")
        c,p1 = words[1:3]
        if c in d2:
          pi = int(p1)
          if (pi in d2[c]):
            d[c][d2[c].index(pi)].dbsnp = words[4]
            d[c][d2[c].index(pi)].dbsnp_val = words[12]
  return d

##################################################### Repeats

def add_repeat_information(d,repeat_bed):
  d_repeats = {}
  with open(repeat_bed,"r") as f:
    for line in f:
      words = line[:-1].split("\t")
      l = d_repeats.get(words[0],[])
      l.append((int(words[1]),int(words[2]),words[3]))
      d_repeats[words[0]] = l
  
  for i in d_repeats.keys():
    d_repeats[i].sort(key=operator.itemgetter(0))
  
  for chr in d_repeats.keys():
    try:
      l_rep = d_repeats[chr]
      l_mut = d[chr]
      mp = 0
      for ep in range(len(l_rep)):
        while mp < len(l_mut):
          if (l_mut[mp].pos >= l_rep[ep][0]) and (l_mut[mp].pos < l_rep[ep][1]): # this mutation in in this repeat, save it and get the next mutation
            l_mut[mp].in_repeat.append(l_rep[ep][2])
            mp += 1
          elif (l_mut[mp].pos < l_rep[ep][0]): # this mutation is before this repeat, get the next mutation
            mp += 1
          else:  # there are no mutaions before this exon or in this repeat, get the next repeat
            break
    except KeyError as e:
      pass
  return d 





##################################################### Expression data

def getMedian(numericValues):
  theValues = sorted(numericValues)
  count = len(theValues)
  try:
    if len(numericValues) == 1:
      return numericValues[0]
    if count % 2 == 1:
      return theValues[(count+1)/2-1]
    else:
      lower = theValues[count/2-1]
      upper = theValues[count/2]
      return (float(lower + upper)) / 2
  except:
    return -1.0

def read_expression_data(infile,sep="\t"):
  expr_table = {}
  if sep == "\\t":
    sep = "\t"
  elif sep == "":
    sep = "\t"
  with open(infile,"r") as f:
    for line in f:
      if line[0:2] != "ID":
        words = line[:-1].split(sep)
        expr_table[words[0]] = getMedian([float(x.replace(",",".")) for x in words[1:]])
  return expr_table


##################################################### Gen symbols


def read_gs(gensymbols,refseq_only):
  col2=4
  rcol1=5
  ucol1=0
  gene_symbol_table = {}
  with open(gensymbols,"r") as f:
    for line in f:
      if line[0] != "#":
        words = line[:-1].split("\t")
        # REFSEQ -> (gene,refseq)
        if (words[rcol1] != ""):
          l = gene_symbol_table.get(words[rcol1],[])
          l.append((words[col2],words[rcol1]))
          gene_symbol_table[words[rcol1]] = l
        # UCSC -> (gene,refseq)
        l = gene_symbol_table.get(words[ucol1],[])
        if refseq_only and (words[rcol1] != ""):
          l.append((words[col2],words[rcol1]))
          gene_symbol_table[words[ucol1]] = l
        elif not refseq_only:
          l.append((words[col2],words[rcol1]))
          gene_symbol_table[words[ucol1]] = l
  return gene_symbol_table


##################################################### BED file handling
def __convert_bed_line_to_exons(words):
  exon_list = []
  (chrom,start,end,name) = words[0:4]
  start = int(start)
  strand = words[5]
  exons = int(words[9])
  exon_starts = words[11].strip(",").split(",")
  exon_sizes = words[10].strip(",").split(",")
  exon_startsI = [start + int(x) for x in exon_starts]
  exon_sizesI = [int(x) for x in exon_sizes]
  r = range(exons)
  if strand == "-":
    r = range(exons-1,-1,-1)
  for ne in r:
    abs_e_end = exon_startsI[ne] + exon_sizesI[ne]
    exon_number_correct_order = (exons - 1 - ne) if strand == "-" else ne
    exon_words = [chrom,exon_startsI[ne],abs_e_end,name + "#exon." + str(exon_number_correct_order),strand,name]
    exon_list.append(exon_words)
  return exon_list

def read_reference_transcripts(bedfile):
  exons = {}
  transcripts = {}
  numt = 0
  
  # read bed file
  with open(bedfile,"r") as f:
    for line in f:
      words = line[:-1].split("\t")
      # store bed entry
      l = transcripts.get(words[3],[])
      l.append(words)
      transcripts[words[3]] = l
      # convert to exons
      es = __convert_bed_line_to_exons(words)
      for e in es:
        l = exons.get(e[0],[])
        l.append(tuple(e[1:]))
        exons[words[0]] = l
      numt += 1
  
  nume = 0
  for i in exons.keys():
      exons[i].sort(key=operator.itemgetter(0))
      nume += len(exons[i])
  
  return transcripts,exons,numt,nume

def no_conversion(id,tab):
  return id

def ucsc_to_refseq_conversion(id,tab):
    return zip(*tab[id])[1][0]
  
##################################################### Muation assignment
def assign_mutations_to_transcipts(d_ref,d_ref_exons,d,gene_expression,exon_expression,expression_sep,gensymbols,splicing_margin = 50):
  d_exon_mutations = {}
  d_mut_props = {}
  d_ref_mutations = {}
  
  # read expression data
  gene_expr_table = {}
  exon_expr_table = {}
  if gene_expression != "":
    gene_expr_table = read_expression_data(gene_expression,expression_sep)
  if exon_expression != "":
    exon_expr_table = read_expression_data(exon_expression,expression_sep)


  #read gensymbols
  gene_symbol_table = read_gs(gensymbols,False)

  # assign mutations to exons and transcripts
  for chr in d_ref_exons.keys():
    try:
      l_exon = d_ref_exons[chr]
      l_mut = d[chr]
      d_mut_props[chr] = [[False,False,False]] * len(l_mut)
      for ep in range(len(l_exon)):
        for mp in range(len(l_mut)):
          if (l_mut[mp].pos >= l_exon[ep][0]) and (l_mut[mp].pos < l_exon[ep][1]):
            exonid = l_exon[ep][2]
            l_mut[mp].exons.append(exonid)
            
            l_mut[mp].exon_expression.append(exon_expr_table.get(exonid,0.0))

            l = d_exon_mutations.get(exonid,[])
            l.append(mp)
            d_exon_mutations[exonid] = l
            
            transcriptid = l_exon[ep][4]
            l_mut[mp].transcripts.append(transcriptid)
            try:
              gsl = gene_symbol_table[transcriptid]
              l_mut[mp].gene_symbols += zip(*gsl)[0]
            except KeyError:
              l_mut[mp].gene_symbols.append("-")
            l_mut[mp].gene_expression.append(gene_expr_table.get(transcriptid,0.0))
            l_mut[mp].transcript_strand.append(l_exon[ep][3])

            l = d_ref_mutations.get(transcriptid,[])
            l.append(mp)
            d_ref_mutations[transcriptid] = l

          elif (l_mut[mp].pos >= (l_exon[ep][0] - splicing_margin)) and (l_mut[mp].pos < (l_exon[ep][1] + splicing_margin)):
            d_mut_props[chr][mp][0] = True
    except KeyError as e:
      pass
  x = 0
  for i in d_exon_mutations.keys():
    x += len(d_exon_mutations[i])
  
  # calculate mutation properties
  l_transcript = d_ref.values()
  for chr in d.keys():
    l_mut = d[chr]
    for ep in range(len(l_transcript)):
      for t in l_transcript[ep]:
        c,s,st,utr_s,utr_st = t[0],int(t[1]),int(t[2]),int(t[6]),int(t[7])
        for mp in range(len(l_mut)):
          if chr == c:
            if (l_mut[mp].pos >= s) and (l_mut[mp].pos < st):
              d_mut_props[chr][mp][1] = True
            if ((l_mut[mp].pos >= s) and (l_mut[mp].pos < utr_s)) or ((l_mut[mp].pos >= utr_st) and (l_mut[mp].pos < st)):
              d_mut_props[chr][mp][2] = True

  for chr in d.keys():
    for mp in range(len(d[chr])):
      if chr in d_mut_props:
        if (d[chr][mp].exons != []):
          d[chr][mp].classification = "transcript"
        elif d_mut_props[chr][mp][0] and d_mut_props[chr][mp][1] and d_mut_props[chr][mp][2]:
          d[chr][mp].classification = "intron_splicing"
        elif d_mut_props[chr][mp][1]:
          d[chr][mp].classification = "intron"
        else:
          d[chr][mp].classification = "intergenic"
      else: # no transcripts on given chromosome (may happen e.g. for chrM)
        d[chr][mp].classification = "intergenic"


  return d_ref_mutations,x,len(d_exon_mutations),gene_symbol_table

  
