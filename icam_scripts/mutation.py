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
import argparse

hetero_snvs =  {"M":("A","C"), "R":("A","G"), "W":("A","T"), "S":("C","G"), "Y":("C","T"), "K":("G","T")}
transition = ["AG","CT"]
transversion = ["AC","AT","CG","GT"]

class mutation(object):
#example input
  """
chr1    67075871        67163158        NM_207014       0       -       67075923        67163102        0       10      196,203,195,156,140,157,113,185,175,226,      0,2868,9883,24546,33769,37180,53553,55628,67600,87061,
chr1    841620  A       G       78      78      35      17      GGGGGGgGgGGGGgGgG       2141-3/:4::/7+7.:
chrX    123456  *       -TTGA/-TTGA     47      111     37      3       -TTGA   *       2       1       0       0       0
  """
  chr = ""
  pos = 0
  ref_base = ""
  is_indel = False
  mut = []
  is_heterozygot = False
  consensus_quality = 0
  snp_quality = 0
  rms = 0
  num_reads = 0
  bases = ""
  base_qualities = ""
  first_allele = ""
  second_allele = ""
  reads_first = 0
  reads_second = 0
  reads_third = 0
  line = ""
  #debug = ""
  def __init__(self,line,base_position=1):
    words = line[:-1].split("\t")
    self.line = line
    self.chr = words[0]
    self.pos = int(words[1]) - base_position
    if (words[2] == "*") and (len(words) > 10):
      self.is_indel = True
      self.reads_first = int(words[10])
      self.reads_second = int(words[11])
      self.reads_third = int(words[12])
    else:
      self.is_indel = False
      self.ref_base = words[2].upper()
      self.bases = words[8]
      self.base_qualities = words[9]

    m = words[3].split("/")
    self.mut = words[3]
    if len(m) == 2:
      self.first_allele = m[0]
      if m[0] != m[1]:
        self.is_heterozygot = True
        self.second_allele = m[1]
    else:
      if hetero_snvs.has_key(m[0]):
        self.is_heterozygot = True
        a = hetero_snvs[m[0]]
        self.first_allele = a[0] if a[0] != self.ref_base else "*"
        self.second_allele = a[1] if a[1] != self.ref_base else "*"
      else:
        self.first_allele = m[0]

    self.set_transition_transversion()
    self.sort_alleles()

    self.consensus_quality = int(words[4])
    self.snp_quality = float(words[5])
    self.rms = int(words[6])
    self.num_reads = int(words[7])

  def set_transition_transversion(self):
    self.is_transition = False
    self.is_transversion = False
    t1 = [self.first_allele, self.ref_base]
    t2 = [self.second_allele, self.ref_base]
    t1.sort()
    t2.sort()
    t1 = "".join(t1).upper()
    t2 = "".join(t2).upper()
    if (t1 in transition) or (t2 in transition):
      self.is_transition = True
    elif (t1 in transversion) or (t2 in transversion):
      self.is_transversion = True

  def sort_alleles(self):
    x = [self.first_allele.upper(),self.second_allele.upper()]
    x.sort()
    self.first_allele = x[0]
    self.second_allele = x[1]
    if (self.first_allele == "*") or (self.first_allele == ""):
      self.first_allele,self.second_allele = self.second_allele,self.first_allele
    

  def get_pos(self):
    return "{0}_{1}".format(self.chr,self.pos)

  def __eq__(self,o):
    if ((self.first_allele.upper() == o.first_allele.upper()) and
        (self.second_allele.upper() == o.second_allele.upper()) and
        (self.pos == o.pos) and
        (self.chr.upper() == o.chr.upper())):
      return True
    return False

  def extended_eq(self,o):
    if self.is_heterozygot and o.is_heterozygot:
      return self == o
    if not ((self.pos == o.pos) and (self.chr.upper() == o.chr.upper())):
      return False
    sfa = hetero_snvs.get(self.first_allele.upper(),self.first_allele.upper())
    ssa = hetero_snvs.get(self.second_allele.upper(),self.second_allele.upper())
    ofa = hetero_snvs.get(o.first_allele.upper(),o.first_allele.upper())
    osa = hetero_snvs.get(o.second_allele.upper(),o.second_allele.upper())
    m1 = False
    m2 = False
    for i in sfa:
      if i in ofa:
        m1 = True
        break
    for i in ssa:
      if i in osa:
        m2 = True
        break
    return (m1 and m2)

  def __ne__(self,o):
    return not self.__eq__(o)

  def __hash__(self):
    return hash("{0}_{1}_{2}_{3}".format(self.chr.upper(),self.pos,self.first_allele.upper(),self.second_allele.upper()))

  def get_alleles(self):
    return [self.first_allele,self.second_allele]
# END class mutation

class jsm_mutation(mutation):
  def __init__(self,line,base_position=1):
    words = line[:-1].split("\t")
    self.line = line
    self.chr = words[0]
    self.pos = int(words[1]) - base_position
    self.ref_base = words[2]
    self.first_allele = words[3]
    self.second_allele = "*"
    s1 = float(words[9]) 
    s2 = float(words[10])
    self.snp_quality = -10 * log10(1-(s1+s2))
    self.num_reads = sum([int(x) for x in words[4:8]])
    self.is_heterozygot = s1 > s2
    if not self.is_heterozygot:
      self.second_allele = ""
    self.sort_alleles()
    self.set_transition_transversion()
# END class jsm_mutation

class cm_mutation(mutation):
  def __init__(self,line,base_position=1):
      words = line[:-1].split("\t")
      self.line = line
      self.chr = words[2]
      self.pos = int(words[3]) - base_position
      self.ref_base = words[20]
      self.first_allele = words[21]
      self.second_allele = "*"
      self.is_heterozygot = True
      self.sort_alleles()
      self.set_transition_transversion()
# END class cm_mutation

class mymut_mutation(mutation):
  def __init__(self,line,refp,altp,base_position=0):
      words = line[:-1].split("\t")
      CHR = "_".join(words[0].split("_")[:-1])
      POS = int(words[0].split("_")[-1]) - base_position
      REF = words[refp]
      ALT = words[altp]
      self.line = line
      self.chr = CHR
      self.pos = POS - base_position
      self.ref_base = REF
      self.first_allele = ALT
      self.second_allele = "*"
      self.is_heterozygot = True
      if (words[2] == words[3]) and (words[2] != "?") and (words[3] != "?"):
        self.second_allele = ""
        self.is_heterozygot = False
      self.sort_alleles()
      self.set_transition_transversion()
# END class mymut_mutation


class vcf_mutation(mutation):
  """
indel
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  b16
chr10   9555543 .       A       AT      .       PASS    AC=11,11;DP=20;MM=0.0,0.22222222;MQ=37.0,37.0;NQSBQ=35.372726,34.95402;NQSMM=0.0,0.0;SC=11,0,9,0     GT       0/1

snp
chr10   3035121 .       A       C       44.36   LowQual DP=3;Dels=0.00;HRun=1;HaplotypeScore=0.00;MQ=37.00;MQ0=0;QD=14.79;SB=-6.99      GT:DP:GL:GQ     0/1:3:-8.51,-0.90,-3.99:30.92
  """
  dbsnp_ref = "."
  filter = ""
  other_alleles = []
  info = {}
  genotype = {}
  def  __init__(self,line,column,base_position=1):
    words = line[:-1].split("\t")
    self.line = line
    self.chr = words[0]
    self.pos = int(words[1]) - base_position
    self.dbsnp_ref = words[2]
    self.ref_base = words[3]
    self.filter = words[6]
    self.snp_quality = float(words[5]) if words[5] != "." else 0.0
    items = words[7].split(";")
    for i in items:
      l = i.split("=")
      if len(l) == 1:
        key,data = l[0],[True]
      else:
        key,data = l
        data = data.split(",")
      self.info[key] = data
    if "DP" in self.info:
      self.num_reads = int(self.info["DP"][0])
    if "MQ" in self.info:
      self.rms = float(self.info["MQ"][0])
    if "QD" in self.info:
      self.qd = float(self.info["QD"][0])
    alt = words[4].split(",")
    allele_list = [self.ref_base] + alt
    if len(words) > 8:
      format = words[8].split(":")
      gt = words[9+column].split(":")
      items = zip(format,gt)
      for i in items:
        key,data = i
        self.genotype[key] = data
      if "GT" in self.genotype:
        gt = self.genotype["GT"].split("/")
        if len(gt) == 1:
          gt = self.genotype["GT"].split("|")
        if gt == ["."] or (len(gt) == 1 and gt[0] == 0):
          gt = [0,0]
        else:
          gt = [int(x) for x in gt]
        #try:
        self.first_allele = allele_list[gt[0]]
        self.second_allele = allele_list[gt[1]]
        #except Exception as e:
        #  print e
        #  print allele_list,gt,self.pos,self.chr
        #  sys.exit()
      else:
        self.first_allele = self.ref_base
        self.second_allele = alt[0]
    else:
      self.first_allele = self.ref_base
      self.second_allele = alt[0]
    if len(self.first_allele) > 1 or len(self.second_allele) > 1 or len(self.ref_base) > 1:
      self.is_indel = True
      # format to match samtools indel format
      self.ref_base = self.ref_base.split("/")[0]
      self.first_allele = self.first_allele.split("/")[0]
      self.second_allele = self.second_allele.split("/")[0]
      if len(self.first_allele) > len(self.ref_base):
        # insertion
        self.first_allele = "+" + self.first_allele[1:]
      elif len(self.first_allele) < len(self.ref_base):
        # deletion
        self.first_allele = "-" + self.ref_base[1:]
      if len(self.second_allele) > len(self.ref_base):
        self.second_allele = "+" + self.second_allele[1:]
      elif len(self.second_allele) < len(self.ref_base):
        self.second_allele = "-" + self.ref_base[1:]
    if self.first_allele != self.second_allele:
      self.is_heterozygot = True
      if self.second_allele == self.ref_base:
        self.second_allele = "*"
      if self.first_allele == self.ref_base:
        self.first_allele = "*"
    else:
      self.second_allele = ""
    if self.is_indel:
      #self.ref_base = ""
      self.ref_base = self.ref_base[0]
    self.set_transition_transversion()
    self.sort_alleles()
# END class


class somatic_mutation(mutation):
  """
chr10   3039319 c       Y       3       1       1       37      3       2
chr10   3039418 A       G       8       30      30      37      1       2
  """
  somatic_score = -1
  tumor_consensus_quality = -1
  tumor_snv_quality = -1
  tumor_rms_mapping_quality = -1
  tumor_depth = -1
  normal_depth = -1
  def __init__(self,line,base_position=1):
    words = line[:-1].split("\t")
    self.line = line
    self.chr = words[0]
    self.pos = int(words[1]) - base_position
    self.ref_base = words[2].upper()
    if hetero_snvs.has_key(words[3]):
      self.is_heterozygot = True
      a = hetero_snvs[words[3]]
      self.first_allele = a[0] if a[0] != self.ref_base else "*"
      self.second_allele = a[1] if a[1] != self.ref_base else "*"
    else:
      self.first_allele = words[3]
    self.somatic_score,self.tumor_consensus_quality,self.tumor_snv_quality,self.tumor_rms_mapping_quality = [float(x) for x in words[4:8]]
    self.tumor_depth,self.normal_depth = [int(x) for x in words[8:10]]
    self.set_transition_transversion()
    self.sort_alleles()
    pass
# END CLASS

class maf_mutation(mutation):
  """A1BG    1       broad.mit.edu   37      19      58858802        58858802        +       Missense_Mutation       SNP     C       C       T """
  """C17orf104       284071  broad.mit.edu   37      17      42746783        42746783        +       Frame_Shift_Del DEL     A       A       -"""
  """GTF2IRD2P1      0       broad.mit.edu   37      7       72685632        72685633        +       RNA     INS     -       -       TAC"""
  def __init__(self,line,base_position=1):
    words = line[:-1].split("\t")
    self.line = line
    self.chr = "chr" + words[4]
    self.pos = int(words[5]) - base_position
    self.first_allele = words[11]
    self.second_allele = words[12]
    self.ref_base = words[10]
    if words[9] == "INS": # insertion
      self.is_indel = True
      self.pos = self.pos - 1 # point to nt before indel
      self.ref_base = ""
      if self.first_allele == "-":
        self.first_allele = ""
      else:
        self.first_allele = "+" + self.first_allele
      if self.second_allele == "-":
        self.second_allele = ""
      else:
        self.second_allele = "+" + self.second_allele
    elif words[9] == "DEL": # deletion
      self.is_indel = True
      self.pos = self.pos - 1 # point to nt before indel
      if self.first_allele == "-":
        self.first_allele = "-" + self.ref_base
      elif self.first_allele == self.ref_base:
        self.first_allele = "*"
      if self.second_allele == "-":
        self.second_allele = "-" + self.ref_base
      elif self.second_allele == self.ref_base:
        self.second_allele = "*"   
      self.ref_base = ""
    if self.first_allele == self.ref_base:
      self.first_allele = "*"
    if self.second_allele == self.ref_base:
      self.second_allele = "*"
    self.set_transition_transversion()
    self.sort_alleles()

class generic_mutation(object):
  pos = 0
  chr = ""
  line = ""
  def __init__(self,line,base_position=1):
    words = line[:-1].split("\t")
    self.line = line
    self.chr = words[0]
    self.pos = int(words[1]) - base_position
  def get_pos(self):
    return "{0}_{1}".format(self.chr,self.pos)
# END class

class tab_mutation(mutation):
  def __init__(self,line,base_position=0):
    words = line[:-1].split("\t")
    self.line = line
    cp = words[0].split("_")
    self.chr = "_".join(cp[:-1])
    self.pos = int(cp[-1]) - base_position
    self.ref_base = words[1]
    self.first_allele = words[2]
    if len(words) > 3:
      self.second_allele = words[3]
    else:
      self.second_allele = "*"
    for s in "+-":
      for x in [self.first_allele,self.second_allele]:
        if s in x:
          self.is_indel = True
          break
    self.set_transition_transversion()
    self.sort_alleles()
# END class tab_mutation


def read_generic_mutations(infile,comment_char="#",base_position=1):
  l = []
  comment = []
  with open(infile,"r") as f:
    for line in f:
      if not (line.startswith(comment_char) or line.startswith("contig") or line.startswith("chrom")):
        if not line[:-1] == "":
          m = generic_mutation(line,base_position)
          l.append(m)
      else:
        comment.append(line[:-1])
  return l,comment

def read_generic_mutations2(infile):
  l,comment = read_generic_mutations(infile)
  return l

def read_jsm_mutations(infile,mint = 0.95):
  l = []
  with open(infile,"r") as f:
   for line in f:
     if line.startswith("chrom"):
       continue
     else:
       w = line[:-1].split("\t")
       if w[3] != "N":
         s = float(w[9]) + float(w[10])
         if s > mint:
           l.append(mutation(line))
  return l

def read_samtools_mutations(infile):
  l = []
  with open(infile,"r") as f:
   for line in f:
     words = line.split("\t")
     if words[2] == words[3]:
       continue
     l.append(mutation(line))
  return l

def read_somatic_sniper(infile):
  l = []
  with open(infile,"r") as f:
   for line in f:
     l.append(somatic_mutation(line))
  return l

def read_vcf(infile,column=0,filter=True):
  l = []
  with open(infile,"r") as f:
   for line in f:
     if not line.startswith('#'):
       m = vcf_mutation(line,column)
       if (m.filter == "PASS") & filter:
         l.append(m)
       if not filter:
         l.append(m)
  return l

def read_vcf1(infile):
  return read_vcf(infile,1,True)
def read_vcf_all(infile):
  return read_vcf(infile,0,False)
def read_vcf1_all(infile):
  return read_vcf(infile,1,False)

def read_table(infile):
  l = []
  with open(infile,"r") as f:
   for line in f:
     if (not line.startswith("position")) and (not line.startswith("Position")):
       l.append(tab_mutation(line))
  return l

def read_cm(infile):
  l = []
  with open(infile,"r") as f:
    r = False
    c = 0
    for line in f:
      c += 1
      if line.startswith("-----") and r:
        break
      if r:
        l.append(cm_mutation(line))
      if (not r) and (line.startswith("mut ID")):  #(c == 5):
        r = True
  return l

def read_mymut(infile):
  l = []
  with open(infile,"r") as f:
    r = False
    c = 0
    refp = -1
    altp = -1
    for line in f:
      c += 1
      if line.startswith("-----") and r:
        break
      if r:
        l.append(mymut_mutation(line,refp,altp))
      if (not r) and (line.startswith("chr_pos")):
        w = line.strip("\n").split("\t")
        refp = w.index("wt nt")
        altp = w.index("mut nt")
        r = True
  return l

def read_maf(infile):
  l = []
  with open(infile,"r") as f:
    for line in f:
      if line.startswith("#"):
        continue
      if line.startswith("Hugo_Symbol"):
        continue
      l.append(maf_mutation(line))
  return l

def get_formats():
  return {"samtools":read_samtools_mutations,"soma":read_somatic_sniper,"vcf":read_vcf,"generic":read_generic_mutations2,"jsm":read_jsm_mutations,"table":read_table,"cm":read_cm,"vcf1":read_vcf1,"maf":read_maf,"vcfa":read_vcf_all,"vcf1a":read_vcf1_all,"mymut":read_mymut}


def main():
  parser = argparse.ArgumentParser(description='Process mutation data in various formats.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('mutation_file', metavar='M', type=str,nargs="?", default=None, help='input file')
  parser.add_argument('format', metavar='F', type=str,nargs="?", default=None, help='file format')
  parser.add_argument('--vcf_genotype_column', dest='col', default=1, type=int, action='store', help='VCF only: which genotype column should be used')
  parser.add_argument('--vcf_passing_only', dest='filter', action='store_true', default=False,help='VCF only: only process entries with filter=PASS')
  parser.add_argument('--list_formats', dest='list', action='store_true', default=False,help='list available formats and exit')
  parser.add_argument('--silent', dest='silent', action='store_true', default=False,help='no output to stdout (overrides --output)')
  parser.add_argument('--output', dest='output', action='store_true', default=False,help='print entries, default=print number of entries')
  parser.add_argument('--simple_tab', dest='tab', action='store', default="",help='print entries as simple tab')
  args = parser.parse_args()

  if args.list:
    for k in get_formats().keys():
      print k
    sys.exit()
  else:
    if (args.mutation_file == None) or (args.format == None):
      parser.print_usage()
      print
      print "Not enough command line arguments!"
      sys.exit()

  
  f = args.mutation_file
  file_format = args.format
  l = []
  if file_format == "vcf":
    l = get_formats()[file_format](f, args.col - 1, args.filter)
  else:
    l = get_formats()[file_format](f)
  if not args.silent:
    if not args.output:
      print len(l)
    else:
      for i in l:
        print i.line[:-1]

  if args.tab != "":
    fo = open(args.tab,"w")
    fo.write("position"+ "\t" + "ref_base" + "\t" + "first_allele"+ "\t" + "second_allele" + "\n")
    for i in l:
      fo.write(i.get_pos() + "\t" + i.ref_base + "\t" + i.first_allele + "\t" + i.second_allele + "\t" + "\n")
    fo.close()


if __name__ == '__main__':
  main()
