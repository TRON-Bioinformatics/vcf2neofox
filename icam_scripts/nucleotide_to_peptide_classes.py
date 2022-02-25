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

import math
import operator

from nucleotide_to_peptide_functions import *
from nucleotide_to_peptide_constants import stop_codon

class mutated_sequence(object):
  """represents a 2x'xd'+1-mer (normally a 27mer) peptide which is centered around a mutation position; also holds the inforamtion for a TROMPApep data line"""
  def __init__(self,seq,pos,mut,seq_u,len_5p,len_3p,strand,xd = 13, xdi = 15, store_mRNA = False, protein_export = False):
    self.mutated_sequence = seq
    self.wt_sequence = seq_u
    self.nucl_pos = pos
    self.mutation = mut
    self.transcript_info = None
    self.x_difference = xd
    if self.mutation.is_indel:
      self.x_difference = xdi
    self.store_mRNA = store_mRNA    
    self.protein_export = protein_export
    self.len_5p = len_5p
    self.len_3p = len_3p
    self.strand = strand

    # _u alwas denotes the unmutated sequence
    self.peptide = ""  # xd- or xdi-mer
    self.peptide_u = ""
    self.mRNA = ""  # includes UTR
    self.mRNA_u = ""
    self.protein = ""  # full length translation
    self.protein_u = ""
    self.CDS = ""  # w/o includes UTRs
    self.CDS_u = ""
    self.peptide_CDS = "" # corresponding RNA to self.peptide
    self.peptide_CDS_u = ""
    
    self.stop_gain = False
    self.stop_loss = False
    self.uniqueID = ""
    self.exon = ""
    self.exon_expression = 0.0
    self.prot_pos = -1
    self.CDS_nucl_pos = -1
    self.aa_change = ""
    self.mhcI = []
    self.mhcII = []
    self.other_peptides_from_this_mutation = 0
    self.next_mutationID = ""
    self.next_mutationDATASET = ""
    self.next_mutationDIST = ""
    self.next_mutationPHASE = ""

    self._init_seqs(len_5p,len_3p,strand)
    self.prot_is_different = self.protein != self.protein_u
    
  def _init_seqs(self,len_5p,len_3p,strand):
    
    self.mRNA = self.mutated_sequence
    self.mRNA_u = self.wt_sequence  
    # handle strand and calc mutation position
    if strand == "+":  
      self.CDS_nucl_pos = self.nucl_pos - len_5p
    elif strand == "-":
      len_5p,len_3p = len_3p,len_5p
      self.mRNA = reverse_complement(self.mRNA)
      self.mRNA_u = reverse_complement(self.mRNA_u)
      self.CDS_nucl_pos = len(self.mRNA) - (self.nucl_pos + 1)
      self.CDS_nucl_pos = self.CDS_nucl_pos - len_5p
      if self.mutation.is_indel:  # indels need special treatment - they modify the length of the transcript, so the "mutation position" has to be adapted accordingly
        m = ""
        if self.mutation.ref_base == self.mutation.first_allele:
          m = self.mutation.second_allele
        else:
          m = self.mutation.first_allele
        if "+" in m:
          self.CDS_nucl_pos = self.CDS_nucl_pos - len(m) # insertion
        else:
          self.CDS_nucl_pos = self.CDS_nucl_pos - 1 # deletion
# NIY --> double nucleotide variant
#      elif self.mutation.is_dnv:
#        self.CDS_nucl_pos = self.CDS_nucl_pos - 1
    self.prot_pos = int(math.floor(self.CDS_nucl_pos/3.0))

   
    # translate
    protein,re = translation(self.mRNA[len_5p:],True)
    protein_u,re_u = translation(self.mRNA_u[len_5p:],True)
 
    # find stop
    # NOTE: genes with internal stop codons in CDS are truncated this way - not sure if this is a wanted behaviour; GPX1 is an example for such a protein
    self.protein = find_stop(protein)
    self.protein_u = find_stop(protein_u)
    
    # calc CDS
    self.CDS = self.mRNA[len_5p:(len_5p + (len(self.protein) * 3))]
    self.CDS_u = self.mRNA_u[len_5p:(len_5p + (len(self.protein_u) * 3))]

    # cut peptides and calculate impact
    if self.mutation.is_indel:
      self.peptide,self.peptide_u = self._cut_indel()
      if self.peptide != self.peptide_u:
        self.mutation.synonymous = False
        if not ((self.CDS_nucl_pos % 3 == 2) and (len(self.mutation.first_allele[1:]) % 3 == 0)):
          self.mutation.frameshift.append(True)
        else:
          self.mutation.frameshift.append(False)   
    else:
      self.stop_loss = len(self.protein) > len(self.protein_u)
      self.peptide,self.peptide_u = self._cut_snv(self.stop_loss)
      if self.peptide != self.peptide_u:
        self.mutation.synonymous = False
        if self.stop_loss:
          self.aa_change = stop_codon + str(self.prot_pos + 1) + self.protein[self.prot_pos]
        elif (len(self.protein) < len(self.protein_u)): # stop gain
            self.aa_change = self.protein_u[self.prot_pos] + str(self.prot_pos + 1) + stop_codon
            self.stop_gain = True
        else:
          self.aa_change = self.protein_u[self.prot_pos] + str(self.prot_pos + 1) + self.protein[self.prot_pos]
        self.mutation.aa_change.append(self.aa_change)
    
    # get peptide_CDS
    x_cds_start = (self.prot_pos * 3) - (self.x_difference * 3)
    x_cds_start = x_cds_start if x_cds_start > 0 else 0
    self.peptide_CDS = self.CDS[x_cds_start:(x_cds_start + (len(self.peptide)*3))]
    self.peptide_CDS_u = self.CDS_u[x_cds_start:(x_cds_start + (len(self.peptide_u)*3))]

    
    # mRNA_env
    if (self.mutation.is_indel and self.peptide != self.peptide_u) or self.store_mRNA:
      p = self.CDS_nucl_pos + len_5p # need UTRs (which is removed in self.CDS) to have full sequence, therefore use sequence from self.mRNA
      mRNA_lower_upper = self._safe_slice_upper_lower(self.mRNA,p)
      mRNA_lower_upper_u = self._safe_slice_upper_lower(self.mRNA_u,p)
      self.mutation.mRNA_environment.append(mRNA_lower_upper)
      self.mutation.mRNA_environment_u.append(mRNA_lower_upper_u)
    if self.protein_export: # hack to export full length mutated/WT protein sequence
      p = self.prot_pos
      protein_lower_upper = self._safe_slice_upper_lower(self.protein,p)
      protein_lower_upper_u= self._safe_slice_upper_lower(self.protein_u,p)
      self.mutation.mRNA_environment.append(protein_lower_upper)
      self.mutation.mRNA_environment_u.append(protein_lower_upper_u)
    
    if not self.mutation.somatic:
      # germline variant - all WT elements are equal to the mutant seqs
      self.peptide_CDS_u = self.peptide_CDS
      self.peptide_u = self.peptide
      self.CDS_u = self.CDS
      self.mRNA_u = self.mRNA
      self.protein_u = self.protein
      self.wt_sequence = self.mutated_sequence

  def _safe_slice_upper_lower(self,seq,p):
     """returns the string 'seq' with all letters in lowercase, except the character at index 'p', which is uppercase"""
     if p > (len(seq)-1):
       return seq.lower()
     x = seq[:p].lower() + seq[p].upper() + seq[(p+1):].lower()
     return x
          
  def _cut_indel(self):
    x_start = self.prot_pos - self.x_difference
    y_end = 0
    # sometimes there are 2*n bases inserted/deleted, - this might result in a resume of the old ORF
    # this is handled here (ther might be cleaner ways, though)
    for i in range(-1,-len(self.protein)-1,-1):
      try:
        if self.protein[i] != self.protein_u[i]:
          y_end = i + 1 + self.x_difference
          break
      except:
        break
    l1 = y_end if y_end < 0 else len(self.protein)
    peptide = self.protein[x_start:l1]
    l2 = y_end if y_end < 0 else len(self.protein_u)
    peptide_u = self.protein_u[x_start:l2]
    return peptide,peptide_u 
  
  def _cut_snv(self,stop_loss):
    x_start = self.prot_pos - self.x_difference
    if stop_loss:
      #y_end = len(self.protein)
      y_end = self.prot_pos + (self.x_difference + 1) # changed this so SNVs always result in only up to 27 AA
    else:
      y_end = self.prot_pos + (self.x_difference + 1)
    start,end = check_sequence_bounds(self.protein,x_start,y_end)
    peptide = self.protein[start:end]
    start,end = check_sequence_bounds(self.protein_u,x_start,y_end)
    peptide_u = self.protein_u[start:end]
    return peptide,peptide_u 
  
  def find_next_mutation(self):
    if len(self.transcript_info.mutid) > 1:
      self.next_mutationDIST = 10000000
      dist = [abs(x-self.prot_pos) for x in self.transcript_info.aapos]
      for i,j in enumerate(dist):
        if j < self.next_mutationDIST:
          if self.transcript_info.mutid[i] != self.mutation.get_pos():
            if (not self.transcript_info.syn[i]) or (j == 0): # both are non-synonymous or in the same codon   --- check this condition ???
              self.next_mutationDIST = j
              self.next_mutationID = self.transcript_info.mutid[i]
              self.next_mutationDATASET = self.transcript_info.dataset[i]
    if self.next_mutationID == "":
      self.next_mutation = ""

  def format_line(self,sep="\t",dec="."):
    """returns representaion of object instance as a TROMPApep data line, with the fields separated by 'sep' and the numbers represented using the descimal separator 'dec'"""
    slist = [ self.uniqueID,
              self.mutation.get_pos(),
              self.transcript_info.gene,
              self.transcript_info.refseq,
              self.transcript_info.id,
              format_float(self.transcript_info.expression,dec),
              self.exon,
              format_float(self.exon_expression,dec),
              self.CDS_nucl_pos,
              self.prot_pos,   # maybe use for CDS length ?
              self.aa_change,
              self.peptide,
              self.peptide_u,
              self.peptide_CDS,
              format_mhc(self.mhcI,sep,dec),
              format_mhc(self.mhcII,sep,dec),
              len(self.transcript_info.mutid),
              self.next_mutationDIST,
              self.next_mutationID,
              self.next_mutationDATASET,
              self.other_peptides_from_this_mutation ]
    if self.transcript_info.phase_check_done:
              slist += [self.next_mutationPHASE]
    s = sep.join([str(x) for x in slist])
    return s
  def header(self,sep):
    """returns the TROMPApep file header line, with the fields separated by 'sep'"""
    # column order is fixed for columns 1(key) - 14(mRNA_for_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL))
    h =      ["key",
              "mutation",
              "gene",
              "RefSeq_transcript",
              "UCSC_transcript",
              "transcript_expression",
              "exon",
              "exon_expression",
              "transcript_position",
              "codon",
              "substitution",
              "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
              "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
              "mRNA_for_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
              "MHC_I_peptide_length_(best_prediction)",
              "MHC_I_allele_(best_prediction)",
              "MHC_I_score_(best_prediction)",
              "MHC_I_epitope_(best_prediction)",
              "MHC_I_epitope_(WT)",
              "MHC_I_score_(WT)",
              "MHC_II_peptide_length_(best_prediction)",
              "MHC_II_allele_(best_prediction)",
              "MHC_II_score_(best_prediction)",
              "MHC_II_epitope_(best_prediction)",
              "MHC_II_epitope_(WT)",
              "MHC_II_score_(WT)",
              "mutations_in_transcript",
              "distance_to_next_mutation(AA_residues)",
              "next_mutation(potential_to_change_27mer)",
              "next_mutation_source",
              "peptide_count_for_this_mutation_in_this_transcript"]
    if self.transcript_info.phase_check_done:
              h += ["phase_of_next_mutation"]
    return sep.join([str(x) for x in h])

class transcript_info(object):
  """represents all mutations in a single transcript"""
  def __init__(self,id,gene):
    """
    Constructor arguments:
    id - transcript ID
    gene - list of tuples, each tuple contains a gene ID and an alternative transcript ID (e.g. Refseq ID)
    """ 
    self.id=id
    g = zip(*gene)
    self.gene=" ".join(g[0])
    self.refseq=" ".join(g[1])
    self.ref=[]
    self.alt1=[]
    self.alt2=[]
    self.cdspos = []
    self.aapos = []
    self.mutid = []
    self.expression = 0.0
    self.dataset = []
    self.syn = []
    self.phase_dict = {}
    self.phase_info = []
    self.phase_check_done = False
    self.mutations = []
  def add_sequence(self,s):
    self.ref.append(s.mutation.ref_base)
    self.alt1.append(s.mutation.first_allele)
    self.alt2.append(s.mutation.second_allele)
    self.cdspos.append(s.CDS_nucl_pos)
    self.aapos.append(s.prot_pos)
    self.mutid.append(s.mutation.get_pos())
    self.dataset.append(s.mutation.filename)
    self.syn.append(s.mutation.synonymous)
    self.mutations.append(s.mutation)
  def order(self):
    if len(self.mutid) > 0:
      l = list(set(zip(self.ref,self.alt1,self.alt2,self.cdspos,self.aapos,self.mutid,self.dataset,self.syn,self.mutations)))
      self.ref,self.alt1,self.alt2,self.cdspos,self.aapos,self.mutid,self.dataset,self.syn,self.mutations = zip(*sorted(l,key=operator.itemgetter(3)))


class mutation_filter(object):
  """represents a filter which can be applied to a TROMPApep line (which has been read into a list)"""
  def __cEq(self,a,b):
    return a == b
  def __cLt(self,a,b):
    return a > b
  def __cSt(self,a,b):
    return a < b
  def __cLet(self,a,b):
    return a >= b
  def __cSet(self,a,b):
    return a<= b
  def __cNe(self,a,b):
    return a != b
  def __cNotIn(self,a,b):
    return b not in a
  def __cNotAtEnd(self,a,b):
    return not a.endswith(b)

  def __init__(self,index,header,converter,method,comparison):
    """ init a filter object, the filter will later (via the method apply) apply the 'method' to the argument of apply at the 'index' and the value of 'converter(comparison)' """
    self.comparisons = {"equals":self.__cEq,
                        "larger_than":self.__cLt,
                        "smaller_than":self.__cSt,
                        "larger_or_equal_than":self.__cLet,
                        "smaller_or_equal_than":self.__cSet,
                        "is_not_equal_to":self.__cNe,
                        "does_not_contain":self.__cNotIn,
                        "does_not_end_with":self.__cNotAtEnd}
    self.compare_function = self.comparisons[method]
    self.method_string = method
    self.converter = converter
    self.target_value = comparison
    self.idx = index
    self.id = header[self.idx]
    self.passed = 0
    
    
  def __test(self,a):
    r = self.compare_function(self.converter(a),self.target_value)
    if r:
      self.passed += 1
    return r

  def apply(self,words):
    if len(words) == 0:
      return words
    if self.__test(words[self.idx]):
      return words
    else:
      return []

  def get_string(self):
    t = str(self.target_value)
    if t == "":
      t = "empty_string"
    return self.id + "_" + self.method_string + "_" + str(self.target_value) + ": " + str(self.passed)





