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

# python imports
import sys
import math
import os
import operator
import time
import optparse

# package imports
import twobit
from mutation import *

# import constant data
from nucleotide_to_peptide_constants import *
# import classes
from nucleotide_to_peptide_classes import *
# import reading functions
from nucleotide_to_peptide_reader import *
# import other functions
from nucleotide_to_peptide_functions import *

class ntp_obj(object):
  """
  class representing an analysis run - this program uses a mutation file (or many) as input and returns the mutation effects
  - goes from nucleotide level information to peptide level information (thus the creative name ...)
  - handles small indels and SNVs
  - germline vs. somatic mutation distinction is built in
  - limited phasing (via short read pairs in BAM files) and correction of peptide effects (two SNVs can effect the same peptide)
  - associates gene and transcript with the peptides
  - associated repeating elements and known SNPs with the mutations/peptides
  - associates expression values with the mutations/peptides
  - interfaces the IEDB tools to do the MHC binding prediction
  - produces TROMPApep/mut files
  NOTE: having this as a class in kind of useless; could have been a function, too ...
  NOTE 2: please excuse the spaghetti code :)
  NOTE 3: here are some very old (and odd) legacy pieces of code - it works fine though and we all know that "Premature optimization is the root of all evil" -- DonaldKnuth
  """
  def __init__(self,cfg_dict,logger=None):
    """
    constructor
    """
    self.configuration = cfg_dict
    self.peptide_output = self.configuration["mutation_set_output"] + ".transcript"
    self.phase_output = self.configuration["mutation_set_output"] + ".phase_report.txt"
    self.logger = logger
 
  def _loginfo(self,msg):
    if self.logger:
      self.logger.info("NTP " + msg)
    #else:
    #  print >> sys.stderr, msg
  def _logerror(self,msg):
    if self.logger:
      self.logger.error(msg)
  def _logdebug(self,msg):
    if self.logger:
      self.logger.debug(msg)

 
  def run(self,infiles):
    # manage configuration and translate to local variables
    mutation_set_output =         self.configuration["mutation_set_output"] # output file basename
    environment_size =            int(self.configuration["environment_size"]) # length of genomic environment of mutations, needed for primer design 
    peptide_len_snv =             int(math.floor((float(self.configuration.get("peptide_len_snv",27)) - 1.0) / 2.0)) # number of additional AA besides mutatt amino acids (N-terminal and C-terminal, for SNVs)
    peptide_len_indel =           int(math.floor((float(self.configuration.get("peptide_len_indel",31)) - 1.0) / 2.0)) # sam as above, but for INDELs
    ref_bed =                     self.configuration["ref_bed"] # bed file with referenc transcriptome coordinates
    repeat_bed =                  self.configuration["repeat_bed"] # bed file with repeating elements (e.g. Repeatmasker track from UCSC genome browser)
    genome_2bit =                 self.configuration["genome_2bit"] # reference genome in 2bit format (as provided by UCSC genome browser)
    gensymbols =                  self.configuration["xref_table"].split(",")[0] # kgxref table from UCSC genome browser
    gene_expression =             self.configuration["gene_expression"] # txt file with gene expression data
    exon_expression =             self.configuration["exon_expression"] # txt fiel with exon expression data
    organism =                    self.configuration["organism"] # organism for MHC prediction (human or mouse)
    mhc_I_selection =             filter(lambda x:len(x)!=0, [x.strip() for x in self.configuration["mhc_I_selection"].split(",")]) # comma separated list of class I HLA alleles
    mhc_II_selection =            filter(lambda x:len(x)!=0, [x.strip() for x in self.configuration["mhc_II_selection"].split(",")]) # comma separated list of class II HLA alleles
    dbsnp_file =                  self.configuration["dbsnp_file"] # dbsnp table from UCSC genome browser
    expression_sep =              self.configuration["expression_sep"] # filed separator used in gene/exon expression files
    data_set_names =              [x.strip() for x in self.configuration["data_set_names"].split(",")] # names of the input mutation files (comma sep. list)

    # number of parallel processes used for IEDB MHC binding prediction
    if self.configuration["nprocs_IEDB"].strip() != "":
      nprocs_IEDB = int(self.configuration["nprocs_IEDB"])
    # the mutation files with the input data can be supplied via the config file or the CLI
    if infiles == []:
      infiles =                   [x.strip() for x in self.configuration["infiles"].split(",")]
    #MHC binding prediction can be skipped - for this the skip_mhc key has to be set to an non-empty value
    skip_mhc = False
    if self.configuration["skip_mhc"].strip() != "":
      skip_mhc = True
    #the next two parameters allow to export full protein sequences or full length mRNA sequences resulting from the mutations
    store_mRNA = False
    if self.configuration.get("store_mRNA","").strip() != "":
      store_mRNA = True
    prot_export = False
    if self.configuration.get("prot_export","").strip() != "":
      prot_export = True
    # field separator for gene/transcript expression data files
    if expression_sep == "":
      expression_sep = "\t"   
    # optional list of BMA files for phasing
    bam_files = []
    if self.configuration.get("bam_files","").strip() != "":
      bam_files = self.configuration["bam_files"].split(",")
    # path to samtools binary, for phasing
    samtools_bin = self.configuration["samtools_bin"]
  
    # read mutation calls
    d,x = read_mutations(infiles,data_set_names)
    self._loginfo("mutation calls " +  str(x))
    d_set = set()
    for i in d: # put all mutation IDs (as provided by get_pos mthod of mutation class) in a set
      d_set.update([mmm.get_pos() for mmm in d[i]])
    if x != len(d_set): # check if mutations are unique (sanity check)
      self._logerror("\n Positions of input mutations not unique ({0} vs {1})! Exiting...".format(x,len(d_set)))
      for i in d:
        xxx = [mmm.get_pos() for mmm in d[i]]
        for ii in xxx:
          c = 0
          for jj in xxx:
            if jj == ii:
              c+=1
          if c > 1:
            self._logdebug(ii)
      sys.exit()    
    del d_set
    
    # check if mutations are in repeating elements, add this information to the mutation objects in dictionary d
    if repeat_bed != "":
      d = add_repeat_information(d,repeat_bed)
    # check if mutations are coverd by entries in dbSNP, add this information to the mutation objects in dictionary d
    if dbsnp_file != "":
      d = add_dbSNP_information(d,dbsnp_file)
    
    # read transcript coordinates from reference bed file, provide dictionary with transcript and exon coordinates
    # the dicts use the chromosome IDs as keys and point to lists of transcript or exon bed lines
    # x1 and x2 are the numbers of transcripts and exons, repsectively
    d_ref,d_ref_exons,x1,x2 = read_reference_transcripts(ref_bed)
    self._loginfo("reference exons " +  str(x2))
    self._loginfo("reference transcripts " +  str(x1))

    # go through exons and assign mutation indices
    # d_ref_mutations is a dictionary, pointing from transcrtipt ID (keys) to a list of muation objects (values), the latter represent the mutations in a given transcript
    # gene_symbol_table is a dictionary, pointing from transcrtipt ID (keys) to gene symbol (values)
    # x is the number of mutations in transcripts
    # x2 is the number of exons with mutations
    # for each mutation in d, the function also assigns if the mutation is in an intergenic region, close to a splice site (currently within 50nt), or otherwise in an intron or in the transcript
    d_ref_mutations,x,x2,gene_symbol_table = assign_mutations_to_transcipts(d_ref,d_ref_exons,d,gene_expression,exon_expression,expression_sep,gensymbols)
    del d_ref_exons # not needed anymore

    # for every transcript, get sequence, mutate, translate both
    self._loginfo("start tr loop")
    dna = twobit.dna_2bit(genome_2bit) # object to get genomic sequences from
    d_mut_peptides = {} # dictionary; keys = window sized (=potential epitope lengths); values are triples of (MPS_ID,mutant_peptide,wt_peptide)  [MPS is for mutant peptide sequence]; needed for MHC prediction
    mutated_sequences = [] # list of all mutated_sequence objects, linkes to d_mut_peptides via the MPS_ID (stored in uniqueID attribute of mutated_sequence class); needed for MHC prediction
    # the next three objects hold info on the phase of mutations; mutation pairs can affect multiple transcripts, but we do need to check the phase pf a pair just once
    phasing_dict = {}  # defined in check_pahse function in n*t*p*_functions.py
    phasing_report = [] # defined in check_pahse function in n*t*p*_functions.py
    phasing_info = {} # defined in check_pahse function in n*t*p*_functions.py
    for i in d_ref.keys(): # loop for each transcript
      if d_ref_mutations.has_key(i): # are there any mutations for this gene ? do nothing if no
        this_transcript_info = transcript_info(i,gene_symbol_table[i]) # init new info object
        # get info on transcrip, as given in the bed file
        words = d_ref[i][0]
        chr = words[0]  #chromosome
        start = int(words[1]) #5' end of transcript, relative to genome
        stop = int(words[2])  #3' end of transcript, relative to genome
        strand = words[5] # direction/strand
        lens = [int(x) for x in words[10][:-1].split(",")] # lengths of blocks (see bed file format)
        begins = [int(x) for x in words[11][:-1].split(",")] # releative starts of blocks (see bed file format)
        _5p_utr = int(words[6]) - start # start of CDS, relative to start position on genome
        _3p_utr = int(words[7]) - start # end of CDS, relative to start position on genome
        # see if transcript in non-coding; by convention the value of field seven is equal to the value of field 1 in a bed line
        non_coding = False
        if _3p_utr == 0:
          _3p_utr = -(stop - start)
          _5p_utr = 0 
          non_coding = True
        # get the whole sequence of the gene - mutations refer to genomic positions
        # using a list makes it easy to track the mutation position in the transcript when doing the virtual splicing
        # always get forward strand, build reverse complement later if needed
        try:
          seq = list(dna.get_seq_range_nucleotides(chr,start,stop-1, "+"))
        except Exception as e:
          self._loginfo("ERROR while accessing 2bit file")
          self._loginfo(e)
          self._loginfo("DNA " + chr + " " + str(start) + " " + str(stop-1))
        # store genomic sequence for later to compare to the mutated ones
        genomic_seq = "".join(seq) # the WT reference genome seq; no mutation in there, so it can remain a string
        # seq is now like [u,u,u,e,e,e,i,i,i,e,e,e,u,u,u] with u=UTR, e=CDS, i=inton
        # genomi_ seq is 'uuueeeiiieeeuuu'

        # get mutation objects
        mut = [] # this list holds all mutations in the trasncript
        for j in d_ref_mutations[i]:
          x = d[chr][j]
          #sanity check -- some refseqs have a non-uniqe name, which may cause trouble here!
          if x.pos >= start and x.pos < stop:
            mut.append(x)
        mut.sort(key=operator.attrgetter('pos')) # sort by genomic position on chromosome
        # another sanity check
        if len(mut) == 0:
          continue
           
        # add expression to transcript_info object; is the same for all mutations in mut
        this_transcript_info.expression = get_mutation_transcript_annotation(mut[0],i)[4]
        mut_utr = []  # list of mutation indices in mut of mutations which are in the UTR
        for mut_ref in range(len(mut)):
          j = mut[mut_ref]
          abs_pos = j.pos - start
          if (j.ref_base.upper() != seq[abs_pos].upper()) and (j.is_indel):
            j.ref_base = seq[abs_pos] # workaround for the old samtools indel format 
          seq[abs_pos] = "({0}|{1}|{2}|{3})".format(j.first_allele,j.second_allele,seq[abs_pos],mut_ref) # insert mutation info in sequence list, needed for "parse_mutations" function
          if (abs_pos < _5p_utr) or (abs_pos >= _3p_utr): # check if mutation is in an UTR
            j.is_in_utr = True
            mut_utr.append(mut_ref)
        # seq is now [u,u,u,e,(T,*,A,0),e,i,i,i,e,e,e,u,u,u] in the example case for a single mutation in the transcript
        # this single mutation is heterozygot, with a single new allele compared to the reference, is A->T and has index 0 in the list mut
        # genomic_seq is unchanged

        # here the introns are removed 
        # also UTR lengths are calculated
        # splice site covering indels are sorted out in the "parse_mutations" function
        seq_spliced = []   # with mutation = tumor
        genomic_seq_spliced = "" # w/o mutation = WT
        for j in range(int(words[9])): # iteratr ofer nubmer of bed blocks
           if _5p_utr >= begins[j] and _5p_utr < (begins[j] + lens[j]):
             _5p_utr = _5p_utr - begins[j] + sum([lens[xxx] for xxx in range(j)]) # position of UTR start on transcript
           if _3p_utr > begins[j] and _3p_utr <= (begins[j] + lens[j]):
             _3p_utr = _3p_utr - begins[j] + sum([lens[xxx] for xxx in range(j)]) # position of UTR start on transcript
           seq_spliced = seq_spliced +  seq[begins[j]:(begins[j] + lens[j])]
           genomic_seq_spliced = genomic_seq_spliced + genomic_seq[begins[j]:(begins[j] + lens[j])]
        rel_3p_utr = len(genomic_seq_spliced) - _3p_utr # length of 3' UTR
        assert len(seq_spliced) == len(genomic_seq_spliced) # same length, remember that a mutation currently is just a single list member, adding one to the length; also true for indels, length changes here are caculated later
        # seq is now [u,u,u,e,(T,*,A,0),e,e,e,e,u,u,u
        # genomic seq is 'uuueeeeeeuuu'

        # create all mutated_sequnece objects for all mutations of this transcript and store in a list
        all_seqs_mutated = []
        for mut_ref,this_mutation in enumerate(mut):
          if (mut_ref in mut_utr): # skip if in UTR in this transcript
            continue
          seqs_mutated_tmp = parse_mutations(mut_ref,seq_spliced,list(genomic_seq_spliced)) # returns tuples with the modified mutated sequence (change from exactly one mutation, and the position of the mutation in seq
          # the first part of the tuple is therefore [u,u,u,e,T,e,e,e,e,u,u,u]
          # the mutated_sequence constructor will remove the UTRs, do the reverse complement if needed and compute the 27mer 
          seqs_mutated = [mutated_sequence("".join(x[0]),x[1],this_mutation,genomic_seq_spliced,_5p_utr,rel_3p_utr,strand,peptide_len_snv,peptide_len_indel,store_mRNA,prot_export) for x in seqs_mutated_tmp] # make the objects which hold the info on the mutant sequence
          if len(seqs_mutated) == 0:
            self._loginfo("no mutated seqs " + this_mutation.get_pos() + " " + this_mutation.first_allele + " " + this_mutation.second_allele)
            continue
          assert len(seqs_mutated) <= 2 # sanity check
          if len(seqs_mutated) == 2:
            #assert not seqs_mutated[1].prot_is_different
            if seqs_mutated[1].prot_is_different:
              seqs_mutated[0].het_non_ref = True
              seqs_mutated[1].het_non_ref = True
              print seqs_mutated[1].mutation.get_pos()

          # cross link the transcript_info and mutated_sequence objects
          # also, add expression values to mutated_sequence object
          for si,s in enumerate(seqs_mutated):
            this_transcript_info.add_sequence(s)
            if len(seqs_mutated) > 1:
              s.other_peptides_from_this_mutation = si
            if s.prot_is_different:
              s.exon,s.exon_expression = get_mutation_transcript_annotation(this_mutation,i)[0:2]
            s.transcript_info = this_transcript_info
          all_seqs_mutated.append(seqs_mutated)
        
        # finalize info object for this transcript
        this_transcript_info.order()  

        # check phase if BAM file(s) is/are available
        # check phase (more than three in phase is an exception) and distance
        if len(bam_files) > 0:
          phasing_dict,phasing_report,phasing_info = check_phase(mut,bam_files,samtools_bin,phasing_dict,phasing_report,phasing_info)
          this_transcript_info.phase_check_done = True
          this_transcript_info.phase_info = phasing_info
          this_transcript_info.phase_dict = phasing_dict
          # merge sequences of mutations which are in phase
          all_seqs_mutated,merge_messages = merge_mutated_seqs_in_phase(this_transcript_info,all_seqs_mutated,mut,peptide_len_snv)
          phasing_report += merge_messages
          

        for seqs_mutated in all_seqs_mutated:
          # prepare and store sequences for MHC binding prediction
          # those sequences are selected in a way that each sliding window of the predictor 
          # will include the variant amino acid(s)
          for si in range(len(seqs_mutated)):
            s = seqs_mutated[si]
            if (s.prot_is_different) and (not "X" in s.protein):
              # store for MHC prediction
              mutated_sequences.append(s)
              s.uniqueID = str(len(mutated_sequences)) # number of stored mutated_sequence onjects so far
              for ws in window_sizes:
                x1 = find_peptides(s.protein,s.protein_u,ws,s.peptide) # find the maximum length tumor unique peptide which allows seuqnece windows with length ws which include the mutation
                # peptides taken from s.protein are not always found in 27mer (stored in s.peptide) in the case of ws=15
                # therfore, the argument s.peptide can be used ti check if the peptides are in the 27mer
                if len(x1) > 1:
                  # in a few cases the mutated peptide can be found in the unmutated protein, too
                  # as there is now "new" peptide, this window size is skipped
                  self._loginfo("mutated peptide in wt seq, " + str(ws)  + " " + s.mutation.get_pos())
                  continue
                elif len(x1) == 0: 
                  if not s.stop_gain: # a stop gain does not provide new sequences
                    self._loginfo("no unique peptides (without stop gain), " + str(ws)  + " " + s.mutation.get_pos())
                  continue
                pep1 = x1[0][0] # mutant peptide for window size ws
                if s.mutation.is_indel: # find WT specific sequences for indels wrt the window size
                  ux1 = find_peptides(s.protein_u,s.protein,ws,"")
                  try:
                    pep2 = ux1[0][0]
                  except:
                    pep2 = ""
                else:
                  start = x1[0][1] # where does the mutant peptide start in the whole 27mer ?
                  pep_len = len(x1[0][0]) # length of mutant peptide
                  pep2 = s.protein_u[start:(start+pep_len)] # WT peptide equivlaent for a SNV
                pl = [find_stop(pep1,ws),find_stop(pep2,ws)] # resolve possible stop codons, i.e. cut sequence if there is any stop codon
                if pl[0] != None: # store them in a dictionary by window size, remember mutated_sequence uniqueID; this allows later to assign the HLA prediction values
                  apl = d_mut_peptides.get(ws,[])
                  apl.append((s.uniqueID,pl[0],pl[1]))
                  d_mut_peptides[ws] = apl
    # TRANSCRIPT LOOP DONE

    if not skip_mhc:
      if len(mhc_I_selection) > 0: # run MHC class 1 prediction
        self._loginfo("mhc prediction")
        self._loginfo("hla class I alleles used: {}".format(" ".join(mhc_I_selection)))
        mutated_sequences = run_mhc_prediction3(mutated_sequences,d_mut_peptides,[8,9,10,11],organism,mhc_I_selection,mhc_class=1,procs=nprocs_IEDB,logger=self.logger)
      if len(mhc_II_selection) > 0: # run MHC class 2 prediction
        self._loginfo("mhcII prediction")
        self._loginfo("hla class II alleles used: {}".format(" ".join(mhc_II_selection)))
        mutated_sequences = run_mhc_prediction3(mutated_sequences,d_mut_peptides,[15],organism,mhc_II_selection,mhc_class=2,procs=nprocs_IEDB,logger=self.logger)
    else:
      self._loginfo("skipping mhc prediction")
 
    # add some data to mutation objects
    for chr in d.keys():
      for i in d[chr]:
        if len(i.gene_expression) != 0:
          i.avg_expression = str(sum(i.gene_expression)/len(i.gene_expression)) # average expression of transcripts affected by mutation
        else:
          i.avg_expression = ""
        if i.genomic_environment == "": # sometimes this is missing, try to add it now
          try:
            env = dna.get_seq_range_nucleotides(i.chr,i.pos - environment_size,i.pos + environment_size, "+")
            i.genomic_environment = env[:environment_size].lower() + env[environment_size].upper() +  env[environment_size + 1:].lower()
          except:
            self._loginfo("no env for " + i.get_pos())
 
    # write the set of mutations in exons/transcripts as TROMPAmut file
    output_header = get_header(human_readable=True)
    with open(mutation_set_output,"w") as f:
      f.write("\t".join(output_header))
      f.write("\n")
      for chr in d.keys():
        for i in d[chr]:
          s = format_line(i,sep="\t",subsep=" ",dec=".",human_readable=True)
          f.write(s)
          f.write("\n") 
    
    # write TROMPApep file
    with open(self.peptide_output,"w") as f:
      if len(mutated_sequences) > 0:
        f.write(mutated_sequences[0].header("\t"))
        f.write("\n")
        for t in mutated_sequences:
          f.write(t.format_line(sep="\t",dec="."))
          f.write("\n")
    
    # write phasing report if BAMs were available, i.e. the phasing was done
    if len(bam_files) > 0:
      with open(self.phase_output,"w") as f:
        for l in phasing_report:
          f.write(l)
          f.write("\n")
    

########### end run method

def main():
  cfg_dict = {}
  no_value_string = "__NO_VALUE__"
  parser = optparse.OptionParser(description="",usage="%prog [options] configFile data_1 format_1 ... data_N format_N \n")
  l=["mutation_set_output","environment_size","repeat_bed","dbsnp_file","ref_bed","genome_2bit","xref_table",
     "gene_expression","exon_expression","expression_sep","organism","mhc_I_selection","mhc_II_selection","skip_mhc","data_set_names","hla_file","store_mRNA","prot_export","samtools_bin","bam_files"]
  for i in l:
    parser.add_option("--" + i,action="store",type="string",dest=i,help="",default=no_value_string)
  (options, args) = parser.parse_args()
  if (len(args) == 0):
    print "no args"
    sys.exit()
  with open(args[0],"r") as cfg:
    for line in cfg:
      line = line[:-1]
      if not (line.startswith("#") or line == ""):
        words = line.split("=")
        words = [x.strip() for x in words]
        if len(words) < 2:
          break;
        cfg_dict[words[0]] = words[1]
  for i in l:
    o = options.ensure_value(i, no_value_string).strip()
    if o == no_value_string:
      pass
    else:
      cfg_dict[i] = o
  indata = args[1:]
  o = ntp_obj(cfg_dict)
  o.run(indata)

if __name__ == '__main__':
  main()
