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
from time import gmtime, strftime
import traceback


# mhc stuff
mhc_tool_version = "20130222"
mhcII_tool_version = ""
window_sizes = (8,9,10,11,15)

mhc_path="/code/IEDB/mhc_i/src/"
mhcII_path="/code/IEDB/mhc_ii"
sys.path.append(mhc_path)
MHC_I_util =  __import__("util",globals(), locals(), [], -1)
MHC_I_seqpredictor =  __import__("seqpredictor",globals(), locals(), [], -1)
sys.path.append(mhcII_path)
mhc_II_binding =  __import__("mhc_II_binding",globals(), locals(), [], -1)
  
def run_mhcII(method, mhc, peptide_length, peptides, return_core=False):
  mhc_exception = ""
  proteins = mhc_II_binding.Proteins()
  for i in peptides:
    if i != None:
      proteins.add_protein(i)
  seq = [('sequence_format', 'auto'), ('sort_output', 'position_in_sequence'), ('cutoff_type', 'none'), ('output_format', 'ascii'), ('allele', mhc), ('sequence_file', ""), ('pred_method', method)]
  form = dict(seq)
  #print proteins.sequences
  mhc_predictor = mhc_II_binding.MHCBindingPredictions(form)
  s = []
  try:
    length,allele,scores = mhc_predictor.predict(proteins)[0]
  except Exception as mhc_call_exception:
    mhc_exception = "prediction not possible due to error in external program\n" + str(mhc_call_exception)
    for i in peptides:
      s.append((peptide_length,mhc,[]))
    return s,mhc_exception
  assert len(scores) == len(peptides)
  for i in scores:
    x = []
    cores = []
    for j in i:
      if method == "consensus3":
        x.append(j[0])
      elif method == "IEDB_recommended":
        x.append(j[0])
        if return_core:
          cores.append(j[4]) # 4 = smm, 7 = nn
      else:
        x.append(j[1])  # j[0] would be the anchor peptide
    s.append((length,allele,x,cores))
  return s,mhc_exception

def run_mhc(method, mhc, peptide_length, peptides, return_core=False):  # return_core=False is for interface compatibility with _run_mhcII
  mhc_exception = ""
  proteins = MHC_I_util.Proteins()
  for i in peptides:
    if i != None:
      proteins.add_protein(i)
  #print MHC_I_util.get_species(mhc).split()
  input = MHC_I_util.InputData(mhc_tool_version,method,[mhc],"",[peptide_length],proteins,MHC_I_util.get_species(mhc).split())
  mhc_predictor = MHC_I_seqpredictor.MHCBindingPredictions(input)
  r = [None] * len(peptides)
  #print mhc_predictor.get_score_unit()
  try:
    mhc_scores = mhc_predictor.predict(input.input_protein.sequences)
    #print mhc_scores
  except Exception as mhc_call_exception:
    mhc_exception = "prediction not possible due to error in external program\n" + str(mhc_call_exception)
    r=[]
    for i in peptides:
      r.append((peptide_length,mhc,[]))
    return r, mhc_exception
  # (results_peptide_length, results_allele, scores) = mhc_scores[0] # v 2.4
  (results_peptide_length, results_allele, scores, method_used) = mhc_scores[0]
  src_cnt = 0
  for i,x in enumerate(peptides):
    assert(x != None)
    if x != None:
      #print >> sys.stderr , i,x,
      if method in ("consensus","IEDB_recommended"):
        s = scores[src_cnt][0]
      else:
        s = scores[src_cnt]
      r[i] = (results_peptide_length, results_allele, s)
      src_cnt += 1
  return r, mhc_exception

def select_mhc_alleles(mhc, selection):
  methods = mhc.keys()
  d = {}
  for i in methods:
    sizes = mhc[i].keys()
    for j in sizes:
      l = mhc[i][j]
      l2 = []
      for k in l:
        for ks in selection:
          if k == ks:
            l2.append(k)
      x = d.get(i,{})
      x[j] = l2
      d[i] = x
  return d

def get_mhc_list(method, s, lengths):
  ms = MHC_I_util.MethodSet()
  method_index = ms.get_method_index(method)
  mhc_list = MHC_I_util.get_mhc_list(method_index)
  d = {}
  d2 = {}
  for (species,mhc, peptide_length) in mhc_list:
    if (s == species) and (peptide_length in lengths):
      l = d2.get(peptide_length,[])
      l.append(mhc)
      d2[peptide_length] = l
  d[method] = d2
  return d

if __name__ == '__main__':
  r = run_mhc("IEDB_recommended","HLA-A*02:01",9,["SYFPEITHI","GSGSGSGSGSGSGS"])
  print r

  r2 = run_mhcII("IEDB_recommended", "HLA-DRB1*04:01", 15, ["SYFPEITHIGSGSGSGSGSGSGS"], True)
  print r2

  sys.exit()

  print strftime("%Y-%m-%d %H:%M:%S", gmtime())
  mhc = _get_mhc_list("netmhcpan","human",[9,10])
  selection = ["HLA-A*02:01","HLA-A*01"]
  mhc = _select_mhc_alleles(mhc,selection)
  print mhc
  print strftime("%Y-%m-%d %H:%M:%S", gmtime())


