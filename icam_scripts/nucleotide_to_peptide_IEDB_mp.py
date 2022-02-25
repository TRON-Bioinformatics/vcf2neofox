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
import multiprocessing

from nucleotide_to_peptide_constants import *
from MHC_prediction_interface_wo_classes import *

def _inner_loop(args):
  mhc_funct,method,mm,ws,lm,idm,lu = args
  res,mhc_excpetion = mhc_funct(method, mm, ws, lm)
  r = []
  best_peptides = []
  unmut_peptides = []
  assert len(res) == len(idm)
  for idx,i in enumerate(res):
    if (len(i[2]) == 0):
      r.append(999.0)
      best_peptides.append("")
      unmut_peptides.append("")
    else:
      r.append(min(i[2]))
      seq_index = i[2].index(r[-1])  # which element is minimal ?
      best_peptides.append(lm[idx][seq_index:(seq_index+ws)])
      if lu[idx] == None:
        unmut_peptides.append("")
      else:
        unmut_peptides.append(lu[idx][seq_index:(seq_index+ws)])
  return (ws,mm,r,best_peptides,unmut_peptides,idm,mhc_excpetion)

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

def run_mhc_prediction3(pep_obj_list,d_mut_peptides,window_sizes,organism,selection,mhc_class=1,procs=20, logger=None):
  if mhc_class == 1:
    method = "IEDB_recommended" #"consensus" #method = "netmhcpan"
    mhc = get_mhc_list(method, organism, window_sizes)
    mhc_funct = run_mhc
  else:  # mhc class 2
    method = "consensus3"
    mhc = mhcII_human
    mhc_funct = run_mhcII
  mhc = select_mhc_alleles(mhc,selection)

  s1 = []
  id1 = []
  jobs = []
  for ws in window_sizes:
    if ws not in d_mut_peptides:
      continue
    idm,lm,lu = zip(*d_mut_peptides[ws])
    assert len(lu) == len(lm)
    m = mhc[method][ws]
    for mm in m:
      jobs.append((mhc_funct,method,mm,ws,lm,idm,lu))
  pool = multiprocessing.Pool(procs)
  r = pool.map(_inner_loop,jobs)
  for i in r:
    s1.append(i[:5])
    id1.append(i[5])
    if i[6] != "":
      logger.info(i[6])
  
  for pep_obj in pep_obj_list:
    id = pep_obj.uniqueID
    x = []
    m1,l1 = get_min_mhc(s1,id1,id) 
    # chr1_1(10, 'H-2-Dd', 69.399999999999991, 'KKKKKKKKKK', 'KKKKKKKKKA')
    m2 = [("","",("",))]
    if len(l1) != 0:
      try:
        if len(m1[4]) >= m1[0]: # in some special cases (e.g. stop loss), the unmutated peptide might be shorter than the window size; catch this here
          m2,mhc_excpetion  = mhc_funct(method,m1[1],m1[0],[m1[4]])
          if mhc_excpetion != "":
            logger.info(mhc_excpetion)
      except MHC_I_util.UnexpectedInputError:
        pass
      except mhc_II_binding.UnexpectedInputError:
        pass
      except IndexError:
        print >> sys.stderr, "INDEX_ERROR in nucleotide_to_peptide_IEDB_mp.py, 94",m1,method,pos,i.get_alleles(),l1
        raise
        #sys.exit()
    m1 = list(m1)
    try:
      m1.append(list(m2[0][2]))
    except:
      print >> sys.stderr, "INDEX_ERROR in nucleotide_to_peptide_IEDB_mp.py, 104",m1,m2,id
      raise
      #sys.exit()
    if mhc_class==1:
      pep_obj.mhcI = m1
    else:
      pep_obj.mhcII = m1
  return pep_obj_list

 


