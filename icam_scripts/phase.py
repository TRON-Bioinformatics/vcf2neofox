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
import mutation
import subprocess


def get_read_len(sambin,bam):
  cmd1 = [sambin,"view",bam]
  cmd2 = ["head","-n1"]
  p1 = subprocess.Popen(cmd1,stdout=subprocess.PIPE)
  p2 = subprocess.Popen(cmd2,stdin=p1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p1.stdout.close()
  out,err = p2.communicate()
  p1.kill()
  return len(out.split("\t")[9])

def get_alignments(pos,bam,sambin):
  if "_" in pos:
    x = pos.split("_")
    chr = "_".join(x[:-1])
    p = x[-1]
    p = int(p) + 1
  elif ":" in pos:
    chr,p = pos.split(":")
    p = int(p)
  else:
    raise Exception("wrong format: " + pos)
  region = chr + ":" + str(p-1) + "-" + str(p+1)
  cmd = sambin + " view " + bam + " " + region
  sys.stderr.flush()
  p1 = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out,err = p1.communicate()
  return out.split("\n")[:-1],p

def parse_cigar(cigar,start1based):
  start = start1based
  # process cigar
  m = []
  r = []
  cs = 0
  read_pos = 0
  for i,c in enumerate(cigar):
    if c == "M":
      l = int(cigar[cs:i])
      m.append((start,l))
      start += l
      r.append((read_pos,l))
      read_pos += l
      cs = i + 1
    elif c in "IS": #move read, not ref
      l = int(cigar[cs:i])
      read_pos += l
      cs = i + 1
    elif c == "H": #do nothing 
      cs = i + 1
    elif c in "DNP=X": #move ref, not read
      l = int(cigar[cs:i])
      start += l
      cs = i + 1
  return m,r

def match_nt(read_seq,read_start,readlen,p,alt,cigar):
  matched_reference_intervals,matched_read_intervals = parse_cigar(cigar,read_start)
  c = -1
  for interval_start,interval_len in matched_reference_intervals:
    c += 1
    if (interval_start > (p - interval_len)) and (interval_start <= p):
      idx = p - interval_start
      read_bases = range(matched_read_intervals[c][0],matched_read_intervals[c][0] + matched_read_intervals[c][1])
      if read_seq[read_bases[idx]].upper() == alt.upper():
        return True
  return False
    
def run(pos1,pos2,alt1,alt2,bam,sambin): # pos1 and pos2 might be one or zero based
  readlen = get_read_len(sambin,bam)
  aln1,p1 = get_alignments(pos1,bam,sambin) # p1 and p2 are one based
  aln2,p2 = get_alignments(pos2,bam,sambin)
  d = {}
  for a in aln1:
    w = a.split("\t")
    st = int(w[3])
    if match_nt(w[9],st,readlen,p1,alt1,w[5]):
      d[w[0]] = w      
  res = []
  nores1 = []
  nores2 = []
  for a in aln2:
    w = a.split("\t")
    st = int(w[3])
    if match_nt(w[9],st,readlen,p2,alt2,w[5]):
      if w[0] in d:
        res.append((w[0],w[9]))
        if w[2] != d[w[0]][2]:
          pass
      else:
        nores1.append((w[0],w[9]))
    else:
      nores2.append((w[0],w[9]))
  return aln1,aln2,res,nores1,nores2

def readable_report(pos1,alt1,pos2,alt2,aln1,aln2,res,nores1,nores2):
  print
  print "First position (mutant allele):", pos1, "(",alt1,")"
  print "Second position (mutant allele):", pos2, "(",alt2,")"
  print len(aln1),"and", len(aln2),"reads were found"
  print "#" * len(str(len(res))) + "##########################################" 
  print "#", len(res),                "(c)DNA fragments carry both mutations #"
  print "#" * len(str(len(res))) + "##########################################"
  print len(nores1), "reads cover position 2 but include only the second mutation"
  print len(nores2), "reads cover position 2 but do not include the second mutation"
  if len(res) > 0:
    x = len(res) / 2
    print "Example read with both mutations\t",res[x][0],"\t",res[x][1]
  if len(nores1) > 0:
    x = len(nores1) / 2
    print "Example read with only second mutation\t",nores1[x][0],"\t",nores1[x][1]
  if len(nores2) > 0:
    x = len(nores2) / 2
    print "Example read without second mutation\t",nores2[x][0],"\t",nores2[x][1]


def run_m(mut1,mut2,bam,sambin):
  alt1 = mut1.first_allele
  alt2 = mut2.first_allele
  pos1 = mut1.get_pos()
  pos2 = mut2.get_pos()
  aln1,aln2,res,nores1,nores2 = run(pos1,pos2,alt1,alt2,bam,sambin)
  return (len(aln1), len(aln2), len(res))


def main():
  sambin = "samtools"
  if len(sys.argv) < 6:
    print >> sys.stderr, "python phase.py mut_1 mut_2 alt_allele_1 alt_allele_2 bam_file [path_to_samtools_binary]"
    print >> sys.stderr, ""
    print >> sys.stderr, "  Uses the reads stored in a BAM file to check the phase of the given mutations"
    print >> sys.stderr, "  mut_1                      first mutation, given as 'chr_zerobasedposition' or 'chr:onebasedposition'"
    print >> sys.stderr, "                             e.g. chr1:12345 or chr1_12344"
    print >> sys.stderr, "  mut_2                      as mut_1, but for the second mutation"
    print >> sys.stderr, "  alt_allele_1               alternative allele of first mutation"
    print >> sys.stderr, "  alt_allele_2               alternative allele of second mutation"
    print >> sys.stderr, "  bam_file                   a BAM file"
    print >> sys.stderr, "  path_to_samtools_binary    optional path to samtools executable (default: samtools)" 
    sys.exit()
  if len(sys.argv) == 7:
    sambin = sys.argv[6]
  aln1,aln2,res,nores1,nores2 = run(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sambin)
  readable_report(sys.argv[1],sys.argv[3],sys.argv[2],sys.argv[4],aln1,aln2,res,nores1,nores2)

if __name__ == '__main__':
  main()


