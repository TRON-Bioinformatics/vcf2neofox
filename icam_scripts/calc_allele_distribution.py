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
import subprocess
import os


import argparse


def main():
  parser = argparse.ArgumentParser(description='Make read count output for iCaM result tables.')
  parser.add_argument('pos_file', type=str, nargs=1, help='position file with positions as chr_pos in the first column (tab separated)')
  parser.add_argument('bam', type=str, nargs='+',help = 'N bam files for processing')
  parser.add_argument('--offset', type=int,default = 0,help = 'default=0')
  parser.add_argument('--add_to_file', action='store_true', default = False,help = 'default=False')
  parser.add_argument('--java', type=str,default = "/kitty/code/jdk1.7.0/bin/java",help = 'default=/kitty/code/jdk1.7.0/bin/java')
  parser.add_argument('--jar', type=str,default = "/kitty/code/iCaM/jar/ReadCounts.jar",help = 'default=/kitty/code/iCaM/jar/ReadCounts.jar')
  parser.add_argument('--tmpfile', type=str,default = "tmpfile123",help = 'default=tmpfile123')
  parser.add_argument('--infile_has_header', action='store_true', default = True,help = 'default=True')
  parser.add_argument('--prefix', nargs="*",type=str,default = "",help = 'default=""')
  parser.add_argument('--other_format',action='store_true', default = False,help = 'default=False')
  parser.add_argument('--samples', type=str,default="",help = 'default=""')
  parser.add_argument('--mem', type=str,default="4",help = 'default="4"')
  args = parser.parse_args()
  run(args.pos_file[0],args.offset,args.add_to_file,args.java,args.jar,args.bam,args.tmpfile,args.prefix,args.infile_has_header,args.other_format,sampletypes=args.samples,mem=args.mem)

def run(pos_file,offset,add_to_file,java,jar,bam,tmpfile,prefix,infile_has_header,other_format,outfile=None,sampletypes=None,variantcols=[1,2],mem="4"):
  pos = []
  inlines = []
  header = []
  outstream = sys.stdout
  if other_format:
    infile_has_header = False
    offset = 1
  if outfile != None:
    outstream = open(outfile,"w")
  with open(pos_file) as f:
    for lc,line in enumerate(f):
      words = line.strip("\n").split("\t")
      if lc == 0:
        if infile_has_header:
          header = words
          continue
        else:
          header = [""] * (len(words))
      inlines.append(line.strip("\n"))
      if other_format:
        if line.startswith("#"):
          continue
        p0 = words[0]
        p1 = words[1]
      else:
        p = words[0].split("_")
        p0 = "_".join(p[0:-1])
        p1 = p[-1]
      pos.append(p0 + "\t" + str(int(p1) - offset))
  if len(pos) == 0:
    raise Exception("No input positions")
  matching_lines = run2(pos,java,jar,tmpfile,bam,mem)
  
  VAF = False
  
  chars_of_interest = "ACGT"
  if add_to_file:
    outstream.write("\t".join(header) + "\t")
  prefixl = [""] * len(bam)
  for i,j in enumerate(prefix):
    prefixl[i] = j
  for i,j in enumerate(prefixl):
    if i > 0:
      outstream.write("\t")
    outstream.write("\t".join([j + "_" + x for x in (["coverage"] + list(chars_of_interest))]))
  if (sampletypes!=None) and (len(sampletypes) >= len(bam)):
    outstream.write("\tVAF_in_tumor\tVAF_in_normal\tVAF_in_RNA")
    VAF = True
  outstream.write("\n")
  for idx,i in enumerate(matching_lines):
    if add_to_file:
      outstream.write(inlines[idx] + "\t")
    if i != "":
      outstream.write("\t".join(i)) 
    else:
      outstream.write("\t".join((["0"]*(1 + len(chars_of_interest)))*len(bam)))
    if VAF:
      outstream.write("\t")
      v = get_VAF_from_line(inlines[idx],i,sampletypes,variantcols)
      outstream.write("\t".join([str(x) for x in v]))
    outstream.write("\n")
  
  if outfile != None:
    outstream.close()

def get_VAF_from_line(inline,i,sampletypes,variantcols):
  try:
    assert len(sampletypes)*5 == len(i)
  except:
    print >> sys.stderr, inline
    print >> sys.stderr, i
    print >> sys.stderr, sampletypes
    sys.exit()

  l = inline.split("\t")
  ref = l[variantcols[0]]
  alt = l[variantcols[1]]
  res = {"r":0,"t":0,"n":0}
  total = {"r":0,"t":0,"n":0}
  # coverage A C G T
  pos = {"A":1, "C":2, "G":3, "T": 4}
  for j,s in enumerate(sampletypes):
    k = i[(j*5):((j*5) + 6)]
    total[s] += int(k[0])
    if alt in pos:
      res[s] += int(k[pos[alt]])
  ll = []
  for j in "tnr":
    ll.append(safediv(res[j],total[j]))
  return ll

def safediv(x,y):
  if y == 0:
    return -1.0
  else:
    return float(x)/float(y) 

def run2(pos,java,jar,tmpfile,bam,mem="4"):
  with open(tmpfile,"w") as f:
    for i in pos:
      f.write(i + "\n")
    f.write("\n")
  """/kitty/code/jdk1.7.0/bin/java -Xmx4G -jar test.jar -interval pos.bed $x1 > test1.out"""
  cmd = java + " -Xmx"+mem+"G -jar " + jar + " -mq 0 -qual 0 -total -interval " + tmpfile + " " +  " ".join(bam)
  sp = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
  (data,err) = sp.communicate()
  if sp.returncode != 0:
    try:
      os.remove(tmpfile)
    except:
      pass
    raise Exception("################ JAVA ERROR ################\n" + err)
  os.remove(tmpfile)
  data = data.split("\n")

  if data[-1] == "":
    data = data[:-1]

  matching_lines = [["0","0","0","0","0"] * len(bam)] * len(pos)
  used = 0
  for line in data:
    words = line.strip("\n").split("\t")
    p = "\t".join([words[0],str(int(words[1])-1)]) #  returned data is one based; input data is zero based
    if p in pos:
      matching_lines[pos.index(p)] = words[2:]
      used +=1

  return matching_lines

if __name__ == '__main__':
  main()
