import sys
import re
from collections import defaultdict
import os

from glob import glob
def extract_tot_reads(stats):
    tot_targ_reads = int(re.findall("([0-9]*).+[0/9].mapped",stats)[0])
    return tot_targ_reads/100e6

def dump(d, fout):#,control):
    keys = sorted(d.keys(), key = lambda k:int(k.split()[-1]))
    for k in keys:
        fout.write("%s\t%s\n"%(k,'\t'.join(d[k])))

dir  = sys.argv[1]
out_dir = sys.argv[2]
chr_list  = ['chr%s'%str(i) for i in ['X','Y']]
open_f = [open(f) for f in glob(dir+'/*cov')]
lopen_f = len(open_f)
stats = map(extract_tot_reads,[open(f).read() for f in glob(dir+'/*report.txt')])
for CHR in chr_list:
    fout = open(out_dir+'/total_%s.res'%CHR,'w')
    d = defaultdict(list)
    for n,file in enumerate(open_f):
        
        START = False
        while 1:
          l = file.readline()
          if not l:
              break
          try:
              chr, pos, cov = l.split()
          except:
              print(l)
          if CHR == chr:
              cov_n = float(cov)/stats[n] ## normalization by number total reads
              if cov_n>10000:
                  pass
              d[chr+'\t'+pos].append(str(cov_n))
              START = True
          elif START:
              file.seek(-len(l),1)
              START = False
              break
    dump(d, fout)
    fout.close()

