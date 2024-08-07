import sys
import re
from collections import defaultdict
import os

from glob import glob
def extract_tot_reads(stats):
    try:
        tot_targ_reads = int(re.findall("([0-9]*).+[0/9].mapped",stats)[0])
        tot_dup_reads =int(re.findall("([0-9]*) \+ 0 *duplicates",stats)[0]) #correct for dups
        return (tot_targ_reads-tot_dup_reads)/100e6
    except:
        return 0

def dump(d, fout):#,control):
    keys = sorted(list(d.keys()), key = lambda k:int(k.split()[-1]))
    for k in keys:
        fout.write("%s\t%s\n"%(k,'\t'.join(d[k])))

dir  = sys.argv[1]
out_dir = sys.argv[2]
chr_list  = ['chr%s'%str(i) for i in list(range(1,23))+['X','Y']]+["MT"]
open_f = [open(f,"rb") for f in glob(dir+'/*cov')[:50]]
lopen_f = len(open_f)
stats = list(map(extract_tot_reads,[open(f).read() for f in glob(dir+'/*report.txt')[:50]]))
for CHR in chr_list:
    fout = open(out_dir+'/total_%s.res'%CHR,'w')
    d = defaultdict(list)
    for n,file in enumerate(open_f):
      if stats[n]>0:
        file.seek(0)
        START = False
        while 1:
          l = file.readline().decode()
          if not l:
              break
          try:
              chr, pos, cov = l.split()
          except:
              print(l)
          if CHR == chr:
              cov_n = float(cov)/stats[n] ## normalized
              #print (cov_n,cov,stats[n])
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

