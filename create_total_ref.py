from numpy import mean,std
import sys

with open(sys.argv[1]) as tot_n:
    for l in tot_n:
        ll = l.strip().split()
        chr, pos = ll[:2]
        mn = mean(map(float,ll[3:]))
        st = std(map(float,ll[3:]))
        print "%s\t%s\t%f\t%f"%(chr,pos,mn,st)
