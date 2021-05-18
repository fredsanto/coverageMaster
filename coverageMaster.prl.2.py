''' 
coverageMaster
 - Federico Santoni

* Updates
  16.09.2016
     Multiquery
     Gene list as input and Gene info retrieval
  1.02.2017
     Automagic detection HMM based
  01.07.2020
     ACGH mode
     stddev control
'''
from optparse import OptionParser
from glob import glob
import sys,os,re
sys.path.append('/home/fsantoni/matplotlib-1.1.0/build/lib.linux-x86_64-2.6')
sys.path.append('.')
import matplotlib
matplotlib.use('PDF')
from pylab import *
from collections import defaultdict
from HMM_CM import *
from ReadCount import *
import time
from sympy import Union, Interval, SympifyError
from multiprocessing import Pool
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from libCoverageMaster import *
from more_itertools import locate




reference = open(wd+"/REFSEQ_hg19.chrM.tab.HUGO").read().strip().split('\n')
#reference = open(wd+"/REFSEQ_hg38_HGMD3").read().strip().split('\n')
exon_reference = defaultdict(list)

for r in open(wd+"/hg19.exons.merged.bed").read().strip().split('\n'):
#for r in open(wd+"/hg38.exons.merged.bed").read().strip().split('\n'):
    rchr,rstart,rend = r.split()
    exon_reference[rchr].append((rstart,rend))


qregion = {}
qregions = []
wid = 1







usage = "usage: %prog [options] <cov_file> <stats_file> <gene list(file or gnames comma separated)|region(chr:start-end)> -r <reference.cov> -o <output_px>"
parser = OptionParser(usage = usage)
parser.add_option("-c",'--control',  dest="control", help = "<optional> txt file with a control file per line (.cov with .report.txt in same folder)")
parser.add_option("-s",'--single-control',  dest="single_control", help = "<optional> single control file (.cov with .report.txt in same folder)")
parser.add_option("-r",'--ref',  dest="ref", help = "reference file", default = "")
parser.add_option("-o",'--out',  dest="output_px", help = "output plot file", default = "")
parser.add_option("-f",'--force',  action="store_true", dest="force", default=False, help = "force output")
parser.add_option("-e",'--exons',  dest="exons", help = " <optional> n. of exons", default = 5)
parser.add_option("-x",'--offset',  dest="offset", help = " <optional> offset to ref", default = 0)
parser.add_option("-d",'--width',  dest="wid", help = " <optional> width of std", default = 1) #choose width of sd
parser.add_option("-l",'--level',  dest="lev", help = " <optional> wavelet level", default = 5) #choose level of compression
parser.add_option("-w",'--wig',  action="store_true", help = " <optional> write wig", default = False)
parser.add_option("-b",'--bed',  action="store_true", help = " <optional> BED input", default = False)
(options, args) = parser.parse_args()
try:
    output_px = options.output_px.split("/")[-1]
    LOGFILE = open(output_px+".CM.log","a")
    offset = float(options.offset)
    cov_region = args[0]
    lev = int(options.lev)
    wid = float(options.wid)
    if '.cov' not in cov_region:
        print("\nMichel, please.. not again... use the .cov!\n")
        raise Exception
    
    stats_file = args[1]
    stats = extract_tot_reads(open(stats_file).read())
    FORCE = options.force
    XIST = gene_reference(['XIST'], reference, LOGFILE, int(options.exons))
    
    if args[2][:-1] == "/":
        raise Exception("GeneList not valid")
    if options.bed: #it's a BED file in input
        for l in open(args[2]):
            qregion = {}
            qregion['chr'],qregion['start'],qregion['end'] = l.strip().split("\t")[0:3]
            qregions.append(qregion)
        qregions += XIST
    elif ":" in args[2]:
        qregions = []
        s1 = args[2].split(':') 
        qregion['chr'] = s1[0]
        qregion['start'], qregion['end'] = s1[1].split('-')[0], s1[1].split('-')[1]
        qregions.append(qregion)
        qregions += XIST
    else:
        try:
            gfilename = args[2]
            glistfile = open(gfilename)
            glist = [x.split()[0] for x in glistfile.read().strip().split(' ')]
            if len(glist)<2:
                glist = [x.split("\n") for x in open(args[2]).read().strip().split(' ')][0]
        except:
                glist = args[2].split(",")
        if ':' not in glist[0]:
            if os.path.exists("%s_coderef"%gfilename):
                qregions = pickle.load(open("%s_coderef"%gfilename,"rb"))
            else:
                qregions = gene_reference(glist, reference, LOGFILE, int(options.exons))+XIST
            ###here > dump if it is non existing
            try:
                if not os.path.exists("%s_coderef"%gfilename) and len(qregions) > 50:
                    pickle.dump(qregions,open("%s_coderef"%gfilename,"wb"))
            except:
                pass #not a glist file
        else:
            for r in glist:
                qregion = {}
                s1 = r.split(':') 
                qregion['chr'] = s1[0]
                qregion['start'], qregion['end'] = s1[1].split('-')[0], s1[1].split('-')[1]
                qregions.append(qregion) 
            #if ":" in r:
            #    qregions += [XIST]
            #else:
            qregions += XIST
    if not options.ref:
        print("\nInclude the reference")
#        raise Exception

#    for f in [stats_file,options.ref,cov_region,options.control,args[2]]:
#        print("%s - %d"%(f, os.path.isfile(f)))
 
    
except:
    print("Unexpected error:", sys.exc_info()[0])
    parser.print_help()
    sys.exit(0)


n_exons = 0
total_bp = 0

#report = open(output_px+".CMreport","w")       
#report2 = open(output_px+".CMexonstruct","w")       
cr = Regions(cov_region)
cr.create_index()
logreport( "Query index created", logfile =LOGFILE)

if options.ref:
    cref = Regions(options.ref)
    cref.create_index()
    logreport( "Reference index created", logfile = LOGFILE)
else:
    cref = None 



##main
def processCoverage(gene):#, cr, cref, stats, ccont, cref, cstats, XIST, pp)

 time.sleep(.5) 
 if gene['chr']!='chrM':   
  
   if 1:#try:
    #reopen for threads purposes
    cref.f = open(cref.f.name)
    cr.f = open(cr.f.name)
    ccont.f = open(ccont.f.name)

    if cref:
        signal, ref_exon_avg, ref_exon_min, enlight, data_n = signalProcessor(gene, cr, cref, stats, red = False)   
        signal = signal + offset
    
    csignal, unused, unused, unused, unused = signalProcessor(gene, ccont, cref, cstats)   
    csignal = csignal + offset
    if gene['chr'] == 'chrX':
        signal = signal/normX
        csignal = csignal/normCX+1e-4
    #data that have zero coverage in regions:
    if sum(signal)==0 and sum(csignal)==0: #no coverage
        return "NO COVERAGE"
    elif sum(csignal)==0:
        logreport("NO Coverage in Control for %s"%gene, logfile=LOGFILE)
        return "CONTROL NO COVERAGE"
        
    sd = array(ref_exon_min)/array(ref_exon_avg)
    genethere = array(enlight)!=0
   
    if sum(abs(csignal))>0 and len(signal)==len(csignal):
        #wid: parameter -d
        stdM = 1+wid*sd
        stdm = 1-wid*sd
        unormsignal = signal-1
        noninfidx = where((abs(unormsignal) <= wid*sd) + (unormsignal<=wid*sd)*(csignal<signal) + (unormsignal>=wid*sd) * (csignal>signal))
        try:
            ratio = signal/(.001 + csignal)
        except:
            logreport("%s aprox failed. Ratio length inconsistency"%gene, logfile=LOGFILE)
            ratio = zeros(len(signal))
            
        ratio[noninfidx] = 1
        approx,alev = HMM(ratio,sd, lev = lev, LOGFILE=LOGFILE, mask=genethere)
        #control,clev = HMM(sd/median(sd),median(sd)*ones(len(sd)), booster = 1000,lev = lev,LOGFILE=LOGFILE)
        maskp = (approx>1)
        maskm = (approx<1)
        if FORCE or ((sum(maskp)+sum(maskm))>0) and (sum(signal*maskp*genethere - maskp*stdM*genethere) > 0) or (-sum(signal*maskm*genethere - maskm*stdm*genethere) > 0):
            
            '''
            cc = control.tolist()
            cc=[1]+cc+[1]
            steps = where(diff(cc)!=0)[0]
            for n,i in enumerate(steps):
                if n%2:
                    art_idx = steps[n-1]-500,steps[n]+500
                    art_idx = (max(0,art_idx[0]),min(len(signal)-1,art_idx[1]))
                    art_idx = list(range(art_idx[0],art_idx[1]+1))
                    control[art_idx] = 0
        
            '''
            #approxe = approx*genethere
            #control = control*genethere
            
            #if FORCE or [x for x in approxe*control if x!=1 and x!=0]:
            #print(FORCE)    
    
            if options.wig:
                wig  = wig_writer(gene["chr"],data_n)
                fwigout = open(output_px+'_'+gene['gene']+'.wig','w')
                fwigout.write('\n'.join(wig))
                fwigout.close()
            #logreport("convertHMM %s"%gene,logfile=LOGFILE)
            call = convertHMM(approx,data_n)
            return([gene,signal,csignal,enlight,stdm,stdM,approx,ref_exon_avg/sum(ref_exon_avg),ratio,call])
 return None

def signalProcessor(gene, cr, cref, stats,LOGFILE=LOGFILE, red=False):
 
  if 0 < (int(gene['end'])-int(gene['start'])) and gene['chr']!='chrM':
    if 'query' in gene:
        qregion = gene['query'] 
        logreport( "Processing Gene: %s - Region %s:%s-%s"%(gene['gene'],gene['chr'],gene['start'],gene['end']), logfile = LOGFILE)
    else:
        qregion  = gene
        qregion['gene']='%s:%s-%s'%(gene['chr'],gene['start'],gene['end'])
        logreport(  "Processing Region %s:%s-%s"%(qregion['chr'],qregion['start'],qregion['end']),logfile = LOGFILE)
    data = []
    if qregion:
        enlight = []
        plot_region = []
        mean_region = []
        mstd_region = []
        pstd_region = []
        pos_region = []
        count = 0
        gpos = 0
        size = (int(qregion['end']) - int(qregion['start'])) # ArryCGH mode
        step = 1 if size < 500e3 else 5
        if cref:
            ref_exon_min = []
            ref_exon_Max = []
            ref_exon_avg = []
            ref = []
            FREEZE = False
            gen_cr = cr.focus_single_prl((qregion['chr'],int(qregion['start']),int(qregion['end'])), step = 1)
            for r in cref.focus_single_prl((qregion['chr'],int(qregion['start']),int(qregion['end'])), step = step):
                rr = r.strip().split()
                chr,pos = rr[:2]
                if len(rr[2:])>2:
                    cov_list = list(map(float,rr[2:]))
                    stdv = std(cov_list)
                else:
                    cov_list = float(rr[2])
                    stdv = float(rr[3])
                mcov = mean(cov_list)
                count+=1
                #print("%d\r"%count, end=' ')
                LOOP = False
                while not LOOP and not FREEZE:
                    try:
                        chr, gpos, cov = next(gen_cr).strip().split()
                        if int(gpos) == (int(pos)+1):
                            FREEZE = True
                        
                    except:
                        chr, pos, cov = qregion['chr'], pos, 0
                        break
                    LOOP = pos == gpos
                if int(gpos) < int(pos):
                    FREEZE = False
                if pos == gpos:
                    FREEZE = False
                    ref_exon_min.append(stdv)
                    ref_exon_avg.append(mcov)
                    ref.append(pos)

                    if red:
                        cov = float(cov)/2
                        
                    enlight = -1.5 + calculate_exon(gene, pos) if  (int(pos)>=int(gene['start'])) and (int(pos)<=int(gene['end'])) else 0                                       
                    #data.append((chr,pos,cov))
                    if isnan(float(cov)):
                        raise "NAN value found! Chr %s Position %s "
                    pos_region.append((pos,float(cov)/stats,enlight))
        
        #pos_region = realign(pos_region, ref) 
        plot_region = [p[1] for p in pos_region]
        enlight = [p[2] for p in pos_region]
        
    if cref:
        plot_region_n = array(plot_region)/array(ref_exon_avg)
        plot_region_n = array([0 if isnan(i) else i for i in plot_region_n])
        data_n = [(i[0],i[1],plot_region_n[n]) for n,i in enumerate(pos_region)]
        return plot_region_n, ref_exon_avg, ref_exon_min, enlight, data_n
    else:
        return plot_region

def print_output(glist,fname,round):
    pass

if __name__ == '__main__':
    XIST = qregions.pop()
    plot_chrX_region, unused, unused, unused, unused = signalProcessor(XIST, cr, cref, stats, LOGFILE)   
    normX = median(plot_chrX_region)
    qregions_failed = []
    if 1:
        if options.control:
            cc = open(options.control).read().strip().split("\n")
        if options.single_control:
            cc = [options.single_control]
        for c in cc:
          if c.split("/")[-1] != cov_region.split("/")[-1]: 
            cstats = c.replace(".cov","")+".report.txt"
            ccont = Regions(c)
            ccont.create_index()
            cstats = extract_tot_reads(open(cstats).read())
            logreport( "Control index created", logfile = LOGFILE)

            control_chrX_region, unused, unused, unused, unused = signalProcessor(XIST, ccont, cref, cstats, LOGFILE)
            normCX = median(control_chrX_region)
 
            logreport("Loop with control %s"%c, logfile = LOGFILE)
            #processCoverage(qregions[0]) #debug mode
            if sys.version_info[0] < 3:
                from contextlib import closing
                with closing(Pool(processes=5)) as p:
                    if len(qregions):
                        repository = p.map(processCoverage, qregions)
                        p.terminate()
                    else:
                        logreport("Problem: empty qregion. try add -b to the command line", logfile = LOGFILE)
            else:
                with Pool(5) as p:
                    repository = p.map(processCoverage, qregions)
            
            if repository:
                qregions_next = [b[0] for b in repository if (b!=None and b!="NO COVERAGE")]
                qregions_failed += [qregions[q] for q in list(locate(repository,lambda x:x=="NO COVERAGE"))]
                logreport("%s has no coverage"%qregions_failed, logfile=LOGFILE)
                qregions = qregions_next
            else:
                break
            logreport("N.calls %d"%len(qregions),logfile=LOGFILE)
            #print_output(qregions,c)
  
    if repository:
        #convertHMM(repository[0][6],repository[0][-1])
        
        plotter(repository, output_px, qregions_failed)
logreport("CoverageMaster is done",logfile = LOGFILE)
