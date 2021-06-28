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
from multiprocessing import Pool,Manager
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from libCoverageMaster import *
from more_itertools import locate
from functools import partial



reference = open(wd+"/REFSEQ_hg19.chrM.tab.HUGO").read().strip().split('\n')
#reference = open(wd+"/REFSEQ_hg38_HGMD3").read().strip().split('\n')
exon_reference = defaultdict(list)

for r in open(wd+"/hg19.exons.merged.bed").read().strip().split('\n'):
#for r in open(wd+"/hg38.exons.merged.bed").read().strip().split('\n'):
    rchr,rstart,rend = r.split()
    exon_reference[rchr].append((rstart,rend))


qregion = {}
qregions = []

minlev = 0
cgd = {}




usage = "usage: %prog [options] <cov_file> <stats_file> <gene list(file or gnames comma separated)|region(chr:start-end)> -r <reference.cov> -o <output_px>"
parser = OptionParser(usage = usage)
parser.add_option("-c",'--control',  dest="control", help = "<optional> txt file with a control file per line (.cov with .report.txt in same folder)")
parser.add_option("-s",'--single-control',  dest="single_control", help = "<optional> single control file (.cov with .report.txt in same folder)")
parser.add_option("-r",'--ref',  dest="ref", help = "reference file", default = "")
parser.add_option("-o",'--out',  dest="output_px", help = "output prefix", default = "")
parser.add_option("-g",'--cgd',  dest="cgd", help = "clinical genomic database", default = "")
parser.add_option("-f",'--force',  action="store_true", dest="force", default=False, help = "force output")
parser.add_option("-e",'--exons',  dest="exons", help = " <optional> n. of exons", default = 5)
parser.add_option("-x",'--offset',  dest="offset", help = " <optional> offset to ref", default = 0)
parser.add_option("-d",'--width',  dest="wid", help = " <optional> width of std", default = 1) #choose width of sd
parser.add_option("-l",'--level',  dest="lev", help = " <optional> wavelet level", default = 5) #choose level of compression
parser.add_option("-m",'--minlevel',  dest="minlev", help = " <optional> min wavelet level", default = 0) #choose level of compression
parser.add_option("-w",'--wig',  action="store_true", help = " <optional> write wig", default = False)
parser.add_option("-b",'--bed',  action="store_true", help = " <optional> BED input", default = False)
(options, args) = parser.parse_args()
try:
    output_px = options.output_px.split("/")[-1]
    LOGFILE = open(output_px+".CM.log","a")
    offset = float(options.offset)
    cov_region = args[0]
    lev = int(options.lev)
    minlev = int(options.minlev)
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
    if options.cgd:
        with open(options.cgd) as f:
            for i in f:
                try:
                    gene,inh = i.strip().split()
                    cgd[gene] = inh
                except:
                    print(i)
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
                if not os.path.exists("%s_coderef"%gfilename) and len(qregions) > 0:
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
logreport( "-"*80+"\n"+"CoverageMaster is warming up", logfile =LOGFILE)  
cr = Regions(cov_region)

logreport( "Query index created", logfile =LOGFILE)

if options.ref:
    cref = Regions(options.ref)
    logreport( "Reference index created", logfile = LOGFILE)
else:
    cref = None 



##main
def processCoverage(terminal,gene,signalBuffer):#, cr, cref, stats, ccont, cref, cstats, XIST, pp)

 time.sleep(.1) 
 #print(signalBuffer.keys())
 if gene['chr']!='chrM':   
   sys.stdout.flush()
   if 1:#try:
    #reopen for threads purposes
    cref.f = open(cref.f.name)
    cr.f = open(cr.f.name)
    ccont.f = open(ccont.f.name)

    
    signal, ref_exon_avg, ref_exon_min, enlight, data_n = signalProcessor(gene, cr, cref, stats, signalBuffer, red = False, store = True)   
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
        logreport("Warning: NO Coverage in Control for %s"%gene, logfile=LOGFILE)
        return "NO COVERAGE"
    if len(signal)<len(csignal):
        logreport("Warning: Signal < Control %s - correction"%gene, logfile=LOGFILE)
        csignal = csignal[:len(signal)]
    if len(signal)>len(csignal): 
        logreport("Warning: Signal > Control %s - correction"%gene, logfile=LOGFILE)
        csignal= csignal.tolist()
        csignal = array(csignal+[1]*(len(signal)-len(csignal)))

    sd = array(ref_exon_min)/array(ref_exon_avg)
    genethere = array(enlight)!=0
   
    if len(signal)==len(csignal):
        genename = gene['gene']
        win = 30
        #wid: parameter -d
        stdM = 1+wid*sd
        stdm = 1-wid*sd
        unormsignal = signal-1
        #noninfidx = where((abs(unormsignal) <= wid*sd) + (unormsignal<=-wid*sd)*(csignal<signal) + (unormsignal>=wid*sd) * (csignal>signal))
        infidx = where((abs(unormsignal)>1/2)*((unormsignal>wid*sd)*(csignal<signal) + (unormsignal<-wid*sd)*(csignal>signal)))
        _tmp  = [arange(i-win,i+win) for i in infidx[0] if win<i<(len(signal)-win)]
        infidx = [i for el in _tmp for i in el]
        infidx = [*{*infidx}]# extend the ROI +/-30
        try:
            ratiofull = signal/(.001 + csignal)
        except:
            logreport("%s aprox failed. Ratio length inconsistency"%gene, logfile=LOGFILE)
            ratio = zeros(len(signal))
            
        #ratio[noninfidx] = 1
        _t = signalBuffer[genename] # workaround for a bug in python for parallel buffer variables
        if len(_t['infidx']) == 0:
            _t['infidx'] = infidx
            #_t['infidx'] =  [*{*_t['infidx']}]#unique(array(_t['infidx'])).tolist()
            #_t['noninfidx'] += noninfidx[0].tolist()
        else:
            _t['infidx'] = list(set(infidx).intersection(_t['infidx']))
        signalBuffer[genename] = _t
        ratio = ones(len(signal))
        ratio[signalBuffer[genename]['infidx']] = ratiofull[signalBuffer[genename]['infidx']]
        #print("sbl %d"%len(signalBuffer[gene['gene']]['noninfidx']))
        try:
            if terminal:
                approx,alev = HMM_long(ratio,sd,gene["gene"], lev = lev, LOGFILE=LOGFILE, mask=genethere, minlev = minlev)
            else:
                approx,alev = HMM_long(ratio,sd,gene["gene"], lev = lev, LOGFILE=LOGFILE, mask=genethere, minlev = minlev)
        except:
            logreport("HMM problem:%s"%gene["gene"],logfile=LOGFILE)
            return None
        #control,clev = HMM(sd/median(sd),median(sd)*ones(len(sd)), booster = 1000,lev = lev,LOGFILE=LOGFILE)
        #maskp = (approx>1)
        #maskm = (approx<1)
        mask = (approx!=1)
        #if FORCE or ((sum(maskp)+sum(maskm))>0 and (sum(signal*maskp*genethere - maskp*stdM*genethere) > 0) or (-sum(signal*maskm*genethere - maskm*stdm*genethere) > 0):
        if FORCE or (sum(mask)>0 and sum(signal*mask*genethere)!= 0):
            
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
            
            call = convertHMM(approx,data_n)
            
            return([gene,signal,csignal,enlight,stdm,stdM,approx,ref_exon_avg/sum(ref_exon_avg),ratio,call])
        else:
            signalBuffer.pop(gene["gene"])
            logreport("Removing %s"%gene["gene"],logfile=LOGFILE)
 return None

def signalProcessor(gene, cr, cref, stats, signalBuffer = None, LOGFILE=LOGFILE, red=False, store = False):
 
  if 0 < (int(gene['end'])-int(gene['start'])) and gene['chr']!='chrM':
    if 'query' in gene:
        qregion = gene['query'] 
        logreport( "Processing Gene: %s - Region %s:%s-%s"%(gene['gene'],gene['chr'],gene['start'],gene['end']), logfile = LOGFILE)
    else:
        qregion  = gene
        qregion['gene']='%s:%s-%s'%(gene['chr'],gene['start'],gene['end'])
        logreport("Processing Region %s:%s-%s"%(qregion['chr'],qregion['start'],qregion['end']),logfile = LOGFILE)
    data = []
    if qregion:
        enlight = []
        plot_region = []
        pos_region = []
        count = 0
        gpos = 0
        size = (int(qregion['end']) - int(qregion['start'])) # ArryCGH mode
        step = 1 if size < 500e3 else 5
        if not store or (gene['gene'] not in signalBuffer.keys()):
            ref_exon_min = []
            ref_exon_Max = []
            ref_exon_avg = []
            ref = []
            FREEZE = False
            FREEZE_B = False
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
                
                LOOP = False
                count=0
                while not LOOP and not FREEZE:
                    count+=1
                    if count>10000:
                        logreport("Error: Loophole in %s %s %s"%(chr,gpos,pos), logfile=LOGFILE)
                        break
                    try:
                        if not FREEZE_B:
                            chr, gpos, cov = next(gen_cr).strip().split()
                        if int(gpos) == (int(pos)+1):
                            FREEZE = True
                        if int(gpos) > (int(pos)+1):
                            FREEZE_B  = True
                            break
                    except:
                        chr, pos, cov = qregion['chr'], pos, 0
                        break
                    LOOP = pos == gpos
                if int(gpos) < int(pos):
                    FREEZE = False
                    FREEZE_B = False
                if pos == gpos:
                    FREEZE = False
                    FREEZE_B = False
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
        
        
            plot_region_n = array(plot_region)/array(ref_exon_avg)
            plot_region_n = array([1 if (isnan(i) or isinf(i)) else i for i in plot_region_n])
            data_n = [(i[0],i[1],plot_region_n[n]) for n,i in enumerate(pos_region)]
            if store:
                signalBuffer[gene["gene"]] = {'plot_region_n':plot_region_n,'ref_exon_avg':ref_exon_avg, "ref_exon_min":ref_exon_min, 'enlight':enlight, 'data_n':data_n, 'noninfidx':[], 'infidx':[]}
            return plot_region_n, ref_exon_avg, ref_exon_min, enlight, data_n
    sB = signalBuffer[gene["gene"]]
    return sB['plot_region_n'], sB['ref_exon_avg'], sB['ref_exon_min'], sB['enlight'], sB['data_n']

def print_output(glist,fname,round):
    pass

if __name__ == '__main__':
    
    XIST = qregions.pop()
    plot_chrX_region, unused, unused, unused, unused = signalProcessor(XIST, cr, cref, stats, LOGFILE)   
    normX = median(plot_chrX_region)
    qregions_failed = []
    try:
        signalBuffer = Manager().dict() 
        if options.control:
            cc = open(options.control).read().strip().split("\n")
        if options.single_control:
            cc = [options.single_control]
        for n,c in enumerate(cc):
            
          if c.split("/")[-1] != cov_region.split("/")[-1]:
            terminal =  n==len(cc)-1  #terminal control loop > Zooming
            cstats = c.replace(".cov","")+".report.txt"
            ccont = Regions(c)
            
            cstats = extract_tot_reads(open(cstats).read())
            logreport( "Control index created", logfile = LOGFILE)

            control_chrX_region, unused, unused, unused, unused = signalProcessor(XIST, ccont, cref, cstats, LOGFILE)
            normCX = median(control_chrX_region)
 
            logreport("Loop with control %s"%c, logfile = LOGFILE)
            processCoverage(terminal, qregions[0],signalBuffer) #debug mode
            pfunc = partial(processCoverage, terminal, signalBuffer=signalBuffer)
            if sys.version_info[0] < 3:
                from contextlib import closing
                with closing(Pool(processes=5)) as p:
                    if len(qregions):
                        if type(qregions)==type("str"):
                            pass
                        
                        repository = p.map(pfunc, qregions)
                        p.terminate()
                    else:
                        logreport("Problem: Invalid qregion. try add -b to the command line", logfile = LOGFILE)
            else:
                with Pool(5) as p:
                    repository = p.map(pfunc, qregions)
            
            if repository:
                qregions_next = [b[0] for b in repository if (b!=None and ("NO COVERAGE" not in b))]
                #passing_genes = [q["gene"] for q in qregions_next]
                qregions_failed += [qregions[q] for q in list(locate(repository,lambda x:"NO COVERAGE" ==  x))]
                logreport("%s has no coverage"%qregions_failed, logfile=LOGFILE)
                qregions = qregions_next
            else:
                break
            logreport("N.calls %d"%len(qregions),logfile=LOGFILE)
            #print_output(qregions,c)
          if FORCE or len(qregions)==0:
              break
    except:
        raise("Error in __main__")
    if repository:
        #convertHMM(repository[0][6],repository[0][-1])
        
        plotter(repository, output_px, qregions_failed, cgd)
logreport("CoverageMaster is done",logfile = LOGFILE)

        