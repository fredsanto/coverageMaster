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
#from scipy import stats
from HMM_CM import *
from ReadCount import *
import time
from sympy import Union, Interval, SympifyError
import pickle
        


def logreport(text, logfile):
    localtime = time.asctime( time.localtime(time.time()) )
    print("%s - %s"%(localtime,text), file=sys.stderr)
    print("%s - %s"%(localtime,text), file=logfile)
    sys.stderr.flush()
    logfile.flush()

wd = os.path.dirname(os.path.abspath(__file__))
reference = open(wd+"/REFSEQ_hg19.chrM.tab.HUGO").read().strip().split('\n')
exon_reference = defaultdict(list)

for r in open(wd+"/hg19.exons.merged.bed").read().strip().split('\n'):
    rchr,rstart,rend = r.split()
    exon_reference[rchr].append((rstart,rend))


qregion = {}
qregions = []


def extract_tot_reads(stats):
    tot_targ_reads =int(re.findall("([0-9]*) \+ 0 *mapped",stats)[0])
    return tot_targ_reads/100e6
    

def wig_writer(chr, data):
    wig = []
    span = 20.0
    head1 = "track type=wiggle_0 name=p_full"
    wig.append(head1)
    
    for snp in data:
      chr, pos, value = chr, snp[0], snp[2]
      chr = "chr"+chr.replace("chr","")
      if 1 or float(value)>2:
        head2 = "variableStep chrom=%s span=%d"
        wig.append(head2%(chr,2*span))
        slope = 0
        for p in range(int(pos)-int(span), int(pos)+int(span)):
            if p!=float(pos):
                #slope += 1/span
                slope = 0
            else:
                slope = 1
            wig.append("%d\t%f"%(int(p), (slope*value)))
    return wig




def exon_boundaries(chr,start,end,exon_reference, nexons=10):
    offset = nexons
    START =  False
    END = False
    for n,(estart,eend) in enumerate(exon_reference[chr]):
        if int(start)<=int(estart) and not START:
            if n>=offset:
                START = exon_reference[chr][n-offset][0]
            elif n<offset:
                START = exon_reference[chr][0][0]
        if int(eend)>=int(end) and not END:
            if n+offset < len(exon_reference[chr]):
                END = exon_reference[chr][n+offset][1]
            else:
                END = exon_reference[chr][-1][1]
    return {'chr':chr,'start':START,'end':END}

def main_trans(transcripts):
    for t in transcripts:
        pass
    return []

def union(data):
    """ Union of a list of intervals e.g. [(1,2),(3,4)] """
    intervals = [Interval(begin, end) for (begin, end) in data]
    u = Union(*intervals)
    return [u] if isinstance(u, Interval) \
       else list(u.args)

def gene_reference(genes, reference, nexons = 10):

 if 0 and 1000<len(genes)<5000 and os.path.isfile(sys.path[0]+"/codedRef"):
     with open(sys.path[0]+"/codedRef",'rb') as fp:
         gs = pickle.load(fp)
         return(gs)
 else:
    gf= []
    genes = sorted(list(set(genes)))
    for g in genes:
        try:
            
            
            #transcripts = map(lambda g: (filter(lambda r:g == r.split()[-1] and 'NR' not in r.split()[0],reference)), [g])[0]
            transcripts = [([r for r in reference if g == r.split()[-1]]) for g in [g]][0]
            #if not transcripts:
            #    transcripts = map(lambda g: (filter(lambda r:g == r.split()[0] ,reference)), [g])[0]
            #gf.append(transcripts[argmax(map(lambda x:int(x.split()[7]),transcripts))])
            #check chromosome
            transcripts_X = [t for t in transcripts if t.split()[1].replace("chr","")=="X"]
            transcripts_Y = [t for t in transcripts if t.split()[1].replace("chr","")=="Y"]
            transcripts = [t for t in transcripts if t not in transcripts_X+transcripts_Y]
            for txs in [transcripts,transcripts_X,transcripts_Y]:
                maxend = 0
                minstart = 9999999999999999
                tmax = []
                exonst = []
                exonen = []
                if txs:
                    for t in txs:
                        tt = t.split()
                        exonst += tt[8].strip(",").split(",")
                        exonen += tt[9].strip(",").split(",")
                        tmax = tt if len(tt)>len(tmax) else tmax 
                        start, end = int(tt[3]),int(tt[4])
                        maxend = max(maxend, end)
                        minstart = min(minstart, start)
                    
                    tmax[3] = str(minstart)
                    tmax[4] = str(maxend)
                    tmax[8] = ','.join(exonst)
                    tmax[9] = ','.join(exonen)
                    tfinal = ['\t'.join(tmax)]
                    if maxend-minstart< 3e6:
                        gf.extend(tfinal)
                        #print(g)
                    else:
                        logreport("the following transcript is too long:\n%s"%str(tfinal), logfile = LOGFILE)
        except:             
            logreport("%s not found in th current reference. Check the gene name"%g,logfile = LOGFILE)
            #sys.exit(1)
    gs = []
    for g in gf:
        
        exfirst = g.split('\t')[3]
        exlast = g.split('\t')[4]
        exstart = list(map(int,g.split('\t')[8].split(',')[:-1])) if len(g.split('\t')[8].split(','))>1 else int(g.split('\t')[8])
        exend = list(map(int,g.split('\t')[9].split(',')[:-1])) if len(g.split('\t')[9].split(','))>1 else int(g.split('\t')[9])
        newex = union(list(zip(exstart,exend))) if type(exend)==list else [Interval(exstart,exend)] #merge and combine exons from different txs
        exstart = [str(x.start) for x in newex]
        exend = [str(x.end) for x in newex]
        #print(g)
        try:
            replace_s = exstart.index([x for x in exstart if float(exfirst)<=float(x)][0])-1
        except:
            replace_s = len(exstart)-1
        try:
            replace_e = exend.index([x for x in exend if float(exlast)>=float(x)][-1])+1
        except:
            replace_e = 0
        replace_s = max(replace_s,0)
        replace_e = min(replace_e,len(exend)-1)
        exstart[replace_s] = exfirst
        exend[replace_e] = exlast
        expairs = list(zip(exstart[replace_s:replace_e+1],exend[replace_s:replace_e+1]))
        gdata = {'NM':g.split()[0],'chr':g.split('\t')[1],'expairs':expairs,'gene':g.split('\t')[-1],'start':exfirst,'end':exlast}
        query = exon_boundaries(gdata['chr'],gdata['start'],gdata['end'],exon_reference, nexons)
        gdata.update({'query':query})
        gs.append(gdata)
        
    return gs



usage = "usage: %prog [options] <cov_file> <stats_file> <gene list(file or gnames comma separated)|region(chr:start-end)> -r <reference.cov> -o <output_px>"
parser = OptionParser(usage = usage)
parser.add_option("-c",'--control',  dest="control", help = "<optional> control file (.cov with .report.txt in same folder)")
parser.add_option("-r",'--ref',  dest="ref", help = "reference file", default = "")
parser.add_option("-o",'--out',  dest="output_px", help = "output plot file", default = "")
parser.add_option("-f",'--force',  action="store_true", dest="force", default=False, help = "force output")
parser.add_option("-e",'--exons',  dest="exons", help = " <optional> n. of exons", default = 10)
parser.add_option("-x",'--offset',  dest="offset", help = " <optional> offset to ref", default = 0)
parser.add_option("-w",'--wig',  action="store_true", help = " <optional> write wig", default = False)
parser.add_option("-b",'--bed',  action="store_true", help = " <optional> BED input", default = False)
(options, args) = parser.parse_args()
try:
    output_px = options.output_px
    LOGFILE = open(output_px+".CM.log","a")
    offset = float(options.offset)
    cov_region = args[0]
    if '.cov' not in cov_region:
        print("\nMichel, please.. not again... use the .cov!\n")
        raise Exception
    
    stats_file = args[1]
    stats = extract_tot_reads(open(stats_file).read())
    WARNING = options.force
    XIST = gene_reference(['XIST'], reference, int(options.exons)).pop()
    if args[2][:-1] == "/":
        raise Exception("GeneList not valid")
    if options.bed: #it's a BED file in input
        for l in open(args[2]):
            qregion['chr'],qregion['start'],qregion['end'] = l.strip().split("\t")
            qregions.append(qregion)
        qregions += [XIST]
    elif ":" in args[2]:
        qregions = []
        s1 = args[2].split(':') 
        qregion['chr'] = s1[0]
        qregion['start'], qregion['end'] = s1[1].split('-')[0], s1[1].split('-')[1]
        qregions += [qregion,XIST]
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
                qregions = gene_reference(glist, reference, int(options.exons))+[XIST]
            ###here > dump if it is non existing
            try:
                if not os.path.exists("%s_coderef"%gfilename):
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
            if ":" in r:
                qregions += [XIST]
            else:
                qregions += XIST
    if not options.ref:
        print("\nInclude the reference")
#        raise Exception

    for f in [stats_file,options.ref,cov_region,options.control,args[2]]:
        print("%s - %d"%(f, os.path.isfile(f)))
 
    
except:
    print("Unexpected error:", sys.exc_info()[0])
    parser.print_help()
    sys.exit(0)


n_exons = 0
total_bp = 0

report = open(output_px+".CMreport","w")       
report2 = open(output_px+".CMexonstruct","w")       
cr = Regions(cov_region)
cr.create_index()
logreport( "Query index created", logfile =LOGFILE)

if options.ref:
    cref = Regions(options.ref)
    cref.create_index()
    logreport( "Reference index created", logfile = LOGFILE)
else:
    cref = None 
if options.control:
    ccov = options.control
    cstats = options.control.replace(".cov","")+".report.txt"
    ccont = Regions(ccov)
    ccont.create_index()
    cstats = extract_tot_reads(open(cstats).read())
    logreport( "Control index created", logfile = LOGFILE)
else:

    parser.print_help()
    sys.exit(0)
  


##main
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



outp_buffer = [] 
scales = [10,20,30,200]#,40,50,100,200]
full_coverage = [0,0,0,0]

def realign ( pos, ref_pos):
  n = 0
  for n,p in enumerate(ref_pos):
    for s,j,f in pos[n:]:
        if int(s)<int(p):
            pos.remove((s,j,f))
        if int(s)==int(p):
            break
        if int(s)>int(p):
            pos.insert(n,(p,0,0))
            break
  pos = pos[:n+1]
  pass
  return pos


def calculate_exon(gene, pos):
    try:
        exn = gene['expairs'].index([x for x in gene['expairs'] if int(x[0])<int(pos)<int(x[1])][0])
        if exn%2:
            return 1
        else:
            return -1
    except:
        return 0

def signalProcessor(gene, cr, cref, stats, red=False):
  
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
        size = (int(qregion['end']) - int(qregion['start'])) # ArryCGH mode
        step = 1 if size < 500e3 else 5
        if cref:
            ref_exon_min = []
            ref_exon_Max = []
            ref_exon_avg = []
            ref = []
            FREEZE = False
            gen_cr = cr.focus_single((qregion['chr'],int(qregion['start']),int(qregion['end'])), step = 1)
            for r in cref.focus_single((qregion['chr'],int(qregion['start']),int(qregion['end'])), step = step):
                rr = r.strip().split()
                chr,pos = rr[:2]
                if pos == "33048718":
                    pass
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
                        
                    except StopIteration:
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


with PdfPages(output_px+'.CMpositives.pdf') as pp:    
 ##XIST normalization for chrX
 XIST = qregions.pop()
 plot_chrX_region, unused, unused, unused, unused = signalProcessor(XIST, cr, cref, stats)   
 normX = median(plot_chrX_region)
 if options.control:
     control_chrX_region, unused, unused, unused, unused = signalProcessor(XIST, ccont, cref, cstats)   
     normCX = median(control_chrX_region)
 for gene in qregions:
  if gene['chr']!='chrM':   
   if not options.force:
       WARNING = False
   try:
    if cref:
        signal, ref_exon_avg, ref_exon_min, enlight, data_n = signalProcessor(gene, cr, cref, stats, red = False)   
        signal = signal + offset
    if options.control:
        csignal, unused, unused, unused, unused = signalProcessor(gene, ccont, cref, cstats)   
        csignal = csignal + offset
    if gene['chr'] == 'chrX':
        signal = signal/normX
        csignal = csignal/normCX+1e-4
    sd = array(ref_exon_min)/array(ref_exon_avg)
    genethere = array(enlight)!=0
    
    if sum(abs(signal))>0 and sum(abs(csignal))>0:
        
        stdM = 1+2*sd
        stdm = 1-sd
        noninfidx = where(abs(signal-1) <= sd)
        ratio = signal/(.001 + csignal)
        ratio[noninfidx] = 1
        approx,alev = HMM(ratio,sd)
        control,clev = HMM(sd/median(sd),median(sd)*ones(len(sd)), booster = 1000)
        
        maskp = (approx>1)
        maskm = (approx<1)
        if ((sum(maskp)+sum(maskm))>0) and (sum(signal*maskp*genethere - maskp*stdM*genethere) > 0) or (-sum(signal*maskm*genethere - maskm*stdm*genethere) > 0):
    
            cc = control.tolist()
            cc=[1]+cc+[1]
            steps = where(diff(cc)!=0)[0]
            for n,i in enumerate(steps):
                if n%2:
                    art_idx = steps[n-1]-500,steps[n]+500
                    art_idx = (max(0,art_idx[0]),min(len(signal)-1,art_idx[1]))
                    art_idx = list(range(art_idx[0],art_idx[1]+1))
                    control[art_idx] = 0

                    
            approxe = approx*genethere
            control = control*genethere

            if 1 or [x for x in approxe*control if x!=1 and x!=0]:
                WARNING = True
    
    
        if WARNING:
        ### wig
            if options.wig:
                wig  = wig_writer(gene["chr"],data_n)
                fwigout = open(output_px+'_'+gene['gene']+'.wig','w')
                fwigout.write('\n'.join(wig))
                fwigout.close()
            
            #plt.hold(True)
            plt.subplot(3,1,1)
            plt.title(gene['gene']+' '+gene['chr'])
            plt.plot(enlight, color = 'g', linewidth=1.0 )
            plt.xticks([], [])
            plt.yticks([], [])
            plt.ylim(-3,1)
            
            plt.subplot(3,1,2)
        
            plt.plot(stdM, 'r')
            plt.plot(stdm, 'r')
            plt.ylim(-.5,2.5)
        
            if options.control:
                plt.plot(csignal, color = '#111111', linewidth=0.8,linestyle ="--" )
            plt.plot(signal, color = 'b', linewidth=1.0 )
            plt.plot(approx, 'k', linewidth=1.5)
            plt.plot(ref_exon_avg/sum(ref_exon_avg),'k',linewidth=1.5)
            plt.xticks([], [])
            plt.yticks([], [])
        
        
        
            if options.control:
                plt.subplot(3,1,3)
                plt.axhline(y=0)
                plt.plot(ratio - 1, 'k', linewidth=1.5)
                axes = plt.gca()
                axes.set_ylim([-1,1])
        
            report.write('%s\n'%gene['gene'])
            try:
                report2.write('%s:%s-%s\n'%(gene['gene'],gene['chr'],gene['expairs']))
            except:
                pass
            plt.savefig(pp, format='pdf')
            plt.close()
    else:
        raise Exception("Gene %s has been skipped"%gene['gene'])
   except Exception as e:
       logreport("Error >%s< or Gene %s is problematic. Continuing"%(e,gene['gene']),logfile = LOGFILE)
logreport("CoverageMaster is done",logfile = LOGFILE)

