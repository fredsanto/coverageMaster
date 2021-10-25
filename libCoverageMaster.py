from sympy import Union, Interval, SympifyError
import sys,os,re
from HMM_CM import *
from ReadCount import *
import time

from multiprocessing import Pool, Manager
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from log import logreport

wd = os.path.dirname(os.path.abspath(__file__))
reference = open(wd+"/REF/REFSEQ_hg19.chrM.tab.HUGO").read().strip().split('\n')
exon_reference = defaultdict(list)

for r in open(wd+"/REF/hg19.exons.merged.bed").read().strip().split('\n'):
    rchr,rstart,rend = r.split()
    exon_reference[rchr].append((rstart,rend))


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


def plotter(repository,output_px, qregions_failed,cgd={}, dgv_xplr=None):

  report = open(output_px+".CMreport","w")       
  report2 = open(output_px+".CMcalls","w")       
  cgd_inh = ""
  
  with PdfPages(output_px+'.CMpositives.pdf') as pp:    
    for box in repository:
      if box and box!="NO COVERAGE":
        gene,signal,csignal,enlight,stdm,stdM,approx,ref_exon_avg,ratio,call = box
        callstr_txt = ""
        if gene["gene"] in cgd.keys():
            cgd_inh = cgd[gene["gene"]]
        else:
            cgd_inh = ""
        plt.subplot(4,1,1)
        if dgv_xplr is not None:
            overlaps = dgv_xplr.get_overlap(chrom=gene['chr'],
                                            alpha_0=int(gene['start']),
                                            beta_0=int(gene['end']))
            if overlaps is not None:
                freq_results = dgv_xplr.get_freq(overlaps)
                overlaps["Gain freq."] = freq_results["freq_g"].tolist()
                #print(overlaps["Gain freq."],overlaps["Name."])
                overlaps["Loss freq."] = freq_results["freq_l"].tolist()
                #print(overlaps["Loss freq."])
                overlaps["Name."] = freq_results["name"].tolist()
        else:
            overlaps = None

        try:
            for c in call:
                callstr_txt = '%s\t%s\t%s\t%s\t%s\t%s'%(gene['gene'], cgd_inh,gene['chr'],c[0],c[1],c[2])
                report2.write("%s" % callstr_txt)
            if overlaps is not None:
                report2.write("\tGains")
                for gain, name in zip(overlaps["Gain freq."].values, overlaps["Name."].values):
                    report2.write("\t%s\t%.3f" % (name,gain))
                report2.write("\tLosses")
                for loss,name in zip(overlaps["Loss freq."].values, overlaps["Name."].values):
                    report2.write("\t%s\t%.3f" % (name,loss))
                report2.write("\n")
            else:
                report2.write("\n")
        except:
            pass
        plt.title(callstr_txt.replace("\t"," "))
        plt.plot(enlight, color = 'g', linewidth=1.0 )
        plt.xticks([], [])
        plt.yticks([], [])
        plt.ylim(-3,1)
            
        plt.subplot(4,1,2)
        
        plt.plot(stdM, 'r')
        plt.plot(stdm, 'r')
        plt.ylim(-.5,2.5)
        
        #if options.control:
        plt.plot(csignal, color = '#111111', linewidth=0.8,linestyle ="--" )
        plt.plot(signal, color = 'b', linewidth=1.0 )
        plt.plot(approx, 'k', linewidth=1.5)
        plt.plot(ref_exon_avg,'k',linewidth=1.5)
        plt.xticks([], [])
        plt.yticks([], [])
        
        
        
        if 1:#options.control:
            plt.subplot(4,1,3)
            plt.axhline(y=0)
            plt.plot(ratio - 1, 'k', linewidth=1.5)
            axes = plt.gca()
            axes.set_ylim([-1,1])
        

        enlight = array(enlight)
        roi = where(((approx-1)*(enlight!=0))!=0)[0]
        if len(roi):
            plt.subplot(4,1,4)
            

            plt.plot(stdM, 'r')
            plt.plot(stdm, 'r')
            plt.ylim(-.5,2.5)
            plt.xlim(min(roi)-100,max(roi)+100)

            plt.plot(csignal, color = '#111111', linewidth=0.8,linestyle ="--" )
            plt.plot(signal, color = 'b', linewidth=1.0 )
            plt.plot(approx, 'k', linewidth=1.5)
            plt.plot(ref_exon_avg,'k',linewidth=1.5)
        #convert_approx_to_genome()        
        report.write('%s\n'%gene['gene'])
        plt.gcf().set_size_inches([7,10])
        plt.savefig(pp, format='pdf')
        plt.close()
  for q in qregions_failed:
      report2.write('%s\t%s\t%s\t%s\n'%(q['chr'], q['start'], q['end'], "NO COVERAGE"))    
  

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

def gene_reference(genes, reference, LOGFILE, nexons = 10):
 
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
        exstart = list(map(int,g.split('\t')[8].split(',')[:])) if len(g.split('\t')[8].split(','))>1 else int(g.split('\t')[8])
        exend = list(map(int,g.split('\t')[9].split(',')[:])) if len(g.split('\t')[9].split(','))>1 else int(g.split('\t')[9])
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

def signalProcessor_old(gene, cr, cref, stats,LOGFILE, red=False):
 
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

def convertHMM(approx, data_n):
        if not len(data_n):
            return("NODATAN")
        call = []
        dd = diff([1]+approx.tolist()+[1])
        transitions = where(dd!=0)[0]
        txp = [transitions[i:i+2] for i in range(0,len(transitions),2)]
        for n,el in enumerate(txp):
            if len(el)==1: ## adiacent del-dup end==start 
                el = (txp[n-1][1],el[0])
            if dd[el[0]]<dd[el[1]]:
                mut="DEL"
            if dd[el[0]]>dd[el[1]]:
                mut="DUP"
            if dd[el[0]]==dd[el[1]]:
                mut="ERROR"
            try:
                call.append((mut,data_n[el[0]][0],data_n[el[1]-1][0]))
            except:
                pass
                print("Problem %s - %s"%(str(el),str(data_n)))
                sys.stdout.flush()
                return None
        return call

