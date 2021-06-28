import sys
from numpy import *
from collections import defaultdict
from scipy.stats import norm
import pywt
from log import logreport


def generate_states(off):
    o = []
    for i in range(1,9):
        i = float(i)
        for comb in range(0,floor(i/2)+1):
            SCNA = log2(i/2*(1-off)+off)
            LOH = ((i-comb)/i*(1-off)+0.5*off-0.5)*100 +(i==1)*4
            o.append((LOH,SCNA))
    return o

def inputread(filename, individual):
    f = open(filename).read().split('\n')
    mg0 = []
    mg1 = []
    posmes = defaultdict(list)
    header = f[0].strip()+'\n'
    #phenotype = array(header.split()[3:],'int')
    oldpos = 0
    for l in f[1:]:
      if individual in l:
        if not l:
            break
        try:
            ll = l.strip().split()
            chr, pos, end = ll[:3]
            meas = abs(float(ll[3]))
            cov = (float(ll[4]))
        #mg0.append(meas)
        #mg1.append(meas[g1])
            if oldpos and int(pos)<oldpos:
             pass
            posmes[chr].append((pos,meas,cov))
            oldpos = int(pos)
        except:
            pass
    return posmes, header

def inputread_old(filename):
    f = open(filename).read().split('\n')
    mg0 = []
    mg1 = []
    pos = []
    header = f[0].strip()+'\n'
    phenotype = array(header.split()[3:],'int')
    g1 = nonzero(phenotype)[0]
    g0 = array(list(set(range(0,len(phenotype))).difference(g1)))
    for l in f[:]:
        if not l:
            break
        ll = l.strip().split()
        p = ll[:3]
        meas = array(ll[3:],int)
        mg0.append(meas)
        #mg1.append(meas[g1])
        pos.append(p)
    return pos, header, mg0, (g0,g1)

def d(i,j):
    return abs(i-j)
    thr = 0.5
    _d = sum(i!=j)/float(len(i))
    return _d#*(_d>=thr)+1*(_d<thr)



def e(i,O,Nstd1 = 7):
    res = []
    os1 = i if i<2 else 2
    if 0.9<i<1.1:
        return(array([.5,0,0]))
    for O_state in O.o_states:
        s1 = O_state
        if s1 == max(O.o_states):
            if s1<2:
                p1 = .5-abs(norm(s1,Nstd1).cdf(os1)-0.5)  
            else:
                p1 = .5
        else:
            p1 = .5-abs(norm(s1,Nstd1).cdf(os1)-0.5) 
            
        res.append(p1)
    return array(res)


def _reduce(vector):
    ret_vector = []
    str_vector = list(map(str,vector))
    for n,v in enumerate(str_vector):
        if str_vector[n:].count(v) == 1:
            ret_vector.append(vector[n])
    return ret_vector
    
class Observations:
    def __init__(self, o_states):
        self.o_states = o_states
        


###Viterbi
def downsample(signal,lev):
    if lev == 0:
        return signal
    baseline = 1e-5
    wvfam = "haar"
    #ds = pywt.wavedec(signal+baseline,wvfam,level = lev)[0]
    ds = pywt.downcoef('a',signal+baseline,wvfam,level = lev)
    #ds = pywt.upcoef('a',coeffs[0],wvfam,level = 1)
    #ds = (ds-min(ds))*(max(signal)-min(signal))/(max(ds)-min(ds)) + min(signal)
    #if signal[0]:
    ds = ds/median(ds)*median(signal)
    return ds

def upsample_trigger(signal,ln):
    signal = around(signal,1)
    signal_pos = where(signal!=1)[0]
    try:
        signal[signal_pos-1]=0.5
        signal[signal_pos+1]=0.5
    except:
        pass #boundaries not amplified
    return repeat(signal.tolist(),ceil(ln/len(signal)))

def upsample(signal,ln):
    try:
        signal = signal.tolist()
    except:
        pass
    return repeat([1]+signal,ceil(ln/len(signal)))
    
def HMM_long(signal,std_o,gene, booster = 1, lev = 5, LOGFILE = "", mask = None, minlev = 0):

    #minlev = 3 HARDWIRED
    
    off = .1
    if len(signal)>1e5:
        lev = 12 #high compression for long regions
    
    o1 = [1,1.5,0.5]#[1,1+2*median(std)+2*off,1-2*median(std)-2*off] ##states
    Obs = Observations(o1) ##to adjust
    v = zeros(len(Obs.o_states))
    pij = zeros((len(Obs.o_states),len(Obs.o_states)))
    
    v = 1*e(Obs.o_states[0],Obs)
    a = 4.5e-6# 1e-7#park se o et al. Nat Gen 2010
    for i,oi in enumerate(Obs.o_states):
        for j,oj in enumerate(Obs.o_states):            
            N = oi
            pij[i,j] = a*booster
            if i==j:
                pij[i,j] = 1-a*booster
    
    fp_list = []
    ZOOM = True
    BEGIN =True

    while(ZOOM): ##### ZOOOMING ENABLED!
        pos_states  = downsample(signal,lev)
        std  = downsample(std_o,lev)
        M = zeros(len(pos_states))
        path_max = []
        _em_v = [e(pos_states[t],Obs,std[t]) for t in range(0,len(pos_states))]
        if BEGIN:
            trigger = 1* array([1.5*(e[0]<e[1]) or 0.5*(e[0]<e[2]) or 1 for e in _em_v])
            #sensitive to del only
            #trigger = 1* array([0.5*(e[0]<e[2]) or 1 for e in _em_v])
            BEGIN = False
        else:
            trigger = downsample(trigger,lev)
        count = 0
        if not sum(abs(trigger-1)):
            trigger = upsample(trigger,len(signal))[:len(signal)]
            logreport( "%s:nothing there"%gene, logfile =LOGFILE)
            return(trigger,o1)
        for t in range(1,len(pos_states)):
            #print("%d\r"%count, end=' ')
            #sys.stdout.flush()
            _v = v*pij
            _m = argmax(_v,1)
            path_max.append(_m)
            _em = _em_v[t]
            v = array([_v[i,_m[i]] *_em[i] for i in range(0,len(_m))])
            v = v /max(v)
            count+=1
        #backtrack
        M = argmax(v)
        i = len(pos_states)-2
        final_path = []
        final_path.append(o1[M])
        while i>-1:
            M = path_max[i][M]
            final_path.append(o1[M])
            i=i-1
        fp_list = final_path[::-1]
        
        approx = upsample(fp_list,len(signal))[:len(signal)]
        if mask is not None:
            trigger = upsample_trigger(trigger,len(signal))[:len(signal)]
            if sum(trigger*mask) and sum((approx-1)*trigger*mask) == 0:
                ZOOM = True
                signal = signal/median(signal) ### signal centering after first inspection (no large aberration detected)
                
                signal = signal*(trigger!=1)+1*(trigger==1)
                lev = lev - 1
                if lev < minlev:
                    break
                logreport( "%s:Zooming %d"%(gene,lev), logfile =LOGFILE)
            else:
                ZOOM = False
        else:
            break #no zooming for regions
        if lev > 6: 
            break #no zooming for big regions
    return(approx,o1) 

def HMM(signal,std,gene, booster = 1, lev = 5, LOGFILE = "", mask = None, minlev = 0):
    #lev is the parameter -l
    off = .1
    if len(signal)>1e6:
        lev = 12
    
    o1 = [1,1.5,0.5]#[1,1+2*median(std)+2*off,1-2*median(std)-2*off] ##states
    Obs = Observations(o1) ##to adjust
    
    pos_states  = downsample(signal,lev)
    
    _em_v = [e(pos_states[t],Obs,std[t]) for t in range(0,len(pos_states))]
    trigger =1* array([1.5*(e[0]<e[1]) or 0.5*(e[0]<e[2]) or 1 for e in _em_v])
    
    trigger = upsample_trigger(trigger,len(signal))[:len(signal)]
    mtrigger = (trigger!=1)*(signal!=1)
    
    return (trigger*mtrigger+trigger*invert(mtrigger),o1)


'''
##save_results
def distance(_fp, pos_meas):
    min_dist = None
    for n,genome in enumerate(_fp):
      dist = 0
      for CHR in list(genome.keys()):
       for (path,hmmlev) in genome:
        res = array([hmmlev[x] for x in path])
        ps = array([x[1:] for x in pos_meas[CHR]])
        dist += sum(sum((ps - res)**2))
      if not min_dist or min_dist[1] > dist:
          min_dist = (n,dist,genome)
    return min_dist
    
##n,d,final = distance(total_fp_list,pos_meas) #final = final_o1,final_path0
'''
