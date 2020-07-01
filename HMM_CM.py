import sys
from numpy import *
from collections import defaultdict
from scipy.stats import norm
import pywt

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
    for O_state in O.o_states:
        s1 = O_state
        if s1 == max(O.o_states):
            if s1<2:
                p1 = .5-abs(norm(s1,Nstd1).cdf(os1)-0.5)  
            else:
                p1 = .5
        else:
            p1 = .5-abs(norm(s1,Nstd1).cdf(os1)-0.5) 
        p = p1
        res.append(p)
    return array(res)


def _reduce(vector):
    ret_vector = []
    str_vector = map(str,vector)
    for n,v in enumerate(str_vector):
        if str_vector[n:].count(v) == 1:
            ret_vector.append(vector[n])
    return ret_vector
    
class Observations:
    def __init__(self, o_states):
        self.o_states = o_states
        


###Viterbi
def downsample(signal,lev):
    coeffs = pywt.wavedec(signal,"haar",level = lev)
    ds = pywt.upcoef('a',coeffs[0],'haar',level = 1)
    ds = ds/mean(ds)*mean(signal)
    return ds

def upsample(signal,ln):
    return repeat(signal,ln/len(signal)+1)

def HMM(signal,std,booster = 1):
    off = .1
    lev = 6
    if len(signal)>1e6:
        lev = 12
    pos_states  = downsample(signal,lev)
    o1 = [1,1.5,0.5]#[1,1+2*median(std)+2*off,1-2*median(std)-2*off] ##states
    Obs = Observations(o1) ##to adjust
    v = zeros(len(Obs.o_states))
    pij = zeros((len(Obs.o_states),len(Obs.o_states)))
    M = zeros(len(pos_states))
    path_max = []
    v = 1*e(Obs.o_states[0],Obs)
    a = 1e-7
    for i,oi in enumerate(Obs.o_states):
        for j,oj in enumerate(Obs.o_states):            
            N = oi
            pij[i,j] = a*booster#0.00000000001#/(1+d(oi,oj))
            if i==j:
                pij[i,j] = 1-a*booster #.99999999998#/(1+d(oi,oj))
    
    fp_list = []
 
    count = 0
    for t in range(1,len(pos_states)):
        print "%d\r"%count,
        sys.stdout.flush()
        _v = v*pij
        _m = argmax(_v,1)
        path_max.append(_m)
        _em = e(pos_states[t],Obs,std[t])
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
    return upsample(fp_list,len(signal))[:len(signal)],o1
 
##save_results
def distance(_fp, pos_meas):
    min_dist = None
    for n,genome in enumerate(_fp):
      dist = 0
      for CHR in genome.keys():
       for (path,hmmlev) in genome:
        res = array(map(lambda x:hmmlev[x],path))
        ps = array(map(lambda x:x[1:],pos_meas[CHR]))
        dist += sum(sum((ps - res)**2))
      if not min_dist or min_dist[1] > dist:
          min_dist = (n,dist,genome)
    return min_dist
    
##n,d,final = distance(total_fp_list,pos_meas) #final = final_o1,final_path0


