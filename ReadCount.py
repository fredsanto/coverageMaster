##!!just requires bed input to be sorted
from optparse import OptionParser  
import sys,os

def index_ucsc(filename):
    print >>sys.stderr, "Indexing %s \n"%filename
    chr_list = []
    f=open(filename)
    for l in f:
        if l[0] == 'c':
            chr = l.split()[0]
        else:
            chr = l.split()[2]
        if chr not in chr_list:
            chr_list.append(chr)
    f.close()
    return chr_list

def chrconv(chr):
    return "chr"+chr.replace("chr","")

class Bisect():
    '''
    Bisection algorithm to find a region (chr:start-end) in a sorted bed file
    init with the filename
    find (chromosome, start) to get a pointer to the region (need further stream)
    '''
    chr_list=('chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY')
    start = len(chr_list)/2
    end = len(chr_list)
    seek_chr = ""
    curr_chr = "chr1"
    size = 0
    fname = ""
    f = None
    pos = 0
    def __init__(self, fname, feeder=None):
        self. size = os.path.getsize(fname)
        self.fname = fname
        self.f = open(self.fname)
        if feeder:
            self.chr_list = index_ucsc(feeder)
        pass
    def create_index(self):
        new_chr_list = []
        self.f.seek(0)
        byte = self.f.read(100)
        while(byte):
            byte = self.f.read(5000)
            l = self.f.readline():
            if l:
                chr = l.split()[0]
                chr = "chr"+chr.replace("chr","")
                if chr in self.chr_list and chr not in new_chr_list:
                    new_chr_list.append(chr)
        self.chr_list = new_chr_list
        self.f.seek(0)
        return
    
    def _which_half(self,chr):
        try:
            if self.chr_list.index(chr)>self.chr_list.index(self.seek_chr):
                return -1
            return 1
        except:
            print >> sys.stderr, "Bad chromosome:%s"%chr
            return 1
    def find(self, chr, start, save_pos = False):  
      
      s = int(start)
      f = self.f
      DIST = 1e0
      SMALL_JUMP = 5e5
      try:
          self.seek_chr = self.chr_list[self.chr_list.index(chr)]
      except:
          print >> sys.stderr, "%s not in my list - skipped"%chr
          return None
      if save_pos and self.seek_chr==self.curr_chr:
          pos = self.pos
          inc = SMALL_JUMP
      else:
          pos = self.size/2
          inc = self.size/4
      oldpos = 0
      
      
      while 1:
        pos = (pos>0)*pos
        f.seek(pos), f.readline()
        _tmp_chr = ""
        while _tmp_chr not in self.chr_list:   #unkn chr - read till a known chr is found
            try:
                _tmp_chr,_tmp_s,_tmp_e = f.readline().split()[:3]
                _tmp_chr = chrconv(_tmp_chr)
            except:
                break
        #pos = f.tell()
        
        self.curr_chr = _tmp_chr
        if _tmp_chr == self.seek_chr:
          if  (int(_tmp_s) +DIST) < s:
            while (int(_tmp_s) +DIST) < s and _tmp_chr == self.seek_chr:
                pos+=SMALL_JUMP
                f.seek(int(pos)), f.readline()
                try:
                    _tmp_chr,_tmp_s,_tmp_e = f.readline().split()[:3]
                    _tmp_chr = chrconv(_tmp_chr)
                except:
                    pos-=SMALL_JUMP
                    f.seek(int(pos)), f.readline()
                    ##close to the end of the file - proceeding line by line
                    #while int(_tmp_s)<s and _tmp_chr == self.seek_chr:
                    #    _tmp_chr,_tmp_s,_tmp_e = f.readline().split()[:3]
                        
                    print >>sys.stderr,"Ouch!! It is the end of file"
                    self.pos = f.tell()
                    return self.f
                    break                
            else:
                pos-=SMALL_JUMP
                pos = (pos>0)*pos
                f.seek(int(pos)), f.readline()
                break
          if (int(_tmp_s) +DIST) >= s:
            while int(_tmp_s) > s and _tmp_chr == self.seek_chr :
                oldpos=pos
                pos-=SMALL_JUMP
                pos = (pos>0)*pos
                f.seek(int(pos)), f.readline()
                if pos<=0:
                    return self.f
                _tmp_chr,_tmp_s,_tmp_e = f.readline().split()[:3]
                _tmp_chr = chrconv(_tmp_chr)
            else:
                #pos = oldpos
                while _tmp_chr != self.seek_chr:   #unkn chr - read till chr is found - assuming that 
                    try:
                        _tmp_chr,_tmp_s,_tmp_e = f.readline().split()[:3]
                        _tmp_chr = chrconv(_tmp_chr)
                    except:
                            break
                break
        oldpos = pos
        pos = pos + (self._which_half(_tmp_chr))*inc
        inc = inc / 2
        if inc == 0:
            break
      self.pos = f.tell()
      return self.f
    def close(self):
      self.f.close()
      
class Regions(Bisect):
    f_pos = None
    mem_file = []

    ###init inherited
    
    def focus(self, (chr,start,end)):
        '''extract all entries in these boundaries''' 
        self.mem_file = []
        chr = chrconv(chr)
        f = self.find(chr,start)
        f.readline()
        for l in f:
          chr_r,e_start,e_end = l.split()[:3]
          chr_r = chrconv(chr_r)
          if chr == chr_r:
            if (start <= int(e_start) and end >=int(e_end)) or (start >= int(e_start) and end <=int(e_end) ):
                self.mem_file.append('%s\t%s\t%s'%(chr_r,e_start,e_end))
            if int(e_start) > end:
                #f.close()
                return

    def focus_single(self, (chr,start,end) , step=1):
      '''extract all entries between given boundaries''' 
      self.mem_file = []
      chr = chrconv(chr)
      f = self.find(chr,int(start))
      count = 0
      if f:
        f.readline()
        for l in f:
          count+=1
          chr_r,e_start = l.split()[:2]
          chr_r = chrconv(chr_r)
          if chr == chr_r:
            if int(start) <= int(e_start) and int(end) >= int(e_start):
                #self.mem_file.append(l.strip())
                if count % step == 0:
                    value = l.strip()
                    yield value
            if int(e_start) > end:
                #f.close()
                return
          else:
              pass

    def into(self, (chr,start,end)):
        start = int(start)
        end = int(end)
        off = (end-start)/2
        for l in self.mem_file :
                chr_r,e_start,e_end = l.split()[:3]
                chr_r = chrconv(chr_r)
                if chr==chr_r:
                    if start>=int(e_start)-off and end<=int(e_end)+off:
                        return 1
                    if int(e_start) > end:
                        self.mem_file = self.mem_file[self.mem_file.index(l)-100:]
                        return 0
        return 0

      
