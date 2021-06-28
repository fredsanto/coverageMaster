

import sys, time


LOGFILE = ""

def logreport(text, logfile):
    localtime = time.asctime( time.localtime(time.time()))
    sys.stderr.write("%s - %s\n"%(localtime,text))
    logfile.write("%s - %s\n"%(localtime,text))
    sys.stderr.flush()
    logfile.flush()

