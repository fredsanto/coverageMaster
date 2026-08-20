import sys
import time


def logreport(text, logfile):
    localtime = time.asctime(time.localtime(time.time()))
    msg = f"{localtime} - {text}\n"
    sys.stderr.write(msg)
    logfile.write(msg)
    sys.stderr.flush()
    logfile.flush()
