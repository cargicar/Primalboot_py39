import datetime
import sys
from queue import Queue # Carlos add
#import Queue

timing_q = Queue()

def t0():
    global timing_q
    timing_q.put(datetime.datetime.now())

def deltat():
    t0 = timing_q.get()
    t1 = datetime.datetime.now()
    return (t1-t0).seconds+((t1-t0).microseconds*1e-6)
