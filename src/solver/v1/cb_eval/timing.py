import datetime
import sys
import Queue

timing_q = Queue.Queue()

def t0():
    global timing_q
    timing_q.put(datetime.datetime.now())

def deltat():
    t0 = timing_q.get()
    t1 = datetime.datetime.now()
    return (t1-t0).seconds+((t1-t0).microseconds*1e-6)
