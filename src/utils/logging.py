# import logging
import datetime
import time
# import config  # Carlos Edit
import os

import socket
import random
import resource

# import stats as st
# import memory as mem

from . import stats as st # Carlos edit
from . import memory as mem # Carlos edit

# can't use this until we setup network aware logging!
def setup_logging():
    logging.basicConfig(format='%(asctime)s %(message)s',level=logging.DEBUG)


PRINT=0
CRITICAL=0
ERROR=1
WARNING=2
INFO=3
PERF=4
DEBUG=5
DEBUG2=6
logfile='../logs/log_ope_lxplus'
logpersig='../logs/log_per_sig'

proc_id='no_proc_id'
# should only be initialized once for each calling process
uniq_id=None

def log(severity, msg):
    global uniq_id
    if uniq_id == None:
        host=socket.gethostname()
        random.seed()
        uniq_id = host + '-' + str(random.randint(0, 100000))

    if severity <= DEBUG2:                                      # Carlos edit
#         dirName='../logs/'
#         try:
#             os.makedirs(dirName)    
#             print("Directory " , dirName ,  " Created ")
#         except FileExistsError:

        fl=logfile + str(datetime.date.today())
        f=open(fl, 'a')
        #f.write('('+str(proc_id)+') '+str(datetime.datetime.now()))
        f.write('('+str(proc_id)+') '+
                time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime()))
        f.write(' ['+uniq_id+']')
        f.write(': ' +msg+'\n')
        f.flush()
        f.close()

def debug(msg):
    log(DEBUG, msg)

# used for temporary debugging
# use debug for more standard debugging
def debug2(msg):
    log(DEBUG2, msg)

def info(msg):
    log(INFO, msg)

def warning(msg):
    log(WARNING, msg)

def error(msg):
    log(ERROR, msg)

def critical(msg):
    log(CRITICAL, msg)


def log_per_sig(msg, sig):
    fl=logpersig + str(datetime.date.today())
    f=open(fl, 'a')
    #f.write('('+str(proc_id)+') '+str(datetime.datetime.now()))
    f.write('per_sig '+
            time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime()))
    f.write(':sigma = '+str(sig) +'::: '+msg+'\n')
    f.flush()
    f.close()


# prints stats of teh current python code 
def stats(msg):
    usage=resource.getrusage(resource.RUSAGE_SELF)
    statmsg= 'STATUS  [%s]: usertime=%s systime=%s mem=%s mb, rss=%s'%(msg,usage[0],usage[1],
                (usage[2]*resource.getpagesize())/1000000.0, mem.rss())
    info(statmsg)
    #debug('STATUS mpfr_array (num, del, diff), (mem, freed, diff):\
    #        (%d,%d,%d) (%d,%d,%d)'%(st.ndarray_init, st.ndarray_del,
    #            st.ndarray_init - st.ndarray_del, st.ndarray_mem,
    #            st.ndarray_free, st.ndarray_mem-st.ndarray_free))
    #debug('STATUS pref_float (num, del, diff), (mem, freed, diff):\
    #        (%d,%d,%d) (%d,%d,%d)'%(st.pf_init,st.pf_del, st.pf_init
    #            - st.pf_del, st.pf_mem, st.pf_free,
    #            st.pf_mem-st.pf_free))
 
