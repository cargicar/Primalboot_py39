import setpath

import numpy as np
import os
import shutil
#import cPickle as pickle
import pickle  # Carlos edit
from subprocess import call,Popen
import time
import datetime

import probspec as ps
import config
import table_funcs as tf

from importlib import reload # Carlos add  
reload(ps)
reload(tf)

if not ps.clusterMode:
    import ope_bound_analyse as analyse
    reload(analyse)

siglist=np.arange(ps.s0, ps.s1, ps.sdelta)

# a dictionary keyed by names of saved batches.  used to avoid resaving batches
# to DB
saved_batches={}

# used to store the results of a single batch (populated in getresults()).
fullresults={}



def newbatch():
    batchdir=ps.workdir+'/'+time.strftime("%Y%m%d-%H%M/", time.localtime())
    tf.checktable()
    print ('new batch with siglist: '+str(siglist) + ' in dir ' + batchdir)
    for i in range(len(siglist)):
        jobdir=batchdir+'/'+str(i) + '/'
        if not os.path.isdir(jobdir):
            os.makedirs(jobdir)
        pickle.dump(siglist[i], open(jobdir + ps.sigfile,'wb'))

def batchlist():
    batches=[]
    cnt=0
    blist= sorted(os.listdir(ps.workdir))[-10:]
    for batch in blist:
        batches+=[batch]
        print ('[%d]: %s'%(cnt,batch))
        cnt+=1
   # var = raw_input("Select batch number (or ENTER for most recent): ")
    var = input("Select batch number (or ENTER for most recent): ") # Carlos edit
    print ("you entered %s" %(var))
    if len(var)==0:
        return batches[-1]
    try:
        jn=int(var)
        return batches[jn]
    except:
        print ('You must enter an integer.')
        raise

def confirm(question):
    var = raw_input(question + " [y/n]")
    if len(var) > 0 and var == 'y':
        return True
    return False


def getresults(return_dir=False):
    global fullresults
    batchdir=ps.workdir+'/'+batchlist()
    print ('\nBatch %s:' %(batchdir) )
    all_done=True
    for jn in sorted(os.listdir(batchdir)):
        jobdir=batchdir+'/'+jn
        resultfile=jobdir+ps.resfile
        sigmafile=jobdir+ps.sigfile
        if os.path.isfile(resultfile):
            sig = pickle.load(open(sigmafile, 'rb'))
            res = pickle.load(open(resultfile, 'rb'))
            fullresults[sig]=res
            print ('job[%s] results stored' %(jn))
        else:
            print ('job[%s] not yet done' %(jn))
            all_done = False
    if return_dir:
        return all_done,batchdir


def savebatch(ignore_saved=False):
    global fullresults, saved_batches
    all_done,batchdir=getresults(return_dir=True)
    if (not ignore_saved) and batchdir in saved_batches.keys():
        print ('Looks like batch already saved.  Call with \
                ignore_saved=True to ignore this warning')
        return
    if all_done:
        print ('Looks like all results are in. Saving...')
        analyse.process(fullresults)
        analyse.save()
        saved_batches[batchdir]=True
        print ('done.')
    else:
        print ('Still some unfinished jobs, not saving.')

# delete a batch folder
def delbatch():
    batchdir=batchlist()
    if confirm("Delete dir [%s]?"%batchdir):
        print ('Deleting %s'%(batchdir))
        shutil.rmtree(ps.workdir + '/' + batchdir)


def status():
    batchdir=ps.workdir+'/'+batchlist()
    for jn in sorted(os.listdir(batchdir)):
        jobdir=batchdir+'/'+jn
        if os.path.isfile(jobdir+ps.resfile):
            print ('Job %s finished.' %(jn))
        elif os.path.isfile(jobdir+ps.runfile):
            print ('Job %s running.' %(jn))
        else:
            print ('Job %s not started.' %(jn))

def runbatch():
    batchdir=ps.workdir+'/'+batchlist()
    for jn in sorted(os.listdir(batchdir)):
        jobdir=batchdir+'/'+jn
        if os.path.isfile(jobdir+ps.resfile):
            print ('Job %s finished.' %(jn))
        elif os.path.isfile(jobdir+ps.runfile):
            print ('Job %s running.' %(jn))
        else:
            print ('Job %s not started.  Starting...' %(jn))
            Popen([ps.pybin, ps.reaper, jobdir],
                    stdout=open(jobdir + '/'+ps.output, 'w'),
                    stderr=open(jobdir + '/'+ps.output, 'w')) 
        #$pybin $reaper $jobdir >& $jobdir/output &
        

def clusterbatch():
    batchdir=ps.workdir+'/'+batchlist()
    for jn in sorted(os.listdir(batchdir)):
        jobdir=batchdir+'/'+jn
        if os.path.isfile(jobdir+ps.resfile):
            print ('Job %s finished.' %(jn))
        elif os.path.isfile(jobdir+ps.runfile):
            print ('Job %s running.' %(jn))
        else:
            jobcmd='"cd %s; %s %s %s >& %s/%s"'%(ps.clusterCWD, ps.clusterpybin,
                    ps.reaper, jobdir, jobdir, ps.output)
            #jobcmd='"ls %s"' %ps.clusterCWD
            print ('Job %s not started.  Starting... (%s)'%(jn,jobcmd))
            #Popen(['bsub', '-q', ps.queue, jobcmd ],
            #        stdout=open(jobdir + '/'+ps.bsuboutput, 'w'),
            #        stderr=open(jobdir + '/'+ps.bsuboutput, 'w')) 
            shellcmd='bsub -q %s -R "type=SLC6" %s'%(ps.queue, jobcmd)
            shellcmd='bsub -q %s %s'%(ps.queue, jobcmd)
            print ('shellcmd: '+shellcmd)
            Popen(shellcmd,
                    stdout=open(jobdir + '/'+ps.bsuboutput, 'w'),
                    stderr=open(jobdir + '/'+ps.bsuboutput, 'w'),
                    shell=True) 
        #$pybin $reaper $jobdir >& $jobdir/output &
        

