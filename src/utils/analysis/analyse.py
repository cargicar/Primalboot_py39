import matplotlib.pyplot as plt
import sqlite3
import numpy as np
from scipy import interpolate

import scan.epsfuncs as epsfuncs
import db.dbcache as dbcache
import config
from utils.ising import spec_exact as ising
reload(dbcache)

#def getspec(results,absmin,absmax):
def getspec(runid,absmin,absmax):
    results,rundata=dbcache.get_results(runid)
    mme=epsfuncs.min_max_eps(results,absmin,absmax)
    mmedict={}
    for i in range(len(mme[0])):
        mmedict[mme[0][i]]=(mme[1][i],mme[2][i])
    speclist={}
    for res in results:
        s=res.keys()[0][0]
        # skip points that maxed out
        if mmedict[s][0] == absmin or mmedict[s][1] == absmax:
            continue
        specstr = res[(s,mmedict[s][0])][2][2]
        solstr = res[(s,mmedict[s][0])][2][3]
        #print specstr
        # spec is returned as a unicode string we have
        # to evaluate to get a spec
        fullspec=eval(specstr)
        sol=[float(opecoeff) for opecoeff in solstr[1:-2].split()]
        #sol = eval(solstr)
        #print len(fullspec), len(sol)
        specdic={}
        for i in range(len(fullspec)):
            if len(fullspec[i]) == 3:
                specdic.setdefault(fullspec[i][1], [])
                specdic[fullspec[i][1]]+=[(fullspec[i][2], sol[i])]
        speclist[s]=specdic
    return speclist

#------------------ INDIVIDUAL RUNS --------------------------

def get_run_minmax_mathematica(runid):
    mmedict=get_run_minmax(runid)
    res='{'
    for k in sorted(mmedict.keys()):
        res+='{'+str(k) + ',' + str(mmedict[k][0]) + '},'
    res=res[:-1] + '}'
    return res


def get_run_minmax(runid):
    results,rundata=dbcache.get_results(runid)
    mme=epsfuncs.min_max_eps(results,0,5)
    mmedict={}
    for i in range(len(mme[0])):
        mmedict[mme[0][i]]=(mme[1][i],mme[2][i], mme[2][i]-mme[1][i])
    return mmedict

def get_run_ops(runid, l, threshold=0.05):
    results,rundata=dbcache.get_results(runid)
    spec=getspec(results, 0, 5)
    ops=condense_ops(get_ops(spec, l))
    return ops,spec

def find_kinks(runids, xlim=None, ylim=None, endpts=30, showplot=True, degree=20):
    """Find kink by minimized second deriv of leading scalar.  Because of
    polynomial interpolation endspoints are often the min so we drop them by
    picking a number of points (endpts) to drop on either end."""
    mins={}
    for runid in runids:
        results,rundata=dbcache.get_results(runid)
        spec=getspec(results, 0, 5)
        ops=condense_ops(get_ops(spec, 0))
        # get a polynomial and plot for the first operator
        siglist=[e[0] for e in ops[0]]
        deltalist=[e[1][0] for e in ops[0]]
        for k in [0,1,2]:
            poly=polyfit(siglist, deltalist, degree=degree, derivs=k)
            if showplot:
                plt.subplot(3, 1, k)
                plt.plot(siglist, poly(siglist))
                set_plot_lim(xlim,ylim)
            sigred=siglist[endpts:-endpts]
            #print 'Min[',k,']: sigma=',sigred[poly(sigred).argmin()]
        # the min from the last plot -- the 2nd deriv min
        mins[runid]=(rundata[1],sigred[poly(sigred).argmin()])
    if showplot:
        plt.show()
    return mins

# interpolate only on range s0-s1 but then find the min on the _smaller_
# interval smin0-smin1
def find_kinks2(runids, xlim=None, ylim=None, s0=-1,s1=-1, smin0=-1, smin1=-1,
        plotderivs=True, plotkinks=False, showplot=True, degree=20):
    """Find kink by minimized second deriv of leading scalar.  Because of
    polynomial interpolation endspoints are often the min so we drop them by
    looking for min only between smin0-smin1 (if != -1).  Also we allow the
    interpolation to run only over a range s0-s1 in case there are some bad
    points."""
    mins={}
    kinksig=[]
    kinkeps=[]
    for runid in runids:
        results,rundata=dbcache.get_results(runid)
        rs0=rundata[2]
        rs1=rundata[3]
        rsdel=rundata[4]
        numpts=int(float(rs1-rs0)/rsdel)
        mme=epsfuncs.min_max_eps(results,0,5)
        # get a polynomial and plot for the first operator
        siglist=[e for e in mme[0]]
        deltalist=[e for e in mme[1]]
        sigdict={siglist[i]: deltalist[i] for i in range(len(siglist))}
        if s0 != -1:
            indl=[i for i in range(len(siglist)) if siglist[i] >= s0 and
                    siglist[i] <= s1]
        else:
            indl=range(len(siglist))
        sigred=[siglist[i] for i in indl]
        delred=[deltalist[i] for i in indl]
        poly=[]
        for k in [0,1,2]:
            poly+=[polyfit(siglist, deltalist, degree=degree, derivs=k)]
            if plotderivs:
                plt.subplot(3, 1, k)
                plt.plot(sigred, poly[k](sigred), label="d=%d"%k)
                set_plot_lim(xlim,ylim)
            #print 'Min[',k,']: sigma=',sigred[poly(sigred).argmin()]
        # we look for min's only between smin0 and smin1
        if smin0 != -1:
            sigminlist=[s for s in sigred if s >= smin0 and s <= smin1]
        elif numpts > 20:
            # drop 10% on either side of the run
            #print 'bounds (%f,%f)'%(rs0+(0.1 * numpts *rsdel),rs1-(0.1 *
            #    numpts * rsdel))
            sigminlist=[s for s in sigred if s >= rs0+(0.1 * numpts *rsdel) and
                    s <= rs1-(0.1 * numpts * rsdel)]
        else:
            sigminlist=sigred
        sigkink=sigminlist[poly[2](sigminlist).argmin()]
        mins[runid]=(rundata[1], (sigkink, sigdict[sigkink]))
        kinksig+=[sigkink]
        kinkeps+=[sigdict[sigkink]]
    if plotkinks:
        plt.plot(kinksig, kinkeps, 'ro')
    if showplot:
        plt.show()
    return mins

# interpolate only on range s0-s1 but then find the min on the _smaller_
# interval smin0-smin1
def find_kinks3(runids, xlim=None, ylim=None, s0=-1,s1=-1, smin0=-1, smin1=-1,
        plotderivs=True, plotkinks=False, showplot=True):
    """Find kink by minimized second deriv of leading scalar.  Because of
    polynomial interpolation endspoints are often the min so we drop them by
    looking for min only between smin0-smin1 (if != -1).  Also we allow the
    interpolation to run only over a range s0-s1 in case there are some bad
    points."""
    mins={}
    kinksig=[]
    kinkeps=[]
    for runid in runids:
        results,rundata=dbcache.get_results(runid)
        rs0=rundata[2]
        rs1=rundata[3]
        rsdel=rundata[4]
        numpts=int(float(rs1-rs0)/rsdel)
        mme=epsfuncs.min_max_eps(results,0,5)
        # get a polynomial and plot for the first operator
        siglist=[e for e in mme[0]]
        deltalist=[e for e in mme[1]]
        sigdict={siglist[i]: deltalist[i] for i in range(len(siglist))}
        if s0 != -1:
            indl=[i for i in range(len(siglist)) if siglist[i] >= s0 and
                    siglist[i] <= s1]
        else:
            indl=range(len(siglist))
        # interpolation fails if two consecutive points have the same value on x
        # or y axis
        sigred=np.array([siglist[i] for i in indl if (i > 0 and deltalist[i] !=
            deltalist[i-1])])
        delred=np.array([deltalist[i] for i in indl
                if (i > 0 and deltalist[i] != deltalist[i-1])])
        tck=interpolate.splrep(sigred, delred, s=0)
        pdata = interpolate.splev(sigred, tck, der=0)
        x=interpolate.splev(np.array([0.12]), tck)
        for k in [0,1,2]:
            if plotderivs:
                plt.subplot(3, 1, k)
                epsdata=interpolate.splev(sigred, tck, der=k)
                #print epsdata
                plt.plot(sigred, epsdata, label="d=%d"%k)
                set_plot_lim(xlim,ylim)
            #print 'Min[',k,']: sigma=',sigred[poly(sigred).argmin()]
        # we look for min's only between smin0 and smin1
        if smin0 != -1:
            sigminlist=[s for s in sigred if s >= smin0 and s <= smin1]
        elif numpts > 20:
            # drop 10% on either side of the run
            #print 'bounds (%f,%f)'%(rs0+(0.1 * numpts *rsdel),rs1-(0.1 *
            #    numpts * rsdel))
            sigminlist=[s for s in sigred if s >= rs0+(0.1 * numpts *rsdel) and
                    s <= rs1-(0.1 * numpts * rsdel)]
        else:
            sigminlist=sigred
        D2epsminlist=interpolate.splev(sigminlist, tck, der=2)
        sigkink=sigminlist[D2epsminlist.argmin()]
        mins[runid]=(rundata[1], (sigkink, sigdict[sigkink]))
        kinksig+=[sigkink]
        kinkeps+=[sigdict[sigkink]]
    if plotkinks:
        plt.plot(kinksig, kinkeps, 'ro')
    if showplot:
        plt.show()
    return mins

#------------------ COMPARING DIFFERENT DIMS --------------------------

def find_sim_ops(runids,l,n,xlim=None, ylim=None):
    res={}
    rundata={}
    spec={}
    oplines={}
    print 'generating operators...',
    for i in runids:
        print i,',',
        res[i],rundata[i]=dbcache.get_results(i)
        spec[i]=getspec(res[i], 0, 5)
        oplines[i]=condense_ops(get_ops(spec[i], l))
    prevrun=oplines[runids[0]][n]  
    lines=[]
    print 'processing...',
    for i in range(len(runids)):
        print i,',',
        dist=1000
        useopline=None
        for opline in oplines[runids[i]]:
            newdist=getdist(opline, prevrun)
            if newdist < dist:
                useopline=opline
                dist=newdist
        if useopline is not None:
            print 'using line for runid[',runids[i],'] w avg: ',opline_avg(useopline)
            lines+=[useopline]
            prevrun=useopline
        else:
            print 'no operator line found for runid',runids[i]
    ptstyle=['r.','b.','g.','c.', 'm.', 'y.', 'k.','r^','b^','g^','c^', 'm^',
            'y^', 'k^']
    for j in range(len(lines)):
        #print line
        plot_opline(lines[j], show=False, scale=False,
                pointspec=ptstyle[j%len(ptstyle)], xlim=xlim, ylim=ylim)
    plt.show()
    return lines,res,spec,oplines

def multispec(runids,l, xlim=None, ylim=None, mins=None, deltamins=0, show=True,
        showexact=True):
    res={}
    rundata={}
    spec={}
    oplines={}
    print 'generating operators...',
    for i in runids:
        print i,',',
        res[i],rundata[i]=dbcache.get_results(i)
        spec[i]=getspec(res[i], 0, 5)
    ptstyle=['r.','b.','g.','c.', 'm.', 'y.', 'k.','r^','b^','g^','c^', 'm^',
            'y^', 'k^']
    for j in range(len(spec.keys())):
        #print line
        plotL(spec[spec.keys()[j]], l, show=False, scale=False,
                pointspec=ptstyle[j%len(ptstyle)], ylim=ylim, xlim=xlim)
        if mins is not None:
            mn=mins[spec.keys()[j]][1]
            plt.plot([mn-deltamins, mn-deltamins], [0, 100], c=ptstyle[j%len(ptstyle)][0])
            plt.plot([mn+deltamins, mn+deltamins], [0, 100], c=ptstyle[j%len(ptstyle)][0])
        if showexact and rundata[runids[j]][1]==2.0:
            Ispec=[s[0] for s in ising[l]]
            Isig=[0.125 for i in Ispec]
            plt.plot(Isig, Ispec,'kx',ms=10.0, mew=2)
    if show:
        plt.show()



#------------------ OPERATOR EXTRACTION --------------------------

def opline_avg(opline):
    return sum([o[1][0] for o in opline])/len(opline)

def opdist(pt1, pt2):
    s1,(delta1,ope1)=pt1
    s2,(delta2,ope2)=pt2
    return ((s1-s2)**2 + (delta1 -delta2)**2)**(0.5)

def getdist(op1, op2):
    return sum([min([opdist(o1, o2) for o1 in op1]) for o2 in op2])/len(op2)

def condense_ops(oplines, cutoff=30):
    newoplines=[]
    for opline in oplines:
        if len(opline) >= cutoff:
            newoplines += [opline]
    for opline in oplines:
        if len(opline) < cutoff:
            dist=100
            addtoline=None
            for goodline in newoplines:
                newdist=getdist(goodline, opline)
                if newdist < dist:
                    addtoline=goodline
            if addtoline is not None:
                addtoline += opline
    newoplines=sorted(newoplines, key=opline_avg)
    return newoplines

def in_opline(opline, s, op, threshold):
    """Checks if the point (s, delta) corresponds to the operator 'op'.  This
    happens when the former are sufficienlty close to some point in the
    latter."""
    for ent in opline:
        #print ent
        if opdist(ent, (s,op)) < threshold:
            return True
    return False


def check_add_op(s, op, oplines, threshold):
    """A utility function of get_ops: check if (s, delta) is part of an existing
    operator or a new operator and adds it approximately to oplist."""
    #print oplist
    for opline in oplines:
        if in_opline(opline, s, op, threshold):
            #print 'found same op'
            opline += [(s, op)]
            return
    opline=[(s, op)]
    oplines+=[opline]
    return


def get_ops(spec, l, threshold=0.05):
    """Return a list of operators "lines" for spin l.  Each entry in oplines is a list of
    (Delta_s, Delta_op, OPE_op) triples that form an approximately continuous
    line.  Threshhold controls the min distance between points for the function to be considered
    continuous."""
    siglist=sorted(spec.keys())
    oplines=[]
    for s in siglist:
        try:
            ops=sorted(spec[s][l])
        except:
            print 'Exception processing (s,l):',(s,l)
            #raise
        # op is (delta, opecoeff)
        for op in ops:
            check_add_op(s, op, oplines, threshold)
    return oplines

#--------------------------------------------------------------------



def plotbounds(runlist, xlim=None, ylim=None, legend=True,showplot=True,
        anom_dim=False):
    pl=[]
    ll=[]
    for runid in runlist:
        results,rundata=dbcache.get_results(runid)
        mme=epsfuncs.min_max_eps(results,0,5)
        dimfree=0
        # plot 
        if anom_dim:
            dim=rundata[1]
            TpDim=float(rundata[7].split('Tp=')[1])
            dimfree=(TpDim-2)/2
            sig = np.array(mme[0]) - dimfree
            eps0 = np.array(mme[1]) - 2*dimfree
            eps1 = np.array(mme[2]) - 2*dimfree
        else:
            sig,eps0,eps1=mme
        p1,p2=plt.plot(sig, eps0, '--', sig,eps1, '-') 
        #p1,p2=plt.plot(mme[0], mme[1], 'ro', mme[0],mme[2], 'go') 
        pl+=[p2]
        ll+=[rundata[7]]
    if legend:
        plt.legend(pl, ll)
    set_plot_lim(xlim,ylim)
    if showplot:
        plt.show() 
        
def plotboundsdelta(runlist, xlim=None, ylim=None, legend=True,showplot=True):
    pl=[]
    ll=[]
    for runid in runlist:
        results,rundata=dbcache.get_results(runid)
        mme=epsfuncs.min_max_eps(results,0,5)
        dlt=[0]
        for i in range(1,len(mme[0])):
            dlt += [mme[1][i]-mme[1][i-1]]
        dlt2=[0]
        for i in range(1,len(dlt)):
            dlt2 += [dlt[i]-dlt[i-1]]
        p1=plt.plot(mme[0], dlt2, 'ro', mme[0], dlt, 'bo')
        #p1,p2=plt.plot(mme[0], mme[1], '--', mme[0],mme[2], '-') 
        #pl+=[p2]
        #ll+=[rundata[7]]
    #if legend:
    #    plt.legend(pl, ll)
    set_plot_lim(xlim,ylim)
    if showplot:
        plt.show() 
        


def polyfit(siglist, deltalist, degree, derivs=0):
        pc=np.polyfit(siglist, deltalist, degree)
        poly=np.poly1d(pc).deriv(derivs)
        return poly


def set_plot_lim(xlim, ylim):
    v=plt.axis()
    lm=[]
    lm+=v
    if xlim is not None:
        lm[0]=xlim[0]
        lm[1]=xlim[1]
    if ylim is not None:
        lm[2]=ylim[0]
        lm[3]=ylim[1]
    plt.axis(lm)

def plot_opline(opline, xlim=None, ylim=None, show=True, scale=True,
        pointspec='ro', polydeg=0):
    siglist=[e[0] for e in opline]
    deltalist=[e[1][0] for e in opline]
    opelist=[np.log(np.abs(e[1][1])) for e in opline]
    mean=sum(opelist)/len(opelist)
    logmean=np.exp(sum([np.log(np.abs(o)) for o in opelist])/len(opelist))
    omax=max(opelist)
    omin=min(opelist)
    factor=100
    opescaled=[factor*(o-omin)/(omax-omin) for o in opelist]
    if scale:
        plt.scatter(siglist,deltalist, c=opescaled, s=opescaled)
    else:
        plt.plot(siglist,deltalist, pointspec)
    poly=None
    if polydeg > 0:
        poly=polyfit(siglist, deltalist, degree=polydeg)
        plt.plot(siglist, poly(siglist), 'r--')
    set_plot_lim(xlim,ylim)
    if show:
        plt.show()
    return plt,siglist,deltalist

def plotL(spec, l, xlim=None, ylim=None, show=True, scale=True, pointspec='ro'):
    '''spec should be returned by getspec'''
    sigs=spec.keys()
    siglist=[]
    deltalist=[]
    opelist=[]
    for s in sigs:
        llist=spec[s][l]
        for dl in llist:
            siglist+=[s]
            deltalist+=[dl[0]]
            opelist+=[dl[1]]
    #plt.plot(xlist,ylist,'ro')
    if scale:
        #plt.scatter(siglist,deltalist, c=opescaled, s=opescaled)
        plt.scatter(siglist,deltalist, c=opelist, s=opelist)
    else:
        plt.plot(siglist,deltalist, pointspec)
    v=plt.axis()
    lm=[]
    lm+=v
    if xlim is not None:
        lm[0]=xlim[0]
        lm[1]=xlim[1]
    if ylim is not None:
        lm[2]=ylim[0]
        lm[3]=ylim[1]
    plt.axis(lm)
    if show:
        plt.show()
    return plt

def plotLdiff(spec, l, xlim=None, ylim=None):
    '''plot difference between operator dimension.
    spec should be returned by getspec'''
    sigs=spec.keys()
    xlist=[]
    ylist=[]
    for s in sigs:
        llist=spec[s][l]
        llist.sort()
        for i in range(1,len(llist)):
            ylist+=[llist[i]-llist[i-1]]
            xlist+=[s]
    plt.plot(xlist,ylist,'ro')
    v=plt.axis()
    lm=[]
    lm+=v
    print lm
    if xlim is not None:
        lm[0]=xlim[0]
        lm[1]=xlim[1]
    if ylim is not None:
        lm[2]=ylim[0]
        lm[3]=ylim[1]
    plt.axis(lm)
    plt.show()
    return plt

