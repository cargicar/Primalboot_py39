import sqlite3

import config


result_cache={}


def del_results(runid):
    if config.dbfile == None:
        print "No db file specified in config!"
        return
    conn = sqlite3.connect(config.dbfile)
    c = conn.cursor()
    bsel='''SELECT rowid FROM BISECTIONS WHERE runid=?'''
    c.execute(bsel, (runid,))
    row = c.fetchone()
    bids=[]
    while row != None:
        #print row[0], row[1] 
        bids+=[row[0]]
        row = c.fetchone()
    results=[]
    print 'bisection ids:',bids
    pdel ='''DELETE FROM points WHERE bisectid=?'''
    for bid in bids:
        points={}
        c.execute(pdel,(bid,))
    bdel='''DELETE FROM bisections WHERE runid=?'''
    c.execute(bdel,(runid,))
    rdel='''DELETE FROM runs WHERE rowid=?'''
    c.execute(rdel,(runid,))
    conn.commit()
    conn.close()


def clear_cache():
    result_cache.clear()

def get_results(runid):
    if config.dbfile == None:
        print "No db file specified in config!"
        return
    if runid in result_cache.keys():
        return result_cache[runid]
    conn = sqlite3.connect(config.dbfile)
    c = conn.cursor()
    #rsel='''SELECT prec, spacetimedim, sigma0, sigma1, sig_delta, bisections,
    #        cblen FROM runs WHERE rowid=?'''
    rsel='''SELECT prec, spacetimedim, sigma0, sigma1, sig_delta, bisections,
            cblen, description FROM runs WHERE rowid=?'''
    bsel='''SELECT rowid,sigma, epsmin_init,epsmax_init,
            epsres,epsmin_fin,epsmax_fin,time FROM bisections WHERE runid=?'''
    psel='''SELECT rowid,sigma,epsilon,status,cost,to_elim,Xb,solution,
            iterations,hotstarted,time FROM points WHERE bisectid=?'''
    c.execute(rsel, (runid,))
    row = c.fetchone()
    rundata=row
    c.execute(bsel, (runid,))
    row = c.fetchone()
    bids=[]
    while row != None:
        #print row[0], row[1] 
        bids+=[row[0]]
        row = c.fetchone()
    results=[]
    for bid in bids:
        points={}
        c.execute(psel,(bid,))
        row=c.fetchone()
        while row != None:
            #print '\t',row[0],row[1],row[2],row[3]
            status=row[3]
            tm=row[10]
            lpdata=[row[4],row[5],row[6],row[7],row[8],row[9]]
            points[(row[1],row[2])]=status,tm,lpdata
            row=c.fetchone()
        results+=[points]
    conn.close()
    result_cache[runid]=results,rundata
    return results,rundata


