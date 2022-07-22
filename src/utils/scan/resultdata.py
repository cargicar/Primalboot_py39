def save_results(scan, sigdata,epsdata,absmin,absmax):
    if config.dbfile == None:
        return
    conn = sqlite3.connect(config.dbfile)
    c = conn.cursor()
    #c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?, -1, 0.0)",
    #             (sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
    c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?,?,?,?, -1, 0.0)",
                (scan.rundata.desc, scan.rundata.prec, scan.rundata.spacetimedim, 
                    sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
    runid=c.lastrowid
    # mme is in different sort order!!
    mme=min_max_eps(scan.results,absmin,absmax)
    mmedict={}
    for i in range(len(mme[0])):
        mmedict[mme[0][i]]=(mme[1][i],mme[2][i])
    resclean=[r for r in scan.results if r is not None and len(r.keys()) > 0]
    #for res in self.results:
    for res in resclean:
        if res is None:
            log.warning("Empty entry in results list")
            continue
        s=res.keys()[0][0]
        epsfinmin,epsfinmax=mmedict[s]
        edata = [e(s) for e in epsdata]
        c.execute("INSERT INTO bisections VALUES (datetime('now'),\
                    ?,?,?,?,?,?,?,?)",(runid, s, edata[0],edata[1],edata[2],
                    epsfinmin,epsfinmax,0.0))
        bisectid=c.lastrowid
        for k in res.keys():
            status,tm,lpdata,cblen=res[k]
            #print lpdata
            c.execute("INSERT INTO points VALUES (datetime('now'),\
                        ?,?,?,?,?,?,?,?,?,?,?,?)",(bisectid, s, k[1], status,
                            lpdata[0], lpdata[1], str(lpdata[2]), str(lpdata[3]),
                            lpdata[4], lpdata[5], cblen, str(tm)))
    # Save (commit) the changes
    conn.commit()
    # We can also close the connection if we are done with it.
    #Just be sure any changes have been committed or they will be lost.
    conn.close()



