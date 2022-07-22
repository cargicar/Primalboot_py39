DROP TABLE points;
DROP TABLE bisections;
DROP TABLE runs;
DROP TABLE operuns;
DROP TABLE opepoints;

CREATE TABLE runs (
	createtime INTEGER,
    description TEXT,
    prec TEXT,
    spacetimedim REAL,
	sigma0 REAL,
	sigma1 REAL,
	sig_delta REAL,
	bisections INTEGER,
	cblen INTEGER,
	time REAL
);

CREATE TABLE bisections (
	createtime INTEGER,
	runid INTEGER,
	sigma REAL,
	epsmin_init REAL,
	epsmax_init REAL,
	epsres REAL,
	epsmin_fin REAL,
	epsmax_fin REAL,
	time REAL,
	FOREIGN KEY(runid) REFERENCES runs(rowid)
);


CREATE TABLE points (
	createtime INTEGER,
	bisectid INTEGER,
	sigma REAL,
	epsilon REAL,
	status TEXT,
	cost REAL,
	to_elim INTEGER,
	Xb TEXT,
	solution TEXT,
	iterations INTEGER,
	hotstarted TEXT,
    cblen INTEGER,
	time REAL,
	FOREIGN KEY(bisectid) REFERENCES bisections(rowid)
);

CREATE TABLE operuns (
	createtime INTEGER,
    description TEXT,
    prec TEXT,
    spacetimedim REAL,
	sigma0 REAL,
	sigma1 REAL,
	sig_delta REAL,
	cblen INTEGER,
	time REAL
);

CREATE TABLE opepoints (
	createtime INTEGER,
	runid INTEGER,
	sigma REAL,
	delta REAL,
	spin INTEGER,
	epsmin REAL,
	bound REAL,
	status TEXT,
	cost REAL,
	to_elim INTEGER,
	Xb TEXT,
	solution TEXT,
	iterations INTEGER,
	hotstarted TEXT,
    cblen INTEGER,
	time REAL,
	FOREIGN KEY(runid) REFERENCES operuns(rowid)
);


CREATE INDEX runid ON bisections (runid);
CREATE INDEX bisectid ON points (bisectid);
CREATE INDEX operunid on opepoints (runid);
