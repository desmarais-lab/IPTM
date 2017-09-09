DROP VIEW IF EXISTS Recipient;
DROP TABLE IF EXISTS RecipientBase;
DROP VIEW IF EXISTS Message;
DROP INDEX IF EXISTS IX_MessageBase_unix_time;
DROP TABLE IF EXISTS MessageBase;
DROP VIEW IF EXISTS EmployeeWithVars;
DROP VIEW IF EXISTS Employee;
DROP TABLE IF EXISTS EmployeeBase;


CREATE TABLE EmployeeBase (
    eid INTEGER,
    name TEXT,
    department TEXT,
    longdepartment TEXT,
    title TEXT,
    gender TEXT,
    seniority TEXT,
    
    PRIMARY KEY(eid ASC)
);

CREATE VIEW Employee AS
SELECT
    eid,
    name,
    longdepartment,
    title,

    gender,
    CASE gender
        WHEN 'Female' THEN 1
        ELSE 0
    END AS genF,

    seniority,
    CASE seniority
        WHEN 'Junior' THEN 1
        ELSE 0
    END AS senJ,

    department,
    CASE department
        WHEN 'Legal' THEN 1
        ELSE 0
    END AS depL,
    CASE department
        WHEN 'Trading' THEN 1
        ELSE 0
    END AS depT
FROM
    EmployeeBase;

CREATE VIEW EmployeeWithVars AS
SELECT
    eid,
    1 AS intercept,
    genF,
    senJ,
    depL,
    depT,
    genF * senJ AS genF_senJ,
    genF * depL AS genF_depL,
    genF * depT AS genF_depT,
    senJ * depL AS senJ_depL,
    senJ * depT AS senJ_depT,
    genF * senJ * depL AS genF_senJ_depL,
    genF * senJ * depT AS genF_senJ_depT
FROM
    Employee;

CREATE TABLE MessageBase (
    mid INTEGER,
    filename TEXT,
    unix_time INTEGER,
    subject TEXT,
    from_eid INTEGER,
    
    PRIMARY KEY(mid ASC),
    FOREIGN KEY(from_eid) REFERENCES Employee(eid)
);

CREATE INDEX IX_MessageBase_unix_time ON MessageBase(unix_time ASC);

CREATE VIEW Message AS
SELECT
    mid,
    filename,
    datetime(unix_time, 'unixepoch') AS time,
    unix_time,
    subject,
    from_eid
FROM
    MessageBase;

CREATE TABLE RecipientBase (
    mid INTEGER,
    rno INTEGER,
    to_eid INTEGER,
    
    PRIMARY KEY(mid ASC, rno ASC)
    FOREIGN KEY(mid) REFERENCES Message(mid)
    FOREIGN KEY(to_eid) REFERENCES Employee(eid)
);

CREATE VIEW Recipient AS
SELECT
    mid,
    rno,
    to_eid
FROM
    RecipientBase;
