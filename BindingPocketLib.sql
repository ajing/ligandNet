DROP DATABASE IF EXISTS ligandNet;

CREATE DATABASE ligandNet;

USE ligandNet;

CREATE TABLE BindingPocketInfo
(
    PDBID           varchar(255),
    ligandName      varchar(255),
    ligandChainID   varchar(1),
    proteinChainID  varchar(1),
    residueName     varchar(3),
    residueNumber   varchar(3),
    AtomName        varchar(3),
    AtomNumber      varchar(1),
    Distance        double(9,8),
    Metal           int             /* Whether it is metal or not */
);

/* BULK INSERT BindingPocketInfo FROM "/users/ajing/ligandNet/pylib/final.txt" WITH (FILEDTERMINATOR = ","); */

LOAD DATA LOCAL INFILE "/users/ajing/ligandNet/pylib/final.txt"
    INTO TABLE BindingPocketInfo
    FIELDS
        TERMINATED BY ","
    LINES
        TERMINATED BY "\n"
