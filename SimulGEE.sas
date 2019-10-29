
libname outd "/home/xinqi/Documents";
    
proc import datafile="/home/xinqi/Downloads/cbdataGEE.csv"
     out=cbdataGEE
     dbms=csv
     replace;
     getnames=yes;
run;

proc genmod data=cbdataGEE descending;
   class dose (ref="0") case;
   model resp = dose / dist=bin link=log;
   repeated  subject=case / type=exch covb corrw;
   by simul;
   ods output GEEEmpPEst=A;
run;

data outd.GEEest;
    set A;
run;
endsas;

/*** format report ***/
options missing = "";
options orientation = portrait;
ods rtf file = "/home/xinqi/Documents/GEE_0719.rtf";

ods rtf off;
