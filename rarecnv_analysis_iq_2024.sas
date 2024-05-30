

options helpbrowser=sas; 
options nocenter nodate nonumber;
footnote;title;
goptions reset=global;

proc import out= data
datafile= "/project/bbl_gur_22q/rarecnv/scripts/revision_apr2024/cnbiq_accuracy_speed_2024.csv"
dbms=csv replace;
 getnames=yes;
guessingrows=1000;
datarow=2; run;




title 'GLM, IQ: Group: Accuracy' ;
proc glm data=data;
 class rarecnv_id Group Del_Dup Locus;
 model imputed_ciq_a = Locus|Del_Dup;
 lsmeans Locus*Del_Dup/stderr pdiff;
 run;


title 'GLM, IQ: Group: Speed' ;
proc glm data=data;
 class rarecnv_id Group Del_Dup Locus;
 model imputed_ciq_s = Locus|Del_Dup;
 lsmeans Locus*Del_Dup/stderr pdiff;
 run;



endsas;
