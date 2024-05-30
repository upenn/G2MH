

options helpbrowser=sas; 
options nocenter nodate nonumber;
footnote;title;
goptions reset=global;

proc import out= data
datafile= "/project/bbl_gur_22q/rarecnv/scripts/revision_apr2024/siteeffects_data.csv"
dbms=csv replace;
 getnames=yes;
guessingrows=10000;
datarow=2; run;



data data_a;
set data;
if measure="s_resid" then delete;



data data_s;
set data;
if measure="a_resid" then delete;

title 'Mixed Model, THE BIG MODEL- Site Effects' ;
proc mixed data=data method=REML;
 class rarecnv_id test measure sitenum;
 model mean_zscore = sitenum|test|measure / DDFM=KenwardRoger;
 repeated test measure / subject=rarecnv_id Type=un@un;
 lsmeans sitenum*test / PDiff; 
 run;

endsas;

title 'Mixed Model, THE BIG MODEL- Site Effects Speed' ;
proc mixed data=data_s method=REML;
 class rarecnv_id Del_Dup Locus test measure sitenum;
 model mean_zscore = sitenum|Locus|Del_Dup|test / DDFM=KenwardRoger;
 repeated test measure / subject=rarecnv_id Type=UN@CS;
 lsmeans sitenum*test*Locus*Del_Dup / PDiff;
 run;


quit;
