OPTIONS COMPRESS=YES REUSE=YES;
proc import out= data
datafile= "/project/bbl_gur_22q/rarecnv/scripts/revision_apr2024/rarecnv_cnb_zscores_winsorized8_qced_imputed_residualized.csv"
dbms=csv replace;
 getnames=yes;
guessingrows=10000;
datarow=2; run;




data data_a;
set data;
if measure="s_resid" then delete;

title '1. Mixed Model Imputed Residualized Data, Accuracy' ;
proc mixed data=data_a method=REML;
 class rarecnv_id Del_Dup Locus test sex sitenum;
 model mean_zscore = Locus|Del_Dup|test sitenum sex test*sitenum / DDFM=KenwardRoger;
 repeated test / subject=rarecnv_id Type=UN;
 lsmeans test*Locus*Del_Dup / PDiff; 
 run;  


data data_s;
set data;
if measure="a_resid" then delete;

title '2. Mixed Model Imputer Residualized Data, Speed' ;
proc mixed data=data_s method=REML;
 class rarecnv_id Del_Dup Locus test sex sitenum;
 model mean_zscore = Locus|Del_Dup|test sex sitenum test*sitenum/ DDFM=KenwardRoger;
 repeated test / subject=rarecnv_id Type=UN;
 lsmeans test*Locus*Del_Dup / PDiff; 
 run;  


title '3. Mixed Model, THE BIG MODEL' ;
proc mixed data=data method=REML;
 class rarecnv_id Del_Dup Locus test measure sex sitenum;
 model mean_zscore = Locus|Del_Dup|test|measure sex sitenum/ DDFM=KenwardRoger;
 repeated test measure / subject=rarecnv_id Type=UN@CS;
 lsmeans test*measure*Locus*Del_Dup / PDiff; 
 run;

quit;
