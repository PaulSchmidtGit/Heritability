/* Data taken from:                                          */
/* John, J. A., and E. R. Williams, 1995 Cyclic and Computer */
/* Generated Designs. Chapman & Hall, London, p.146          */
/* yield trial with oats laid out as an a-design.            */
/*                                                           */
/* The trial had 24 genotypes, three complete replications,  */
/* and six incomplete blocks within each replication. The    */
/* block size was four. The data were analyzed by a linear   */
/* mixed model with effects for genotypes, replicates, and   */
/* incomplete blocks. Blocks were modeled as independent     */
/* random effects to recover interblock information          */

data a;
input
rep   block     gen     y;
datalines;
 1       1       11    4.1172
 1       1        4    4.4461
 1       1        5    5.8757
 1       1       22    4.5784
 1       2       21    4.6540
 1       2       10    4.1736
 1       2       20    4.0141
 1       2        2    4.3350
 1       3       23    4.2323
 1       3       14    4.7572
 1       3       16    4.4906
 1       3       18    3.9737
 1       4       13    4.2530
 1       4        3    3.3420
 1       4       19    4.7269
 1       4        8    4.9989
 1       5       17    4.7876
 1       5       15    5.0902
 1       5        7    4.1505
 1       5        1    5.1202
 1       6        6    4.7085
 1       6       12    5.2560
 1       6       24    4.9577
 1       6        9    3.3986
 2       1        8    3.9926
 2       1       20    3.6056
 2       1       14    4.5294
 2       1        4    4.3599
 2       2       24    3.9039
 2       2       15    4.9114
 2       2        3    3.7999
 2       2       23    4.3042
 2       3       12    5.3127
 2       3       11    5.1163
 2       3       21    5.3802
 2       3       17    5.0744
 2       4        5    5.1202
 2       4        9    4.2955
 2       4       10    4.9057
 2       4        1    5.7161
 2       5        2    5.1566
 2       5       18    5.0988
 2       5       13    5.4840
 2       5       22    5.0969
 2       6       19    5.3148
 2       6        7    4.6297
 2       6        6    5.1751
 2       6       16    5.3024
 3       1       11    3.9205
 3       1        1    4.6512
 3       1       14    4.3887
 3       1       19    4.5552
 3       2        2    4.0510
 3       2       15    4.6783
 3       2        9    3.1407
 3       2        8    3.9821
 3       3       17    4.3234
 3       3       18    4.2486
 3       3        4    4.3960
 3       3        6    4.2474
 3       4       12    4.1746
 3       4       13    4.7512
 3       4       10    4.0875
 3       4       23    3.8721
 3       5       21    4.4130
 3       5       22    4.2397
 3       5       16    4.3852
 3       5       24    3.5655
 3       6        3    2.8873
 3       6        5    4.1972
 3       6       20    3.7349
 3       6        7    3.6096
; RUN;

/**************************************/
/* include macro directly from github */
/**************************************/

/* Macro %H2_piepho */
filename _inbox "%sysfunc(getoption(work))/Macro H2_Piepho.sas";
proc http method="get" 
url="https://raw.githubusercontent.com/PaulSchmidtGit/Heritability/master/Alternative%20Heritability%20Measures/SAS/MACRO%20H2_Piepho.sas" out=_inbox;
run; %Include _inbox; filename _inbox clear;

ODS HTML CLOSE; *Turn html results viewer off;

/**************/
/* fit models */
/**************/

/* Genotype as random effect */
proc mixed data=a;
class rep block gen;
model y = rep /ddfm=kr;
random gen rep*block;
ods output Covparms=Covparms;        * obtain estimated variance components;
run;

data Covparms_no_gen; 
set  Covparms;
if Covparm = "gen" then delete;
run;

/* Genotype as fixed effect */
proc mixed data=a;
class rep block gen;
model y = gen rep /ddfm=kr;
random rep*block;
lsmeans gen /pdiff;
parms /pdata=Covparms_no_gen noiter; * fixed variance components estimated in model above for better comparability;
ods output Diffs=Diffs;              * obtain pairwise differences table;
run;

/*****************/
/* H2 estimation */
/*****************/
%H2_piepho(ENTRY_NAME=gen, COVPARMS=Covparms, DIFFS=Diffs, OUTPUT=H2_piepho);

ods html; *Turn html results viewer on;

/* Show results */
title "Ad hoc H2 'Piepho'"; 
proc print data=H2_piepho label; 
run;





