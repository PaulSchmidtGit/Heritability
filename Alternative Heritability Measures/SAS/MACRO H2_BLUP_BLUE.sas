/************************************************************************************************/
/*                     _  _ ___   ___ _   _   _ ___    ___ _   _   _ ___                        */
/*                    | || |_  ) | _ ) | | | | | __|__| _ ) | | | | | _ \                       */
/*                    | __ |/ /  | _ \ |_| |_| | _|___| _ \ |_| |_| |  _/                       */
/*                    |_||_/___| |___/____\___/|___|  |___/____\___/|_|                         */
/*                                                                                              */
/* This macro computes two measures of heritability (H² Reg & H² SumDiv) based on genotypic     */
/* BLUEs and BLUPs.			                                                                    */
/*																								*/
/* Example application code can be found on https://github.com/PaulSchmidtGit/Heritability      */
/*																								*/
/* The methods are based on unpublished work of the author of this macro.						*/
/*																								*/
/*	Requirements/Input:																			*/
/*		The model that is used to analyze the data beforehand should have a genotype main 		*/
/*		effect. This model should then be fitted in two versions: once with the genotype main  	*/
/*		effect random to obtain genotype BLUPs and once with the genotype main effect           */
/*		fixed to obtain LS-means for the genotypes (i.e. genotype BLUEs). Furthermore, the      */
/*      model with the random genotype main effect should have a pseudo-intercept effect in     */
/*      the model statement and the /noint option. See example code for this macro.             */
/*																								*/
/*		SAS/STAT																				*/
/*			Name of genetic effect																*/
/*				ENTRY_NAME=	specifies the genotypic main effect (e.g. var, entry, gen, g).	    */
/*			Dataset 'LSM_MU'																	*/
/*				LSM_MU= specifies the MIXED / GLIMMIX ODS output with the least squares mean	*/
/*				for the intercept "Mu" from a model where the genotype main effect is random.	*/
/*			Dataset 'SOLUTIONR'																	*/
/*				SOLUTIONR= specifies the MIXED / GLIMMIX ODS output with the random-effects 	*/
/*				solution vector	from the same model as SOLUTIONF.								*/
/*			Dataset 'LSM_G'																	    */
/*				LSM_G= specifies the MIXED / GLIMMIX ODS output with the least squares means    */
/*				of the genotypes. It must be produced from a model where the genotypic effect 	*/
/*				is fixed.																		*/
/*			Name for output file																*/
/*				OUTPUT= specifies the name for the output dataset.								*/
/*																								*/
/*	Note that in order to prevent complications due to overwritten data, one should not use 	*/
/*	dataset or variable names starting with "xm_" as some are used in this macro.				*/
/*																								*/
/*	Version 02 October 2018																		*/
/*																								*/
/*	Written by: Paul Schmidt (Paul.Schmidt@uni-hohenheim.de)									*/
/*																								*/
/************************************************************************************************/

%MACRO H2_BLUE_BLUP(ENTRY_NAME=, LSM_MU=, SOLUTIONR=, LSM_G=, OUTPUT=);

/* Save overall mean from in macro variable "xm_Mu_ran" */
DATA &LSM_Mu.; SET &LSM_Mu.;
CALL SYMPUT("xm_Mu_ran", Estimate);
RUN;

/* Format genotype BLUE Table */
DATA xm_blues;
LENGTH Effect $ 32; 
SET &LSM_G.;
WHERE Effect="&ENTRY_NAME.";
KEEP Effect &ENTRY_NAME. Estimate;
RENAME Estimate=BLUE;
RUN;

/* Format genotype BLUP Table */
DATA xm_blups;
LENGTH Effect $ 32; 
SET &SolutionR.;
WHERE Effect="&ENTRY_NAME.";
KEEP Effect &ENTRY_NAME. Estimate;
RENAME Estimate=BLUP;
RUN;

/* Merge genotype BLUEs and BLUPs into BLUE-BLUP table */
PROC SORT DATA=xm_blues; BY &ENTRY_NAME.; RUN;
PROC SORT DATA=xm_blups; BY &ENTRY_NAME.; RUN;
DATA xm_full;
MERGE xm_blues xm_blups;
BY &ENTRY_NAME.; 
scaledGBLUE    = BLUE-&xm_Mu_ran.;
absBLUP        = abs(BLUP);
absscaledGBLUE = abs(scaledGBLUE);
RUN;

/* H2 Reg */
PROC MIXED DATA=xm_full; 
MODEL BLUP=scaledGBLUE /S NOINT;
ODS OUTPUT SolutionF=xm_reg;
RUN;

DATA xm_reg; SET xm_reg;
KEEP h2_reg;
H2_reg=Estimate;
RUN;

/* H2 SumDiv */
PROC MEANS DATA=xm_full SUM NOPRINT; VAR absBLUP absscaledGBLUE;
OUTPUT OUT=xm_sum;
RUN;

DATA xm_sum; 
SET xm_sum; 
WHERE _stat_='MEAN';
H2_SumDiv=(absblup/absscaledGBLUE);
KEEP H2_SumDiv;
RUN; 

/* Final formatting */
DATA &OUTPUT.;
RETAIN H2_reg H2_SumDiv;
MERGE xm_reg xm_sum;
FORMAT _numeric_ 10.3;
LABEL 	H2_reg="H² Reg"
		H2_SumDiv="H² SumDiv";
RUN;

/* Delete temporary files */
PROC DATASETS LIBRARY=work;
   DELETE xm_blues xm_blups xm_full xm_reg xm_sum;
RUN;

%MEND H2_BLUE_BLUP;
