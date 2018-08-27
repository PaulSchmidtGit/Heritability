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
/*		fixed to obtain LS-means for the genotypes (i.e. genotype BLUEs).				        */
/*																								*/
/*		SAS/STAT																				*/
/*			Name of genetic effect																*/
/*				ENTRY_NAME=	specifies the genotypic treatment factor (e.g. var, entry, gen, g).	*/	
/*			Dataset 'SolutionF'																	*/
/*				SOLUTIONF= specifies the MIXED / GLIMMIX ODS output with the fixed-effects 		*/
/*				solution vector	from a model where the genotype main effect is random.			*/
/*			Dataset 'SolutionR'																	*/
/*				SOLUTIONR= specifies the MIXED / GLIMMIX ODS output with the random-effects 	*/
/*				solution vector	from the same model as SOLUTIONF.								*/
/*			Dataset 'LSMeans'																	*/
/*				LSMEANS= specifies the MIXED / GLIMMIX ODS output with the least squares means  */
/*				of the genotypes. It must be produced from a model where the genotypic effect 	*/
/*				is fixed.																		*/
/*			Name for output file																*/
/*				OUTPUT= specifies the name for the output dataset.								*/
/*																								*/
/*	Note that in order to prevent complications due to overwritten data, one should not use 	*/
/*	dataset names starting with "xm_" as some are used in this macro.							*/
/*																								*/
/*	First version 27 August 2018																*/
/*																								*/
/*	Written by: Paul Schmidt (Paul.Schmidt@uni-hohenheim.de)									*/
/*																								*/
/************************************************************************************************/

%MACRO H2_BLUE_BLUP(ENTRY_NAME=, LSM_Mu=, SolutionR=, LSM_G=, OUTPUT=);

/* Format genotype BLUE Table */
DATA xm_blues; SET &LSM_G.;
WHERE Effect="&ENTRY_NAME.";
KEEP Effect &ENTRY_NAME. Estimate;
RENAME Estimate=BLUE;
RUN;

/* Format genotype BLUP Table */
DATA xm_blups; SET &SolutionR.;
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
RUN;

/* Save overall mean from in macro variable "Mu_ran" */
DATA &LSM_Mu.; SET &LSM_Mu.;
CALL SYMPUT("Mu_ran", Estimate);
RUN;

/* Expanding BLUE-BLUP table */
DATA xm_full; SET xm_full;
scaledGBLUE=BLUE-&Mu_ran.;
absBLUP=abs(BLUP);
absscaledGBLUE=abs(scaledGBLUE);
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
PROC MEANS DATA=xm_full SUM NOPRINT; VAR absBLUP;
OUTPUT OUT=xm_sumabsblup SUM=sum_absblup;
RUN;
PROC MEANS DATA=xm_full SUM NOPRINT; VAR absscaledGBLUE;
OUTPUT OUT=xm_sumabsscaledGBLUE SUM=sum_absscaledGBLUE;
RUN;

DATA xm_sum;
MERGE xm_sumabsblup xm_sumabsscaledGBLUE;
H2_SumDiv=(sum_absblup/sum_absscaledGBLUE);
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
   DELETE xm_blues xm_blups xm_full xm_reg xm_sumabsblup xm_sumabsscaledGBLUE xm_sum;
RUN;

%MEND H2_BLUE_BLUP;
