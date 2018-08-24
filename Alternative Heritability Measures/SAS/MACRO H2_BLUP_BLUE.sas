/************************************************************************************************/
/*    _  _ __   ___    _   _            _             ___ _   _   _ ___    ___ _   _   _ ___    */
/*   | || |_ ) | __|__| |_(_)_ __  __ _| |_ ___ ___  | _ ) | | | | | __|__| _ ) | | | | | _ \   */
/*   | __ /__| | _|(_-<  _| | '  \/ _` |  _/ -_|_-<  | _ \ |_| |_| | _|___| _ \ |_| |_| |  _/   */
/*   |_||_|    |___/__/\__|_|_|_|_\__,_|\__\___/__/  |___/____\___/|___|  |___/____\___/|_|     */
/*                                                                                              */
/* This macro computes two measures of heritability based on genotypic BLUEs and BLUPs.			*/
/*																								*/
/* Example code can be found below at the end of this file.										*/
/*																								*/
/* The methods are based on unpublished work of the author of this macro.						*/
/*																								*/
/*	Requirements/Input:																			*/
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
/*	First version 07 November 2016																*/
/*																								*/
/*	Written by: Paul Schmidt (Paul.Schmidt@uni-hohenheim.de)									*/
/*																								*/
/************************************************************************************************/

%MACRO H2_BLUE_BLUP(ENTRY_NAME=, LSM_Mu=, SolutionR=, LSM_G=, OUTPUT=);

/* Prepare BLUE-BLUP Table */
DATA m_blues; SET &LSMEANS.;
WHERE Effect="&ENTRY_NAME.";
KEEP Effect &ENTRY_NAME. Estimate;
RENAME Estimate=BLUE;
RUN;

DATA m_blups; SET &SolutionR.;
WHERE Effect="&ENTRY_NAME.";
KEEP Effect &ENTRY_NAME. Estimate;
RENAME Estimate=BLUP;
RUN;

PROC SORT DATA=m_blues; BY &ENTRY_NAME.; RUN;
PROC SORT DATA=m_blups; BY &ENTRY_NAME.; RUN;
DATA m_full;
MERGE m_blues m_blups;
BY &ENTRY_NAME.; 
RUN;

/* Overall mean in g_ran */
DATA &LSM_Mu.; SET &LSM_Mu.;
CALL SYMPUT("Mu_ran", Estimate);
RUN;

/* Expanding BLUE-BLUP Table */
DATA m_full; SET m_full;
scaledGBLUE=BLUE-&Mu_ran.;
absBLUP=abs(BLUP);
absscaledGBLUE=abs(scaledGBLUE);
RUN;

/* H2 Regression */
PROC MIXED DATA=m_full; 
MODEL BLUP=scaledGBLUE /S NOINT;
ODS OUTPUT SolutionF=m_reg;
RUN;
DATA m_reg; SET m_reg;
KEEP h2_reg;
h2_reg=Estimate*100;
RUN;


/* H2 sumdiv */
PROC MEANS DATA=m_full SUM NOPRINT; VAR absBLUP;
OUTPUT OUT=m_sumabsblup SUM=sum_absblup;
RUN;
PROC MEANS DATA=m_full SUM NOPRINT; VAR absscaledGBLUE;
OUTPUT OUT=m_sumabsscaledGBLUE SUM=sum_absscaledGBLUE;
RUN;

DATA m_sum;
MERGE m_sumabsblup m_sumabsscaledGBLUE;
H2_sumdiv=(sum_absblup/sum_absscaledGBLUE)*100;
KEEP H2_sumdiv;
RUN; 

/* Output */
DATA &OUTPUT.;
RETAIN h2_reg H2_sumdiv;
MERGE m_reg m_sum;
FORMAT _numeric_ 8.2;
LABEL 	h2_reg="H² BLUP~BLUE"
		h2_sum="H² Sumdiv";
RUN;

%MEND H2_BLUE_BLUP;
