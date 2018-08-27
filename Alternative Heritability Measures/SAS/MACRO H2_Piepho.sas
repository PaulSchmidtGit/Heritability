/************************************************************************************************/
/*                              _  _ __   ___ _          _                                      */
/*                             | || |_ ) | _ (_)___ _ __| |_  ___                               */
/*                             | __ /__| |  _/ / -_) '_ \ ' \/ _ \                              */
/*                             |_||_|    |_| |_\___| .__/_||_\___/                              */
/*                                                 |_|                                          */
/*																								*/
/*	This macro computes the entry-based heritability based on the mean variance of a 		    */
/*	difference of two adjusted genotype means (i.e. genotypic BLUEs).							*/				
/*																								*/
/*	Example application code can be found on https://github.com/PaulSchmidtGit/Heritability     */
/*																								*/
/*	This method is based on																		*/
/*		Piepho, H.-P., and J. Möhring. 2007. Computing heritability and selection response from */
/*		unbalanced plant breeding trials. Genetics 177(3):1881–1888. 							*/
/*		[p. 1884, eq.(19)]																	    */
/*																								*/
/*	Requirements/Input:																			*/
/*		The model that is used to analyze the data beforehand should have a genotype main 		*/
/*		effect. This model should then be fitted in two versions: once with the genotype main  	*/
/*		effect random to obtain its variance component and once with the genotype main effect   */
/*		fixed to obtain (the variance of) pairwise differences between LS-means.				*/
/*																								*/
/*		SAS/STAT																				*/
/*			Name of genetic effect																*/
/*				ENTRY_NAME=	specifies the genotypic treatment factor (e.g. var, entry, gen, g).	*/	
/*			Dataset 'covparms'																	*/
/*				COVPARMS= specifies the MIXED / GLIMMIX ODS output with variance components.	*/
/*				It must be produced from a model where the genotypic effect is random.			*/
/*			Dataset 'diffs'																		*/
/*				DIFFS= specifies the MIXED / GLIMMIX ODS output with pairwise comparisons of	*/
/*				adjusted genotype means (BLUE). It must be produced from a model where the 		*/
/*				genotypic effect is fixed.														*/
/*			Name for output file																*/
/*				OUTPUT= specifies the name for the output dataset.								*/
/*																								*/
/*	Note that in order to prevent complications due to overwritten data, one should not use 	*/
/*	dataset names starting with "xm_" as some are used in this macro.							*/
/*																								*/
/*	Version 27 August 2018  																	*/
/*																								*/
/*	Written by: Paul Schmidt (Paul.Schmidt@uni-hohenheim.de)									*/
/*																								*/
/************************************************************************************************/
;

%MACRO H2_piepho(ENTRY_NAME=, COVPARMS=, DIFFS=, OUTPUT=);

/* Extract genotypic variance component from COVPARM output and save it in macro variable "xm_gen_var" */
DATA xm_cp;
	SET  &COVPARMS.;
	WHERE CovParm="&ENTRY_NAME.";
	CALL SYMPUT("xm_gen_var", Estimate);
RUN;

/* Extract standard errors for all pairwise comparisons of genotypes from DIFFS output */
DATA xm_diffs; 
	SET &DIFFS.; 
	var_diffs=stderr**2;	* Obtain variances of a difference "var_diffs" by squaring standard errors of a difference;
RUN;	

/* Obtain mean variance of a difference "xm_avdBLUE_g" */
PROC MEANS DATA=xm_diffs MEAN NOPRINT;											
	VAR var_diffs; 
	OUTPUT  OUT=xm_mean_var MEAN=xm_avdBLUE_g; 
	RUN;

/* Final calculation & formatting */
DATA &OUTPUT.;
	KEEP   	H2_Piepho xm_gen_var xm_avdBLUE_g;
	RETAIN 	H2_Piepho xm_gen_var xm_avdBLUE_g;
	SET xm_mean_var;
	xm_gen_var 	= &xm_gen_var.;
	H2_Piepho	=(&xm_gen_var. / (&xm_gen_var. + 0.5 * xm_avdBLUE_g)); *See Piepho & Möhring (2007) p.1884 eq.(19);
	FORMAT 	H2_Piepho 10.3;
	LABEL 	H2_Piepho	 ="H² Piepho"
			xm_gen_var	 ="Genotypic variance component"
			xm_avdBLUE_g ="Mean variance of a difference of two genotypic BLUEs";
RUN;

/* Delete temporary files */
PROC DATASETS LIBRARY=work;
   DELETE xm_cp xm_diffs xm_mean_var;
RUN;

%MEND H2_piepho;

