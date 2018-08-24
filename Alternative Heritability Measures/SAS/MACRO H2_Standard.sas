/************************************************************************************************/
/* 						 _  _ __      _                _             _ 							*/
/* 						| || |_ )  __| |_ __ _ _ _  __| |__ _ _ _ __| |							*/
/* 						| __ /__| (_-<  _/ _` | ' \/ _` / _` | '_/ _` |							*/
/* 						|_||_|    /__/\__\__,_|_||_\__,_\__,_|_| \__,_|							*/
/*                                                												*/
/*																								*/
/*	This macro computes the standard, broad-sense heritability for multi-environment trials		*/
/*	given a dataset with genotype means per environment (=year-location-combination).			*/	
/*																								*/
/*	Example data and code can be found below at the end of this file.							*/
/*																								*/
/*	Requirements/Input:																			*/	
/*		The model that is used to analyze the MET data beforehand should have random effects 	*/
/*		for year, location, genotype, all their two-way interaction effects and their three-way */
/*		interaction effect (see example).														*/
/*																								*/
/*		SAS/STAT																				*/
/*			Dataset 'data'																		*/
/*				DATA= specifies the input dataset for your model. For this macro it is assumed  */
/*  			that it is a dataset of genotype means per environment (=year-location-			*/
/*				combination).																	*/
/*			Dataset 'covparms'																	*/
/*				COVPARMS= specifies the MIXED / GLIMMIX ODS output with variance components.	*/
/*				It must be produced from a model where the genotypic effect is random.			*/
/*			Number of replicates 																*/
/*				NUMBER_REP= specifies the number of replicates that were used in the trial*     */
/*				[e.g. 2, 3, 4]																	*/
/*			Name of genetic effect																*/
/*				ENTRY_NAME=	specifies the genotypic treatment factor as used in the model		*/	
/*				[e.g. G, gen, var, entry, Sorte]												*/
/*			Name of year effect																	*/
/*				YEAR_NAME= specifies the year factor as used in the model						*/
/*				[e.g. Y, year, season, Jahr]													*/
/*			Name of location effect																*/
/*				LOCATION_NAME= specifies the location factor as used in the model				*/
/*				[e.g. L, loc, site, Ort]														*/
/*			Name for output file																*/
/*				OUTPUT= specifies the name for the two output datasets.							*/
/*																								*/
/*	Note that in order to prevent complications due to overwritten data, one should not use 	*/
/*	dataset names starting with "xm_" as they are used in this macro.							*/
/*																								*/
/*	Version 13 September 2017																	*/
/*																								*/
/*	Written by: Paul Schmidt (Paul.Schmidt@uni-hohenheim.de)									*/
/*																								*/
/************************************************************************************************/

%MACRO H2_standard(DATA=, COVPARMS=, NUMBER_REP=, ENTRY_NAME=, YEAR_NAME=, LOCATION_NAME=, OUTPUT=);

/* Get number of years and number of locations */
PROC FREQ DATA=&DATA. NLEVELS; 
	TABLES &YEAR_NAME. &LOCATION_NAME. ; 
	ODS OUTPUT NLevels=xm_nlevel; 
	RUN;	
DATA xm_nlevel;
	SET  xm_nlevel;
	shortname="xxxxx";
	IF TableVar="&YEAR_NAME." 	  THEN shortname="num_y";
	IF TableVar="&LOCATION_NAME." THEN shortname="num_l";
	CALL SYMPUT(shortname, NLevels);
	RUN;

/* Combine VC and n */
DATA xm_comb;
	SET &COVPARMS.;
	number_&YEAR_NAME.=1; 
	number_&LOCATION_NAME.=1;  
	number_rep=1; 			  
	IF FIND(CovParm,"&YEAR_NAME.")>0 	 THEN DO; number_&YEAR_NAME.	=&num_y.; END;
	IF FIND(CovParm,"&LOCATION_NAME.")>0 THEN DO; number_&LOCATION_NAME.=&num_l.; END;
	IF FIND(CovParm,"Residual")>0 	     THEN DO; number_&YEAR_NAME.	=&num_y.; 
												  number_&LOCATION_NAME.=&num_l.; 
												  number_rep			=&NUMBER_REP.;  
											  END;
	divisor=number_&YEAR_NAME.*number_&LOCATION_NAME.*number_rep;
	sigma_quotient=estimate/divisor;
	IF (FIND(CovParm,"Residual")=0 AND FIND(CovParm,"&ENTRY_NAME.")=0) OR CovParm = "&ENTRY_NAME." THEN DO; sigma_quotient=.; 
																											divisor=.; 
																											number_&YEAR_NAME.=.; 	
																											number_&LOCATION_NAME.=.; 
																											number_rep=.; 
																										END;
	IF number_&YEAR_NAME.=1 	THEN number_&YEAR_NAME.=.;
	IF number_&LOCATION_NAME.=1 THEN number_&LOCATION_NAME.=.;
	IF number_rep=1 			THEN number_rep=.;
	RUN;

/* Obtain genotypic variance as macro variable "gen_var" */
DATA xm_cp;
	SET  &COVPARMS.;
	WHERE CovParm="&ENTRY_NAME.";
	CALL SYMPUT("gen_var", estimate);
	RUN;

/* Obtain phenotypic variance */
PROC MEANS DATA=xm_comb SUM NOPRINT;
	WHERE Covparm ne "&ENTRY_NAME.";
	VAR sigma_quotient;
	OUTPUT OUT=xm_phen_var(DROP=_TYPE_ _FREQ_) SUM=phen_var; 
	RUN;

/* Final calculation & formatting of two outputs*/
DATA &OUTPUT._1;
	KEEP   	H2_Standard gen_var phen_var;
	RETAIN 	H2_Standard gen_var phen_var;
	SET xm_phen_var;
	gen_var 	= &gen_var.;
	H2_Standard	=(&gen_var. / (&gen_var. + phen_var));
	FORMAT 	H2_Standard 8.2 
			gen_var phen_var 10.3;
	LABEL 	H2_Piepho ="Standard H²"
			gen_var	  ="Genetic variance component"
			phen_var  ="Phenotypic variance";
	RUN;

DATA &OUTPUT._2;
	SET xm_comb;
	FORMAT estimate sigma_quotient 8.2;
	LABEL number_&YEAR_NAME.	 ="Number &YEAR_NAME."
		  number_&LOCATION_NAME. ="Number &LOCATION_NAME."
		  number_rep			 ="Number of replicates"
		  estimate				 ="Variance component = dividends for phenotypic variance summands"
		  divisor				 ="Number product = Divisors for phenotypic variance summands"
		  sigma_quotient		 ="Quotients = phenotypic variance summands"; 
	RUN;

%MEND;
