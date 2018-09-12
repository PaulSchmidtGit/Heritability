/************************************************************************************************/
/*                             _  _ ___    ___     _ _ _                                        */
/*                            | || |_  )  / __|  _| | (_)___                                    */
/*                            | __ |/ /  | (_| || | | | (_-<                                    */
/*                            |_||_/___|  \___\_,_|_|_|_/__/                                    */
/*																								*/
/*	This macro computes the entry-based heritability based on the mean variance of a 		    */
/*	difference of two genotype BLUPs.                                 							*/				
/*																								*/
/*	Example application code can be found on https://github.com/PaulSchmidtGit/Heritability     */
/*																								*/
/*	This method is based on																		*/
/*		Cullis, B. R., A. B. Smith, and N. E. Coombes. 2006. On the design of early generation  */
/*      variety trials with correlated data. Journal of Agricultural, Biological, and           */
/*      Environmental Statistics 11(4):381–393. 							                    */
/*		[p. 385]																	            */
/*																								*/
/*	Requirements/Input:																			*/
/*		The model that is used to analyze the data beforehand should have a random genotype     */
/*      main effect in order to obtain its variance component and the variance of pairwise      */
/*      differences between genotype BLUPs. Furthermore, the genotype main effect must be       */
/*		the first random effect written in the model for this macro to work.                    */
/*																								*/
/*		SAS/STAT																				*/
/*			Name of genetic effect																*/
/*				ENTRY_NAME=	specifies the genotypic treatment factor (e.g. var, entry, gen, g).	*/	
/*			Dataset 'covparms'																	*/
/*				COVPARMS= specifies the MIXED / GLIMMIX ODS output with variance components.	*/
/*				It must be produced from a model where the genotypic effect is random.			*/
/*			Dataset 'm_c22g'																	*/
/*				This dataset should contain the estimated variance-covariance matrix of the     */
/*              genotype BLUPs. Note that it is not straightforward to obtain this matrix in    */
/*              SAS, as it needs to be extracted from the MMEqSol= MIXED / GLIMMIX ODS output.  */
/*              We do, however, provide another MACRO named "getC22g" that does exactly that.   */
/*              It is also available on: https://github.com/PaulSchmidtGit/Heritability         */
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

%MACRO H2_cullis(ENTRY_NAME=, COVPARMS=, MMEQSOL=, OUTPUT=);

	/* Extract genotypic variance component from COVPARM output and save it in macro variable "xm_gen_var" */
	DATA xm_cp; SET &COVPARMS.;
		WHERE CovParm="&ENTRY_NAME.";
		CALL SYMPUT("xm_gen_var", Estimate);
		RUN;

	/* Extract C22g Matrix "m_c22g" from MMEQSOL via getC22g Macro from GitHub */
	filename _inbox "%sysfunc(getoption(work))/MACROS getC22g getGFD getGamma.sas";
		proc http method="get" 
		url="https://raw.githubusercontent.com/PaulSchmidtGit/Heritability/master/Alternative%20Heritability%20Measures/SAS/MACROS%20getC22g%20getGFD%20getGamma.sas" out=_inbox;
		run; %Include _inbox; filename _inbox clear;
	%getC22g(ENTRY_NAME=&ENTRY_NAME., MMEQSOL=&MMEQSOL.);

	/* Obtain average variance of a difference between genotype BLUPs "xm_avdBLUP_g" and calculate H2_cullis */
	PROC IML;
		USE m_c22g; READ ALL INTO xm_c22g;
		n_g          = nrow(xm_c22g);
		xm_avdBLUP_g = 2/n_g*(trace(xm_c22g)-(sum(xm_c22g)-trace(xm_c22g))/(n_g-1));
		H2_cullis    = 100*(1-xm_avdBLUP_g/(2*&xm_gen_var.));
		CREATE xm_1 FROM xm_avdBLUP_g; APPEND FROM xm_avdBLUP_g;
		CREATE xm_2 FROM H2_cullis; APPEND FROM H2_cullis;
		QUIT; RUN;

	/* Final calculation & formatting */
	DATA &OUTPUT.;
		KEEP   	H2_Cullis xm_gen_var xm_avdBLUP_g;
		RETAIN 	H2_Cullis xm_gen_var xm_avdBLUP_g;
		MERGE xm_2(RENAME=(COL1=H2_Cullis)) xm_1(RENAME=(COL1=xm_avdBLUP_g));
		xm_gen_var=&xm_gen_var.;
		FORMAT 	H2_Cullis 8.2 
				xm_gen_var xm_avdBLUP_g 10.3;
		LABEL 	H2_Cullis	 ="H² Cullis"
				xm_gen_var	 ="Genotypic variance component"
				xm_avdBLUP_g ="Average variance of a difference of two genotypic BLUPs";
		RUN;

	/* Delete temporary files */
	PROC DATASETS LIBRARY=work;
   		DELETE xm_1 xm_2;
	RUN;

%MEND H2_cullis;
