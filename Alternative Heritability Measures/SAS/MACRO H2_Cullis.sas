%MACRO H2_cullis(ENTRY_NAME=, COVPARMS=, m_c22g=, OUTPUT=);

	/* Extract genotypic variance component from COVPARM output and save it in macro variable "xm_gen_var" */
	DATA xm_cp; SET &COVPARMS.;
		WHERE CovParm="&ENTRY_NAME.";
		CALL SYMPUT("xm_gen_var", Estimate);
		RUN;

	/* Obtain average variance of a difference between genotype BLUPs "xm_avdBLUP_g" and calculate H2_cullis */
	PROC IML;
		USE &m_c22g.; READ ALL INTO xm_c22g;
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

%MEND H2_cullis;
