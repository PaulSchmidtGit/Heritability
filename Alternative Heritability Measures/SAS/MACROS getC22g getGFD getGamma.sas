/************************************************************************************************
           _      ___ ___ ___                _    ___ ___ ___      
  __ _ ___| |_   / __|_  )_  )  __ _ _ _  __| |  / __|_  )_  )__ _ 
 / _` / -_)  _| | (__ / / / /  / _` | ' \/ _` | | (__ / / / // _` |
 \__, \___|\__|  \___/___/___| \__,_|_||_\__,_|  \___/___/___\__, |
 |___/                                                       |___/ 
************************************************************************************************/
%MACRO getC22g(ENTRY_NAME=, MMEQSOL=, SOLUTIONF=, OUT_C22=xm_c22, OUT_C22g=xm_C22g);
	/* Reduce ODS Output to numeric matrix */
	/***************************************/

	/* Obtain starting row/column of C22 with respect to C */
	DATA xm_max;												
		SET &SOLUTIONF.;
		xm_n=_n_;
		KEEP xm_n;
		RUN;
	PROC MEANS DATA=xm_max NOPRINT;						
		VAR xm_n; 
		OUTPUT max=xm_max OUT=xm_max;
		RUN;
	DATA xm_max;
		SET xm_max;
		CALL SYMPUT ("X",xm_max);
		RUN;
	%LET x=%SYSEVALF(&x+1);										

	/* C22.g */
	PROC MEANS DATA=&MMEQSOL. NOPRINT;							
		WHERE Effect="&ENTRY_NAME."; 	
		VAR row; 
		OUTPUT OUT=xm_temp;
		RUN;
	DATA xm_temp; SET xm_temp;									/* Create variables &min., &max., &colmax. and &colmin. */			
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 
	DATA &OUT_C22g.; 											/* Drop most unwanted rows and columns of C */
		SET &MMEQSOL.; 
		if row<&min. then delete;
		if row>&max. then delete;
		KEEP &colmin.-&colmax.; 
		RUN;

	/* C22 */
	PROC MEANS DATA=&MMEQSOL. NOPRINT;							
		VAR row; 
		OUTPUT OUT=xm_temp2;
		RUN;
	DATA xm_temp2; SET xm_temp2;								/* Create variables &min., &max., &colmax. and &colmin. for MMEQSOL */			
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 
	DATA &OUT_C22.; 											/* Drop most unwanted columns of C */
		SET &MMEQSOL.; 
		if row<&x. then delete;
		KEEP col&x.-&colmax.; 
		RUN;

	PROC DATASETS LIBRARY=work;								
   		DELETE xm_temp xm_temp2 xm_max ;
		RUN; QUIT;
%Mend getc22g;

/************************************************************************************************
           _      ___     ___                _   ___  
  __ _ ___| |_   / __|   | __|  __ _ _ _  __| | |   \ 
 / _` / -_)  _| | (_ |_  | _|  / _` | ' \/ _` | | |) |
 \__, \___|\__|  \___( ) |_|   \__,_|_||_\__,_| |___/ 
 |___/               |/                               
************************************************************************************************/
%MACRO getGFD (ENTRY_NAME=, G=, OUT_m_g=m_g, OUT_m_f=m_f, OUT_m_d=m_d);
	/* Obtain and save dimensions of G in variables &min., &max., &comax. and &colmin. */
	PROC MEANS DATA=&G. NOPRINT;							
		VAR row; 
		OUTPUT OUT=xm_temp;
		RUN;
	DATA xm_temp; SET xm_temp;											
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 

	/* Reduce ODS output to numeric matrix G */
	DATA m_G; 											
		SET &G.;
		KEEP &colmin.-&colmax.; 
		RUN;

	DATA m_f; 											
		SET &G.;
		row=_n_;
		KEEP &colmin.-&colmax.; 
		RUN;

	/* Obtain and save dimensions of D (i.e. part of G referring to entry main effect, also referred to as G.g) in variables &min., &max., &comax. and &colmin. */
	PROC MEANS DATA=&G. NOPRINT;						
		WHERE Effect="&ENTRY_NAME."; 
		VAR row; OUTPUT OUT=xm_temp;
		RUN;
	DATA xm_temp; SET xm_temp;							/* Save dimensions in variables &min., &max. and &colmin. */			
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 

	DATA m_f; 											
		SET m_f;
		row=_n_;
		IF row<&min. THEN DELETE; 
		IF row>&max. THEN DELETE; 
		DROP row;
		RUN;

	DATA m_d; 											
		SET &G.;
		row=_n_;
		IF row<&min. THEN DELETE; 
		IF row>&max. THEN DELETE; 
		KEEP &colmin.-&colmax.; 
		RUN;

%MEND getGFD;

/************************************************************************************************
           _      ___                       
  __ _ ___| |_   / __|__ _ _ __  _ __  __ _ 
 / _` / -_)  _| | (_ / _` | '  \| '  \/ _` |
 \__, \___|\__|  \___\__,_|_|_|_|_|_|_\__,_|
 |___/                                      
************************************************************************************************/
%MACRO getGamma(m_C22=xm_C22, m_G=m_g, m_F=m_f, m_D= m_d);
	PROC IML;
		USE &m_C22.; READ ALL INTO C22;
		USE &m_G.;   READ ALL INTO G;
		USE &m_F.;   READ ALL INTO F;
		USE &m_D.;   READ ALL INTO D;

		M		 = G-C22; 					/* (12) See Mclean, Sanders, Stroup (1991)      */
		inv_G	 = inv(G);					/* Inverse of G 								*/
		Q		 = F*inv_G*M*inv_G*t(F);	/* (10) Q = F G^-1 M G^-1 F' 					*/
		omega	 = (D||Q)//(Q||Q);  		/* (10) Create Omega 							*/
		CALL SVD(u, lambda, v, omega); 		/* (13) Cholesky Decompostion of omega part I  	*/
		m_gamma = u*DIAG(SQRT(lambda)); 	/* (13) Cholesky Decompostion of omega part II 	*/

		CREATE m_gamma FROM m_gamma; APPEND FROM m_gamma;
	QUIT; RUN;
%MEND getGamma;
