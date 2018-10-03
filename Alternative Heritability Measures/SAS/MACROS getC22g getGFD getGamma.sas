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
	DATA xm_max;SET &SOLUTIONF.;xm_n=_n_;keep xm_n;run;
	PROC MEANS DATA=xm_max NOPRINT;						/* Obtain starting row/column of C22 with respect to C */
		VAR xm_n; OUTPUT max=xm_max OUT=xm_max;
		RUN;
	DATA xm_max;set xm_max;
		CALL SYMPUT ("X",xm_max);
		RUN;
		%let x=%sysevalf(&x+1);
	PROC MEANS DATA=&MMEQSOL. NOPRINT;						/* Obtain starting row/column of C22 with respect to C */
		WHERE Effect="&ENTRY_NAME."; 
		VAR row; OUTPUT OUT=xm_temp;
		RUN;
	PROC MEANS DATA=&MMEQSOL. NOPRINT;						/* Obtain starting row/column of C22 with respect to C */
		VAR row; OUTPUT OUT=xm_temp2;
		RUN;
	DATA xm_temp; SET xm_temp;								/* Create variables &min., &max., &colmax. and &colmin. */			
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 
	DATA &OUT_C22g.; 											/* Drop most unwanted rows and columns of C */
		SET &MMEQSOL.; 
		if row<&min then delete;
		if row>&max then delete;
		KEEP &colmin.-&colmax; 
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

%MACRO getGFD(G=, ENTRY_NAME=, EXCLUDE_ZEROS=);
	/* Reduce ODS Output to numeric matrix */
	/***************************************/
	PROC MEANS DATA=&G. NOPRINT;						/* Obtain dimensions of D (=part of G referring to entry main effect) */
		WHERE Effect="&ENTRY_NAME."; 
		VAR row; OUTPUT OUT=xm_temp;
		RUN;
	DATA xm_temp; SET xm_temp;								/* Save dimensions in variables &min., &max. and &colmin. */			
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 
	DATA xm_Gx; 											/* Reduce ODS output G to purely numeric G-matrix */
		KEEP &colmin.-Col999999999; 
		SET &G.; 
		RUN;

	PROC IML;												/* Create G-, F- and D-Matrix */
		USE xm_Gx; READ ALL INTO xm_G;
		m_G  = xm_G;
		m_F = xm_G[&min.:&max.,            ];
		m_D = xm_G[&min.:&max., &min.:&max.];
		CREATE m_G FROM m_G; APPEND FROM m_G; 
		CREATE m_F FROM m_F; APPEND FROM m_F; 
		CREATE m_D FROM m_D; APPEND FROM m_D; 
		QUIT; RUN;

	/* Clean up: delete temporary file */
	/***********************************/
	PROC DATASETS LIBRARY=work;
	   	DELETE xm_temp xm_temp2 xm_temp3 xm_Gx;
		RUN; QUIT;

%MEND getGFD;

/************************************************************************************************
           _      ___                       
  __ _ ___| |_   / __|__ _ _ __  _ __  __ _ 
 / _` / -_)  _| | (_ / _` | '  \| '  \/ _` |
 \__, \___|\__|  \___\__,_|_|_|_|_|_|_\__,_|
 |___/                                      
************************************************************************************************/
%MACRO getGamma(m_C22=, m_G=, m_F=, m_D= );
	PROC IML;
		USE &m_C22.; READ ALL INTO C22;
		USE &m_G.;   READ ALL INTO G;
		USE &m_F.;   READ ALL INTO F;
		USE &m_D.;   READ ALL INTO D;

		M		 = G-C22; 					/* (12) Still unclear why, but definitely true. See Mclean, Sanders, Stroup (1991) */
		inv_G	 = inv(G);					/* Inverse of G 								*/
		Q		 = F*inv_G*M*inv_G*t(F);	/* (10) Q = F G^-1 M G^-1 F' 					*/
		omega	 = (D||Q)//(Q||Q);  		/* (10) Create Omega 							*/
		CALL SVD(u, lambda, v, omega); 		/* (13) Cholesky Decompostion of omega part I  	*/
		m_gamma = u*DIAG(SQRT(lambda)); 	/* (13) Cholesky Decompostion of omega part II 	*/

		CREATE m_gamma FROM m_gamma; APPEND FROM m_gamma;
		QUIT; RUN;
%MEND getGamma;
