/************************************************************************************************
           _      ___ ___ ___                _    ___ ___ ___      
  __ _ ___| |_   / __|_  )_  )  __ _ _ _  __| |  / __|_  )_  )__ _ 
 / _` / -_)  _| | (__ / / / /  / _` | ' \/ _` | | (__ / / / // _` |
 \__, \___|\__|  \___/___/___| \__,_|_||_\__,_|  \___/___/___\__, |
 |___/                                                       |___/ 
************************************************************************************************/

%MACRO getC22g(ENTRY_NAME=, MMEQSOL=, EXCLUDE_ZEROS=); 	/* Genotype main effect must be first random effect in model */
	/* Reduce ODS Output to numeric matrix */
	/***************************************/
	PROC MEANS DATA=&MMEQSOL. NOPRINT;						/* Obtain starting row/column of C22 with respect to C */
		WHERE Effect="&ENTRY_NAME."; 
		VAR row; OUTPUT OUT=xm_temp;
		RUN;
	DATA xm_temp; SET xm_temp;								/* Save variables &min. and &colmin. */			
		CALL SYMPUT(CATS("Col",_STAT_),CATS("Col",Row));
		CALL SYMPUT(           _STAT_ ,           Row);
		RUN; 
	DATA xm_temp; 											/* Drop most unwanted columns of C */
		KEEP &colmin.-Col999999999; 
		SET &MMEQSOL.; 
		RUN;
	PROC IML;
		USE xm_temp; READ ALL INTO xm_C22;						/* Drop all unwanted colums and rows of C to obtain C22 */
		m_C22 = xm_C22[&min.:nrow(xm_C22), 1:ncol(xm_C22)-1];		
		CREATE m_C22 FROM m_C22; APPEND FROM m_C22;
		QUIT; RUN;

	%IF EXCLUDE_ZEROS="TRUE" %THEN *** This step is still under construction ;
		%DO;
			/* Eliminate 0 rows and columns of matrix       */
			/* in case at leas one VC was estimated to be 0 */
			/************************************************/
			DATA xm_temp;									/* Calculate sum of every row */
				SET m_C22;
			 	ARRAY n {*} _NUMERIC_;
			 	sum=SUM(of n[*]);
				RUN;
			DATA xm_temp; 									/* Add running number N */			
				RETAIN N sum;
				SET xm_temp;
				N=_N_;
			RUN;
			DATA xm_temp2;									/* Delete 0-rows */
				DROP N sum;
				SET xm_temp;
				WHERE sum ne 0;
			RUN;
			DATA xm_temp3;									/* Obtain names of 0-columns */
				KEEP N nr dropcols;
				SET xm_temp;
				WHERE sum=0;
				nr=_N_;
				dropcols=CATS("Col",N);
			RUN;
				/* DO loop only executed if 0-lines are present */
				%LET check = %SYSFUNC(OPEN(work.xm_temp3,is));	
				%IF %SYSFUNC(ATTRN(&check,NOBS))>0 %THEN %DO;	
					PROC SQL NOPRINT;								/* Save 0-column names into list */
						SELECT strip(dropcols) INTO :droplist 
						separated BY ' ' FROM xm_temp3;
						QUIT;
					DATA xm_temp;									/* Delete 0-columns */
					  	SET xm_temp2 (DROP=&droplist
						RUN;
					PROC IML;										/* Rename Columns to Col1-Coln */
						USE xm_temp; READ ALL INTO m_C22;
						CREATE m_C22 FROM m_C22; APPEND FROM m_C22; 
						QUIT; RUN
				%END;
				%LET rc=%SYSFUNC(CLOSE(&check));
	%END;

	/* C22 is obtained. Now C22g */
	/*****************************/
	PROC MEANS DATA=&MMEQSOL. NOPRINT;
		WHERE Effect="&ENTRY_NAME."; 
		VAR row;
		OUTPUT OUT=xm_dim(DROP= _FREQ_ _TYPE_);
		RUN;
	DATA xm_dim; SET xm_dim;
		CALL SYMPUT(_STAT_,Row);
		RUN;
	PROC IML;
		USE m_C22; READ ALL INTO xm_C22;
		m_C22g			= xm_C22[1:&N.,1:&N.]; 
		CREATE m_C22g FROM m_C22g; APPEND FROM m_C22g; 
		QUIT; RUN;

	/* Clean up: delete temporary file */
	/***********************************/
	PROC DATASETS LIBRARY=work;								
   		DELETE xm_temp xm_temp2 xm_temp3 xm_dim;
		RUN; QUIT;

%MEND getC22g;

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

	%IF EXCLUDE_ZEROS="TRUE" %THEN 
		%DO;
			/* Eliminate 0 rows and columns of matrix       */
			/* in case at leas one VC was estimated to be 0 */
			/************************************************/
			DATA xm_temp;									/* Calculate sum of every row */
				SET xm_Gx;
			 	ARRAY n {*} _NUMERIC_;
			 	sum=SUM(of n[*]);
				RUN;
			DATA xm_temp; 									/* Add running number N */			
				RETAIN N sum;
				SET xm_temp;
				N=_N_;
			RUN;
			DATA xm_temp2;									/* Delete 0-rows */
				DROP N sum;
				SET xm_temp;
				WHERE sum ne 0;
			RUN;
			DATA xm_temp3;									/* Obtain names of 0-columns */
				KEEP N nr dropcols;
				SET xm_temp;
				WHERE sum=0;
				nr=_N_;
				dropcols=CATS("Col",N);
			RUN;

				/* DO loop only executed if 0-lines are present */
				%LET check = %SYSFUNC(OPEN(work.xm_temp3,is));	
				%IF %SYSFUNC(ATTRN(&check,NOBS))>0 %THEN %DO;
					PROC SQL NOPRINT;						/* Save 0-column names into list */
						SELECT strip(dropcols) INTO :droplist 
						separated BY ' ' FROM xm_temp3;
						QUIT;
					DATA xm_Gx;								/* Delete 0-columns */
					  	SET xm_temp2 (DROP=&droplist);
						RUN;
				%END;
				%LET rc=%SYSFUNC(CLOSE(&check));

	%END;

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
