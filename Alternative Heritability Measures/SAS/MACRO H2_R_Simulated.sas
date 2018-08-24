/************************************************************************************************/
/*			  ___ _           _      _          _    _  _ __                _   ___ 			*/
/*			 / __(_)_ __ _  _| |__ _| |_ ___ __| |  | || |_ )  __ _ _ _  __| | | _ \			*/
/*			 \__ \ | '  \ || | / _` |  _/ -_) _` |  | __ /__| / _` | ' \/ _` | |   /			*/
/*			 |___/_|_|_|_\_,_|_\__,_|\__\___\__,_|  |_||_|    \__,_|_||_\__,_| |_|_\			*/
/*                                                                        						*/
/*	This macro computes (i) the measure of heritability as the simulated expected correlation  	*/
/*	of predicted and true genotypic value and (ii) a simulated expected response to selection.  */
/*																								*/
/*	Example code can be found below at the end of this file.									*/
/*																								*/
/*	This method is based on																		*/
/*		Piepho, H.-P., and J. Möhring. 2007. Computing heritability and selection response from */
/*		unbalanced plant breeding trials. Genetics 177(3):1881–1888.							*/
/*	Comments in the code refer to respective equations from this article [i.e. (1), (2) etc.]   */
/*																								*/
/*	This macro actually consists of 4 macros:													*/
/*	   (i) %GET_C22   extracts the C22 matrix from ods output of SAS							*/
/*	  (ii) %GET_GFD   extracts the G, F and D matrices from ods output of SAS					*/
/*	 (iii) %GET_GAMMA calculates Gamma from C22, G, F and D				 						*/
/*	  (iv) %GET_HRsim simulates H² and R based on Gamma											*/
/*																								*/
/*	Requirements/Input:																			*/
/*		Preliminary fitted mixed model with random genotypic treatment factor					*/
/*			A genotypic treatment factor should be defined as the >>>first<<< effect in the  	*/
/*			RANDOM statement of PROC MIXED.														*/
/*		SAS/STAT SAS/IML																		*/
/*			Name of genetic effect 'ENTRY_NAME'													*/
/*				ENTRY_NAME=	specifies the genotypic treatment factor (e.g. var, entry, gen, g).	*/
/*			Dataset 'MMEQSOL'																	*/
/*				MMEQSOL= specifies the MIXED / GLIMMIX ODS output with the mixed model			*/
/*				equations solution. 															*/
/*			Dataset 'GMATRIX'																	*/
/*				G= specifies the MIXED / GLIMMIX ODS output with the estimated G matrix.		*/
/*																								*/
/*	Output:																						*/
/*		C22, G, F, D and Gamma matrices as datasets xm_c22, xm_g, xm_f, xm_d and xm_gamma		*/
/*		Heritability estimate as xm_h2_gg														*/
/*		Response to selection depending on the number of selected genotypes as xm_R      		*/
/*	Note that in order to prevent complications due to overwritten data, one should not use 	*/
/*	dataset names starting with "xm_" as some are used in this macro.							*/
/*																								*/
/*	First version 11 October 2017																*/
/*																								*/
/*	Written by: Paul Schmidt (Paul.Schmidt@uni-hohenheim.de)									*/
/*																								*/
/************************************************************************************************/

/************************************************************************************************
  ___ _           _      _         _  _ __                _   ___ 
 / __(_)_ __ _  _| |__ _| |_ ___  | || |_ )  __ _ _ _  __| | | _ \
 \__ \ | '  \ || | / _` |  _/ -_) | __ /__| / _` | ' \/ _` | |   /
 |___/_|_|_|_\_,_|_\__,_|\__\___| |_||_|    \__,_|_||_\__,_| |_|_\
                                                                  
************************************************************************************************/
%MACRO H2RSim(m_Gamma=, n_sim=, H_OUT=, R_OUT=);
	PROC IML;
		USE &m_Gamma.; READ ALL INTO gamma;
		n = nrow(gamma);				/* n = number of rows in gamma */
		z = j(n,1,0);					/* z = vector with n rows, full of 0s */
		r2=0; a=0; b=0; c=0;			/* Predefine r2[correlation g-g_hat], a[Cov g-g_hat], b[Var_g] and c[Var g_hat] as 0 */
		sim_max  = &n_sim.;				/* sim_max 	= number of simulations */
		sel_mean = j(n/2,1,0);			/* sel_mean	= vector half the size of gamma full of 0s */

		DO sim=1 TO sim_max;
		  DO i=1 TO n;
		    z[i]=normal(sim);			/*  replace 0s in z with standard normally distributed random numbers. Seed=-1 */
		  END;
		  w 	  = gamma*z;			/* (14) Multiply decomposed Omega (i.e. gamma) with standard normally distributed random numbers (i.e. z) to create simulated vector w */
		  g 	  = w[1:n/2];       	/* g     =      true genetic value = upper half of w */
		  g_hat   = w[n/2+1:n];			/* g_hat = estimated genetic value = lower half of w */
		  r_g_hat = rank(g_hat);		/* Give ranks to all g_hats 						 */
		  v 	  = g_hat||r_g_hat||g; 	/* Create v by putting estimated genetic value, rank of estimated genetic value and true genetic value next to each other */
		  CALL sort(v, {1}, {1}); 		/* Sort v so that biggest g_hat is on top 			 */
			  DO j=1 TO n/2;
			    sel_mean[j] = sel_mean[j]+sum(v[1:j,3])/j;	/* Selection index for all possible numbers of selected genotypes from best 1 to best n */
			  END;
		  g		= g-g[:];						
		  g_hat	= g_hat-g_hat[:];			
		  r2	= r2+ ( t(g_hat)*g )**2/( t(g)*g*t(g_hat)*g_hat ); 	/* Sum up all Correlations g-g_hat from each simulation */
		  a		= a+t(g_hat)*g;										/* Sum up all Covariances  g-g_hat 	*/
		  b		= b+t(g)*g;											/* Sum up all Variances g 			*/
		  c		= c+t(g_hat)*g_hat;									/* Sum up all Variances g_hat 		*/
		END;

		H2_gg	=(r2/sim_max)*100; 									/* Get average Correlations g-g_hat */
		/*H2_gg_b =(a*a/(b*c))*100*/								/* Get average Correlations g-g_hat (alternative calculation) */
		sel_mean=sel_mean/sim_max;									/* Get average Selection Gain/Response to Selection R 		  */

		CREATE xm_H2_gg VAR {H2_gg /*H2_gg_b*/}; APPEND;
		CREATE xm_R     VAR {sel_mean};      	 APPEND;
	QUIT; RUN; 

	DATA &H_OUT.; 
		SET xm_H2_gg;
		LABEL  H2_gg  ="H² as r² of (g-g^)"
			 /*H2_gg_b="H² as r² of (g-g^) [alternative]"*/;
		FORMAT H2_gg 
 			 /*H2_gg_b*/ 8.2;
		RUN;

	DATA &R_OUT.;
		RETAIN n_selected SEL_MEAN;
		KEEP   n_selected SEL_MEAN;
		SET xm_R;
		n_selected=_N_;
		LABEL 	n_selected="Number of selected genotypes"
				SEL_MEAN="Simulated R";
		FORMAT	SEL_MEAN 8.3;
		RUN;
%MEND H2RSim;