/* H2 Oakey */

%MACRO H2Oakey(m_C22g=, m_D= ,OUTPUT=H2Oakey);

	PROC IML;
		USE &m_C22g.; READ ALL INTO C22g;
		USE &m_D.;    READ ALL INTO D;

		n_g		 = nrow(D);						  /* number of genotypes */
		inv_D	 = inv(D);						  /* inverse of D */
		M        = I(nrow(D)) - (inv_D * C22g);
		CALL SVD(u, eM, v, M); 					  /* get eigenvalues */
		xm_H2Oak  = sum(eM)/(nrow(eM)-1);         /* mean of eigenvalues */

		CREATE xm_H2Oak FROM xm_H2Oak; APPEND FROM xm_H2Oak;
	QUIT; RUN;

	/* Final formatting */
	DATA &OUTPUT.;
		SET xm_H2Oak;
		LABEL  COL1  ="H² Oakey";
		FORMAT COL1 8.2;
		RUN;

%MEND H2Oakey;

