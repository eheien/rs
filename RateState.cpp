#include "RateState.h"

realtype *interaction;

int run_rate_state_sim(std::vector<std::vector<realtype> > &results, RSParams &params) {
	realtype		long_term_reltol, event_reltol, t, tout, tbase=0;
	N_Vector		y, long_term_abstol, event_abstol;
	unsigned int	i, n;
	int				flag, err_code;
	void			*long_term_cvode, *event_cvode, *current_cvode;
	
	// Create serial vector of length NEQ for I.C. and abstol
	y = N_VNew_Serial(params.num_eqs()*params.num_blocks());
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
	long_term_abstol = N_VNew_Serial(params.num_eqs()*params.num_blocks()); 
	if (check_flag((void *)long_term_abstol, "N_VNew_Serial", 0)) return(1);
	event_abstol = N_VNew_Serial(params.num_eqs()*params.num_blocks()); 
	if (check_flag((void *)event_abstol, "N_VNew_Serial", 0)) return(1);
	
	// Initialize y
	for (i=0;i<params.num_blocks();++i) {
		NV_Ith_S(y,i*params.num_eqs()+EQ_X) = params.init_val(i, EQ_X);
		NV_Ith_S(y,i*params.num_eqs()+EQ_V) = params.init_val(i, EQ_V);
		NV_Ith_S(y,i*params.num_eqs()+EQ_H) = params.init_val(i, EQ_H);
	}
	
	/* Initialize interactions */
	/*interaction = new realtype[NBLOCKS*NBLOCKS];
	double int_level = 1e-2;
	double dropoff = 1.1;
	for (i=0;i<NBLOCKS;++i) {
		for (n=0;n<NBLOCKS;++n) {
			interaction[i*NBLOCKS+n] = (i==n?(1.0-int_level):int_level);
		}
	}*/
	
	/* Set the scalar relative tolerance */
	long_term_reltol = RCONST(1.0e-12);
	event_reltol = RCONST(1.0e-12);
	/* Set the vector absolute tolerance */
	for (i=0;i<params.num_blocks();++i) {
		Xth(long_term_abstol,i) = RCONST(1.0e-12);
		Vth(long_term_abstol,i) = RCONST(1.0e-12);
		Hth(long_term_abstol,i) = RCONST(1.0e-12);
		Xth(event_abstol,i) = RCONST(1.0e-12);
		Vth(event_abstol,i) = RCONST(1.0e-12);
		Hth(event_abstol,i) = RCONST(1.0e-12);
	}
	
	/* Call CVodeCreate to create the solver memory and specify the 
	 * Backward Differentiation Formula and the use of a Newton iteration */
	long_term_cvode = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)long_term_cvode, "CVodeCreate", 0)) return(1);
	event_cvode = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)event_cvode, "CVodeCreate", 0)) return(1);
	
	// Turn off error messages
	//CVodeSetErrFile(long_term_cvode, NULL);
	//CVodeSetErrFile(event_cvode, NULL);
	
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function in y'=f(t,y), the inital time T0, and
	 * the initial dependent variable vector y. */
	flag = CVodeInit(long_term_cvode, func, T0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	flag = CVodeInit(event_cvode, func, T0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(long_term_cvode, long_term_reltol, long_term_abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);
	flag = CVodeSVtolerances(event_cvode, event_reltol, event_abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);
	
	/* Set the root finding function */
	//flag = CVodeRootInit(long_term_cvode, params.num_blocks(), vel_switch_finder);
	//flag = CVodeRootInit(event_cvode, params.num_blocks(), vel_switch_finder);
	//if (check_flag(&flag, "CVodeRootInit", 1)) return(1);
	
	/* Call CVDense to specify the CVDENSE dense linear solver */
	//flag = CVSpbcg(cvode_mem, PREC_NONE, 0);
	//if (check_flag(&flag, "CVSpbcg", 1)) return(1);
	//flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
	//if (check_flag(&flag, "CVSpgmr", 1)) return(1);
	flag = CVDense(long_term_cvode, params.num_eqs()*params.num_blocks());
	if (check_flag(&flag, "CVDense", 1)) return(1);
	flag = CVDense(event_cvode, params.num_eqs()*params.num_blocks());
	if (check_flag(&flag, "CVDense", 1)) return(1);
	
	flag = CVodeSetUserData(long_term_cvode, &params);
	if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);
	flag = CVodeSetUserData(event_cvode, &params);
	if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);
	
	CVodeSetMaxNumSteps(long_term_cvode, 100000);
	CVodeSetMaxNumSteps(event_cvode, 100000);
	
	/* Set the Jacobian routine to Jac (user-supplied) */
	flag = CVDlsSetDenseJacFn(long_term_cvode, Jac);
	if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
	flag = CVDlsSetDenseJacFn(event_cvode, Jac);
	if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
	
	/* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
	
	tout = T0+params.time_step();
	err_code = 0;
	int mode = 0;
	int	num_res_vals = params.num_eqs()*params.num_blocks();
	while(t+tbase < params.end_time()) {
		switch (mode) {
			case 0:
				current_cvode = long_term_cvode;
				break;
			case 1:
				current_cvode = event_cvode;
				break;
		}
		flag = CVode(current_cvode, tout, y, &t, CV_NORMAL);
		record_results(results, y, t, tbase, num_res_vals);
		
		if (check_flag(&flag, "CVode", 1)) {
			err_code = flag;
			break;
		}
		if (flag == CV_ROOT_RETURN) {
			int flagr;
			int rootsfound[params.num_blocks()];
			flagr = CVodeGetRootInfo(current_cvode, rootsfound);
			if (check_flag(&flagr, "CVodeGetRootInfo", 1)) return(1);
			//mode = !mode;
			//std::cerr << rootsfound[0] << std::endl;
		}
		if (flag == CV_SUCCESS) tout += params.time_step();
		if (t > 10) {
			t -= 10;
			tout -= 10;
			Xth(y,0) -= 10;
			tbase += 10;
			CVodeReInit(current_cvode, t, y);
		}
	}
	
	std::cerr << err_code <<
		" X:" << Xth(y,0) <<
		" V:" << Vth(y,0) <<
		" H:" << Hth(y,0) <<
		" F:" << F(0,Vth(y,0),Hth(y,0),params) <<
		" dV:" << (Xth(y,0)-t-params.param(0, K_PARAM)*F(0,Vth(y,0),Hth(y,0),params))/params.param(0, R_PARAM) <<
		" dH:" << -Hth(y,0)*Vth(y,0)*log(Hth(y,0)*Vth(y,0)) << std::endl;
	/* Print some final statistics */
	//PrintFinalStats(cvode_mem);
	
	/* Free y and abstol vectors */
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(long_term_abstol);
	N_VDestroy_Serial(event_abstol);
	
	/* Free integrator memory */
	CVodeFree(&long_term_cvode);
	
	return err_code;
}

void record_results(std::vector<std::vector<realtype> > &results, N_Vector y, realtype t, realtype tbase, int num_res_vals) {
	std::vector<realtype>	cur_step_results;
	int						i;
	
	cur_step_results = std::vector<realtype>(num_res_vals+1);
	cur_step_results[0] = t+tbase;
	for (i=0;i<num_res_vals;++i) {
		if (i==0) cur_step_results[i+1] = NV_Ith_S(y, i)+tbase;
		else cur_step_results[i+1] = NV_Ith_S(y, i);
	}
	
	results.push_back(cur_step_results);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

realtype F(int bnum, realtype v, realtype h, RSParams &params) {
	return 1 + params.param(bnum, A_PARAM)*log(v) + params.param(bnum, B_PARAM)*log(h);
}

realtype FdV(int bnum, realtype v, realtype h, RSParams &params) {
	return params.param(bnum, A_PARAM)/v;
}

realtype FdH(int bnum, realtype v, realtype h, RSParams &params) {
	return params.param(bnum, B_PARAM)/h;
}

//#define SLOWNESS_LAW

static int func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype	fval;
	int			i;
	RSParams	*params = (RSParams*)(user_data);

	for (i=0;i<params->num_blocks();++i) {
		if (Vth(y,i) <= 0 || Hth(y,i) <= 0) return 1;
	}
	
	for (i=0;i<params->num_blocks();++i) {
		fval = params->param(i, K_PARAM)*F(i,Vth(y,i),Hth(y,i),*params);
		
		Xth(ydot,i) = Vth(y,i);
		Vth(ydot,i) = (t-Xth(y,i)-fval)/params->param(i, R_PARAM) - pow(params->param(i, R_PARAM), 1.5)*params->param(i, W_PARAM)*Vth(y,i);
#ifdef SLOWNESS_LAW
		Hth(ydot,i) = RCONST(1) - Hth(y,i)*Vth(y,i);
#else
		Hth(ydot,i) = -Hth(y,i)*Vth(y,i)*log(Hth(y,i)*Vth(y,i));
#endif
	}
	
	return 0;
}

int vel_switch_finder(realtype t, N_Vector y, realtype *gout, void *user_data)
{
	RSParams	*params = (RSParams*)(user_data);
	int			i;
	
	// Check whether each block has hit V=0.1, in which case we move to fast mode
	for (i=0;i<params->num_blocks();++i) {
		if (Vth(y,i) < 0) return -1;
		gout[i] = log10(Vth(y,i))+1;
	}
	
	return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

int Jac(long int N, realtype t,
		N_Vector y, N_Vector fy, DlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	RSParams	*params = (RSParams*)(user_data);
	
	realtype x, v, h;
	
	x = Xth(y,0); v = Vth(y,0); h = Hth(y,0);
	
	if (v <= 0 || h <= 0) return 1;
	
	XXth(J,0,0) = RCONST(0);
	XVth(J,0,0) = RCONST(1);
	XHth(J,0,0) = RCONST(0);
	VXth(J,0,0) = -RCONST(1)/params->param(0, R_PARAM);
	VVth(J,0,0) = -params->param(0, K_PARAM)*params->param(0, A_PARAM)/(params->param(0, R_PARAM)*v);
	VHth(J,0,0) = -params->param(0, K_PARAM)*params->param(0, B_PARAM)/(params->param(0, R_PARAM)*h);
#ifdef SLOWNESS_LAW
	HXth(J,0,0) = RCONST(0); 
	HVth(J,0,0) = -h;
	HHth(J,0,0) = -v;
#else
	HXth(J,0,0) = RCONST(0); 
	HVth(J,0,0) = -h*(log(h*v)+1);
	HHth(J,0,0) = -v*(log(v*h)+1);
#endif
	
	return 0;
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

void PrintOutput(char *filename, std::vector<realtype*> &results, RSParams &params) {
	realtype	xi, vi, hi, *data;
	int			i;
	FILE		*fp;
	std::vector<realtype*>::const_iterator		it;
	
	fp = fopen(filename, "w");
	fprintf(fp, "t ");
	for (i=0;i<params.num_blocks();++i) fprintf(fp, "x%d v%d h%d F%d ", i, i, i, i);
	fprintf(fp, "\n");
	
	for (it=results.begin();it!=results.end();++it) {
		data = *it;
		fprintf(fp, "%0.7e ", data[0]);
		for (i=0;i<params.num_blocks();++i) {
			xi = data[i*params.num_eqs()+1];
			vi = data[i*params.num_eqs()+2];
			hi = data[i*params.num_eqs()+3];
			fprintf(fp, "%14.6e %14.6e %14.6e %14.6e ", xi, vi, hi, F(i, vi, hi, params));
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
	long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
	int flag;
	
	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	check_flag(&flag, "CVodeGetNumRhsEvals", 1);
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	check_flag(&flag, "CVodeGetNumErrTestFails", 1);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
	
	flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
	check_flag(&flag, "CVDlsGetNumJacEvals", 1);
	flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
	check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
	
	flag = CVodeGetNumGEvals(cvode_mem, &nge);
	check_flag(&flag, "CVodeGetNumGEvals", 1);
	
	printf("\nFinal Statistics:\n");
	printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
		   nst, nfe, nsetups, nfeLS, nje);
	printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
		   nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
	int *errflag;
	
	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		//fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
		return(1); }
	
	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			//fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
			return(1); }}
	
	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		//fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
		return(1); }
	
	return(0);
}
