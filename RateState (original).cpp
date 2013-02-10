/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/29 22:21:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <iostream>

/* Header files with a description of contents used */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
 the mathematical problem description given above.
 
 Ith(v,i) references the ith component of the vector v, where i is in
 the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
 using the N_VIth macro in nvector.h. N_VIth numbers the components of
 a vector starting from 0.
 
 IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
 i and j are in the range [1..NEQ]. The IJth macro is defined using the
 DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
 dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   3                /* number of equations  */
#define X0    RCONST(-10.0)
#define V0    RCONST(1.0)	// V0 must be greater than 0
#define H0    RCONST(1.0)	// H0 must be greater than 0
#define RTOL  RCONST(1.0e-7)	// scalar relative tolerance
#define ATOLX RCONST(1.0e-12)	// vector absolute tolerance components
#define ATOLV RCONST(1.0e-12)
#define ATOLH RCONST(1.0e-12)
#define T0    RCONST(0.0)		// initial time
#define TSTEP RCONST(0.25)		// output time factor
#define TMAX  500

#define A_PARAM	0.0625
#define B_PARAM	0.125
#define R_PARAM	1e-5
#define K_PARAM	18

//#define A_PARAM(i)		params[i*4+0]
//#define B_PARAM(i)		params[i*4+1]
//#define R_PARAM(i)		params[i*4+2]
//#define K_PARAM(i)		params[i*4+3]

#define NUM_BLOCKS	1

/* Functions Called by the Solver */

static double F(realtype v, realtype h);
static int func(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(FILE *fp, realtype t, realtype y1, realtype y2, realtype y3);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

double *params;		// size 4*num_blocks, every 4 entries are params [a, b, r, k]

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
	realtype reltol, t, tout;
	N_Vector y, abstol;
	void *cvode_mem;
	int flag, flagr;
	FILE *fp;
	
	y = abstol = NULL;
	cvode_mem = NULL;
	
	// Set up params
	params = new double[NUM_BLOCKS];
	//A_PARAM(0) = 0.0625;
	//B_PARAM(0) = 0.125;
	//R_PARAM(0) = 1e-5;
	//K_PARAM(0) = 20;
	
	/* Create serial vector of length NEQ for I.C. and abstol */
	y = N_VNew_Serial(NEQ);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
	abstol = N_VNew_Serial(NEQ); 
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);
	
	/* Initialize y */
	Ith(y,1) = X0;
	Ith(y,2) = V0;
	Ith(y,3) = H0;
	
	/* Set the scalar relative tolerance */
	reltol = RTOL;
	/* Set the vector absolute tolerance */
	Ith(abstol,1) = ATOLX;
	Ith(abstol,2) = ATOLV;
	Ith(abstol,3) = ATOLH;
	
	/* Call CVodeCreate to create the solver memory and specify the 
	 * Backward Differentiation Formula and the use of a Newton iteration */
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function in y'=f(t,y), the inital time T0, and
	 * the initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, func, T0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);
	
	/* Call CVDense to specify the CVDENSE dense linear solver */
	//flag = CVSpbcg(cvode_mem, PREC_NONE, 0);
	//if (check_flag(&flag, "CVSpbcg", 1)) return(1);
	//flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
	//if (check_flag(&flag, "CVSpgmr", 1)) return(1);
	flag = CVDense(cvode_mem, NEQ);
	if (check_flag(&flag, "CVDense", 1)) return(1);
	
	CVodeSetMaxNumSteps(cvode_mem, 100000);
	/* Set the Jacobian routine to Jac (user-supplied) */
	flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
	if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
	
	/* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
	
	fp = fopen("out.txt", "w");
	fprintf(fp, "t x v theta F\n");
	
	tout = T0+TSTEP;
	while(t < TMAX) {
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		PrintOutput(fp, t, Ith(y,1), Ith(y,2), Ith(y,3));
		
		if (check_flag(&flag, "CVode", 1)) {
			std::cerr << Ith(y,1) << " " << Ith(y,2) << " " << Ith(y,3) << std::endl;
			break;
		}
		if (flag == CV_SUCCESS) {
			tout += TSTEP;
		}
	}
	
	fclose(fp);
	/* Print some final statistics */
	PrintFinalStats(cvode_mem);
	
	/* Free y and abstol vectors */
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(abstol);
	
	/* Free integrator memory */
	CVodeFree(&cvode_mem);
	
	return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

static double F(realtype v, realtype h) {
	return 1 + A_PARAM*log(v) + B_PARAM*log(h);
}

static int func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype x, v, h, xd, vd, hd, fval;
	
	x = Ith(y,1); v = Ith(y,2); h = Ith(y,3);
	
	if (v <= 0 || h <= 0) return 1;
	
	fval = F(v,h);
	
	xd = Ith(ydot,1) = v;
	vd = Ith(ydot,2) = (t-x-K_PARAM*fval)/R_PARAM;
	hd = Ith(ydot,3) = RCONST(1) - h*v;
	
	return 0;
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	realtype x, v, h;
	
	x = Ith(y,1); v = Ith(y,2); h = Ith(y,3);
	
	if (v <= 0 || h <= 0) return 1;
	
	IJth(J,1,1) = RCONST(0);
	IJth(J,1,2) = RCONST(1);
	IJth(J,1,3) = RCONST(0);
	IJth(J,2,1) = -RCONST(1)/R_PARAM;
	IJth(J,2,2) = -K_PARAM*A_PARAM/(R_PARAM*v);
	IJth(J,2,3) = -K_PARAM*B_PARAM/(R_PARAM*h);
	IJth(J,3,1) = RCONST(0); 
	IJth(J,3,2) = -h;
	IJth(J,3,3) = -v;
	
	return 0;
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(FILE *fp, realtype t, realtype y1, realtype y2, realtype y3)
{
	fprintf(fp, "%0.7e %14.6e %14.6e %14.6e %14.6e\n", t, y1, y2, y3, F(y2, y3));
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
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }
	
	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
					funcname, *errflag);
			return(1); }}
	
	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }
	
	return(0);
}
