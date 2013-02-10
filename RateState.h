#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>

// Format of result vector: [x0, v0, h0, x1, v1, h1, ..., xn, vn, hn]

#define Xth(y,block)	NV_Ith_S(y,block*NEQ+0)
#define Vth(y,block)	NV_Ith_S(y,block*NEQ+1)
#define Hth(y,block)	NV_Ith_S(y,block*NEQ+2)

#define XXth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+0,bj*NEQ+0)
#define XVth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+0,bj*NEQ+1)
#define XHth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+0,bj*NEQ+2)
#define VXth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+1,bj*NEQ+0)
#define VVth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+1,bj*NEQ+1)
#define VHth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+1,bj*NEQ+2)
#define HXth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+2,bj*NEQ+0)
#define HVth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+2,bj*NEQ+1)
#define HHth(A,bi,bj)	DENSE_ELEM(A,bi*NEQ+2,bj*NEQ+2)


/* Problem Constants */

#define NBLOCKS	1
#define T0    RCONST(0.0)		// initial time

#define NPARAMS			5
#define A_PARAM			0
#define B_PARAM			1
#define R_PARAM			2
#define K_PARAM			3
#define W_PARAM			4

#define NEQ		3         /* number of equations per block  */
#define EQ_X			0
#define EQ_V			1
#define EQ_H			2

class RSParams {
private:
	int			nblocks;
	int			neqs;
	int			nparams;
	realtype	*init_vals;
	realtype	*params;
	
	realtype	tmax, tstep;
	
public:
	RSParams(int numblocks, int numeqs, int num_params, realtype time_step, realtype end_time) :
			nblocks(numblocks), neqs(numeqs), nparams(num_params), tstep(time_step), tmax(end_time) {
		init_vals = new realtype[nblocks*neqs];
		params = new realtype[nblocks*nparams];
	};
	
	~RSParams(void) { delete init_vals; delete params; };
	
	int	num_blocks(void) const { return nblocks; };
	int	num_eqs(void) const { return neqs; };
	int	num_params(void) const { return nparams; };
	realtype time_step(void) const { return tstep; };
	realtype end_time(void) const { return tmax; };
	
	realtype& init_val(int bnum, int eqnum) {
		return init_vals[bnum*neqs+eqnum];
	}
	
	realtype& param(int bnum, int pnum) {
		return params[bnum*nparams+pnum];
	}
};

/* Functions Called by the Solver */
void record_results(std::vector<std::vector<realtype> > &results, N_Vector y, realtype t, realtype tbase, int num_res_vals);

realtype F(int bnum, realtype v, realtype h, RSParams &params);
realtype FdH(int bnum, realtype v, realtype h, RSParams &params);
realtype FdV(int bnum, realtype v, realtype h, RSParams &params);

static int func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int vel_switch_finder(realtype t, N_Vector y, realtype *gout, void *user_data);

int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

void PrintOutput(char *filename, std::vector<realtype*> &results, RSParams &params);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

extern realtype *interaction;

int run_rate_state_sim(std::vector<std::vector<realtype> > &results, RSParams &params);

