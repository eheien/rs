#include "RateState.h"
#include "RateStateSimWindow.h"
#include <iostream>

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

//#define PARAM_SWEEP

int main(int argc, char **argv)
{
#ifdef PARAM_SWEEP
	std::vector<std::vector<realtype> >	results;
	RSParams				params(NBLOCKS, NEQ, NPARAMS, 0.01, 500);
	unsigned int			i, npoints;
	int						res, rlog;
	double					xi, vi, hi;
	double					param_a, param_b, param_k, param_r;
	double					a_step = 0.01, b_step = 0.01;
	
	for (param_a=a_step;param_a<0.4+a_step;param_a+=a_step) {
		for (param_b=param_a/0.52553300962692;param_b<1.0+b_step;param_b+=b_step) {
			for (rlog=-5;rlog<-4;rlog++) {
				for (param_k=20;param_k<24;param_k+=4) {
					param_r = pow(10, rlog);
					for (i=0;i<params.num_blocks();++i) {
						results.clear();
						params.param(i, A_PARAM) = RCONST(param_a);
						params.param(i, B_PARAM) = RCONST(param_b);
						params.param(i, K_PARAM) = RCONST(param_k);
						params.param(i, R_PARAM) = RCONST(param_r);
						params.init_val(i, EQ_X) = RCONST(-15.0);
						params.init_val(i, EQ_V) = RCONST(1.0);
						params.init_val(i, EQ_H) = RCONST(1.0);
						
					}
					res = run_rate_state_sim(results, params);
					std::cerr << res << " " << param_a << " " << param_b << " " << param_k << " " << rlog << std::endl;
				}
			}
		}
	}
	
#else
    QApplication			a(argc, argv);
	RateStateSimWindow		rs_win;
	
	rs_win.show();
	
	/*
	results.clear();
	run_rate_state_sim(results, params);
	
	PrintOutput("out.txt", results, params);
	 */	
	return a.exec(); 
#endif
}
