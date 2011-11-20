#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NR_LIMIT pow(10,-10) 
#define INIT_GUESS 0.0
#define PI 3.141592654

#define MICRO pow(10,-6)
#define MILLI pow(10,-3)
#define EXPV1V2 (  exp( (v[0] - v[1])/(25*MILLI) )  )
#define EXPV2 (  exp( v[1]/(25*MILLI) )  )

#define PRINT_NR


void solveNR(double init_v1, double init_v2);

int main(void){
#ifdef PRINT_NR
	printf("%4s,%13s,%13s,%13s,%13s,%10s\n","itr","V_DA","V_DB","f[0]","f[1]","residue");
	solveNR(INIT_GUESS, INIT_GUESS);
#endif

	return 0;
}


/**
	solves the non-linear eq. of question 2 using Newton Raphson 
	returns the final flux value after reaching NR_LIMIT
*/
void solveNR(double init_v1, double init_v2) {
	double v[2] = {init_v1, init_v2};
		
	double f[2], fprime[2][2], fprime_inv[2][2], fprime_det, residue;
	int iteration = 0;
	do{	
		/* iteration k */
		f[0] = 2*MILLI*v[0] + 0.6*MICRO*(EXPV1V2-1) - 440*MICRO;
		f[1] = 0.6*MICRO*(EXPV1V2-1) - 1.2*MICRO*(EXPV2-1);
		fprime[0][0] = 2*MILLI + 24*MICRO*EXPV1V2;
		fprime[0][1] = -24*MICRO*EXPV1V2;
		fprime[1][0] = 24*MICRO*EXPV1V2;
		fprime[1][1] = -24*MICRO*EXPV1V2 -48*MICRO*EXPV2;
		fprime_det = fprime[0][0]*fprime[1][1] - fprime[0][1]*fprime[1][0]; /* ad - bc */
		fprime_inv[0][0] = 1/fprime_det * fprime[1][1];
		fprime_inv[0][1] = 1/fprime_det * -fprime[0][1];
		fprime_inv[1][0] = 1/fprime_det * -fprime[1][0];
		fprime_inv[1][1] = 1/fprime_det * fprime[0][0];
		
		residue = (residue = (f[0]+f[1])/2) < 0 ? -residue:residue;
		printf("%4d,%13f,%13f,%13e,%13e,%10e\n",iteration,v[0]-v[1],v[1],f[0],f[1],residue);

		/* iteration k+1 */
		iteration++; 				/* increment iteration for next loop */
		v[0] = v[0] - fprime_inv[0][0]*f[0] - fprime_inv[0][1]*f[1];
		v[1] = v[1] - fprime_inv[1][0]*f[0] - fprime_inv[1][1]*f[1];

	}while(residue > NR_LIMIT);
	
}