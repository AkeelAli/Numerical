#include <stdio.h>
#include <stdlib.h>

#define NR_LIMIT 0.0000000001 /* run algorithm until |f/f(0)| < 10^-6 */
#define SUB_LIMIT 0.0000000001
#define INIT_GUESS 0.0 /* starting with zero flux */
#define PI 3.141592654

#define SUB

/* global table holding the BH function from table 1 in the assignment */
double BH_table[][2] = { {0.0, 0.0}, {0.2, 14.7}, {0.4, 36.5}, {0.6, 71.7},
	{0.8, 121.4}, {1.0, 197.4}, {1.1, 256.2}, {1.2, 348.7}, {1.3, 540.6},
	{1.4, 1062.8}, {1.5, 2318.0}, {1.6, 4781.9}, {1.7, 8687.4}, {1.8, 13924.3},
	{1.9, 22650.2} };
	
#define B_SIZE sizeof(BH_table)/sizeof(BH_table[0])

double solveNR(double init_flux);
double solveSub(double init_flux);

int main(void){
#ifdef PRINT_NR
	printf("%4s,%15s,%15s,%15s,%15s\n","itr","flux","f","f'","residue");
	solveNR(INIT_GUESS);
#endif
	
#ifdef SUB
	printf("%4s,%15s,%15s\n","itr","flux","f");
	solveSub(INIT_GUESS);
#endif

	return 0;
}

double getH(double flux);

/**
	Solves the non linear equation of the form f(flux) = 0
	using successive substitution
*/
double solveSub(double init_flux){
	double flux = init_flux;
	double expr = 125000000 / PI; /* expression appearing twice in eqn */
	
	double f, residue;
	int iteration = 0;
	
	do{
		f = ((0.3 * getH(flux)) + (expr * flux) - 8000);

		printf("%4d, %15e, %15e\n",iteration,flux,f);

		iteration++;
		flux = flux - f;
	
		residue = f < 0 ? -f:f;
	}while(residue > SUB_LIMIT);
	
	return flux;
}

double slope(double flux);

/**
	solves the non-linear eq. of question 2 using Newton Raphson 
	returns the final flux value after reaching NR_LIMIT
*/
double solveNR(double init_flux){
	double flux = init_flux;
	double expr = 125000000 / PI; /* expression appearing twice in eqn */

	double fzero = -8000; /* f(0) */
	
	double f, fprime, residue;
	int iteration = 0;
	do{	
		/* iteration k */
		f = ((0.3 * getH(flux)) + (expr * flux) - 8000);
		fprime = 3000 * slope(flux) + expr;
#ifdef PRINT_NR
		printf("%4d,%15e, %15e, %15e, %15e\n",iteration,flux,f,fprime,residue);
#endif
		/* iteration k+1 */
		iteration++; 				/* increment iteration for next loop */
		flux = flux - (f / fprime); 	/* computer flux(k+1) here */
		residue = (residue = f / fzero) < 0 ? -residue:residue;
	}while(residue > NR_LIMIT);
	
	return flux;
}

/**
	derives the H corresponding to the given flux using 
	a piecewise linear interpolation from the BH_table
	returns the H
*/
double getH(double flux){
	double B = 10000 * flux;
	
	/* B outside range to the left */
	if (B < BH_table[0][0])
		return slope(flux) * (B - BH_table[0][0]) + BH_table[0][1];
	
	/* search for 2-points range containing B */
	int i;
	for (i = 0; i < B_SIZE - 1; i++){
		if (B >= BH_table[i][0] && B <= BH_table[i+1][0])
			return slope(flux) * (B - BH_table[i][0]) + BH_table[i][1];
	}
	
	/* B outside BH_table to the right */ 
	return slope(flux) * (B - BH_table[B_SIZE - 2][0]) + BH_table[B_SIZE - 2][1];
}

/**
	returns the slope of the nearest or containing range
	in the BH relation
*/
double slope(double flux){
	double B = 10000 * flux;
	
	if (B < BH_table[0][0])
		return (BH_table[1][1] - BH_table[0][1]) / 
			(BH_table[1][0] - BH_table[0][0]);
	
	int i;
	for (i = 0; i < B_SIZE - 1; i++){
		if (B >= BH_table[i][0] && B <= BH_table[i+1][0])
			return (BH_table[i+1][1] - BH_table[i][1]) / (
				BH_table[i+1][0] - BH_table[i][0]);
	}
	
	return (BH_table[B_SIZE - 1][1] - BH_table[B_SIZE - 2][1]) / 
		((BH_table[B_SIZE - 1][0] - BH_table[B_SIZE - 2][0]));
	
}	