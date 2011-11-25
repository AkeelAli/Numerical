#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 1 point gauss legendre for question 4 (c) */
double GaussLegendreIntegral(double start, double end, int N_start, int N_end){
	double integral = 0;
	
	double segments = N_end - N_start + 1;
	double interval = (end - start) / segments;
	double argument_i = interval;
	double argument_o = interval / 2 + start;

	int i;
	for (i = 0; i < segments; i++){
		integral += interval * log(argument_i*i + argument_o);
	}
	
	return integral;	
}  

int main(void) {
	double integral, error;
	
	int N;
	int i;
	
	/* Question 4 (a) */
	/*
	for (N = 1; N <= 20; N++){
		integral = 0;
		error = 0;
		
		for(i = 0; i < N; i++) {
			integral += 1./N * sin((2*i+1)/(2.*N));
		}
		
		error = fabs(1 - cos(1) - integral);	
		printf("%d,%.10f,%.10f\n",N,integral,error);
	}*/
	
	/* Question 4 (b) */
	/*
	for (N = 10; N <= 200; N+=10){
		integral = 0;
		error = 0;
		
		for(i = 0; i < N; i++) {
			integral += 1./N * log((2*i+1)/(2.*N));
		}
		
		error = fabs(-1 - integral);	
		printf("%d,%.10f,%.10f\n",N,integral,error);
	}*/
	
	/* Quesiton 4 (c) */
	int N_div = 6;
	double divider = 0.1;
	
	integral = 0;
	error = 0;
	for (N_div = 6; N_div < 10; N_div++) {
		for (divider = 0.05; divider < 0.65; divider+=0.05) {
			integral = GaussLegendreIntegral(0, divider, 1, N_div) 
				+ GaussLegendreIntegral(divider, 1, N_div + 1, 10);
			
			error = fabs(-1 - integral);	
			printf("%d,%.2f,%.10f,%.10f\n",N_div,divider,integral,error);
		}
	}	

	
	return 0;
}

