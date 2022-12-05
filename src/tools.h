/*
 * tools.h
 *
 *  Created on: 24 nov. 2022
 *      Author: Kévin Carrier
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <complex.h>



/*
 * @brief Compute the Euclidean norm by ignoring the punctured positions
 */
double norm(int *v, unsigned int q, int length){
    double res = 0.0, tmp;
    for (int i = 0 ; i < length ; ++i){
    	if (v[i] >= 0){
    		tmp = (v[i] > q/2.0 ? (double)(v[i]) - (double)q : (double)(v[i]));
    		res += ( tmp * tmp );
    	}
    }
    return sqrt(res);
}


/*
 * @brief Compute the Euclidean distance by ignoring the punctured positions
 */
double dist(int *u, int *v, unsigned int q, int length){
    double res = 0.0, tmp;
    for (int i = 0 ; i < length ; ++i){
    	if ((v[i] >= 0) && (u[i] >= 0)){
    		tmp = fabs((double)(v[i]) - (double)(u[i]));
    		if (tmp > q/2.0){
    			tmp -= (double)q;
    		}
    		res += ( tmp * tmp );
    	}
    }
    return sqrt(res);
}



int partition(double *arr, int l, int u) {
	int i, j;
	double v, tmp;
	v = arr[l];
	i = l;
	j = u + 1;
	do {
		do {
			++i;
		} while ( (arr[i] < v) && (i <= u) );
		do {
			--j;
		} while (v < arr[j]);
		if (i < j) {
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
		}
	} while (i < j);
	arr[l] = arr[j];
	arr[j] = v;
	return j;
}

void quick_sort(double *arr, int l, int u){
	int j;
	if (l < u) {
		j = partition(arr, l, u);
		quick_sort(arr, l, j-1);
		quick_sort(arr, j+1, u);
	}
}


void fft(complex double * f, int n){
	if(n != 1){
		int n_over_2 = n/2;
		complex double * tmp1 = (complex double*)malloc(sizeof(complex double)*n_over_2);
		complex double * tmp2 = (complex double*)malloc(sizeof(complex double)*n_over_2);
		for (int k = 0 ; k < n_over_2 ; ++k){
			tmp1[k] = f[2*k];
			tmp2[k] = f[2*k + 1];
		}
		fft(tmp1, n_over_2);
		fft(tmp2, n_over_2);

		complex double W=cexp(-2.0*I*M_PI/n), W_tmp=1;
		for( int k = 0; k<n_over_2 ; ++k){
			f[k+n_over_2] = tmp1[k] - W_tmp * tmp2[k];
			f[k] = tmp1[k] + W_tmp * tmp2[k];
			W_tmp *= W;
		}
		free(tmp1);
		free(tmp2);
	}
}


void ifft_(complex double * f, int n){
	if(n != 1){
		int n_over_2 = n/2;
		complex double * tmp1 = (complex double*)malloc(sizeof(complex double)*n_over_2);
		complex double * tmp2 = (complex double*)malloc(sizeof(complex double)*n_over_2);
		for (int k = 0 ; k < n_over_2 ; ++k){
			tmp1[k] = f[2*k];
			tmp2[k] = f[2*k + 1];
		}
		ifft_(tmp1, n_over_2);
		ifft_(tmp2, n_over_2);

		complex double W=cexp(2.0*I*M_PI/n), W_tmp=1;
		for( int k = 0; k<n_over_2 ; ++k){
			f[k+n_over_2] = tmp1[k] - W_tmp * tmp2[k];
			f[k] = tmp1[k] + W_tmp * tmp2[k];
			W_tmp *= W;
		}
		free(tmp1);
		free(tmp2);
	}
}

void ifft(complex double * f, int n){
	ifft_(f, n);
	for (int k = 0 ; k < n ; ++k){
		f[k] /= (double)n;
	}
}

#endif /* TOOLS_H_ */
