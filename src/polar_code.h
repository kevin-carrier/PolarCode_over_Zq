/*
 * polar_code.h
 *
 *  Created on: 24 nov. 2022
 *      Author: Kévin Carrier
 *
 * Description: Non Binary Polar Code construction based on:
 *                --> E. Sasoglu, E. Telatar, and E. Arikan, “Polarization for arbitrary discrete memoryless channels,” in IEEE Information Theory Workshop, 2009.
 *                --> M.-C. Chiu, “Non-binary polar codes with channel symbol permutations,” in Int. Symp. on Info. Theory and its Applications (ISITA), 2014.
 */

#ifndef POLAR_CODE_H_
#define POLAR_CODE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ring.h"
#include "gaussian_channel.h"
#include "tools.h"

typedef struct PolarCode {
	ring_t *ring;
	channel_t *channel;
	unsigned int m;
	unsigned int n; // length of the NON punctured code. The real length is n - nb_punc
	unsigned int k;
	unsigned int nb_punc;
	unsigned int **coeffs;
	double *frozen_probs;
	int *frozen_bits;
} polar_t;


void free_polar(polar_t *code){
	free_ring(code->ring);
	free_channel(code->channel);
	free(code->frozen_probs);
	free(code->frozen_bits);
	for( int i = 0 ; i < code->m ; ++i){
		free(code->coeffs[i]);
	}
	free(code->coeffs);
}


/*
 * @brief Save a description of the polar code in a file.
 * @param fd output file descriptor.
 * @param code the polar code to output
 */
void print_polar(FILE *fd, polar_t *code){
	fprintf(fd, "s=%d\n", code->ring->s);
	fprintf(fd, "m=%d\n", code->m);
	fprintf(fd, "k=%d\n", code->k);
	fprintf(fd, "nb_punc=%d\n\n", code->nb_punc);

	fprintf(fd, "coeffs=\n");
	for (int i = 0 ; i < code->m - 1 ; ++i){
		for (int j = 0 ; j < (unsigned int)(1 << i) - 1 ; ++j){
			fprintf(fd, "%d\t", code->coeffs[i][j]);
		}
		fprintf(fd, "%d\n", code->coeffs[i][(unsigned int)(1 << i) - 1]);
	}
	for (int j = 0 ; j < (code->n / 2) - 1 ; ++j){
		fprintf(fd, "%d\t", code->coeffs[code->m - 1][j]);
	}
	fprintf(fd, "%d\n\n", code->coeffs[code->m - 1][(code->n / 2) - 1]);

	fprintf(fd, "frozen_probs=");
	for (int i = 0 ; i < code->n - 1 ; ++i){
		fprintf(fd, "%.12f\t", code->frozen_probs[i]);
	}
	fprintf(fd, "%.12f\n", code->frozen_probs[code->n - 1]);

	fprintf(fd, "frozen_bits=");
	for (int i = 0 ; i < code->n - 1 ; ++i){
		fprintf(fd, "%d\t", code->frozen_bits[i]);
	}
	fprintf(fd, "%d\n", code->frozen_bits[code->n - 1]);
}

polar_t *copy_polar(polar_t *code_src){
	polar_t *code = (polar_t *)malloc(sizeof(polar_t));

	code->ring = gen_ring(code_src->ring->s);
	code->channel = gen_channel(code_src->ring->q, code_src->channel->std_dev);
	code->m = code_src->m;
	code->n = code_src->n;
	code->k = code_src->k;
	code->nb_punc = code_src->nb_punc;

	code->coeffs = (unsigned int **)malloc(sizeof(unsigned int *) * code_src->m);
	for (int i = 0 ; i < code_src->m ; ++i){
		code->coeffs[i] = (unsigned int *)malloc(sizeof(unsigned int)*((int)(code_src->n >> 1)));
		for (int j = 0 ; j < (unsigned int)(1 << i) ; ++j){
			code->coeffs[i][j] = code_src->coeffs[i][j];
		}
	}

	code->frozen_probs = (double *)malloc(sizeof(double) * code_src->n);
	for (int i = 0 ; i < code_src->n ; ++i){
		code->frozen_probs[i] = code_src->frozen_probs[i];
	}

	code->frozen_bits = (int *)malloc(sizeof(int) * code_src->n);
	for (int i = 0 ; i < code_src->n ; ++i){
		code->frozen_bits[i] = code_src->frozen_bits[i];
	}

	return code;
}

/*
 * @brief Build a non binary polar code from a file descriptor.
 * @param fd must be a file output by print_polar whose contains the description of a non binary polar code.
 *
 * @return a polar code whose description is given in the file fd.
 */
polar_t *polar_from_file(FILE *fd){
	// see format of print_polar
	int line_size = 131072;
	rewind(fd);

	polar_t *code = (polar_t *)malloc(sizeof(polar_t));
	code->coeffs = NULL;
	code->frozen_bits = NULL;

	char str[line_size];
	char *ptr;

	while (fgets(str, line_size - 1, fd) != NULL){
		ptr = strtok(str, "=");

		if (strcmp(ptr, "s") == 0){
			ptr = strtok(NULL, "=");
			code->ring = gen_ring(atoi(ptr));
		} else if (strcmp(ptr, "m") == 0){
			ptr = strtok(NULL, "=");
			code->m = atoi(ptr);
			code->n = (int)(1 << code->m);
		} else if (strcmp(ptr, "nb_punc") == 0){
			ptr = strtok(NULL, "=");
			code->nb_punc = atoi(ptr);
		} else if (strcmp(ptr, "k") == 0){
			ptr = strtok(NULL, "=");
			code->k = atoi(ptr);
		} else if (strcmp(ptr, "coeffs") == 0){
			code->coeffs = (unsigned int **)malloc(sizeof(unsigned int *) * code->m);
			for (int i = 0 ; i < code->m ; ++i){
				code->coeffs[i] = (unsigned int *)malloc(sizeof(unsigned int)*((int)(code->n >> 1)));
				if (fgets(str, line_size - 1, fd) == NULL){
					fprintf(stderr, "fgets error: file hasn't the good format \n");
					exit(EXIT_FAILURE);
				}
				ptr = strtok(str, "\t");
				for (int j = 0 ; j < (int)(1 << i) ; ++j){
					code->coeffs[i][j] = atoi(ptr);
					ptr = strtok(NULL, "\t");
				}
			}
		} else if (strcmp(ptr, "frozen_probs") == 0){
			code->frozen_probs = (double *)malloc(sizeof(double) * code->n);
			for (int i = 0 ; i < code->n ; ++i){
				ptr = strtok(NULL, "\t");
				code->frozen_probs[i] = strtod(ptr, NULL);
			}
		} else if (strcmp(ptr, "frozen_bits") == 0){
			code->frozen_bits = (int *)malloc(sizeof(int) * code->n);
			for (int i = 0 ; i < code->n ; ++i){
				ptr = strtok(NULL, "\t");
				code->frozen_bits[i] = atoi(ptr);
			}
		}
	}
	double std_dev = pow(code->ring->q, 1.0 - (double)(code->k)/((double)(code->n - code->nb_punc))) / sqrt(2.0 * M_PI * M_E);
	code->channel = gen_channel(code->ring->q, std_dev);

	return code;
}




/*
 * @brief Build a (punctured) polar code over ZZ/(2^s)ZZ of length 2^m and dimension k.
 * @param s log_2 of the size of the ring
 * @param n length of the punctured polar code
 * @param k dimension of the plar code
 * @param nb_test number of simulation to generate the polarized channels
 * @param fd read the parameters from fd (can be partial)
 *
 * @return a random non binary polar code over ZZ/(2^s)ZZ of length 2^m and dimension k.
 */
polar_t *gen_polar(unsigned int s, unsigned int n, unsigned int k, double nb_tests, FILE *fd){
	int i, j, jj, u;

	polar_t *code;
	if (fd != NULL){
		code = polar_from_file(fd);
	} else {
		code = (polar_t *)malloc(sizeof(polar_t));
		code->ring = gen_ring(s);
		code->m = floor(log2(n-1)) + 1;
		code->n = (int)(1 << (code->m));
		code->nb_punc = code->n - n;
		code->k = k;
		double std_dev = pow(code->ring->q, 1.0 - (double)k/((double)(n))) / sqrt(2.0 * M_PI * M_E);
		code->channel = gen_channel(code->ring->q, std_dev);
		code->coeffs = NULL;
		code->frozen_bits = NULL;
	}

	if (code->coeffs == NULL){
		// random choice of the coefficient of the generator (they must be invertible <==> odd)
		code->coeffs = (unsigned int **)malloc(sizeof(unsigned int *) * code->m);
		for( i = 0 ; i < code->m ; ++i){
			code->coeffs[i] = (unsigned int *)malloc(sizeof(unsigned int)*((int)(code->n >> 1)));
			for( j = 0 ; j < (int)(1 << i) ; ++j){
				code->coeffs[i][j] = ( rand() % (code->ring->q / 2) )*2 + 1; // random odd int between 1 and q-1
			}
		}
	}

	if (code->frozen_bits == NULL){
		// initialization of the probabilities channels
		code->frozen_probs = (double *)malloc(sizeof(double) * code->n);
		for ( i = 0 ; i < code->n ; ++i){
			code->frozen_probs[i] = 0.0;
		}

		int *word = (int *)malloc(sizeof(int) * code->n);
		int *codeword = (int *)malloc(sizeof(int) * code->n);
		int *received = (int *)malloc(sizeof(int) * code->n);
		double **probs = (double **)malloc(sizeof(double *) * code->n);
		complex double *probs_tmp1 = (complex double *)malloc(sizeof(complex double) * code->ring->q);// complex is needed for fft in convolutional product
		complex double *probs_tmp2 = (complex double *)malloc(sizeof(complex double) * code->ring->q);// complex is needed for fft in convolutional product
		int *decoded = (int *)malloc(sizeof(int) * code->n);

		for (int i_test = 0 ; i_test < nb_tests ; ++i_test){

			// generate a random word then applying the U|U+V tree
			for ( i = 0 ; i < code->n ; ++i){
				word[i] = rand() % code->ring->q;
			}
			memcpy((void *)codeword, (void *)word, sizeof(int) * code->n);
			for ( i = code->m - 1 ; i >= 0 ; --i){
				for ( j = 0 ; j < (int)(1 << i) ; ++j){
					unsigned int n_tmp = (int)(1 << (code->m - i - 1));
					for ( jj = j * 2 * n_tmp ; jj < (j * 2 * n_tmp + n_tmp) ; ++jj){
						codeword[jj] = ( ( codeword[jj] + codeword[jj + n_tmp] ) % code->ring->q );
						codeword[jj + n_tmp] = ( ( code->coeffs[i][j] * codeword[jj + n_tmp] ) % code->ring->q );
					}
				}
			}

			// Gaussian noise
			for ( i = 0 ; i < code->nb_punc ; ++i){
				received[i] = -1;
			}
			for ( i = code->nb_punc ; i < code->n ; ++i){
				received[i] = ( (codeword[i] + gen_error(code->channel)) % code->ring->q );
			}

			// Probs received
			for ( i = 0 ; i < code->nb_punc ; ++i){
				probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
				for ( u = 0 ; u < code->ring->q ; ++u){
					probs[i][u] = 1.0 / (double)(code->ring->q);
				}
			}
			for ( i = code->nb_punc ; i < code->n ; ++i){
				probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
				for ( u = 0 ; u < code->ring->q ; ++u){
					probs[i][u] = code->channel->distrib[(received[i] - u) % code->ring->q];
				}
			}

			// Genie-aided Success Cancelation decoding
			memcpy((void *)decoded, (void *)codeword, sizeof(int) * code->n);
			for ( i = 0 ; i < code->m ; ++i){
				for ( j = 0 ; j < (int)(1 << i) ; ++j){
					unsigned int n_tmp = (int)(1 << (code->m - i - 1));
					for ( jj = j * 2 * n_tmp ; jj < (j * 2 * n_tmp + n_tmp) ; ++jj){

						decoded[jj + n_tmp] = ( ( invert(code->ring, code->coeffs[i][j]) * decoded[jj + n_tmp] ) % code->ring->q );
						decoded[jj] = ( ( decoded[jj] - decoded[jj + n_tmp] ) % code->ring->q );
						for ( u = 0 ; u < code->ring->q ; ++u){
							probs_tmp1[u] = probs[jj][u];
							probs_tmp2[u] = probs[jj + n_tmp][( -code->coeffs[i][j] * u ) % code->ring->q];
						}

						double norm_fact = 0;
						for ( u = 0 ; u < code->ring->q ; ++u){
							/*
							probs[jj][u] = 0;
							for ( int uu = 0 ; uu < code->ring->q ; ++uu){
								probs[jj][u] += creal(probs_tmp2[-uu % code->ring->q]) * creal(probs_tmp1[(u + uu) % code->ring->q]);
							}
							*/
							probs[jj + n_tmp][u] = creal(probs_tmp2[-u % code->ring->q]) * creal(probs_tmp1[(decoded[jj] + u) % code->ring->q]);
							norm_fact += creal(probs[jj + n_tmp][u]);
						}
						for ( u = 0 ; u < code->ring->q ; ++u){
							probs[jj + n_tmp][u] /= norm_fact;
						}

						// convolutional product
						fft(probs_tmp1, code->ring->q);
						fft(probs_tmp2, code->ring->q);

						for ( u = 0 ; u < code->ring->q ; ++u){
							probs_tmp1[u] *= probs_tmp2[u];
						}
						ifft(probs_tmp1, code->ring->q);

						for ( u = 0 ; u < code->ring->q ; ++u){
							probs[jj][u] = creal(probs_tmp1[u]);
						}
					}
				}
			}

			/*
			// uncomment if you want to print the vector probabilities
			for ( u = 0 ; u < code.ring.q ; ++u){
				printf("probs[%d] = ", u);
				for ( i = 0 ; i < code.n ; ++i){
					printf("\t%.4f", probs[i][u]);
				}
				printf("\n");
			}
			printf("\n\n");
			*/


			for ( i = 0 ; i < code->n ; ++i){
				code->frozen_probs[i] += (1.0 - probs[i][word[i]]);
			}

		}

		// average
		for ( i = 0 ; i < code->n ; ++i){
			code->frozen_probs[i] /= nb_tests;
		}

		// froze the n - k worst channels
		double *sorted_frozen_probs = (double *)malloc(sizeof(double) * code->n);
		memcpy((void *)sorted_frozen_probs, (void *)(code->frozen_probs), sizeof(double) * code->n);
		quick_sort(sorted_frozen_probs, 0, code->n - 1);

		code->frozen_bits = (int *)malloc(sizeof(int) * code->n);
		for ( i = 0 ; i < code->n ; ++i){
			code->frozen_bits[i] = ( code->frozen_probs[i] < sorted_frozen_probs[code->k] ? -1 : 0);
			//code->frozen_bits[i] = ( code->frozen_probs[i] < sorted_frozen_probs[code->k] ? -1 : rand() % code->ring->q); // coset of the polar code
		}

		// free memory
		free(word);
		free(codeword);
		free(received);
		for ( i = 0 ; i < code->n ; ++i){
			free(probs[i]);
		}
		free(probs);
		free(probs_tmp1);
		free(probs_tmp2);
		free(decoded);
		free(sorted_frozen_probs);
	}

	return code;
}



/*
 * @brief recursive function of the probabilistic or deterministic Success Cancelation decoding.
 * @param code the polar code.
 * @param probs the vector probabilities for each coordinate of the received word
 *              --> for all i in {0..n-1} and j in {0..q-1}, probs[i][j] = probability that the i-th coordinate is j.
 * @param n length of the current word to decode.
 * @param depth level in the global polar code.
 * @param index the coordinate to considere in the global polar code.
 * @param flag_probabilistic =1 if probabilistic decoding or =0 if deterministic.
 *
 * @return A codeword with a high probability.
 */
int *sc_decoder(polar_t *code, double **probs, int n, int depth, int index, int flag_probabilistic){
	unsigned int i, u;
	int *decoded = (int *)malloc(sizeof(int) * n);

	if (n == 1){ // leaf treatment
		if (code->frozen_bits[index] < 0){ // no frozen bit

			if (flag_probabilistic){ // probabilistic decoder
				double rand_tmp = (double)(rand())/(double)RAND_MAX;
				double cumul_prob = 0;
				for ( u = 0 ; u < code->ring->q ; ++u){
					cumul_prob += probs[0][u];
					if (rand_tmp < cumul_prob){
						break;
					}
				}
				if (u == code->ring->q){
					u = code->ring->q - 1;
				}
				decoded[0] = u;
			} else { // deterministic decoder (taking the best probabilities)
				unsigned int u_opt = 0;
				double prob_opt = probs[0][0];
				for (unsigned int u = 1 ; u < code->ring->q ; ++u){
					if (prob_opt < probs[0][u]){
						prob_opt = probs[0][u];
						u_opt = u;
					}
				}
				decoded[0] = u_opt;
			}
		} else { // frozen bit (the value of the bit is given by the code)

			decoded[0] = code->frozen_bits[index];
		}
	} else { // recursive call

		double **probs_V = (double **)malloc(sizeof(double *) * n/2);
		complex double *probs_tmp1 = (complex double *)malloc(sizeof(complex double) * code->ring->q);// complex is needed for fft in convolutional product
		complex double *probs_tmp2 = (complex double *)malloc(sizeof(complex double) * code->ring->q);// complex is needed for fft in convolutional product

		for ( i = 0 ; i < n/2 ; ++i){

			probs_V[i] = (double *)malloc(sizeof(double) * code->ring->q);


			for ( u = 0 ; u < code->ring->q ; ++u){
				probs_tmp1[u] = probs[i][u];
				probs_tmp2[u] = probs[i + n/2][( -code->coeffs[depth][index] * u ) % code->ring->q];
			}

			// convolutional product
			fft(probs_tmp1, code->ring->q);
			fft(probs_tmp2, code->ring->q);

			for ( u = 0 ; u < code->ring->q ; ++u){
				probs_tmp1[u] *= probs_tmp2[u];
			}
			ifft(probs_tmp1, code->ring->q);

			for ( u = 0 ; u < code->ring->q ; ++u){
				probs_V[i][u] = creal(probs_tmp1[u]);
			}

		}
		int *decoded_V = sc_decoder(code, probs_V, n/2, depth + 1, 2*index, flag_probabilistic);

		double **probs_U = (double **)malloc(sizeof(double *) * n/2);
		for ( i = 0 ; i < n/2 ; ++i){
			probs_U[i] = (double *)malloc(sizeof(double) * code->ring->q);
			double norm_fact = 0;
			for ( u = 0 ; u < code->ring->q ; ++u){
				probs_U[i][u] = probs[i + n/2][( code->coeffs[depth][index] * u ) % code->ring->q] * probs[i][(u + decoded_V[i]) % code->ring->q];
				norm_fact += probs_U[i][u];
			}
			for ( u = 0 ; u < code->ring->q ; ++u){
				probs_U[i][u] /= norm_fact;
			}
		}
		int *decoded_U = sc_decoder(code, probs_U, n/2, depth + 1, 2*index + 1, flag_probabilistic);

		for ( i = 0 ; i < n/2 ; ++i){
			decoded[i] = ( ( decoded_V[i] + decoded_U[i] ) % code->ring->q );
			decoded[i + n/2] = ( ( code->coeffs[depth][index] * decoded_U[i] ) % code->ring->q );
		}

		// free memory
		for ( i = 0 ; i < n/2 ; ++i){
			free(probs_V[i]);
			free(probs_U[i]);
		}
		free(probs_V);
		free(probs_U);
		free(decoded_V);
		free(decoded_U);
		free(probs_tmp1);
		free(probs_tmp2);
	}

	return decoded;
}


/*
 * @brief List decoding of the polar code using probabilistic SC decoder.
 * @param code the polar code.
 * @param probs the vector probabilities for each coordinate of the received word
 *              --> for all i in {0..n-1} and j in {0..q-1}, probs[i][j] = probability that the i-th coordinate is j.
 * @param received the word to decode.
 * @param list_size Size of the list decoding.
 *
 * Note that the negative coordinates in @received correspond to some punctured positions...
 * In that case, the probability vectors corresponding to those positions must be initialized with a maximal entropy
 * and those positions must be ignored in the calculation of distance.
 *
 * @return A codeword with a high probability.
 */
int *psc_list_decoder(polar_t *code, double **probs, int *received, unsigned int list_size){
	int i;

	// detection of punctered positions
	for ( i = 0 ; i < code->n ; ++i){
		if (received[i] < 0){
			for (int j = 0 ; j < code->ring->q ; ++j){
				probs[i][j] = 1.0/(double)(code->ring->q);
			}
		}
	}

	int *decoded_min = sc_decoder(code, probs, code->n, 0, 0, 0); // at least one decoding taking the best probabilities
	double dist_min = dist(decoded_min, received, code->ring->q, code->n);

	for ( i = 1 ; i < list_size ; ++i){
		int *decoded_cur = sc_decoder(code, probs, code->n, 0, 0, 1);
		double dist_cur = dist(decoded_cur, received, code->ring->q, code->n);
		if (dist_cur < dist_min){
			dist_min = dist_cur;
			memcpy((void *)decoded_min, (void *)decoded_cur, sizeof(int) * code->n);
		}
		free(decoded_cur);
	}

	return decoded_min;
}

#endif /* POLAR_CODE_H_ */
