/*
 ============================================================================
 Name        : PolarCode_over_Zq.c
 Author      : Kévin Carrier
 Version     :
 Copyright   : open source
 Description : Distortion of a Polar Code over ZZ/qZZ where q=2^s.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "polar_code.h"


void psc_list_decoder_test(unsigned int s, unsigned int n, unsigned int k, unsigned int list_size){
	polar_t *code = gen_polar(s, n, k, 100, NULL);
	//polar_t *code = gen_polar(s, n, k, 100);
	print_polar(stdout, code);
	printf("\n\n");

	printf("number of punctured positions: %d\n\n", code->nb_punc);

	int *word = (int *)malloc(sizeof(int) * code->n);
	int *codeword = (int *)malloc(sizeof(int) * code->n);
	int *received = (int *)malloc(sizeof(int) * code->n); // negative integer for the punctured positions
	double **probs = (double **)malloc(sizeof(double *) * code->n);

	/*
	 * generate a random word then applying the U|U+V tree
	 */
	for (int i = 0 ; i < code->n ; ++i){
		if (code->frozen_bits[i] == -1){
			word[i] = rand() % code->ring->q;
		} else {
			word[i] = code->frozen_bits[i];
		}
	}
	printf("    word =\t");
	for (int i = 0 ; i < code->n ; ++i){
		printf("%d\t", word[i]);
	}
	printf("\n");

	memcpy((void *)codeword, (void *)word, sizeof(int) * code->n);
	for (int i = code->m - 1 ; i >= 0 ; --i){
		for (int j = 0 ; j < (int)(1 << i) ; ++j){
			unsigned int n_tmp = (int)(1 << (code->m - i - 1));
			for (int jj = j * 2 * n_tmp ; jj < (j * 2 * n_tmp + n_tmp) ; ++jj){
				codeword[jj] = ( ( codeword[jj] + codeword[jj + n_tmp] ) % code->ring->q );
				codeword[jj + n_tmp] = ( ( code->coeffs[i][j] * codeword[jj + n_tmp] ) % code->ring->q );
			}
		}
	}
	printf("codeword =\t");
	for (int i = 0 ; i < code->n ; ++i){
		printf("%d\t", codeword[i]);
	}
	printf("\n");

	/*
	 * Gaussian noise + puncturing
	 */
	for (int i = 0 ; i < code->nb_punc ; ++i){
		received[i] = -1;
	}
	for (int i = code->nb_punc ; i < code->n ; ++i){
		received[i] = ( (codeword[i] + gen_error(code->channel)) % code->ring->q );
	}


	printf("received =\t");
	for (int i = 0 ; i < code->n ; ++i){
		printf("%d\t", received[i]);
	}
	printf("\n");


	/*
	 * Probs received
	 */
	for (int i = 0 ; i < code->n ; ++i){
		probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
		if (received[i] < 0){
			for (int u = 0 ; u < code->ring->q ; ++u){
				probs[i][u] = 1.0/(double)(code->ring->q);
			}
		} else {
			for (int u = 0 ; u < code->ring->q ; ++u){
				probs[i][u] = code->channel->distrib[(received[i] - u) % code->ring->q];
			}
		}
	}


	/*
	 * List decoding using Probabilistic Success Cancelation
	 */
	int *decoded = psc_list_decoder(code, probs, received, list_size);

	printf(" decoded =\t");
	for (int i = 0 ; i < code->n ; ++i){
		printf("%d\t", decoded[i]);
	}
	printf("\n\n");

	printf("dist(codeword, received) = %.12f\n", dist(codeword, received, code->ring->q, code->n));
	printf("dist(decoded,  received) = %.12f\n", dist(decoded, received, code->ring->q, code->n));



	for (int i = 0 ; i < code->n ; ++i){
		free(probs[i]);
	}
	free(probs);
	free(word);
	free(codeword);
	free(received);
	free(decoded);
}




double distortion(polar_t *code, unsigned int list_size, unsigned int nb_tests){
	int *randword = (int *)malloc(sizeof(int) * code->n); // negative integers for the punctured positions
	double **probs = (double **)malloc(sizeof(double *) * code->n);

	double sum_dist = 0;
	for (int i_test = 0 ; i_test < nb_tests ; ++i_test){

		/*
		 * generate a punctured random word
		 */
		for (int i = 0 ; i < code->nb_punc ; ++i){
			randword[i] = -1;
		}
		for (int i = code->nb_punc ; i < code->n ; ++i){
			randword[i] = rand() % code->ring->q;
		}
		/*
		printf("random word =\t");
		for (int i = 0 ; i < code->n ; ++i){
			printf("%d\t", randword[i]);
		}
		printf("\n");
		*/

		/*
		 * Probs
		 */
		for (int i = 0 ; i < code->nb_punc ; ++i){
			probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
			for (int u = 0 ; u < code->ring->q ; ++u){
				probs[i][u] = 1.0 / (double)(code->ring->q);
			}
		}
		for (int i = code->nb_punc ; i < code->n ; ++i){
			probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
			for (int u = 0 ; u < code->ring->q ; ++u){
				probs[i][u] = code->channel->distrib[(randword[i] - u) % code->ring->q];
			}
		}

		/*
		 * List decoding using Probabilistic Success Cancelation
		 */
		int *decoded = psc_list_decoder(code, probs, randword, list_size);
		/*
			printf("decoded =\t");
			for (int i = 0 ; i < code->n ; ++i){
				printf("%d\t", decoded[i]);
			}
			printf("\n\n");
		 */

		sum_dist += dist(decoded, randword, code->ring->q, code->n);
		//printf("dist(decoded, randword) = %.12f\n", sum_dist / (double)(i_test + 1));
		free(decoded);

	}

	// free memory
	for (int i = 0 ; i < code->n ; ++i){
		free(probs[i]);
	}
	free(probs);
	free(randword);
	return sum_dist/(double)nb_tests;
}



unsigned int s = 0, n = 0, k = 0 ;

unsigned int nb_samples = 100;
unsigned int nb_codes = 100;
unsigned int nb_randwords = 100;

int parse_cmdline(int argc, char *argv[]){
	int i;

	for (i = 1 ; i < argc ; ++i){
		if (argv[i][0] == '-') {
			switch (argv[i][1]){
			case 's':
				if (!argv[i][2]){
					++i;
					s = atoi(argv[i]);
				}else{
					s = atoi(argv[i]+2);
				}
				break;
			case 'n':
				if (!argv[i][2]){
					++i;
					n = atoi(argv[i]);
				}else{
					n = atoi(argv[i]+2);
				}
				break;
			case 'k':
				if (!argv[i][2]){
					++i;
					k = atoi(argv[i]);
				}else{
					k = atoi(argv[i]+2);
				}
				break;
			case 'S':
				if (!argv[i][2]){
					++i;
					nb_samples = atoi(argv[i]);
				}else{
					nb_samples = atoi(argv[i]+2);
				}
				break;
			case 'N':
				if (!argv[i][2]){
					++i;
					nb_codes = atoi(argv[i]);
				}else{
					nb_codes = atoi(argv[i]+2);
				}
				break;
			case 'R':
				if (!argv[i][2]){
					++i;
					nb_randwords = atoi(argv[i]);
				}else{
					nb_randwords = atoi(argv[i]+2);
				}
				break;
			case 'h' :
				printf("Usage : %s <arguments>\n", argv[0]);
				printf("\t\"-h\" : help\n");
				printf("\t\"-s\" : q = 2^s (unsigned int REQUIRED)\n");
				printf("\t\"-n\" : lengths of the polar code (unsigned int REQUIRED)\n");
				printf("\t\"-k\" : dimension of the polar code (unsigned int REQUIRED)\n");
				printf("\t\"-S\" : number of samples for estimating the frozen positions (unsigned int OPTIONAL : DEFAULT=100)\n");
				printf("\t\"-N\" : number of codes to test (unsigned int OPTIONAL : DEFAULT=100)\n");
				printf("\t\"-R\" : number of random words to decode for each distortion computation (unsigned int OPTIONAL : DEFAULT=100)Ón");
				return 2;
				break;
			default :
				fprintf(stderr, "Arguments parse error : argument \"-%c\" unknown !\n", argv[i][1]);
				fprintf(stderr, " -h for help !\n\n");
				return 1;
			}
		}
	}
	if (!s || !n || !k){
		fprintf(stderr, "Arguments parse error : \n");
		if (!s){
			fprintf(stderr, "\t argument -s is required !\n");

		}
		if (!n){
			fprintf(stderr, "\t argument -n is required !\n");

		}
		if (-!k){
			fprintf(stderr, "\t argument -k is required !\n");

		}
		fprintf(stderr, " -h for help !\n\n");
		return 1;
	}

	return 0;
}




int main(int argc, char * argv[]){
	srand(time(NULL));


	// parsing arguments
	switch(parse_cmdline(argc, argv)){
	case 0 :
		break;
	case 1 :
		exit(EXIT_FAILURE);
	case 2 :
		exit(EXIT_SUCCESS);
	}


	//psc_list_decoder_test(2, 14, 4, 20);

	double d, d_min = sqrt(n) * pow(2, s - 1);
	polar_t *best_code = NULL;
	polar_t *code = NULL;
	for ( int i = 0 ; i < nb_codes ; ++i){
		code = gen_polar(s, n, k, nb_samples, NULL);
		d = distortion(code, 1, nb_randwords);
		printf("distortion of the code %d: %.5f\n",i, d);
		if (d < d_min){
			d_min = d;
			if (best_code != NULL){
				free_polar(best_code);
			}
			best_code = copy_polar(code);
		}
		free_polar(code);
	}

	printf("\n\n ************************ BEGIN: The polar code parameters ************************");
	print_polar(stdout, best_code);
	printf("\n ************************  END: The polar code parameters  ************************\n\n");


	printf("dgv = %.12f\n", sqrt(best_code->n - best_code->nb_punc) * best_code->channel->std_dev);
	printf("decoding list of size 1 \n     --> average distortion = %.5f\n\n", distortion(best_code, 1, nb_randwords));
	printf("decoding list of size 8 \n     --> average distortion = %.5f\n\n", distortion(best_code, 8, nb_randwords));
	printf("decoding list of size 16 \n     --> average distortion = %.5f\n\n", distortion(best_code, 16, nb_randwords));


	free_polar(best_code);

	return EXIT_SUCCESS;
}
