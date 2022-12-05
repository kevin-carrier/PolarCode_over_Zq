/*
 ============================================================================
 Name        : BothMay_simulations.c
 Author      : Kévin Carrier
 Version     :
 Copyright   : open source
 Description : Simulation for the Both-May algorithm in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// practical measurement tools

typedef struct Vector
{
    uint64_t *v;
    unsigned int len;
    unsigned int nb_blocs;
} Vector;

typedef struct Matrix
{
	Vector **cols;
    unsigned int nrows;
    unsigned int ncols;
} Matrix;

typedef struct SparseVector
{
    unsigned int *pos;
    unsigned int len;
    unsigned int weight;
} SparseVector;

typedef struct s_List_data List_data;
struct s_List_data
{
    List_data *next;
    SparseVector *info;
    Vector *red1;
    Vector *red2;
    Vector *red3;
};

typedef struct s_List List;
struct s_List
{
    List_data *list_data;
    unsigned long long size;
};

/*
 *  @brief Free the memory occupied by a vector
 *  @param vect the vector to free
 */
void vector_clear(Vector *vect)
{
	free(vect->v);
	free(vect);
}

/*
 *  @brief Free the memory occupied by a matrix
 *  @param M the matrix to free
 */
void matrix_clear(Matrix *M)
{
	for (int i = 0 ; i < M->ncols ; ++i)
	{
		vector_clear(M->cols[i]);
	}
	free(M->cols);
	free(M);
}

/*
 *  @brief Free the memory occupied by a sparse vector
 *  @param sv the sparse vector to free
 */
void sv_clear(SparseVector *sv)
{
	free(sv->pos);
	free(sv);
}

/*
 *  @brief Generate a first sparse vector of given weight
 *  @param len the length of the sparse vector
 *  @param weight the weight of the sparse vector
 *  @return the vector of length len and weight weight where all the 1 are placed to the left: 11111...000000000...
 */
SparseVector *sv_init(unsigned int len, unsigned int weight)
{
	SparseVector *sv = (SparseVector *)malloc(sizeof(SparseVector));

	sv->len = len;
	sv->weight = weight;
	sv->pos = (unsigned int *)malloc(sizeof(unsigned int)*weight);

	for (int i = 0 ; i < weight ; ++i)
	{
		sv->pos[i] = i;
	}

	return sv;
}

/*
 *  @brief Xor two sparse vectors
 *  @param sv1, sv2 the sparse vectors to xor
 *  @return sv1 + sv2 where + stands for the xor operation
 */
SparseVector *sv_add(SparseVector *sv1, SparseVector *sv2)
{
	if (sv1->len != sv2->len)
	{
		fprintf(stderr, "sv_add error: sv1 and sv2 cannot be added!\n");
		exit(EXIT_FAILURE);
	}
	SparseVector *sv = sv_init(sv1->len, sv1->weight + sv2->weight);
	int i = 0, j = 0;
	sv->weight = 0;
	while ((i < sv1->weight) && (j < sv2->weight))
	{
		if (sv1->pos[i] < sv2->pos[j])
		{
			sv->pos[sv->weight] = sv1->pos[i];
			++(sv->weight);
			++i;
		}
		else if (sv1->pos[i] > sv2->pos[j])
		{
			sv->pos[sv->weight] = sv2->pos[j];
			++(sv->weight);
			++j;
		}
		else
		{
			++i;
			++j;
		}
	}
	while (i < sv1->weight)
	{
		sv->pos[sv->weight] = sv1->pos[i];
		++(sv->weight);
		++i;
	}
	while (j < sv2->weight)
	{
		sv->pos[sv->weight] = sv2->pos[j];
		++(sv->weight);
		++j;
	}
	return sv;
}

/*
 *  @brief Concatenate two sparse vectors
 *  @params sv1, sv2 the sparse vectors to concatenate
 *  #return sv1 | sv2 where | stands for the concatenation
 */
SparseVector *sv_concat(SparseVector *sv1, SparseVector *sv2)
{
	SparseVector *sv = sv_init(sv1->len + sv2->len, sv1->weight + sv2->weight);
	for (int i = 0 ; i < sv1->weight ; ++i){
		sv->pos[i] = sv1->pos[i];
	}
	for (int i = 0 ; i < sv2->weight ; ++i){
		sv->pos[sv1->weight + i] = (sv1->len + sv2->pos[i]);
	}
	return sv;
}


/*
 *  @brief For enumerating all the sparse vectors of weight sv->weight. Start with sv_int() and end with NULL pointer.
 *  @params sv the sparse vector to "increment"
 */
void sv_next(SparseVector *sv)
{
	if (sv->pos[0] == sv->len - sv->weight)
	{
		free(sv->pos);
		sv->pos = 0;
	}
	else
	{
		for (int i = sv->weight - 1 ; i >= 0 ; --i){
			++(sv->pos[i]);
			if (sv->pos[i] != sv->len - sv->weight + 1 + i){
				for (int j = i + 1 ; j < sv->weight ; ++j){
					sv->pos[j] = sv->pos[j-1] + 1;
				}
				break;
			}
		}
	}
	return;
}

/*
 *  @brief A comparison function for sparse vectors
 *  @params sv1, sv2 the sparse vectors to compare
 *  #return -1 if sv1 < sv2, 0 if sv1 == sv2 and 1 if sv1 > sv2
 */
int sv_compare(SparseVector *sv1, SparseVector *sv2){
	if (sv1->len != sv2->len){
		fprintf(stderr, "sv_compare error: sv1 and sv2 cannot be compared!\n");
		exit(EXIT_FAILURE);
	}
	unsigned int min_weight = fmin(sv1->weight, sv2->weight);
	for (unsigned int i = 0 ; i < min_weight ; ++i){
		if(sv1->pos[i] < sv2->pos[i]){
			return -1;
		} else if (sv1->pos[i] > sv2->pos[i]) {
			return 1;
		}
	}
	if (sv1->weight < sv2->weight){
		return -1;
	} else if (sv1->weight > sv2->weight) {
		return 1;
	}
	return 0;
}

/*
 *  @brief Convert a vector into a string of binary values
 *  @param vect the vector to convert
 *  #return the binary representation of the vector vect
 */
char *vect_to_str(Vector *vect)
{
	char *vect_str = (char *)malloc(sizeof(char)*vect->len);
	int j = 0;
	for (int i = 0 ; i < vect->nb_blocs ; ++i)
	{
		uint64_t tmp = vect->v[i];
		do
		{
			vect_str[j] = (tmp&1 ? '1' : '0');
			tmp >>= 1;
			j++;
		} while (j%64);
	}

	vect_str[vect->len] = '\0';
	return vect_str;
}

/*
 *  @brief Generate the all zero vector
 *  @param len the length of the generated vector
 *  #return the all zero vector of length len
 */
Vector *zero_vector(unsigned int len){
	Vector *vect = (Vector *)malloc(sizeof(Vector));
	vect->len = len;
	vect->nb_blocs = (vect->len/64 + (vect->len%64 ? 1 : 0));
	vect->v = (uint64_t*)calloc(vect->nb_blocs, sizeof(uint64_t));
	return vect;
}


/*
 *  @brief Generate a random vector
 *  @param len the length of the generated vector
 *  #return a random vector of length len
 */
Vector *random_vector(unsigned int len){
	Vector *vect = (Vector *)malloc(sizeof(Vector));
	vect->len = len;
	vect->nb_blocs = (vect->len/64 + (vect->len%64 ? 1 : 0));
	vect->v = (uint64_t*)malloc(vect->nb_blocs*sizeof(uint64_t)*8);
	for (int i = 0 ; i < vect->nb_blocs ; ++i)
	{
		vect->v[i] = (((uint64_t)rand()) << 60) ^ (((uint64_t)rand()) << 30) ^ rand();
	}
	return vect;
}

/*
 *  @brief Compute the Hamming weight of a vector
 *  @param vect a vector
 *  #return the Hamming weight of vect
 */
unsigned int hamming_weight(Vector *vect)
{
	unsigned int hw = 0;
	if ((vect->len)%64)
	{
		for (int i = 0 ; i < vect->nb_blocs - 1 ; ++i){
			hw += __builtin_popcount(vect->v[i] & 0xffffffff);
			hw += __builtin_popcount((vect->v[i] >> 32) & 0xffffffff);
		}
		uint64_t mask = (((uint64_t)1) << ((vect->len)%64)) - 1;
		hw += __builtin_popcount((vect->v[vect->nb_blocs - 1]) & mask & 0xffffffff);
		hw += __builtin_popcount((((vect->v[vect->nb_blocs - 1]) & mask) >> 32) & 0xffffffff);
	}
	else
	{
		for (int i = 0 ; i < vect->nb_blocs ; ++i){
			hw += __builtin_popcount(vect->v[i] & 0xffffffff);
			hw += __builtin_popcount((vect->v[i] >> 32) & 0xffffffff);
		}
	}
	return hw;
}

/*
 *  @brief Generate a random matrix
 *  @param nrows the number of rows
 *  @param ncols the number of columns
 *  #return a random matrix of size nrows x ncols
 */
Matrix *random_matrix(unsigned int nrows, unsigned int ncols)
{
	Matrix *M = (Matrix *)malloc(sizeof(Matrix));
	M->cols = (Vector **)malloc(sizeof(Vector *) * ncols);
	M->ncols = ncols;
	M->nrows = nrows;
	for (int i = 0 ; i < ncols ; ++i)
	{
		M->cols[i] = random_vector(nrows);
	}
	return M;
}

/*
 *  @brief Vector / matrix product
 *  @param sv a sparse vector
 *  @param M a matrix having as much columns as the length of sv
 *  #return the product sv x M
 */
Vector *prod_sv_matrix(SparseVector *sv, Matrix *M)
{
	if (sv->len != M->ncols){
		fprintf(stderr, "prod_sv_matrix error: the length of sv must be equal to the number of columns in M!\n");
		exit(EXIT_FAILURE);
	}
	Vector *res = zero_vector(M->nrows);
	for (int i = 0 ; i < sv->weight ; ++i)
	{
		for (int j = 0 ; j < res->nb_blocs ; ++j)
		{
			res->v[j] ^= M->cols[sv->pos[i]]->v[j];
		}
	}
	return res;
}

/*
 *  @brief Xor the vector v to the vector dest : dest = dest + v
 *  @param dest the destination vector
 *  @param v the vector to add
 */
void vector_add(Vector *dest, Vector *v)
{
	if (dest->len != v->len){
		fprintf(stderr, "vector_add error: dest and v must have the same length!\n");
		exit(EXIT_FAILURE);
	}
	for (int j = 0 ; j < dest->nb_blocs ; ++j)
	{
		dest->v[j] ^= v->v[j];
	}
}

/*
 *  @brief Create an empty list
 *  @return an empty list
 */
List *list_create()
{
	List *list = (List *)malloc(sizeof(List));
	list->list_data = NULL;
	list->size = 0;
	return list;
}

/*
 *  @brief Free the memory occupied by a list
 *  @param list the list to free
 */
void list_clear(List *list)
{
	List_data *next;
	List_data *list_data = list->list_data;
	while (list_data != NULL){
		next = list_data->next;
		sv_clear(list_data->info);
		vector_clear(list_data->red1);
		vector_clear(list_data->red2);
		vector_clear(list_data->red3);
		free(list_data);
		list_data = next;
	}
	free(list);
}

/*
 *  @brief Append an element in a list if the element is not already in the list
 *  @param list the list to modify
 *  @params info, red1, red2, red3 the data to append in the list
 */
void list_append(List *list, SparseVector *info, Vector *red1, Vector *red2, Vector *red3)
{
	if (list->list_data == NULL)
	{
		List_data *list_data = (List_data *)malloc(sizeof(List_data));
		list_data->next = NULL;
		list_data->info = info;
		list_data->red1 = red1;
		list_data->red2 = red2;
		list_data->red3 = red3;
		list->list_data = list_data;
		++(list->size);
		return;
	}

	int cmp = sv_compare(list->list_data->info, info);
	if (cmp < 0)
	{
		List_data *list_data = (List_data *)malloc(sizeof(List_data));
		list_data->next = list->list_data;
		list_data->info = info;
		list_data->red1 = red1;
		list_data->red2 = red2;
		list_data->red3 = red3;
		list->list_data = list_data;
		++(list->size);
		return;
	}
	if (cmp > 0)
	{
		List_data *list_data = list->list_data;
		while (list_data->next != NULL){
			cmp = sv_compare(list_data->next->info, info);
			if (cmp <= 0){
				break;
			}
			list_data = list_data->next;
		}
		if (list_data->next == NULL)
		{
			List_data *list_data_tmp = (List_data *)malloc(sizeof(List_data));
			list_data_tmp->next = NULL;
			list_data_tmp->info = info;
			list_data_tmp->red1 = red1;
			list_data_tmp->red2 = red2;
			list_data_tmp->red3 = red3;
			list_data->next = list_data_tmp;
			++(list->size);
			return;
		}
		if (cmp < 0){
			List_data *list_data_tmp = (List_data *)malloc(sizeof(List_data));
			list_data_tmp->next = list_data->next;
			list_data_tmp->info = info;
			list_data_tmp->red1 = red1;
			list_data_tmp->red2 = red2;
			list_data_tmp->red3 = red3;
			list_data->next = list_data_tmp;
			++(list->size);
			return;
		}

	}
}


// theoretical measurement tools

/*
 *  @brief The binomial function
 *  @param n size of the vector
 *  @param k destination weight
 *  @return the number of n length vectors of weight k
 */
double binomial(unsigned int n, unsigned int k)
{
	if ((k < 0) || (k > n))
	{
		return 0.0;
	}

	double res = 1.0;
	if (n > 2*k)
	{
		for (unsigned int i = n ; i > n - k ; --i)
		{
			res *= i;
		}
		for (unsigned int i = k ; i > 1 ; --i)
		{
			res /= i;
		}
	}
	else
	{
		for (unsigned int i = n ; i > k ; --i)
		{
			res *= i;
		}
		for (unsigned int i = n - k ; i > 1 ; --i)
		{
			res /= i;
		}
	}
	return res;
}


/*
 *  @brief Number of pairs of vectors both of weight ww such that their sum is a given vector of weight w
 *  @param n size of the vectors
 *  @param w destination weight
 *  @param ww source weight
 *  @return the number of pairs of vectors both of weight ww such that their sum is a given vector of weight w
 */
double nb_repr(unsigned int n, unsigned int w, unsigned int ww)
{
    return binomial(n - w, ww - (w/2))*binomial(w, (w/2));
}


/*
 *  @brief Probability that the sum of two vectors both of weight ww is a vector of weight w
 *  @param n size of the vectors
 *  @param w destination weight
 *  @param ww source weight
 *  @return the probability that the sum of two vectors both of weight ww is a vector of weight w
 */
double prob_repr(unsigned int n, unsigned int w, unsigned int ww)
{
    return binomial(n, w) * nb_repr(n, w, ww) / (binomial(n, ww)*binomial(n, ww));
}


// main
int main(void)
{
	srand( time( NULL ) );

    // parameters to test
	unsigned int n = 120, k = 30;
	unsigned int r1 = 32, r2 = 13, r3 = n - k - r1 - r2;
	unsigned int p1 = 6, p0 = p1/2, p2 = 8, p3 = 10;
	unsigned int w11 = 12, w12 = 24, w13 = 8;
	unsigned int w22 = 8, w23 = 10;
	unsigned int w33 = 14;

	printf("*** Recalling parameters ***\n\n");
	printf("\tDecoding problem parameters: n = %u,  k = %u,  w = %u\n", n, k, p3+w13+w23+w33);
	printf("\tSize of the redundancy subsets: r1 = %u,  r2 = %u,  r3 = %u\n", r1, r2, r3);
	printf("\tStage 1 parameters: w11 = %u,  w12 = %u,  w13 = %u\n", w11, w12, w13);
	printf("\tStage 2 parameters: w22 = %u,  w23 = %u\n", w22, w23);
	printf("\tStage 3 parameters: w33 = %u\n", w33);


	// theoretical values
	printf("\n*** Theoretical values ***\n\n");
	double E2 = nb_repr(k, p2, p1) * nb_repr(r1, w12, w11) / (double)(((unsigned long long)1) << r1);
	double E3 = nb_repr(k, p3, p2) * (nb_repr(r1, w13, w12) / (double)(((unsigned long long)1) << r1)) * (nb_repr(r2, w23, w22) / (double)(((unsigned long long)1) << r2));

	double S0 =  binomial(k/2, p0);

	//double S1 = S0 * S0 * binomial(r1, w11) / (double)(((unsigned long long)1) << r1);
	double S1 = binomial(k, p1) * binomial(r1, w11) / (double)(((unsigned long long)1) << r1);

	double S2_BM18 = binomial(k, p2) * (binomial(r2, w22) / (double)(((unsigned long long)1) << r2)) * prob_repr(r1, w12, w11);
	double S2_bound = binomial(k, p2) * (binomial(r2, w22) / (double)(((unsigned long long)1) << r2)) * (binomial(r1, w12) / (double)(((unsigned long long)1) << r1));
	double S2 = S1 * S1 * prob_repr(k, p2, p1) * (binomial(r2, w22) / (double)(((unsigned long long)1) << r2)) * prob_repr(r1, w12, w11);
	S2 = fmin(S2_bound, S2);

	double S3_BM18 = binomial(k, p3) * (binomial(r3, w33)/(double)(((unsigned long long)1) << r3)) * prob_repr(r1, w13, w12) * prob_repr(r2, w23, w22);
	double S3_bound = binomial(k, p3) * (binomial(r1, w13)/(double)(((unsigned long long)1) << r1)) * (binomial(r2, w23)/(double)(((unsigned long long)1) << r2)) * (binomial(r3, w33)/(double)(((unsigned long long)1) << r3));
	double S3 = S2 * S2 * (binomial(r3, w33)/(double)(((unsigned long long)1) << r3)) * prob_repr(k, p3, p2) * prob_repr(r1, w13, w12) * prob_repr(r2, w23, w22);
	S3= fmin(S3_bound, S3);

	printf("Correctness lemma: \n\tE2 = %.12f\n\tE3 = %.12f\n", E2, E3);
	printf("Size of the lists (bound from [BM18]): \n\tS0 <= %.12f\n\tS1 <= %.12f\n\tS2 <= %.12f\n\tS3 <= %.12f\n", S0, S1, S2_BM18, S3_BM18);
	printf("Size of the lists (corrected): \n\tS0 = %.12f\n\tS1 = %.12f\n\tS2 = %.12f\n\tS3 = %.12f\n", S0, S1, S2, S3);


	// practical values

	printf("\n*** Practical size of the lists ***\n\n");


	unsigned int nb_tests = 10;
	unsigned long long sum_L1_size = 0;
	unsigned long long sum_L2_size = 0;
	unsigned long long sum_L3_size = 0;

	for (unsigned int i = 0 ; i < nb_tests ; ++i)
	{
/*
		// random source (information part of the parity-check matrix)
		Matrix *H1_left = random_matrix(r1, k/2);
		Matrix *H1_right = random_matrix(r1, k/2);
		Matrix *H2 = random_matrix(r2, k);
		Matrix *H3 = random_matrix(r3, k);


		// computing L1
		List *L1 = list_create();
		SparseVector *info1 = sv_init(k/2, p0);
		while (info1->pos){
			Vector *y1_tmp = prod_sv_matrix(info1, H1_left);
			SparseVector *info2 = sv_init(k/2, p0);
			while (info2->pos){
				Vector *y1 = prod_sv_matrix(info2, H1_right);
				vector_add(y1, y1_tmp);

				if (hamming_weight(y1) == w11)
				{
					SparseVector *info = sv_concat(info1, info2);
					Vector *y2 = prod_sv_matrix(info, H2);
					Vector *y3 = prod_sv_matrix(info, H3);
					list_append(L1, info, y1, y2, y3);
				}
				else
				{
					vector_clear(y1);
				}

				sv_next(info2);
			}
			sv_next(info1);
			sv_clear(info2);
		}
		sv_clear(info1);
		matrix_clear(H1_left);
		matrix_clear(H1_right);
		matrix_clear(H2);
		matrix_clear(H3);

		sum_L1_size += L1->size;
*/

		// random source (information part of the parity-check matrix)
		Matrix *H1 = random_matrix(r1, k);
		Matrix *H2 = random_matrix(r2, k);
		Matrix *H3 = random_matrix(r3, k);


		// computing S1 without separating the information set in two distinct parts (asymptotically equivalent)
		List *L1 = list_create();
		SparseVector *info = sv_init(k, p1);
		while (info->pos){
			Vector *y1 = prod_sv_matrix(info, H1);
			if (hamming_weight(y1) == w11)
			{
				SparseVector *info_copy = sv_init(k, p1);
				for (int i = 0 ; i < p1 ; ++i)
				{
					info_copy->pos[i] = info->pos[i];
				}
				Vector *y2 = prod_sv_matrix(info, H2);
				Vector *y3 = prod_sv_matrix(info, H3);
				list_append(L1, info_copy, y1, y2, y3);

			}
			else
			{
				vector_clear(y1);
			}
			sv_next(info);
		}
		sv_clear(info);
		matrix_clear(H1);
		matrix_clear(H2);
		matrix_clear(H3);

		sum_L1_size += L1->size;
		printf("size of the lists L1: %llu\n", L1->size);

		// computing S2
		List *L2 = list_create();
		List_data *data1;
		List_data *data2;
		data1 = L1->list_data;
		while(data1)
		{
			data2 = L1->list_data;
			while(data2)
			{
				SparseVector *info = sv_add(data1->info, data2->info);
				if (info->weight == p2)
				{
					Vector *y2 = zero_vector(r2);
					vector_add(y2, data1->red2);
					vector_add(y2, data2->red2);
					if (hamming_weight(y2) == w22)
					{
						Vector *y1 = zero_vector(r1);
						vector_add(y1, data1->red1);
						vector_add(y1, data2->red1);

						if (hamming_weight(y1) == w12)
						{
							Vector *y3 = zero_vector(r3);
							vector_add(y3, data1->red3);
							vector_add(y3, data2->red3);
							list_append(L2, info, y1, y2, y3);
						}
						else
						{
							sv_clear(info);
							vector_clear(y1);
							vector_clear(y2);
						}
					}
					else
					{
						sv_clear(info);
						vector_clear(y2);
					}
				}
				else
				{
					sv_clear(info);
				}
				data2 = data2->next;
			}
			data1 = data1->next;
		}
		list_clear(L1);

		sum_L2_size += L2->size;
		printf("size of the lists L2: %llu\n", L2->size);

		//computing S3
		List *L3 = list_create();
		data1 = L2->list_data;
		while(data1)
		{
			data2 = L2->list_data;
			while(data2)
			{
				SparseVector *info = sv_add(data1->info, data2->info);
				if (info->weight == p3)
				{
					Vector *y3 = zero_vector(r3);
					vector_add(y3, data1->red3);
					vector_add(y3, data2->red3);
					if (hamming_weight(y3) == w33)
					{
						Vector *y2 = zero_vector(r2);
						vector_add(y2, data1->red2);
						vector_add(y2, data2->red2);
						if (hamming_weight(y2) == w23)
						{
							Vector *y1 = zero_vector(r1);
							vector_add(y1, data1->red1);
							vector_add(y1, data2->red1);
							if (hamming_weight(y1) == w13)
							{
								list_append(L3, info, y1, y2, y3);
							}
							else
							{
								sv_clear(info);
								vector_clear(y1);
								vector_clear(y2);
								vector_clear(y3);
							}
						}
						else
						{
							sv_clear(info);
							vector_clear(y2);
							vector_clear(y3);
						}
					}
					else
					{
						sv_clear(info);
						vector_clear(y3);
					}
				}
				else
				{
					sv_clear(info);
				}
				data2 = data2->next;
			}
			data1 = data1->next;
		}
		list_clear(L2);

		sum_L3_size += L3->size;
		printf("size of the lists L3: %f\n", (double)(L3->size));
		list_clear(L3);
		printf("\n");
	}

	printf("Average size of the lists L1: %f\n", ((double)sum_L1_size)/((double)nb_tests));
	printf("Average size of the lists L2: %f\n", ((double)sum_L2_size)/((double)nb_tests));
	printf("Average size of the lists L3: %f\n", ((double)sum_L3_size)/((double)nb_tests));

	return EXIT_SUCCESS;
}
