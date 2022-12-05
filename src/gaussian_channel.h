/*
 * gaussian_channel.h
 *
 *  Created on: 24 nov. 2022
 *      Author: Kévin Carrier
 *
 * Description: Simulation of a Gaussian channel
 */

#ifndef GAUSSIAN_CHANNEL_H_
#define GAUSSIAN_CHANNEL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RAND_PRECISION 1000000
#define LAYOVER 10

typedef struct Channel {
	double std_dev;
	unsigned int q;
	double *distrib;
	unsigned int *rand_set;
} channel_t;

void free_channel(channel_t *channel){
	free(channel->distrib);
	free(channel->rand_set);
}

channel_t *gen_channel(unsigned int q, double std_dev){
	channel_t *channel = (channel_t *)malloc(sizeof(channel_t));
	channel->q = q;
	channel->std_dev = std_dev;
    channel->distrib = (double *)malloc(sizeof(double) * q);
    double prob;

    double variance = std_dev * std_dev;
    double C = 0;
    for (int x = (int)(- LAYOVER * q - ((q - 2)>>1)) ; x <= (int)(LAYOVER * q + (q >> 1)) ; ++x){
    	C += exp(- (double)(x * x)/(double)(2 * variance));
    }

    for ( int i = (int)(-((q - 2)>>1)) ; i <= (int)(q >> 1) ; ++i){
    	prob = 0;
    	for ( int x = -LAYOVER ; x <= LAYOVER ; ++x){
    		prob += exp(-(double)((i + q * x) * (i + q * x))/(double)(2 * variance));
    	}
    	prob /= C;
    	channel->distrib[i % q] = prob;
    }

    channel->rand_set = (unsigned int *)malloc(sizeof(unsigned int)*RAND_PRECISION);

    unsigned int value = 0;
    double bound = RAND_PRECISION * channel->distrib[value];

    for (int i = 0 ; i < RAND_PRECISION ; ){
    	if (i > bound){
    		value += 1;
    		if (value >= q){
    			break;
    		}
    		bound += RAND_PRECISION * channel->distrib[value];
    	} else {
    		channel->rand_set[i++] = value;
    	}
    }

    return channel;
}

void print_channel(FILE *fd, channel_t *channel){
	double sum_probs = 0;

	fprintf(fd, "Distribution:\n");
	for (int i = 0 ; i < channel->q ; ++i){
		sum_probs += channel->distrib[i];
		printf("[%d , %.12f],\n", i, channel->distrib[i]);
	}

	fprintf(fd, "\n\nSum of probabilities: %.12f\n\n", sum_probs);
}

unsigned int gen_error(channel_t *channel){
	return channel->rand_set[rand() % RAND_PRECISION];
}

#endif /* GAUSSIAN_CHANNEL_H_ */
