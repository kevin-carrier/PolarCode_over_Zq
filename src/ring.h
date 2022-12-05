/*
 * ring.h
 *
 *  Created on: 24 nov. 2022
 *      Author: Kévin Carrier
 *
 * Description: Ring ZZ/qZZ where q=2^s
 */

#ifndef RING_H_
#define RING_H_

#include <stdio.h>
#include <stdlib.h>

typedef struct Ring {
	unsigned int s;
	unsigned int q;
	unsigned int *inverts;
} ring_t;

void free_ring(ring_t *ring){
	free(ring->inverts);
}

ring_t *gen_ring(unsigned int s){
	ring_t *ring = (ring_t *)malloc(sizeof(ring_t));
	ring->s = s;
	ring->q = (int)(1 << s);
	ring->inverts = (unsigned int *)malloc(sizeof(unsigned int) * ring->q);
	ring->inverts[0] = -1;
	for (unsigned int x = 0; x < ring->q; ++x){
		ring->inverts[x] = -1;
		for (unsigned int y = 0; y < ring->q; ++y){
			if (((x * y) % ring->q) == 1){
				ring->inverts[x] = y;
				break;
			}
		}
	}
	return ring;
}

unsigned int invert(ring_t *ring, unsigned int x){
	if ((x % 2) == 0){
		fprintf(stderr, "invert error: %d cannot be inverted!\n", x);
		exit(EXIT_FAILURE);
	}
	return ring->inverts[x];
}

#endif /* RING_H_ */
