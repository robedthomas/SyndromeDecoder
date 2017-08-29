/**
@file SyndromeDecoding.c
@author Rob Thomas
@brief This file contains code for correcting errors in messages using syndrome 
decoding. All codes are presumed to be modulo 2.
*/

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Matrix.h"
#include "Queue.h"

#include "SyndromeDecoding.h"



int main (int argc, char *argv[])
{
	/*** Prompt the user for the values of n (the size in bytes of each message)
	     and k (the number of vectors making up the code space's basis). ***/
	int n, k;
	while ( !getNK(&n, &k) )
	{}

	/*** Receive from the user the k vectors of size n that make up the basis of
	     the codeword space. ***/
	Matrix *basis;

	initBasis(&basis, k, n);
	while ( !getBasis(basis, k, n) )
	{}

	puts("\nGenerator Matrix G:");
	Mx_print(*basis);

	/*** Determine the parity check matrix of the codeword space from its
	     generator matrix. ***/
	Matrix *H;
	getParityCheckMatrix(&H, basis, k, n);

	puts("\nParity Check Matrix H:");
	Mx_print(*H);

	/*** The codeword space C has 2^(n-k) cosets, each of which can be 
	     represented by its corresponding syndrome. Beginning with vectors of 
	     the lowest weight, find the syndrome of each vector that is NOT in C.
	     If it is a syndrome that has not been seen yet, set that vector to be
	     the coset leader for that syndrome. Continue until all 2^(n-k) 
	     cosets have been assigned coset leaders. ***/

	/* First, assemble a list of all of the (non-zero) vectors in (F_2)^n. */
	Queue *AllVectors = Q_Init(1, sizeof(Matrix *), false);
	getAllVectors(AllVectors, n);

	/* Now, find the syndrome of each of these vectors. A syndrome of 0 means 
	   the vector is a codeword. Whenever a new non-zero syndrome is found, set
	   the vector that produced it to be the coset leader associated with that
	   syndrome. */
	int numCosets = getNumCosets(k, n);

	Matrix *CosetLeaders[numCosets];

	findCosetLeaders(CosetLeaders, numCosets, AllVectors, *H);

	/* @DEBUG: print coset leaders and their syndromes */
	/*Matrix *syn;
	printf("NumCosets: %d\n", numCosets);
	for ( int i = 0; i < numCosets; i++ )
	{
		printf("Coset leader:\n");
		Mx_print(*CosetLeaders[i]);

		syn = syndrome(*CosetLeaders[i], *H);
		printf("Syndrome:\n");
		Mx_print(*syn);

		Mx_FreeMatrix(syn);
	} */

	/*** Enter the decoding loop: first, prompt the user for a vector. Find the
	     syndrome of that vector. The coset leader corresponding to that 
	     syndrome is the error vector e. The decoded message is the sum of the
	     user's vector and the error vector. Repeat until the user chooses to
	     stop. ***/
	puts("\t---TO EXIT: enter a zero vector of length n---");

	handleInput(CosetLeaders, *H, n);

	/* Once the user is done, free all matrices and exit. */
	puts("Goodbye...");

	cleanup(basis, H, AllVectors, CosetLeaders, numCosets);
}

/**
@fn getNK
@brief Prompts the user for the values of n (the size of each message) and k
(the number of vectors in the codeword subspace's basis).
@param n Pointer to where n should be stored.
@param k Pointer to where k should be stored.
@return True if n and k were successfully submitted and n >= k. False otherwise.
*/
bool getNK (int *n, int *k)
{
	char buffer[MAX_VECTOR_LENGTH];

	printf("Value of n (the size of each message): ");
	fgets(buffer, MAX_VECTOR_LENGTH, stdin);
	sscanf(buffer, "%d", n);

	printf("Value of k (number of vectors in code space basis): ");
	fgets(buffer, MAX_VECTOR_LENGTH, stdin);
	sscanf(buffer, "%d", k);

	if ( *n < *k )
	{
		fprintf(stderr, "ERROR: k must be less than n.\n");

		return false;
	}

	return true;
}

/**
@fn initBasis
@brief Allocates new matrices to fill the basis of the codeword space.
@param basis Double pointer which will become an array of Matrix pointers.
@param k The number of vectors in the codeword space's basis.
@param n The length of each vector in the codeword space's basis. 
*/
void initBasis (Matrix **basis, int k, int n)
{
	*basis = Mx_NewMatrix(NULL, k, n);
}

/**
@fn getBasis
@brief Prompts the user to enter k vectors that will fill the basis of the 
codeword space. These vectors must be in order such that the first k columns of
the basis yield the k x k identity matrix.
@param basis The list of Matrix pointers that represents the basis.
@param k The number of vectors in the codeword subspace's basis.
@param n The size of each vector in the codeword subspace's basis.
@return True if the user successfully enters a generator matrix in standard 
form. False otherwise.
*/
bool getBasis (Matrix *basis, int k, int n)
{
	/* Notify the user what kind of vectors to submit. */
	printf("Please enter %d vectors of length %d which make up the generator matrix.\n", k, n);
	printf("(Elements of vectors must be space-delimited)\n");

	char nextVector[MAX_VECTOR_LENGTH];
	char *nextToken;
	bool tokenized = false;

	/* Prompt for the next vector in the basis. */
	for ( int i = 1; i <= k; i++ )
	{
		printf("%d: ", i);

		/* Read the user's next vector. */
		fgets(nextVector, MAX_VECTOR_LENGTH, stdin);

		tokenized = false;

		/* Parse the vector from string form into individual elements. */
		for ( int j = 0; j < n; j++ )
		{
			/* Get the next token, which will be the next element in the 
			   current vector. */
			if ( !tokenized )
			{
				nextToken = strtok(nextVector, DELIMITERS);

				tokenized = true;
			}
			else
			{
				nextToken = strtok(NULL, DELIMITERS);
			}

			/* Save this element into the corresponding basis vector. */
			sscanf(nextToken, "%d", &(basis->elems[i - 1][j]));

			/* Make sure the vector has only elements of 0 and 1. */
			if ( basis->elems[i - 1][j] < 0 || basis->elems[i - 1][j] >= MODULUS )
			{
				fprintf(stderr, "ERROR: Vectors must be modulo %d.\n\n", MODULUS);

				return false;
			}
		}
	}

	/* Check the first k columns of the basis. If they make up the k x k 
	   identity matrix, then the user has successfully entered a generator
	   matrix. */
	Matrix *m = Mx_submatrix(*basis, 0, k - 1, 0, k - 1);

	if ( !Mx_isIdentityMatrix(*m) )
	{
		fprintf(stderr, "ERROR: The first k = %d columns of the basis must be the identity matrix.\n\n", k);

		return false;
	}

	return true;
}

/**
@fn getParityCheckMatrix
@brief Determines the parity check matrix of a codeword space from its generator
matrix.
@param H Double pointer indicating where H will be stored. H should NOT have 
been allocated for at the time of calling this function.
@param G Pointer to the generator matrix. Must be in standard form.
@param k The number of vectors in G.
@param n The size of each vector in G.
*/
void getParityCheckMatrix (Matrix **H, Matrix *G, int k, int n)
{
	/* Obtain the last n-k columns of G, which contain the submatrix P. */
	Matrix *P = Mx_submatrix(*G, 0, k - 1, k, n - 1);

	/* Obtain the transverse of P. */
	Matrix *P_T = Mx_transverse(*P);

	/* To obtain H, append the (n-k) x (n-k) identity matrix to the end of 
	   the transverse of P. */
	Matrix *I = Mx_getIdentityMatrix(n - k);
	*H = Mx_append(*P_T, *I);

	/* Free the matrices produced here (except H). */
	Mx_FreeMatrix(P);
	Mx_FreeMatrix(P_T);
	Mx_FreeMatrix(I);
}

/**
@fn getAllVectors
@brief Fills the given queue with all of the non-zero vectors in (F_2)^n in 
order of increasing weight.
@param list Pointer to the queue which will be filled with vectors.
@param n The size of each vector.
*/
void getAllVectors (Queue *list, int n)
{
	/* Initialize another queue which will be used for preserving the order of
	   the vectors. */
	Queue *queue = Q_Init(1, sizeof(Pair), false);

	/* First, add the zero vector to the queue. */
	Matrix *zeroVector = Mx_NewMatrix(NULL, 1, n);

	Pair *zeroPair = (Pair *)malloc( sizeof(Pair) );
	zeroPair->matrix = zeroVector;
	zeroPair->activeIndex = 0;

	Q_enqueue(queue, zeroPair);

	free(zeroPair);

	/* As long as the queue is not empty, remove the next pair from it and add
	   each of that vector's children to the queue. */
	Pair p;

	while ( queue->NumItems > 0 )
	{
		/* Dequeue the next pair. */
		Q_dequeue(queue, &p);

		/* Add this pair's children to the queue and list of vectors. */
		addNextWeight(list, queue, p.matrix, p.activeIndex);
	}

	/* Finally, free the queue of pairs. */
	Q_Free(queue);
}

/**
@fn addNextWeight
@brief Adds all children of the given vector/index pair to the list of vectors
and queue of pairs.
@param list The list to be filled with all vectors in the linear space in order 
of increasing weight.
@param queue The queue of vector/index pairs used for maintaining the order of
the vector list.
@param V The vector whose children will be found.
@param activeIndex The index at which adding of 1's to V will begin.
*/
void addNextWeight (Queue *list, Queue *queue, Matrix *V, int activeIndex)
{
	Matrix *nextV;
	Pair nextPair;

	for ( int i = activeIndex; i < V->cols; i++ )
	{
		/* Create a copy of the current matrix. */
		nextV = Mx_CopyMatrix(V);

		/* Add a 1 to the copy at the active index. */
		nextV->elems[0][i] = 1;

		/* Add the copy to both the vector list and the queue of pairs. */
		Q_enqueue(list, &nextV);

		nextPair.matrix = nextV;
		nextPair.activeIndex = i + 1;

		Q_enqueue(queue, &nextPair);
	}
}

/**
@fn getNumCosets
@brief Determines the number of cosets present for the given values of n and k.
There are 2^(n-k) - 1 cosets of any [n, k] code.
@param k The number of vectors in the basis of the codeword space.
@param n The size of each vectors in the basis.
@return The number of cosets of the codeword space.
*/
int getNumCosets (int k, int n)
{
	int c = 1;

	for ( int i = 0; i < n - k; i++ )
	{
		c *= 2;
	}

	return c;
}

/**
@fn findCosetLeaders
@brief Finds the coset leader of each coset of the codeword space.
@param CosetLeaders An array of Matrix pointers to be filled with coset leaders.
@param numCosets The number of cosets of the codeword space.
@param AllVectors A queue filled with all of the vectors in the linear subspace
(F_2)^n (excluding the 0 vector since it must be a codeword).
@param H The parity check matrix of the codeword space.
*/
void findCosetLeaders (Matrix **CosetLeaders, int numCosets, Queue *AllVectors, 
					   Matrix H)
{
	/* First, fill the CosetLeaders list with empty Matrices to indicate that
	   no coset leaders have been found. */
	for ( int i = 0; i < numCosets; i++ )
	{
		CosetLeaders[i] = NULL;
	}

	Matrix *nextVector;
	Matrix *nextSyndrome;
	int index;

	while ( AllVectors->NumItems > 0 )
	{
		/* Get the next vector from the queue. */
		Q_dequeue(AllVectors, &nextVector);

		/* Find the syndrome of this vector. */
		nextSyndrome = syndrome(*nextVector, H);

		/* Get the index corresponding to this vector. */
		index = indexOfSyndrome(*nextSyndrome);

		/* If there is not yet a coset leader for this syndrome, set this vector 
		   to be the coset leader. */
		if ( !CosetLeaders[index] )
		{
			CosetLeaders[index] = Mx_CopyMatrix(nextVector);
		}

		/* Free the vector and syndrome. */
		Mx_FreeMatrix(nextVector);
		Mx_FreeMatrix(nextSyndrome);
	}	
}

/**
@fn syndrome
@brief Finds the syndrome of a given vector.
@param vector The vector to find the syndrome of.
@param H The parity check matrix of the codeword space.
@return A pointer to a newly allocated Matrix struct containing the syndrome of
vector.
*/
Matrix *syndrome (Matrix vector, Matrix H)
{
	Matrix *vector_T = Mx_transverse(vector);

	Matrix *syn = Mx_multiplyMod(H, *vector_T, MODULUS);

	Mx_FreeMatrix(vector_T);

	/* @DEBUG */
	/*printf("Syndrome of vector:\n");
	Mx_print(vector);
	Mx_print(*syn);
	puts("");*/

	return syn;
}

/**
@fn indexOfSyndrome
@brief Finds the index corresponding to the given syndrome.
@param s The syndrome vector.
@return The index in the list of coset leaders corresponding to the given
syndrome.
*/
int indexOfSyndrome (Matrix s)
{
	int index = 0;

	for ( int i = 0; i < s.rows; i++ )
	{
		if ( s.elems[i][0] == 1 )
		{
			int x = 1;
			for ( int j = 0; j < i; j++ )
			{
				x *= 2;
			}

			index += x;
		}
	}

	return index;
}

/**
@fn handleInput
@brief Handles the user's input vectors and outputs these vectors' closest
codewords.
@param CosetLeaders The list of all coset leaders.
@param H The parity check matrix of the codeword space.
@param n The size of each codeword.
*/
void handleInput ( Matrix **CosetLeaders, Matrix H, int n )
{
	Matrix *userVector, *codeword, *syn, *cosetLeader;
	bool userWantsToExit = false;

	do
	{
		/* Prompt the user to enter an erroneous vector which will be decoded. */
		userVector = getVectorFromUser(n);

		/* If the vector is 0, then the user wants to exit. */
		if ( Mx_isZeroMatrix(*userVector) )
		{
			userWantsToExit = true;
		}
		/* Otherwise, decode the vector by finding its syndrome. */
		else
		{
			/* Get the vector's syndrome. */
			syn = syndrome(*userVector, H);

			/* If the syndrome is 0, then the vector is already a codeword. */
			if ( Mx_isZeroMatrix(*syn) )
			{
				puts("This vector is already a codeword!\n");
			}
			else
			{
				/* Add to the user's vector the coset leader corresponding to that 
		       syndrome. The result is the the decoded vector. */
				cosetLeader = CosetLeaders[indexOfSyndrome(*syn)];

				codeword = Mx_addMod(*userVector, *cosetLeader, MODULUS);

				puts("Codeword:");
				Mx_print(*codeword);
				puts("");

				/* Free the codeword. */
				Mx_FreeMatrix(codeword);
			}

			/* Free the syndrome. */
			Mx_FreeMatrix(syn);
			
		}

		/* Free the user's vector. */
		Mx_FreeMatrix(userVector);

	} while ( !userWantsToExit );
}

/**
@fn getVectorFromUser
@brief Prompts the user to enter a (space-delimited) binary vector to be 
decoded. Reads this vector and returns it as a newly allocated Matrix struct.
@param n The length of each vector in the codeword space.
@return A pointer to a newly allocated matrix struct containing the vector 
entered by the user.
*/
Matrix *getVectorFromUser (int n)
{
	char VectorString[MAX_VECTOR_LENGTH];

	/* Prompt the user to enter a vector of length n. */
	printf("Enter a vector of length %d to be decoded: ", n);

	/* Read the vector from stdin. */
	fgets(VectorString, MAX_VECTOR_LENGTH, stdin);

	char *nextToken;
	Matrix *Vector = Mx_NewMatrix(NULL, 1, n);
	bool tokenized = false;

	/* Parse the vector from string form into individual elements. */
	for ( int j = 0; j < n; j++ )
	{
		/* Get the next token, which will be the next element in the 
		   current vector. */
		if ( !tokenized )
		{
			nextToken = strtok(VectorString, DELIMITERS);

			tokenized = true;
		}
		else
		{
			nextToken = strtok(NULL, DELIMITERS);
		}

		/* Save this element into the corresponding basis vector. */
		sscanf(nextToken, "%d", &(Vector->elems[0][j]));

		/* Make sure the vector has only elements of 0 and 1. */
		if ( Vector->elems[0][j] < 0 || Vector->elems[0][j] >= MODULUS )
		{
			fprintf(stderr, "ERROR: Vectors must be modulo %d.\n\n", MODULUS);

			return NULL;
		}
	}

	return Vector;
}

/**
@fn cleanup
@brief Frees all dynamically allocated memory.
@param basis Pointer to the generator matrix.
@param H Pointer to the parity check matrix.
@param AllVectors The queue of all vectors in the linear space.
@param CosetLeaders Array of pointers to coset leader vectors.
@param numCosets The size of the CosetLeaders array.
*/
void cleanup (Matrix *basis, Matrix *H, Queue *AllVectors, 
			  Matrix **CosetLeaders, int numCosets)
{
	Mx_FreeMatrix(basis);
	Mx_FreeMatrix(H);

	Q_Free(AllVectors);

	for ( int i = 0; i < numCosets; i++ )
	{
		free(CosetLeaders[i]);
	}
}