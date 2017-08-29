/**
@file SyndromeDecoding.h
@author Rob Thomas
@brief This file contains code for constructing tables of coset leaders and 
their associated syndromes. Whenever a table is constructed for a code, it is 
saved into a .txt file for later use. All codes are assumed to be in modulo 2.
*/

#ifndef SYNDROMEDECODING_H
#define SYNDROMEDECODING_H

#define MODULUS 2
#define MAX_VECTOR_LENGTH 256
#define DELIMITERS " \t\n\0"


/**
@typedef Pair
@brief A struct representing a single vector and where it will be modified next.
This struct is only used for assembling a list of all vectors in a linear space.
*/
typedef struct 
{
	Matrix *matrix;
	int activeIndex;
} Pair;



/**
@fn getNK
@brief Prompts the user for the values of n (the size of each message) and k
(the number of vectors in the codeword subspace's basis).
@param n Pointer to where n should be stored.
@param k Pointer to where k should be stored.
@return True if n and k were successfully submitted and n >= k. False otherwise.
*/
bool getNK (int *n, int *k);

/**
@fn initBasis
@brief Allocates new matrices to fill the basis of the codeword space.
@param basis Double pointer which will become an array of Matrix pointers.
@param k The number of vectors in the codeword space's basis.
@param n The length of each vector in the codeword space's basis. 
*/
void initBasis (Matrix **basis, int k, int n);

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
bool getBasis (Matrix *basis, int k, int n);

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
void getParityCheckMatrix (Matrix **H, Matrix *G, int k, int n);

/**
@fn getAllVectors
@brief Fills the given queue with all of the non-zero vectors in (F_2)^n in 
order of increasing weight.
@param list Pointer to the queue which will be filled with vectors.
@param n The size of each vector.
*/
void getAllVectors (Queue *list, int n);

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
void addNextWeight (Queue *list, Queue *queue, Matrix *V, int activeIndex);

/**
@fn getNumCosets
@brief Determines the number of cosets present for the given values of n and k.
There are 2^(n-k) - 1 cosets of any [n, k] code.
@param k The number of vectors in the basis of the codeword space.
@param n The size of each vectors in the basis.
@return The number of cosets of the codeword space.
*/
int getNumCosets (int k, int n);

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
					   Matrix H);

/**
@fn syndrome
@brief Finds the syndrome of a given vector.
@param vector The vector to find the syndrome of.
@param H The parity check matrix of the codeword space.
@return A pointer to a newly allocated Matrix struct containing the syndrome of
vector.
*/
Matrix *syndrome (Matrix vector, Matrix H);

/**
@fn indexOfSyndrome
@brief Finds the index corresponding to the given syndrome.
@param s The syndrome vector.
@return The index in the list of coset leaders corresponding to the given
syndrome.
*/
int indexOfSyndrome (Matrix s);

/**
@fn handleInput
@brief Handles the user's input vectors and outputs these vectors' closest
codewords.
@param CosetLeaders The list of all coset leaders.
@param H The parity check matrix of the codeword space.
@param n The size of each codeword.
*/
void handleInput ( Matrix **CosetLeaders, Matrix H, int n );

/**
@fn getVectorFromUser
@brief Prompts the user to enter a (space-delimited) binary vector to be 
decoded. Reads this vector and returns it as a newly allocated Matrix struct.
@param n The length of each vector in the codeword space.
@return A pointer to a newly allocated matrix struct containing the vector 
entered by the user.
*/
Matrix *getVectorFromUser (int n);

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
			  Matrix **CosetLeaders, int numCosets);


#endif /* SYNDROMEDECODING_H */