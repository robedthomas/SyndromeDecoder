/**
@file Matrix.h
@author Rob Thomas
@brief This file contains implementations of basic matrix functions and 
properties.
*/

#ifndef MATRIX_H
#define MATRIX_H



/**
@typedef Matrix
@brief A struct representing a single matrix of elements.
*/
typedef struct 
{
	int rows;
	int cols;

	int **elems;
} Matrix;



/**
@fn Mx_NewMatrix
@brief Initializes a matrix's dimensions to those given and elements to zero.
If a pointer to an already existing Matrix struct was passed in, simply 
resets this Matrix's data members. Otherwise, allocates a new Matrix structs
and initializes its data members.
@param m A pointer to the Matrix struct to be initialized. If NULL, a Matrix
struct will be allocated for.
@param rows The number of rows this matrix will have.
@param cols The number of columns this matrix will have.
@return If a pointer to a Matrix was passed in, then that pointer is returned.
Otherwise, a pointer to the newly allocated Matrix is returned.
*/
Matrix *Mx_NewMatrix (Matrix *m, int rows, int cols);

/**
@fn Mx_CopyMatrix
@brief Copies an existing matrix into a newly allocated Matrix struct.
@param m Pointer to the Matrix struct to be copied.
@return A pointer to a dynamically allocated Matrix struct that is a copy of m.
*/
Matrix *Mx_CopyMatrix (Matrix *m);

/**
@fn Mx_FreeMatrix
@brief Frees a Matrix struct.
@param m Pointer to the Matrix struct to be free.
*/
void Mx_FreeMatrix (Matrix *m);

/**
@fn Mx_haveEqualDimensions
@brief Reports whether or not two matrices have the same number of rows and 
columns.
@param a One of the two matrices.
@param b The other of the two matrices.
@return True if these matrices have the same number of rows and the same number
of columns. Otherwise, returns false.
*/
bool Mx_haveEqualDimensions (Matrix a, Matrix b);

/**
@fn Mx_applyModulus
@brief Applies a modulus to each element of a matrix.
@param m Pointer to the matrix which the modulus will be applied to.
@param modulus The modulus to be applied to matrix m.
*/
void Mx_applyModulus (Matrix *m, int modulus);

/**
@fn Mx_add
@brief Adds two matrices of the same dimensions.
@param a One of the two matrices to add.
@param b The other matrix to add.
@return If a and b have the same dimensions, a new Matrix struct will be 
allocated with the dimensions of a and b whose elements will be the sum of a's
elements and b's elements. A pointer to this allocated Matrix will be returned.
If a and b do NOT have the same dimensions, a NULL pointer will be returned.
*/
Matrix *Mx_add (Matrix a, Matrix b);

/**
@fn Mx_addMod
@brief Adds two matrices together, then applies a modulus to each element in the
sum matrix.
@param a One of the matrices to be summed.
@param b The other matrix to be summed.
@param modulus The modulus to be applied to the sum matrix.
@return If a and b have the same dimensions, a new Matrix struct will be 
allocated with the dimensions of a and b whose elements will be the sum of a's
elements and b's elements with the given modulus applied. A pointer to this 
allocated Matrix will be returned. If a and b do NOT have the same dimensions, a 
NULL pointer will be returned.
*/
Matrix *Mx_addMod (Matrix a, Matrix b, int modulus);

/**
@fn Mx_multiply
@brief Multiplies two matrices together.
@param a The first of the two matrices to multiply together (ORDER MATTERS!).
@param b The second of the two matrices to multiply together.
@return If a and b have compatible dimensions, a new matrix will be allocated and
returned containing the matrix product of a and b. Otherwise, a NULL pointer 
will be returned.
*/
Matrix *Mx_multiply (Matrix a, Matrix b);

/**
@fn Mx_multiplyMod
@brief Multiplies two matrices together, then applies a modulus to each element 
in the sum matrix.
@param a One of the matrices to be multiplied.
@param b The other matrix to be multiplied.
@param modulus The modulus to be applied to the product matrix.
@return If a and b have the same dimensions, a new Matrix struct will be 
allocated with the dimensions of a and b. This matrix will be filled with the 
product of a and b. A pointer to this allocated Matrix will be returned. If a 
and b do NOT have compatible dimensions, a NULL pointer will be returned.
*/
Matrix *Mx_multiplyMod (Matrix a, Matrix b, int modulus);

/**
@fn Mx_transverse
@brief Finds the transverse of a matrix.
@param m The matrix whose transverse will be found.
@return A newly allocated matrix that is the transverse of matrix m. If m had
improper dimensions, a NULL pointer will be returned instead.
*/
Matrix *Mx_transverse (Matrix m);

/**
@fn Mx_submatrix
@brief Copies a portion of an existing matrix into a new matrix.
@param m The matrix from which a portion will be copied.
@param startRow The first row to be copied (inclusive).
@param endRow The last row to be copied (inclusive).
@param startCol The first column to be copied (inclusive).
@param endCol The last column to be copied (inclusive).
@return A new matrix containing the elements of matrix m that are encapsulated
by the given rows and columns.
*/
Matrix *Mx_submatrix (Matrix m, int startRow, int endRow, int startCol, 
					  int endCol);

/**
@fn Mx_append
@brief Appends one matrix onto another.
@param a The left portion of the new matrix.
@param b The right portion of the new matrix.
@return A pointer to a newly allocated Matrix struct whose left side is the
elements of the matrix a and whose right side is the elements of the matrix b.
A NULL pointer if a and b do not have the same number of rows.
*/
Matrix *Mx_append (Matrix a, Matrix b);

/**
@fn Mx_isZeroMatrix
@brief Determines whether or not a matrix is a zero matrix.
@param m The matrix to be tested.
@return True if each of m's elements is 0. False otherwise.
*/
bool Mx_isZeroMatrix ( Matrix m );

/**
@fn Mx_isIdentityMatrix
@brief Determines whether or not a matrix is an identity matrix.
@param m The matrix to be examined.
@return True if m is an identity matrix (of any size). False otherwise.
*/
bool Mx_isIdentityMatrix (Matrix m);

/**
@fn Mx_getIdentityMatrix
@brief Yields a newly allocated identity matrix of the given size.
@param size The dimension of the identity matrix to be produced.
@return A pointer to a newly allocated identity matrix of dimensions 
size x size.
*/
Matrix *Mx_getIdentityMatrix (int size);

/**
@fn Mx_print
@brief Prints a matrix to stdout.
@param m The Matrix struct to be printed.
*/
void Mx_print (Matrix m);


#endif /* MATRIX_H */