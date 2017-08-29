/**
@file Matrix.c
@author Rob Thomas
@brief This file contains implementations of basic matrix functions and 
properties.
*/

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

#include "Matrix.h"



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
Matrix *Mx_NewMatrix (Matrix *m, int rows, int cols)
{
	/* Make sure that the rows and columns given are non-negative. */
	if ( rows < 0 || cols < 0 )
	{
		return NULL;
	}

	/* If an existing Matrix was passed in, free its elements list and create
	   a new one of the new dimensions. */
	if (m)
	{
		if ( m->elems )
		{
			for ( int i = 0; i < m->rows; i++ )
			{
				free (m->elems[i]);
			}

			free(m->elems);
		}

		/* Set the Matrix's rows and columns to those passed in. */
		m->rows = rows;
		m->cols = cols;

		/* Allocate a new list of elements and initialize each element to zero. */
		m->elems = (int **)malloc( sizeof(int *) * m->rows );

		for ( int i = 0; i < m->rows; i++)
		{
			m->elems[i] = (int *)malloc( sizeof(int) * m->cols );

			for ( int j = 0; j < m->cols; j++ )
			{
				m->elems[i][j] = 0;
			}
		}

		return m;
	}

	/* Otherwise, allocate a new matrix. */
	Matrix *M = (Matrix *)malloc( sizeof(Matrix) );

	M->rows = rows;
	M->cols = cols;

	M->elems = (int **)malloc( sizeof(int *) * M->rows );

	for ( int i = 0; i < M->rows; i++)
	{
		M->elems[i] = (int *)malloc( sizeof(int) * M->cols );

		for ( int j = 0; j < M->cols; j++ )
		{
			M->elems[i][j] = 0;
		}
	}

	return M;
}

/**
@fn Mx_CopyMatrix
@brief Copies an existing matrix into a newly allocated Matrix struct.
@param m Pointer to the Matrix struct to be copied.
@return A pointer to a dynamically allocated Matrix struct that is a copy of m.
*/
Matrix *Mx_CopyMatrix (Matrix *m)
{
	Matrix *c = Mx_NewMatrix(NULL, m->rows, m->cols);

	for ( int i = 0; i < c->rows; i++ )
	{
		for ( int j = 0; j < c->cols; j++ )
		{
			c->elems[i][j] = m->elems[i][j];
		}
	}

	return c;
}

/**
@fn Mx_FreeMatrix
@brief Frees a Matrix struct.
@param m Pointer to the Matrix struct to be free.
*/
void Mx_FreeMatrix (Matrix *m)
{
	for ( int i = 0; i < m->rows; i++ )
	{
		free(m->elems[i]);
	}

	free(m->elems);
	free(m);
}

/**
@fn Mx_haveEqualDimensions
@brief Reports whether or not two matrices have the same number of rows and 
columns.
@param a One of the two matrices.
@param b The other of the two matrices.
@return True if these matrices have the same number of rows and the same number
of columns. Otherwise, returns false.
*/
bool Mx_haveEqualDimensions (Matrix a, Matrix b)
{
	if ( a.rows == b.rows && a.cols == b.cols )
	{
		return true;
	}

	return false;
}

/**
@fn Mx_applyModulus
@brief Applies a modulus to each element of a matrix.
@param m Pointer to the matrix which the modulus will be applied to.
@param modulus The modulus to be applied to matrix m.
*/
void Mx_applyModulus (Matrix *m, int modulus)
{
	/* Apply modulus to each of m's elements. */
	for ( int i = 0; i < m->rows; i++ )
	{
		for ( int j = 0; j < m->cols; j++ )
		{
			m->elems[i][j] = m->elems[i][j] % modulus;
		}
	}
}

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
Matrix *Mx_add (Matrix a, Matrix b)
{
	/* Ensure that a and b have the same dimensions. */
	if ( !Mx_haveEqualDimensions(a, b) )
	{
		return NULL;
	}

	/* Allocate a new matrix to hold the sum of a and b. */
	Matrix *sumMatrix = Mx_NewMatrix(NULL, a.rows, a.cols);

	/* Sum a and b together into the new matrix. */
	for ( int i = 0; i < a.rows; i++ )
	{
		for ( int j = 0; j < a.cols; j++ )
		{
			sumMatrix->elems[i][j] = a.elems[i][j] + b.elems[i][j];
		}
	}

	return sumMatrix;
}

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
Matrix *Mx_addMod (Matrix a, Matrix b, int modulus)
{
	/* Add the given matrices together. */
	Matrix *sumMatrix = Mx_add(a, b);

	/* If the two matrices had different dimensions, return NULL. */
	if ( !sumMatrix )
	{
		return NULL;
	}

	/* Otherwise, apply the modulus to the sum matrix and return it. */
	Mx_applyModulus(sumMatrix, modulus);

	return sumMatrix;
}

/**
@fn Mx_multiply
@brief Multiplies two matrices together.
@param a The first of the two matrices to multiply together (ORDER MATTERS!).
@param b The second of the two matrices to multiply together.
@return If a and b have compatible dimensions, a new matrix will be allocated and
returned containing the matrix product of a and b. Otherwise, a NULL pointer 
will be returned.
*/
Matrix *Mx_multiply (Matrix a, Matrix b)
{
	/* Make sure a's number of columns is equal to b's number of rows. */
	if ( a.cols != b.rows )
	{
		return NULL;
	}

	/* Allocate a new CxD matrix, where C is a's number of rows and D is b's 
	   number of columns. */
	Matrix *product = Mx_NewMatrix(NULL, a.rows, b.cols);

	/* Fill the product matrices' elements using row x column multiplication. */
	int sum;

	for ( int i = 0; i < product->rows; i++ )
	{
		for ( int j = 0; j < product->cols; j++ )
		{
			sum = 0;

			for ( int x = 0; x < a.cols; x++ )
			{
				sum += a.elems[i][x] * b.elems[x][j];
			}

			product->elems[i][j] = sum;
		}
	}

	return product;
}

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
Matrix *Mx_multiplyMod (Matrix a, Matrix b, int modulus)
{
	/* Multiply the two matrices together. */
	Matrix *product = Mx_multiply(a, b);

	/* Check if a and b have compatible rows. */
	if ( !product )
	{
		return NULL;
	}

	/* Apply the modulus to the matrix product. */
	Mx_applyModulus(product, modulus);

	return product;
}

/**
@fn Mx_transverse
@brief Finds the transverse of a matrix.
@param m The matrix whose transverse will be found.
@return A newly allocated matrix that is the transverse of matrix m. If m had
improper dimensions, a NULL pointer will be returned instead.
*/
Matrix *Mx_transverse (Matrix m)
{
	/* Make sure m doesn't have improper (non-positive) dimensions. */
	if ( m.rows <= 0 || m.cols <= 0 )
	{
		return NULL;
	}

	/* Allocate a new matrix. */
	Matrix *t = Mx_NewMatrix(NULL, m.cols, m.rows);

	/* Set t's rows to be m's columns and t's columns to be m's rows. */
	for ( int i = 0; i < t->rows; i++ )
	{
		for ( int j = 0; j < t->cols; j++ )
		{
			t->elems[i][j] = m.elems[j][i];
		}
	}

	return t;
}

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
					  int endCol)
{
	if ( startRow > endRow || startCol > endCol )
	{
		return NULL;
	}

	if ( endRow >= m.rows || endCol >= m.cols )
	{
		return NULL;
	}

	/* Determine the dimensions of the submatrix. */
	int numRows = endRow - startRow + 1;
	int numCols = endCol - startCol + 1;

	/* Allocate a new matrix with these dimensions. */
	Matrix *sub = Mx_NewMatrix(NULL, numRows, numCols);

	/* Fill this matrix's elements. */
	for ( int i = 0; i < sub->rows; i++ )
	{
		for ( int j = 0; j < sub->cols; j++ )
		{
			sub->elems[i][j] = m.elems[i + startRow][j + startCol];
		}
	}

	return sub;
}

/**
@fn Mx_append
@brief Appends one matrix onto another.
@param a The left portion of the new matrix.
@param b The right portion of the new matrix.
@return A pointer to a newly allocated Matrix struct whose left side is the
elements of the matrix a and whose right side is the elements of the matrix b.
A NULL pointer if a and b do not have the same number of rows.
*/
Matrix *Mx_append (Matrix a, Matrix b)
{
	if ( a.rows <= 0 || a.cols <= 0 || b.rows <= 0 || b.cols <= 0 )
	{
		return NULL;
	}

	if ( a.rows != b.rows )
	{
		return NULL;
	}

	/* Determine the dimensions of the new matrix. */
	int numRows = a.rows;
	int numCols = a.cols + b.cols;

	/* Allocate a new matrix with these dimensions. */
	Matrix *m = Mx_NewMatrix(NULL, numRows, numCols);

	/* Fill this matrix's elements. The left side of the matrix will have a's
	   elements, while the right side of the matrix will have b's elements. */
	for ( int i = 0; i < m->rows; i++ )
	{
		for ( int j = 0; j < a.cols; j++ )
		{
			m->elems[i][j] = a.elems[i][j];
		}

		for ( int j = 0; j < b.cols; j++ )
		{
			m->elems[i][j + a.cols] = b.elems[i][j];
		}
	}

	return m;
}

/**
@fn Mx_isZeroMatrix
@brief Determines whether or not a matrix is a zero matrix.
@param m The matrix to be tested.
@return True if each of m's elements is 0. False otherwise.
*/
bool Mx_isZeroMatrix ( Matrix m )
{
	for ( int i = 0; i < m.rows; i++ )
	{
		for ( int j = 0; j < m.cols; j++ )
		{
			if ( m.elems[i][j] != 0 )
			{
				return false;
			}
		}
	} 

	return true;
}

/**
@fn Mx_isIdentityMatrix
@brief Determines whether or not a matrix is an identity matrix.
@param m The matrix to be examined.
@return True if m is an identity matrix (of any size). False otherwise.
*/
bool Mx_isIdentityMatrix (Matrix m)
{
	/* Check if m has equal numbers of rows and columns. If not, it cannot be an
	   identity matrix. */
	if ( m.rows != m.cols )
	{
		return false;
	}

	/* Make sure diagonal elements are 1 and non-diagonal elements are 0. */
	for ( int i = 0; i < m.rows; i++ )
	{
		for ( int j = 0; j < m.cols; j++ )
		{
			if ( j == i && m.elems[i][j] != 1 )
			{
				/* @DEBUG: */
				//printf("non-1 diagonal element found at (%d, %d).\n", i, j);
				return false;
			}
			else if ( j != i && m.elems[i][j] != 0 )
			{
				/* @DEBUG: */
				//printf("non-0 element found at (%d, %d).\n", i, j);
				return false;
			}
		}
	}

	return true;
}

/**
@fn Mx_getIdentityMatrix
@brief Yields a newly allocated identity matrix of the given size.
@param size The dimension of the identity matrix to be produced.
@return A pointer to a newly allocated identity matrix of dimensions 
size x size.
*/
Matrix *Mx_getIdentityMatrix (int size)
{
	Matrix *I = Mx_NewMatrix(NULL, size, size);

	for ( int i = 0; i < size; i++ )
	{
		I->elems[i][i] = 1;
	}

	return I;
}

/**
@fn Mx_print
@brief Prints a matrix to stdout.
@param m The Matri struct to be printed.
*/
void Mx_print (Matrix m)
{
	for ( int i = 0; i < m.rows; i++ )
	{
		/* Print the left cap. If this is the first and/or last row, print {. */
		if ( i == 0 || i == m.rows - 1 )
		{
			printf("{ ");
		}
		else
		{
			printf("( ");
		}

		for ( int j = 0 ; j < m.cols; j++ )
		{
			/* Print the next element. */
			printf("%d ", m.elems[i][j]);
		}

		/* Print the right cap. If this is the first and/or last row, print }. */
		if ( i == 0 || i == m.rows - 1 )
		{
			printf("}\n");
		}
		else
		{
			printf(")\n");
		}
	}
}