// matrix.h	 Created by Burak BORHAN,   2002
// Definition of class Matrix
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

using std::ostream;

class Matrix {
	
	friend ostream &operator<<( ostream &, const Matrix & );						   
	friend Matrix operator+( const double, const Matrix & ); // Matrix as rvalue;		
	friend Matrix operator-( const double, const Matrix & ); // Matrix as rvalue;		 
	friend Matrix operator*( const double, const Matrix & ); // Matrix as rvalue;		 
																					 
public:
	Matrix( int = 1, int = 1 );	// constructor
	Matrix( const Matrix & );	// copy constructor
	~Matrix();					// destructor

	// copy a single pointer array in to the matrix
	const Matrix &copyArray( int r, int c, const double *array );  

	int getSize() const;		// get size of the matrix	
	int rowsCount() const;		// get number of rows		
	int colsCount() const;		// get number of coloumns	

	double getMaxVal() const;	// return the max value of the matrix	
	double getMaxColVal( const int c ) const; // return the max value of the cth col.	 
    double getMaxRowVal( const int r ) const; // return the max value of the rth row.	
	
	double getMinVal() const; // return the min value of the matrix		
	double getMinColVal( const int c ) const; // return the min value of the cth col   
	double getMinRowVal( const int r ) const; // return the min value of the rth row   
	
	double getMean() const;		// return the mean of the values of the matrix	  
	double getColMean( const int c ) const; // return the mean of the values of the cth col.  
	double getRowMean( const int r ) const; // return the mean of the values of the rth row	  

	// find the value in the matrix and return as a nx2 matrix
	const Matrix findVal( const double v ) const;				 
	int getValCount( const double v ) const;	// return the number of value in the matrix	   

	bool isSquare() const;		// return true if the matrix is square, false o/w. 

	// operator overleading.........

	const Matrix &operator=( const Matrix & );	// assign matrix's		  
	const Matrix &operator=( const char );		// I icin				  
	
	// Overloaded comparison operators 
	bool operator==( const Matrix & ) const;	// compare equal		  
	
	// Determine if two matrix's are equal and
	// return true, otherwise return false (uses operator==).
	bool operator!=( const Matrix &right ) const					 
	{ return ! ( *this == right ); }
	
	// Overloaded subscript operators
	double &operator()( int, int );				// subscript operator		 
	const double &operator()( int, int ) const;	// subscript operator		 
	
	// Overloaded math operators
	
	Matrix operator+( const Matrix & );	// add two Matrices					
	Matrix operator+( const double );	// add a value to a matrix			

	// add two Matrices modify this											 
	void operator+=( const Matrix & );
	// add a value to a matrix modify this									 
	void operator+=( const double );  
	
	Matrix operator-( const Matrix & ); // substract one Matrix from an other  
	Matrix operator-( const double );	// substract a value from a matrix	   

    // substract one Matrix from an other modify this matrix				  
	void operator-=( const Matrix & ); 										  
	// substract a value from a matrix modify this matrix					  
	void operator-=( const double );  
	
	Matrix operator*( const Matrix & );	// multiply a matrix with an other	   
	Matrix operator*( const double );										 
																			  
	// multiply a matrix with an other modify this matrix
	void operator*=( const Matrix & );
	// multiply a matrix with a value modify this matrix					   
	void operator*=( const double );

	Matrix operator/( const double );	// divide the matrix by a value;	   

	// divide the matrix by a value modify the this matrix					   
	void operator/=( const double );

	Matrix operator^( const int );		// take a power of the matrix

	// methods of the Matrix class
	
	const Matrix &convToIdentity(); // convert matrix to identity matrix of the same size	  
	const Matrix &convToInvDiagonal(); // convert matrix to inverse diagonal				  

	// get the transpose of a matrix
	Matrix &getTranspose();				// 
    // calculate the LUD of the matrix
	void LUBKS( int phy_size, int *outvect, double *output );
	// calculate the LUD of the matrix
	int LUD( int phy_size, int *outvect, int output );
	// get the inverse of the matrix and return the matrix.
	Matrix &getInverse();										 
	// extract a submatrix from a matrix
	Matrix getSubMatrix( const int startRow, const int startCol,
							const int rowSize, const int colSize  ) const;	   

	Matrix &deleteRow( const int r );	// Delete a row from the matrix		    
	Matrix &deleteCol( const int c );	// Delete a coloumn from the matrix	    
	// inserts the row matrix to the given row shifting other rows down
	Matrix &insertRow( const int r, const Matrix &rVector );  					
	// inserts the coloumn matrix to the given coloumn shifting other coloumns right
	Matrix &insertCol( const int c, const Matrix &cVector ); 				
	// interchange rows in the matrix ( r1 <-> r2 )
	Matrix &interchangeRows( const int r1, const int r2 );	 				
	// interchange columns in the matrix ( c1 <-> c2 )
	Matrix &interchangeCols( const int c1, const int c2 );	 					

	// add the right matrix on top of this
	Matrix &attachTop( const Matrix & );					 
	// add the right matrix under this
	Matrix &attachBottom( const Matrix & );				 
	// add the right matrix to Left of this
	Matrix &attachLeft( const Matrix & );					 
	// add the right matrix to Right of this
	Matrix &attachRight( const Matrix & );				


	static int getMatrixCount();		// matrix's instantiated  
		 
private:
	int rows;	// number of rows
	int cols;	// number of coloumns
	double *ptr;		// pointer to the first element of matrix
	static int matrixCount;	// # of matrix's instantiated
	
};

#endif	// MATRIX_H
