/*
 * MathMatrix.h
 *
 *      Author: Theodore Faust based upon course material provided in UCLA Math 280 (Prof. Chris Anderson).
 */

#include <vector>
#include <iostream>
#include <initializer_list>
#include <stdexcept> // For standard exceptions
#include <sstream>   // For stringstream
#include <cmath>
#include <limits>
#include "MathVector.h"

#ifndef MATHMATRIX_H_
#define MATHMATRIX_H_

class MathMatrix
{

public:

    /// Null constructor
    MathMatrix()                                : matrixData()
    {
    	this->rows = 0;
    	this->cols = 0;
    }

    /// Copy constructor
    MathMatrix(const MathMatrix& V)             : matrixData(V.matrixData)
    {
    	this->rows = V.rows;
    	this->cols = V.cols;
    }

    /// Creates matrix of size rows x cols with entries set to val (defaults to 0.0)
    MathMatrix(size_t rows, size_t cols, double val = 0.0)               : matrixData(rows*cols, val)
    {
    	this->rows = rows;
    	this->cols = cols;
    }

    /*!
     Values list constructor. Values specified by { { row 0 data } {row 1 data} ... {row M data} }
     The number of columns is determined by the number of elements in the
     first row specification. If subsequent row lists do not have the
     identical number of elements than an exeption is thrown.
    */

    MathMatrix(std::initializer_list< std::initializer_list<double> > init)  : matrixData()
    {
        this->rows = init.size();
        this->cols = init.begin()->size();

        matrixData.resize(rows*cols,0.0);

        size_t iIndex = 0;
        size_t jIndex = 0;

        for (auto i =init.begin(); i != init.end(); ++i)
        {
            jIndex = 0;
            for (auto j=i->begin(); j != i->end(); ++j)
            {
                if(jIndex >= this->cols)
                {
                throw std::invalid_argument("Unacceptable list specification for MathMatrix");
                }
                this->operator()(iIndex,jIndex) = *j;
                jIndex++;
            }
            if(jIndex < this->cols)
            {
            throw std::invalid_argument("Unacceptable list specification for MathMatrix");
            }
            iIndex++;
        }
    };

    /// Destructor
    virtual ~MathMatrix()
    {};


    /// Returns number of rows in the mathMatrix
    size_t getRowSize() const
    {
        return this->rows;
    }

    /// Returns number of columns in the mathMatrix
    size_t getColSize() const
    {
    	return this->cols;
    }

    /// Returns element value at (rowIndex, colIndex) (indexing starts at 0).
    /// The matrix data is assumed to be stored by columns    (Fortran style)

    inline double& operator()(long rowIndex, long colIndex)
    {
//        std::string errorString;
//        if(indexError(rowIndex,colIndex,errorString))
//        {throw std::out_of_range(errorString);}

        return matrixData[rowIndex + (rows)*colIndex];
    }

    /// Returns element value at (rowIndex, colIndex) (indexing starts at 0).
    /// The matrix data is assumed to be stored by columns    (Fortran style)

    inline const double& operator()(long rowIndex, long colIndex) const
    {
//        std::string errorString;
//        if(indexError(rowIndex,colIndex,errorString))
//        {throw std::out_of_range(errorString);}

        return matrixData[rowIndex + (rows)*colIndex];
    }


    ///  (*this) = V and returns (*this).
    MathMatrix& operator=(const MathMatrix& V)
    {
    	std::string  errorString;
    	if (not(this->rows == 0 && this->cols == 0) && sizeError(this->rows,this->cols,V.rows,V.cols,errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
    		this->initialize(V);
    		return *this;
    	}

    }

    ///  Returns   (*this) + V
    MathMatrix operator+(const MathMatrix& V) const
    {
    	std::string  errorString;
    	if (sizeError(this->rows,this->cols,V.rows,V.cols,errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
        MathMatrix R(*this);
        for(size_t i = 0; i  < this->rows*this->cols; ++i)
        {
            R.matrixData[i] += V.matrixData[i];
        }
        return R;
    	}
    }

    ///  Returns  (*this) - V
    MathMatrix operator-(const MathMatrix& V) const
    {
    	std::string  errorString;
    	if (sizeError(this->rows,this->cols,V.rows,V.cols,errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
        MathMatrix R(*this);
        for(size_t i = 0; i  < this->rows*this->cols; ++i)
        {
            R.matrixData[i] -= V.matrixData[i];
        }
        return R;
    	}
    }

    ///  Incremental matrix addition : (*this) = (*this) +  V

    void operator += (const MathMatrix& V)
    {
    	std::string  errorString;
    	if (sizeError(this->rows,this->cols,V.rows,V.cols,errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
        for(size_t i = 0; i  < this->rows*this->cols; ++i)
        {
            this->matrixData[i] += V.matrixData[i];
        }
    	}
    }

    ///  Incremental matrix subtraction : (*this) = (*this) -  V
    void operator -= (const MathMatrix& V)
    {
    	std::string  errorString;
    	if (sizeError(this->rows,this->cols,V.rows,V.cols,errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
        for(size_t i = 0; i  < this->rows*this->cols; ++i)
        {
            this->matrixData[i] -= V.matrixData[i];
        }
    	}
    }

    ///  Incremental matrix multiplication : (*this) = (*this)*alpha
    void operator *= (double alpha)
    {
        for(size_t i = 0; i  < this->rows*this->cols; ++i)
        {
            this->matrixData[i] *= alpha;
        }
    }

    ///  Incremental matrix division : (*this) = (*this)/alpha
    void operator /= (double alpha)
    {
    	for(size_t i = 0; i  < this->rows*this->cols; ++i)
    	{
    		this->matrixData[i] /= alpha;
    	}
    }


    /// Returns V*alpha
    MathMatrix operator*(double alpha) const
    {
        MathMatrix R(*this);
        for(size_t i = 0; i  < this->rows*this->cols; ++i)
        {
            R.matrixData[i] *= alpha;
        }
        return R;
    }

    /// Returns V*v
    MathVector operator*(const MathVector& v) const
    {
    	std::string  errorString;
    	if (multArgError(this->cols, v.size(), errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
    		MathVector w(this->rows);
    		for(size_t i = 0; i  < this->rows; ++i)
    		{
    			for (size_t j = 0; j < this->cols; ++j)
    			{
    				w(i) += this->operator()(i, j) * v(j);
    			}
    		}
    		return w;
    	}
    }

    /// Returns V*T
    MathMatrix operator*(const MathMatrix& T) const
    {
    	std::string errorString;
    	if (multArgError(this->cols, T.rows, errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
    		MathMatrix W(this->rows, T.cols);
    		for(size_t k = 0; k < T.cols; ++k)
    		{
    			for(size_t i = 0; i  < this->rows; ++i)
    			{
    				for (size_t j = 0; j < this->cols; ++j)
    				{
    					W(i,k) += this->operator()(i, j) * T(j,k);
    				}
    			}
    		}
    	return W;
    	}
    }

    /// Returns  alpha*V
    friend MathMatrix operator*(double alpha, const MathMatrix& B)
    {
        MathMatrix R(B);
        for(size_t i = 0; i  < B.rows*B.cols; ++i)
        {
            R.matrixData[i] *= alpha;
        }
        return R;
    }

    /// Returns V/alpha
    MathMatrix operator/(double alpha) const
    {
    	MathMatrix R(*this);
    	for(size_t i = 0; i  < this->rows*this->cols; ++i)
    	{
    		R.matrixData[i] /= alpha;
    	}
    	return R;
    }

    /// Overload of << for output to ostream
    friend std::ostream& operator<<(std::ostream& outStream, const MathMatrix& V)
    {
    	for(size_t i = 0; i <  V.rows; ++i)
    	{
    		for(size_t j = 0; j<V.cols; ++j)
    		{
    			outStream <<  V(i,j) << " ";
    		}
    		outStream << "\n";
    	}
    	return outStream;
    }

    /// Null initializer
    void initialize()
    {
    	this->matrixData.clear();
    	this->rows = 0;
    	this->cols = 0;
    }

    /// Copy initializer
    void initialize(const MathMatrix& V)
    {
    	this->matrixData.clear();
    	this->matrixData = V.matrixData;
    	this->rows = V.rows;
    	this->cols = V.cols;
    }

    /// Initialize matrix of size rows x cols with entries set to val (defaults to 0.0)
    void initialize(size_t rows, size_t cols, double val = 0.0)
    {
    	this->matrixData.clear();
    	this->matrixData.resize(rows*cols, val);
    	this->rows = rows;
    	this->cols = cols;
    }

    /// Values list initializer
    void initialize(std::initializer_list<std::initializer_list<double> > init)
        {
            matrixData.clear();

            this->rows = init.size();
            this->cols = init.begin()->size();

            matrixData.resize(rows*cols,0.0);

            size_t iIndex = 0;
            size_t jIndex = 0;

            for (auto i =init.begin(); i != init.end(); ++i)
            {
                jIndex = 0;
                for (auto j=i->begin(); j != i->end(); ++j)
                {
                    if(jIndex >= this->cols)
                    {
                    throw std::invalid_argument("Unacceptable list specification for MathMatrix");
                    }
                    this->operator()(iIndex,jIndex) = *j;
                    jIndex++;
                }
                if(jIndex < this->cols)
                {
                throw std::invalid_argument("Unacceptable list specification for MathMatrix");
                }
                iIndex++;
            }
        }

    /// Returns a pointer to the first element of the array of doubles storing the
    /// matrix data.

    double* data()
    {
        return matrixData.data();
    }

    /// Returns a pointer to the first element of the array of doubles storing the
    /// matrix data.

    const double* data() const
    {
        return matrixData.data();
    }


    /// Conversion of 1 X 1 matrix to double

    operator double()
    {
        if((this->rows != 1)||(this->cols != 1))
        {throw std::invalid_argument("Invalid conversion of MathMatrix to double");}

        return matrixData[0];
    }


    /// Returns the 1-norm of the matrix (maximal column sum)

    double norm1() const
    {
        double colSum;
        double maxColSum = 0.0;

        for(size_t j = 0; j < cols; ++j)
        {
            colSum = 0.0;

            for(size_t i = 0; i < rows; ++i)
            {
                colSum += std::abs((*this)(i,j));
            }

            maxColSum = (maxColSum > colSum) ? maxColSum : colSum;
        }

        return maxColSum;
    }

    /// Returns the inf-norm of the matrix (maximal row sum)

    double normInf() const
    {
        double rowSum;
        double maxRowSum = 0.0;

        for(size_t i = 0; i < rows; ++i)
        {
            rowSum = 0.0;

            for(size_t j = 0; j  < cols; ++j)
            {
                rowSum += std::abs((*this)(i,j));
            }

            maxRowSum = (maxRowSum > rowSum) ? maxRowSum : rowSum;
        }

        return maxRowSum;
    }

    /// vector*matrix operation only created if the matrix has only one row

    friend MathMatrix operator*(const MathVector& v, const MathMatrix& M)
    {
        if(M.rows != 1){throw std::invalid_argument("vector*matrix operation only allowed if matrix row size is 1");}

        MathMatrix R(v.size(),M.cols);
        for(size_t i = 0; i < v.size(); ++i)
        {
            for(size_t j = 0; j < M.cols; ++j)
            {
                R(i,j) = v(i)*M(0,j);
            }
        }

        return R;
    }

/*
 *  Returns the solution to (*this) x = b or least squares solution in the case
 *  when the number of rows is greater than the number of columns.
 *
 *  If the number of rows is less than the number of columns an invalid_argument
 *  exception is thrown.
 *
 *  When invoking the inverse operator, a runtime exception  is thrown if
  * the matrix is singular or nearly singular. In such cases, the exception
 *  is caught and a warning message is output -- the program execution is
 *  not terminated.
 *
 *
*/

    MathVector applyInverse(const MathVector& b) const
    {
    std::string errorString;
    if(invArgError(rows, b.size(),errorString)) {throw std::invalid_argument(errorString);}

    long M        =   this->rows;
    long N        =   this->cols;
    long NRHS     =   1;

    // Create COPIES of input matrix and right hand side
    // to pass to the solver.

    MathMatrix Atemp(*this);
    MathVector Btemp(b);

    // Access the data pointers for each

    double* AdataPtr = Atemp.data();
    double* BdataPtr = Btemp.data();

    // Create return argument and capture data pointer

    MathVector X(N);
    double* XdataPtr = X.data();

    // Call the solver. Catch ill-condition system as a runtime error.
    // Print out a warning message, rather than terminate the program.

    try
    {
    applyInverse(AdataPtr, M, N, BdataPtr, NRHS ,XdataPtr);
    }
    catch(const std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
    }

    return X;
    }


    /// Returns inverse operator (assumes full rank)

    MathMatrix inverse() const
    {
        if(rows != cols)
        {throw std::invalid_argument("MathMatrix inverse not defined for non-square matrices ");}

        MathMatrix Id(rows,rows);
        for(size_t i = 0; i < rows; ++i) {Id(i,i) = 1.0;}


        // Call the solver. Catch ill-condition system as a runtime error.
        // Print out a warning message, rather than terminate the program.

        try
        {
            return applyInverse(Id);
        }
        catch(const std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
        }

        return MathMatrix(); // Return null instance to suppress warnings
    }

    /// Returns pseudo-inverse operator (assumes full column rank)

    MathMatrix pseudoInverse() const
    {
        MathMatrix Id(rows,rows);

        for(size_t i = 0; i < rows; ++i) {Id(i,i) = 1.0;}

        // Call the solver. Catch ill-condition system as a runtime error.
        // Print out a warning message, rather than terminate the program.

        try
        {
            return applyInverse(Id);
        }
        catch(const std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
        }

        return MathMatrix(); // Return null instance to suppress warnings
    }


    MathMatrix applyInverse(const MathMatrix& B) const
    {
        std::string errorString;
        if(invArgError(this->rows, B.rows,errorString)) {throw std::invalid_argument(errorString);}

        long M        =   this->rows;
        long N        =   this->cols;
        long NRHS     =   B.cols;

        // Create COPIES of input matrix and right hand side
        // to pass to the solver.

        MathMatrix Atemp(*this);
        MathMatrix Btemp(B);

        // Access the data pointers for each

        double* AdataPtr = Atemp.data();
        double* BdataPtr = Btemp.data();

        // Create return argument and capture data pointer

        MathMatrix X(N,NRHS);

        double* XdataPtr = X.data();

        // Call the solver. Catch ill-condition system as a runtime error.
        // Print out a warning message, rather than terminate the program.

        try
        {
            applyInverse(AdataPtr, M, N, BdataPtr, NRHS ,XdataPtr);
        }
        catch(const std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
        }


        return X;
    }

    private:

    std::vector<double> matrixData;
    size_t rows;
    size_t cols;

    //
    // This is a "low level" interface to a procedure for computing
    // solutions, or least squares solutions, to A*X = B.
    //
    // double* Adata = pointer to first element of data of an M x N matrix.
    // It is assumed that the data is stored by column (column major ordering)
    // in a single contiguous array
    //
    // !!! This routine overwrites the data in Adata
    //
    // double* Bdata =  pointer to first element of data of an M x NRHS matrix of right hand sides.
    //                  It is assumed that the data is stored by column (column major ordering)
    //                  in a single contiguous array
    //
    // !!! This routine overwrites the data in Bdata
    //
    //
    // double* Xdata = pointer to the first element of an N x NRHS matrix whose columns contain the
    //                 solution (or least squares solution) to A*X = B.
    //                 It is assumed that the data is stored by column (column major ordering)
    //                 in a single contiguous array
    //
    // No bounds checking is done.
    //
    // If the inverse procedure is successful and A is not detected as being singular
    // or nearly singular, then the return value is 0, otherwise the
    // the return value is non-zero.
    //
    // The solver assumes that M >= N and A is a full rank matrix. If A is singular
    // or nearly singular, then a runtime_error exception is thrown. Generally this
    // will be caught and the error message printed, rather than terminating the
    // program.
    //

    int applyInverse(double* AdataPtr, long M, long N, double* BdataPtr, long NRHS, double* XdataPtr) const
    {

        std::numeric_limits<double>   numLimits;
        double local_machine_eps = 10.0*numLimits.epsilon();

        int returnValue = 0;

        long m     =   M;
        long n     =   N;
        long p     =   NRHS;

        // Create required temporaries

        std::vector<double> S(m,0.0);
        std::vector<double> C(m,0.0);
        std::vector<double> T(m,0.0);

        // Capture data pointers

        double *S_data = S.data();
        double *C_data = C.data();
        double *T_data = T.data();

        //
        // The following code determines solutions and least squares solutions
        // to linear systems of equations using a QR factorization.  It was created
        // circa 1996 to supply functionality to member functions that determine solutions
        // (and least squares solutions) of linear systems without the need to use
        // an externally defined linear algebra routines, i.e. LAPACK. The heavy use of
        // pointers in the implementation was done to improve computational efficiency.
        // In retrospect, this was not a good idea since it renders the implementation
        // unintelligible. In all likelihood, were the original implementation available
        // and compiled using today's optimizing compilers, the result would be faster
        // than this implementation.
        //
        // Since the procedure is based upon QR factorization, the method is
        // robust and can create both solutions of square systems of equations and
        // least square solutions of non-square systems. For square systems, there is a
        // corresponding increase in computational cost over creatig and using an LU factorization,
        // but the goal was reliability rather than ultimate computational efficiency.
        // (For the latter one should use LAPACK routines).
        //
        //  The procedure assumes the matrix data is stored by columns.
        //
        //  The QR procedure and estimation of the condition number is based
        //  upon material from an early version of "Matrix Computations"
        //  by G. Golub and C. Van Loan. On ones to-do list would be a
        //  new version based upon QR with column pivoting, such as that
        //  described in section 5.4  of "Matrix Computations"
        //  by G. Golub and C. Van Loan 4th edition.
        //

        double tau;

        double* Aptr; double* Sptr;
        double* Cptr; double* Tptr;
        double* Xptr;

        register double* Top;
        register double* piv_elem;

        for(long k = 1;  k <= n; k++)
        {
          piv_elem = AdataPtr + (k-1)*m + (k-1);
          for( Aptr  = piv_elem + 1, Sptr = S_data, Cptr = C_data;
               Aptr <= piv_elem + (m-k); Aptr++, Sptr++, Cptr++)
          {
            if( *Aptr == 0.0 ){ *Cptr = 1.0; *Sptr =0.0;}
            else
            { if(std::abs(*Aptr) > std::abs(*piv_elem) )
              { tau   = -(*piv_elem)/(*Aptr);
               *Sptr  = 1.0/std::sqrt(1.0 + tau*tau);
                *Cptr = (*Sptr)*tau;
              }
              else
              { tau   = -(*Aptr)/(*piv_elem);
                *Cptr = 1.0/std::sqrt(1.0 + tau*tau);
                *Sptr = (*Cptr)*tau;
              }
            }
           *piv_elem = ((*Cptr) * (*piv_elem)) - ((*Sptr) * (*Aptr));
          }

          for(long j = k+1; j <= n; j++)
          {
           Top = AdataPtr + (j-1)*m + (k-1);
           for(Aptr = Top + 1, Sptr = S_data, Cptr = C_data, Tptr = T_data;
               Aptr <= Top + (m - k); Aptr++, Sptr++, Cptr++, Tptr++)
            {  *Tptr = (*Sptr) * (*Top);
               *Top  = ((*Cptr) * (*Top)) - ((*Sptr) * (*Aptr));
               *Aptr *= (*Cptr);
            }

           for(Aptr = Top + 1, Tptr = T_data; Aptr <= Top + (m - k);  Aptr++, Tptr++)
               { *Aptr += *Tptr;}
          }
        //
        //    Transform the right hand side
        //
          for(long j = 1; j <= p ; j++)
          {
             Top = BdataPtr + (j-1)*m + (k-1);
             for( Xptr = Top + 1, Sptr = S_data, Cptr = C_data, Tptr = T_data;
              Xptr <= Top + (m - k); Xptr++, Sptr++, Cptr++, Tptr++)
              { *Tptr = (*Sptr) * (*Top);
                *Top  = ((*Cptr) * (*Top)) - ((*Sptr) * (*Xptr));
                *Xptr *= (*Cptr);
              }

             for( Xptr  = Top + 1, Tptr = T_data; Xptr <= Top + (m - k);
              Xptr++, Tptr++)
              { *Xptr += *Tptr; }
          }

        }
        //
        //  Estimate the condition number
        //
        double R_norm = 0.0;
        double R_col_norm;
        for(long k = 1; k <= n; k++)
        {  R_col_norm = 0.0;
           for( Aptr = AdataPtr + (k-1)*m; Aptr < AdataPtr + (k-1)*m + k; Aptr++ )
           {R_col_norm += std::abs(*Aptr);}
           if(R_norm < R_col_norm ) R_norm = R_col_norm;
         }

        long singular_flag = 0;
        double       RCOND = 1.0;

        for(long j=1; j <= n ; j++)
        {
            if(std::abs(*(AdataPtr + (j-1) + (j-1)*m)) <= local_machine_eps*R_norm )
            {
             singular_flag = 1;
             if(RCOND > std::abs(*(AdataPtr + (j-1) + (j-1)*m))/R_norm)
             {
                 RCOND = std::abs(*(AdataPtr + (j-1) + (j-1)*m))/R_norm;
             }
            }
        }

        //
        //  Back Substitute
        //
        register double XJ;
        for(long k = 1; k <= p; k++)
        {
          for(long j=  n; j >= 2; j--)
          {
           XJ = (*(BdataPtr +(k-1)*m + (j-1))) / (*(AdataPtr + (j-1)*m + (j-1)));
           *(BdataPtr +(k-1)*m + (j-1)) = XJ;
            for( Xptr = BdataPtr + (k-1)*m, Aptr = AdataPtr + (j-1)*m;
             Xptr < BdataPtr + (k-1)*m + (j-1); Xptr++, Aptr++)
            *Xptr -= XJ*(*Aptr);
          }
        *(BdataPtr + (k-1)*m) = (*(BdataPtr + (k-1)*m))/(*(AdataPtr));
        }

        for(long i = 0; i < n; i++)
        {
        for(long j = 0; j < p; j++)
        {
         *(XdataPtr + i + j*n) = *(BdataPtr + i + j*m);
        }}

        // Throw a runtime_error if during the solution process the
        // condition number estimate of the system indicates a singular
        // or nearly singular system of equations
        //
        // Catch this error to output a warning message instead
        // of terminating the program

        std::stringstream errorStream;
        if(singular_flag == 1)
        {
            returnValue = 1;
            errorStream  << "XXX Matrix system singular or close to singular XXX "   << std::endl;
            errorStream  << "XXX      Computed results may be inaccurate     XXX"    << std::endl;
            errorStream  << "XXX Estimated reciprocal condition number (RCOND) : "   << RCOND << std::endl;
            throw std::runtime_error(errorStream.str());
        }

        return returnValue;
    }

    /*! Checks for equivalence of left and right matrix operand dimensions
        If not equal : An errorString is set to an error message
                       and the function returns true.

        If equal      : returns false.
    */
    bool sizeError(long leftRowSize, long leftColSize,
                   long rghtRowSize, long rghtColSize, std::string& errorString) const
    {
        bool returnFlag = false;

        std::stringstream errorStream;    // Use a std::string stream to compose error message
           if((leftRowSize != rghtRowSize) || (leftColSize != rghtColSize))
        {
         errorStream   << "MathMatrix Operand Error" << std::endl;
         errorStream   << "Left operand size : "      << leftRowSize << " X " << leftColSize  << std::endl;
         errorStream   << "Rght operand size : "      << rghtRowSize << " X " << rghtColSize  << std::endl;

         errorString = errorStream.str();
         returnFlag = true;
            }

        return returnFlag;
    }

    /*! Checks for out of range index error.
        If out of range  : An errorString is set to an error message
                           and the function returns true.

        Not out of range : returns false.
    */

//    bool indexError(long i, long j, std::string& errorString) const
//    {
//        bool returnFlag = false;
//        std::stringstream errorStream;    // Use a std::string stream to compose error message
//
//        if(((i < 0) || (i >= (long)this->rows)) || ((j < 0) || (j >= (long)this->cols)))
//        {
//             errorStream   << "MathMatrix Index Error" << std::endl;
//             errorStream   << "Offending index       : ("    << i << ", " << j << ") "<< std::endl;
//             errorStream   << "Acceptable row index range : [0, "  << this->rows-1 << "] "<< std::endl;
//             errorStream   << "Acceptable col index range : [0, "  << this->cols-1 << "] "<< std::endl;
//
//             errorString = errorStream.str();
//             returnFlag = true;
//            }
//         return returnFlag;
//    }

    /*! Checks for appropriate matrix multiplication operands left columns = right rows
        If not equal : An errorString is set to an error message
                       and the function returns true.

        If equal      : returns false.
    */

    bool multArgError(long leftColSize, long rghtRowSize,  std::string& errorString) const
    {
        bool returnFlag = false;

        std::stringstream errorStream;    // Use a std::string stream to compose error message
           if(leftColSize != rghtRowSize)
        {
         errorStream   << "MathMatrix Multiplication Error" << std::endl;
         errorStream   << "Left operand cols : "      << leftColSize <<  std::endl;
         errorStream   << "Rght operand rows : "      << rghtRowSize <<  std::endl;

         errorString = errorStream.str();
         returnFlag = true;
            }

        return returnFlag;
    }


    /*! Checks for appropriate matrix inverse/least squares solution operands
        left rows >=  right rows
        If not satisfied : An errorString is set to an error message
                           and the function returns true.

        If equal          : returns false.
    */

    bool invArgError(long leftRowSize, long rghtRowSize,  std::string& errorString) const
    {
        bool returnFlag = false;

        std::stringstream errorStream;    // Use a std::string stream to compose error message
        if(leftRowSize < rghtRowSize)
        {
         errorStream   << "MathMatrix Inverse/Least Squares Error" << std::endl;
         errorStream   << "Left operand rows : "      << leftRowSize <<  std::endl;
         errorStream   << "Rght operand rows : "      << rghtRowSize <<  std::endl;

         errorString = errorStream.str();
         returnFlag = true;
         }

        return returnFlag;
    }
};



#endif /* MATHMATRIX_H_ */
