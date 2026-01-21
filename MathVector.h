/*
 * MathVector.h
 *
 *      Author: Theodore Faust based upon course material provided in UCLA Math 280 (Prof. Chris Anderson).
 *
 *
 *
 */

#include <vector>
#include <iostream>
#include <initializer_list>
#include <stdexcept> // For standard exceptions
#include <sstream>   // For stringstream
#include <cmath>

#ifndef MATHVECTOR_H_
#define MATHVECTOR_H_

class MathVector
{

public:

    /// Null constructor
    MathVector()                                : vData()
    {
    }

    /// Copy constructor
    MathVector(const MathVector& V)             : vData(V.vData)
    {
    }

    /// Creates vector of size vectorSize with all elements initialized to val (which defaults to 0 if not specified)
    MathVector(size_t vectorSize, double val = 0.0)   : vData(vectorSize,val)
    {
    }

    /// Values list initializer
    MathVector(std::initializer_list<double>  init)  : vData(init)
    {
    };

    /// Destructor
    virtual ~MathVector()
    {};


    /// Returns size of the MathVector
    size_t size() const
    {
        return vData.size();
    }


    /// Resizes to the vector to newSize.
    void resize(size_t newSize)
    {
    	(this->vData).resize(newSize);
    }

    /// Resizes the vector to newSize and if newSize > current size, initializes new values with val
    void resize(size_t newSize, double val)
    {
    	(this->vData).resize(newSize,val);
    }

    /// Resizes the vector to size 0
    void clear()
    {
    	(this->vData).resize(0);
    }

    /// Returns reference to element value at location index (indexing starts at 0).
    inline double& operator()(long index)
    {
//    	std::string  errorString;
//    	if (indexError(index,errorString))
//    	{
//    		throw std::out_of_range(errorString);
//    	}
//    	else
//    	{
        return vData[index];
//    	}
    }

    /// Returns reference to element value at location index (indexing starts at 0).
    inline const double& operator()(long index) const
    {
//    	std::string  errorString;
//    	if (indexError(index,errorString))
//    	{
//    		throw std::out_of_range(errorString);
//    	}
//    	else
//    	{
        return vData[index];
//    	}
    }


    ///  (*this) = V and returns (*this). The size of (*this) returned is adjusted to that of V.
    MathVector& operator=(const MathVector& V)
    {
    	std::string  errorString;
    	if (this->size() != 0 && sizeError(this->size(),V.size(),errorString))
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
    MathVector operator+(const MathVector& V) const
    {
    	std::string  errorString;
    	if (sizeError(this->size(),V.size(),errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
        MathVector R(*this);
        for(size_t i = 0; i  < this->size(); ++i)
        {
            R.vData[i] += V.vData[i];
        }
        return R;
    	}
    }

    ///  Returns  (*this) - V
    MathVector operator-(const MathVector& V) const
    {
    	std::string  errorString;
    	if (sizeError(this->size(),V.size(),errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
    	MathVector R(*this);
    	for(size_t i = 0; i  < this->size(); ++i)
    	{
    		R.vData[i] -= V.vData[i];
    	}
    	return R;
    	}
    }

    ///  Incremental vector addition : (*this) = (*this) +  V
    void operator += (const MathVector& V)
    {
    	std::string  errorString;
    	if (sizeError(this->size(),V.size(),errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
            for(size_t i = 0; i  < this->size(); ++i)
            {
                this->vData[i] += V.vData[i];
            }
    	}

    }

    ///  Incremental vector subtraction : (*this) = (*this) -  V
    void operator -= (const MathVector& V)
    {
    	std::string  errorString;
    	if (sizeError(this->size(),V.size(),errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
        	for(size_t i = 0; i  < this->size(); ++i)
        	{
        		this->vData[i] -= V.vData[i];
        	}
    	}

    }

    ///  Incremental scalar multiplication : (*this) = (*this)*alpha
    void operator *= (double alpha)
    {
        for(size_t i = 0; i  < this->size(); ++i)
        {
            this->vData[i] *= alpha;
        }
    }

    ///  Incremental scalar division : (*this) = (*this)/alpha
    void operator /= (double alpha)
    {
    	for(size_t i = 0; i  < this->size(); ++i)
    	{
    		this->vData[i] /= alpha;
    	}
    }


    /// Returns V*alpha
    MathVector operator*(double alpha) const
    {
        MathVector R(*this);
        for(size_t i = 0; i  < this->size(); ++i)
        {
            R.vData[i] *= alpha;
        }
        return R;
    }

    /// Returns  alpha*V
    friend MathVector operator*(double alpha, const MathVector& B)
    {
        MathVector R(B);
        for(size_t i = 0; i  < B.size(); ++i)
        {
            R.vData[i] *= alpha;
        }
        return R;
    }

    /// Returns V/alpha
    MathVector operator/(double alpha) const
    {
    	MathVector R(*this);
    	for(size_t i = 0; i  < this->size(); ++i)
    	{
    		R.vData[i] /= alpha;
    	}
    	return R;
    }


    /// Sets row display flag (true = horizontal output, false = vertical output)

    static void setRowDisplay(bool flag)
    {
    	rowDisplay = flag;
    }

    /// Overload of << for output to ostream
    friend std::ostream& operator<<(std::ostream& outStream, const MathVector& V)
    {
    	if (V.rowDisplay)
		{
    		 for(size_t i = 0; i <  V.size(); ++i)
    		 {
    			 outStream <<  V(i) << " ";
    		 }
    		 return outStream;
		}
    	else
    	{
    		for(size_t i = 0; i <  V.size()-1; ++i)
    		{
    			outStream << V(i) << "\n";
    		}
    		outStream << V(V.size()-1);
    		return outStream;
    	}
    }

    /// Computes the (standard) dot product of two MathVectors of the same size
    double dot(const MathVector& V) const
    {
    	std::string  errorString;
    	if (sizeError(this->size(),V.size(),errorString))
    	{
    		throw std::invalid_argument(errorString);
    	}
    	else
    	{
    		double total = 0.0;
    		for (size_t i = 0; i < V.size(); i++)
    		{
    			total += (V.vData)[i]*(this->vData)[i];
    		}
    		return total;
    	}
    }

    /// Null initializer
    void initialize()
    {
    	this->vData.clear();
    }

    /// Copy initializer
    void initialize(const MathVector& V)
    {
    	this->vData.clear();
    	this->vData = V.vData;
    }

    /// Initializer of vector of size vectorSize with all elements initialized to val (defaults to 0 if not specified)
    void initialize(size_t vectorSize, double val = 0.0)
    {
    	this->vData.clear();
    	this->vData.resize(vectorSize, val);
    }


    /// Values list initializer
    void initialize(std::initializer_list<double>  init)
    {
        this->vData.clear();
        this->vData.insert(vData.end(), init.begin(), init.end());
    };

    /// Returns a pointer to the first element of the array of doubles storing the
        /// vector data.

        double* data()
        {
            return vData.data();
        }

        /// Returns a pointer to the first element of the array of doubles storing the
        /// vector data.

        const double* data() const
        {
            return vData.data();
        }


        operator double()
        {
            if(this->size() != 1)
            {throw std::invalid_argument("Invalid conversion of MathVector to double ");}
            return vData[0];
        }

        double normInf()
        {
            double maxAbs = 0.0;
            for(size_t i = 0; i < this->size(); ++i)
            {
                maxAbs = (maxAbs > std::abs(vData[i])) ? maxAbs : vData[i];
            }
            return maxAbs;
        }

        double norm2()
        {
            double norm2sum = 0.0;
            for(size_t i = 0; i < this->size(); ++i)
            {
                norm2sum  += vData[i]*vData[i];
            }
            return std::sqrt(std::abs(norm2sum));
        }

        double norm1()
        {
            double norm1sum = 0.0;
            for(size_t i = 0; i < this->size(); ++i)
            {
                norm1sum  += std::abs(vData[i]);
            }
            return norm1sum;
        }

    private:

    std::vector<double> vData;
    static bool           rowDisplay;

    /*! Checks for equivalence of left and right operands.
            If not equal : An errorString is set to an error message
                           and the function returns true.

            If equal      : returns false.
        */
        bool sizeError(long leftSize, long rghtSize, std::string& errorString) const
        {
            bool returnFlag = false;

            std::stringstream errorStream;    // Use a string stream to compose error message
               if(leftSize != rghtSize)
            {
             errorStream   << "MathVector Operand Error"  << std::endl;
             errorStream   << "Left operand size : "      << leftSize  << std::endl;
             errorStream   << "Rght operand size : "      << rghtSize  << std::endl;

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

//        bool indexError(long index, std::string& errorString) const
//        {
//            bool returnFlag = false;
//            std::stringstream errorStream;    // Use a string stream to compose error message
//            if((index < 0) || (index >= (long)this->size()))
//            {
//                 errorStream   << "MathVector Index Error" << std::endl;
//                 errorStream   << "Offending index   : "      << index << std::endl;
//                 errorStream   << "Acceptable range  : [0, "  << this->size()-1 << "] "<< std::endl;
//
//                 errorString = errorStream.str();
//                 returnFlag = true;
//                }
//             return returnFlag;
//        }
};

bool MathVector::rowDisplay = true;

#endif /* MATHVECTOR_H_ */
