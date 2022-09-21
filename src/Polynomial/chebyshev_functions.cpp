/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#include "../../include/Polynomial/chebyshev_functions.h"

namespace smartuq::polynomial {

    // Import math functions into namespace to fix overloading issues
    using std::sin;
    using std::cos;
    using std::atan2;
    using std::sqrt;
    using std::pow;

    /************************************************/
    /*                 ATAN2                        */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> atan2(const chebyshev_polynomial<T> &y, const chebyshev_polynomial<T> &x){

        T tix = x.get_coeffs()[0];
        T tiy = y.get_coeffs()[0];

        T titheta= atan2(tiy,tix);

        T tixy = sqrt(tix*tix+tiy*tiy);
        T sinti = tiy/tixy;
        T costi = tix/tixy;

        chebyshev_polynomial<T> xx = costi*x+sinti*y;
        chebyshev_polynomial<T> yy = -sinti*x+costi*y;

        return titheta+atan(yy/xx);

    }
    template class chebyshev_polynomial<double>
    atan2(const chebyshev_polynomial<double> &, const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    atan2(const chebyshev_polynomial<float> &, const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    atan2(const chebyshev_polynomial<long double> &, const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  SIN                         */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> sin(const chebyshev_polynomial<T> &other){

        return chebyshev_polynomial<T>::approximation(std::sin,other);

    }
    template class chebyshev_polynomial<double>
    sin(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    sin(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    sin(const chebyshev_polynomial<long double> &);


    /************************************************/
    /*                  COS                         */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> cos(const chebyshev_polynomial<T> &other){

        return chebyshev_polynomial<T>::approximation(std::cos,other);

    }
    template class chebyshev_polynomial<double>
    cos(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    cos(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    cos(const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  TAN                         */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> tan(const chebyshev_polynomial<T> &other){

        return chebyshev_polynomial<T>::approximation(std::tan,other);

    }
    template class chebyshev_polynomial<double>
    tan(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    tan(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    tan(const chebyshev_polynomial<long double> &);


    /************************************************/
    /*                  ASIN                        */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> asin(const chebyshev_polynomial<T> &other){
        std::vector<T> range = other.get_range();
        range[0]=std::max(range[0],(T) -1.0);
        range[1]=std::min(range[1],(T) 1.0);

        return chebyshev_polynomial<T>::approximation(std::asin,other);

    }
    template class chebyshev_polynomial<double>
    asin(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    asin(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    asin(const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  ACOS                        */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> acos(const chebyshev_polynomial<T> &other){
        std::vector<T> range = other.get_range();
        range[0]=std::max(range[0],(T) -1.0);
        range[1]=std::min(range[1],(T) 1.0);

        return chebyshev_polynomial<T>::approximation(std::acos,other);

    }
    template class chebyshev_polynomial<double>
    acos(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    acos(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    acos(const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  ATAN                        */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> atan(const chebyshev_polynomial<T> &other){

        return chebyshev_polynomial<T>::approximation(std::atan,other);

    }
    template class chebyshev_polynomial<double>
    atan(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    atan(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    atan(const chebyshev_polynomial<long double> &);

    // OTHERS
    //EXPONENTIAL FUNCTION
    template <class T>
    /************************************************/
    /*                  EXP                         */
    /************************************************/
    chebyshev_polynomial<T> exp(const chebyshev_polynomial<T> &other){

        return chebyshev_polynomial<T>::approximation(std::exp,other);

    }
    template class chebyshev_polynomial<double>
    exp(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    exp(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    exp(const chebyshev_polynomial<long double> &);

    //TANGENT HYPERBOLIC FUNCTION
    template <class T>
    /************************************************/
    /*                  TANH                         */
    /************************************************/
    chebyshev_polynomial<T> tanh(const chebyshev_polynomial<T> &other){

        return chebyshev_polynomial<T>::approximation(std::tanh,other);

    }
    template class chebyshev_polynomial<double>
    tanh(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    tanh(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    tanh(const chebyshev_polynomial<long double> &);


    /************************************************/
    /*                  SQRT                        */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> sqrt(const chebyshev_polynomial<T> &other){

        std::vector<T> range = other.get_range();
        range[0]=std::max(range[0],(T) (ZERO*ZERO));

        return chebyshev_polynomial<T>::approximation(std::sqrt, other, range);

    }
    template class chebyshev_polynomial<double>
    sqrt(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    sqrt(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    sqrt(const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  LOG                         */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> log(const chebyshev_polynomial<T> &other){

        std::vector<T> range = other.get_range();
        range[0]=std::max(range[0],(T) (ZERO*ZERO));

        return chebyshev_polynomial<T>::approximation(std::log, other, range);

    }
    template class chebyshev_polynomial<double>
    log(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    log(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    log(const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  LOG10                       */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> log10(const chebyshev_polynomial<T> &other){

        std::vector<T> range = other.get_range();
        range[0]=std::max(range[0],(T) (ZERO*ZERO));

        return chebyshev_polynomial<T>::approximation(std::log10, other, range);

    }
    template class chebyshev_polynomial<double>
    log10(const chebyshev_polynomial<double> &);
    template class chebyshev_polynomial<float>
    log10(const chebyshev_polynomial<float> &);
    template class chebyshev_polynomial<long double>
    log10(const chebyshev_polynomial<long double> &);

    /************************************************/
    /*                  POW                         */
    /************************************************/
    template <class T>
    chebyshev_polynomial<T> pow(const chebyshev_polynomial<T> &other, const double &exponent){

        // NOTE: right now when doing pow(p(x),-k) it approximates x^(-k) and composes with p. Maybe it should approximate x^k and compose with 1/p?

        int nvar =  other.get_nvar();
        int degree = other.get_degree();
        chebyshev_polynomial<T> res(nvar,degree,other.is_monomial_base());

        // trivial case
        if (exponent == 0) {
            res = 1.0;
        // natural exponent
        } else if ( exponent == floor(fabs(exponent)) ){//fancier way of checking if integer without doing exact comparison of two doubles?
            res = other;
            for(int i=1; i<exponent; i++)
                res *= other;
        } else {
            std::vector<T> range = other.get_range();
            range[0]=std::max(range[0],(T) (ZERO*ZERO));

            static auto f = [exponent](const T x)->T {return std::pow(x, exponent);};
            T (*ff)(T x) = [](const T x)->T { return f(x); };

            res = chebyshev_polynomial<T>::approximation(ff, other, range);
        }

        return res;
    }

    template class chebyshev_polynomial<double>
    pow(const chebyshev_polynomial<double> &, const double &);
    template class chebyshev_polynomial<float>
    pow(const chebyshev_polynomial<float> &, const double &);
    template class chebyshev_polynomial<long double>
    pow(const chebyshev_polynomial<long double> &, const double &);

}