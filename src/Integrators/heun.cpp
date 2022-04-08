/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Integrators/heun.h"


using namespace smartuq;
using namespace smartuq::integrator;

template < class T >
heun<T>::heun(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Heun's method of order 2", dyn)
{
}

template < class T >
heun<T>::~heun(){

}

template < class T >
int heun<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

        // sanity checks
        if(ti<0 || tend<0)
                smart_throw(m_name+": initial time and final time must be greater or equal to 0");
        if(tend<ti)
                smart_throw(m_name+": final time must be greater than initial time");

	xfinal.clear();

	std::size_t size = x0.size();
	std::vector<T> dx, dx_temp;
	std::vector<T> x(x0), x_temp(x0);

	double h = (tend-ti)/nsteps;

    for(int i=0; i<nsteps; i++){
		m_dyn->evaluate(ti+i*h, x, dx);

		for(std::size_t j=0; j<size; j++){
		    x_temp[j] = x[j] + h*dx[j];
		}

		m_dyn->evaluate(ti+(i+1)*h, x_temp, dx_temp);
		for(std::size_t j=0; j<size; j++){
		    x[j] += h/2.0*(dx[j] + dx_temp[j]);
		}

	}

	for(std::size_t i=0; i<x0.size(); i++)
	    xfinal.push_back(x[i]);

	return 0;
}


template class heun<double>;
template class heun<float>;
template class heun<long double>;
template class heun<polynomial::chebyshev_polynomial<double> >;
template class heun<polynomial::chebyshev_polynomial<float> >;
template class heun<polynomial::chebyshev_polynomial<long double> >;
template class heun<polynomial::taylor_polynomial<double> >;
template class heun<polynomial::taylor_polynomial<float> >;
template class heun<polynomial::taylor_polynomial<long double> >;


