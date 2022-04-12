/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Integrators/euler.h"


namespace smartuq::integrator {

	template < class T >
	euler<T>::euler(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Forward Euler integration scheme", dyn)
	{
	}

	template < class T >
	euler<T>::~euler(){

	}

	template < class T >
	int euler<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

		// sanity checks
			if(ti<0 || tend<0)
					smart_throw(m_name+": initial time and final time must be greater or equal to 0");
			if(tend<ti)
					smart_throw(m_name+": final time must be greater than initial time");

		xfinal.clear();

		std::vector<T> dx;
		std::vector<T> x(x0);

		double h = (tend-ti)/nsteps;

		for(int i=0; i<nsteps; i++){
			m_dyn->evaluate(ti+i*h, x, dx);
			for(std::size_t j=0; j<x.size(); j++){
				x[j] += h*dx[j];
			}
		}

		for(std::size_t i=0; i<x0.size(); i++)
			xfinal.push_back(x[i]);

		return 0;
	}


	template class euler<double>;
	template class euler<float>;
	template class euler<long double>;
	template class euler<polynomial::chebyshev_polynomial<double> >;
	template class euler<polynomial::chebyshev_polynomial<float> >;
	template class euler<polynomial::chebyshev_polynomial<long double> >;
	template class euler<polynomial::taylor_polynomial<double> >;
	template class euler<polynomial::taylor_polynomial<float> >;
	template class euler<polynomial::taylor_polynomial<long double> >;

}