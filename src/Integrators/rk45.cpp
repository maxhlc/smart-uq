/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
---------------- Copyright (C) 2022 Imperial College London ----------------
-------------- E-mail: m.hallgarten-la-casta21@imperial.ac.uk --------------
-------------------- Author: Max I. Hallgarten La Casta --------------------
*/

#include "../../include/Integrators/rk45.h"

namespace smartuq::integrator {

    template<class T>
    rk45<T>::rk45(const dynamics::base_dynamics<T>* dyn, const double atol, const double rtol) : base_integrator<T>("Runge Kutta Fehlberg 4(5) variable step time", dyn), m_atol(atol), m_rtol(rtol) {

    }

    template<class T>
    rk45<T>::~rk45() {

    }

    template <class T>
    int rk45<T>::integrate(const double& ti, const double& tend, const int& nsteps, const std::vector<T>& x0, std::vector<T>& xfinal) const {
        // Number of state variables
        std::size_t n = x0.size();

        // Declare temporary vectors
        std::vector<T> x(x0), xtemp(x0), k1, k2, k3, k4, k5, k6;

        // Set current time to initial time
        double t = ti;

        // Calculate initial step size, and declare new step size variable
        double h = (tend-ti)/nsteps;
        double hnew;

        // Declare truncation error variable
        double TE;

        // Declare integration status flags
        bool isIntegrationComplete = false, isToleranceSatisfied = true, isStepComplete = false;

        while(!isIntegrationComplete){
            // Ensure step remains within bounds
            if((t + h) > tend)
                h = tend - t;

            // Check if integration has finished (h -> 0)
            if(h < 1e-15)
                break;

            while(!isStepComplete){
                // Evaluate times
                double t1 = t;
                double t2 = t + 2.0/9.0 * h;
                double t3 = t + 1.0/3.0 * h;
                double t4 = t + 3.0/4.0 * h;
                double t5 = t + 1.0/1.0 * h;
                double t6 = t + 5.0/6.0 * h;

                // Evaluate k1
                for(std::size_t ii=0; ii<n; ii++)
                    xtemp[ii] = x[ii];
                m_dyn->evaluate(t1, xtemp, k1);

                // Evaluate k2
                for(std::size_t ii=0; ii<n; ii++)
                    xtemp[ii] = x[ii] + 2.0/9.0*k1[ii]*h;
                m_dyn->evaluate(t2, xtemp, k2);

                // Evaluate k3
                for(std::size_t ii=0; ii<n; ii++)
                    xtemp[ii] = x[ii] + (1.0/12.0*k1[ii] + 1.0/4.0*k2[ii])*h;
                m_dyn->evaluate(t3, xtemp, k3);

                // Evaluate k4
                for(std::size_t ii=0; ii<n; ii++)
                    xtemp[ii] = x[ii] + (69.0/128.0*k1[ii] - 243.0/128.0*k2[ii] + 135.0/64.0*k3[ii])*h;
                m_dyn->evaluate(t4, xtemp, k4);

                // Evaluate k5
                for(std::size_t ii=0; ii<n; ii++)
                    xtemp[ii] = x[ii] + (-17.0/12.0*k1[ii] + 27.0/4.0*k2[ii] - 27.0/5.0*k3[ii] + 16.0/15.0*k4[ii])*h;
                m_dyn->evaluate(t5, xtemp, k5);

                // Evaluate k6
                for(std::size_t ii=0; ii<n; ii++)
                    xtemp[ii] = x[ii] + (65.0/432.0*k1[ii] - 5.0/16.0*k2[ii] + 13.0/16.0*k3[ii] + 4.0/27.0*k4[ii] + 5.0/144.0*k5[ii])*h;
                m_dyn->evaluate(t6, xtemp, k6);

                // Evaluate truncation error
                TE = 0.0;
                for(std::size_t ii=0; ii<n; ii++){
                    TE += std::pow((
                        -  1.0/150.0 * k1[ii].get_coeffs()[0]
                     // +  0.0/0.0   * k2[ii].get_coeffs()[0]
                        +  3.0/100.0 * k3[ii].get_coeffs()[0]
                        - 16.0/75.0  * k4[ii].get_coeffs()[0]
                        -  1.0/20.0  * k5[ii].get_coeffs()[0]
                        +  6.0/25.0  * k6[ii].get_coeffs()[0]
                        )*h,2);
                }
                TE = std::sqrt(TE);

                // Calculate new step size
                hnew = 0.9*h*std::pow(m_atol/TE, 0.2);

                // Check convergence
                isToleranceSatisfied = (TE <= m_atol);

                // Terminate step if converged
                if(isToleranceSatisfied){
                    isStepComplete = true;
                    isToleranceSatisfied = true;
                } else {
                    h = hnew;
                }
            }

            // Evaluate weighted average
            for(std::size_t ii=0; ii<n; ii++)
                x[ii] = x[ii] + (
                        + 47.0/450.0 * k1[ii]
                     // +  0.0/0.0   * k2[ii]
                        + 12.0/25.0  * k3[ii]
                        + 32.0/225.0 * k4[ii]
                        +  1.0/30.0  * k5[ii]
                        +  6.0/25.0  * k6[ii]
                        )*h;

            // Update time
            t += h;

            // Update next step size
            h = hnew;

            // Set flag for next step
            isStepComplete = false;     
        }

        // Update final state
        for(std::size_t ii=0; ii<n; ii++)
            xfinal[ii] = x[ii];

        // Return zero
        return 0;
    }

    template class rk45<polynomial::chebyshev_polynomial<double>>;
    template class rk45<polynomial::chebyshev_polynomial<float>>;
    template class rk45<polynomial::chebyshev_polynomial<long double>>;
    template class rk45<polynomial::taylor_polynomial<double>>;
    template class rk45<polynomial::taylor_polynomial<float>>;
    template class rk45<polynomial::taylor_polynomial<long double>>;

}
