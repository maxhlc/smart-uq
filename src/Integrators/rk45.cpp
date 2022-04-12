/*
MIT License

Copyright (c) 2021-2022 Max Hallgarten La Casta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
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

        std::size_t n = x0.size();
        std::vector<T> x(x0), xtemp(x0), k1, k2, k3, k4, k5, k6;

        double t = ti;

        double h = (tend-ti)/nsteps;
        double hnew;
        std::vector<double> TEvec(n);
        double TE;

        bool isIntegrationComplete = false, isStepComplete = false;

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
                for(std::size_t ii=0; ii<n; ii++)
                    TEvec[ii] = (
                                -  1.0/150.0 * k1[ii].get_coeffs()[0]
                              //+  0.0/0.0   * k2[ii].get_coeffs()[0]
                                +  3.0/100.0 * k3[ii].get_coeffs()[0]
                                - 16.0/75.0  * k4[ii].get_coeffs()[0]
                                -  1.0/20.0  * k5[ii].get_coeffs()[0]
                                +  6.0/25.0  * k6[ii].get_coeffs()[0]
                                )*h;
                TE = 0.0;
                for(std::size_t ii=0; ii<n; ii++)
                    TE += pow(TEvec[ii], 2);
                TE = sqrt(TE);

                // Calculate new step size
                hnew = 0.9*h*pow(m_atol/TE, 0.2);

                // Check convergence
                if(TE > m_atol){
                    h = hnew;
                } else {
                    isStepComplete = true;
                }
            }

            // Evaluate weighted average
            for(std::size_t ii=0; ii<n; ii++)
                x[ii] = x[ii] + (
                        + 47.0/450.0 * k1[ii]
                      //+  0.0/0.0   * k2[ii]
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
