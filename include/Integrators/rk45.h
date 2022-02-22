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

#ifndef SMARTUQ_RK45_H
#define SMARTUQ_RK45_H

#include "base_integrator.h"

namespace smartuq::integrator {
    
    template<class T> class base_integrator;
    /**
     * @brief The Runge-Kutta-Fehlberg 4(5) integrator scheme.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-22
     * 
     */
    template<class T>
    class rk45 : public base_integrator<T> {

        private:

            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

            const double m_tol;

        public:

            /**
             * @brief rk45 constructor.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-02-22
             * 
             * @param[in] dyn Dynamics model
             * @param[in] tol Solver tolerance
             */
            rk45(const dynamics::base_dynamics<T>* dyn, const double tol);

            /**
             * @brief rk45 desctructor.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-02-22
             * 
             */
            ~rk45();

            /**
             * @brief Integrate method to integrate between two given times with an initial step length and initial state.
             * 
             * The method implements the RK4(5) method. The number of steps is used to calculate the initial step, and maintained
             * in this format for compatibility with the other integrators.
             * 
             * @param[in] ti Initial time.
             * @param[in] tend End time.
             * @param[in] nsteps Number of integration steps (used to calculate initial step size).
             * @param[in] x0 Initial state.
             * @param[out] xfinal Final state.
             * @return int 
             */
            int integrate(const double& ti, const double& tend, const int& nsteps, const std::vector<T>& x0, std::vector<T>& xfinal) const;

    };

}

#endif