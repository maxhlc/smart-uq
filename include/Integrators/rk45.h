/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
---------------- Copyright (C) 2022 Imperial College London ----------------
-------------- E-mail: m.hallgarten-la-casta21@imperial.ac.uk --------------
-------------------- Author: Max I. Hallgarten La Casta --------------------
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

            const double m_atol;
            const double m_rtol;

        public:

            /**
             * @brief rk45 constructor.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-02-22
             * 
             * @param[in] dyn Dynamics model
             * @param[in] atol Solver absolute tolerance
             * @param[in] rtol Solver relative tolerance
             */
            rk45(const dynamics::base_dynamics<T>* dyn, const double atol, const double rtol);

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