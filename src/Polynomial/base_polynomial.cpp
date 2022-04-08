/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Polynomial/base_polynomial.h"

using namespace smartuq;
using namespace polynomial;



/**************/
/*CONSTRUCTORS*/
/**************/
template < class T >
base_polynomial<T>::base_polynomial(const int &vars, const int &order, const std::vector<T> &a, const std::vector<T> &b): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_monomial_base(false), m_J(0), m_N(0), m_a(a), m_b(b){

    if(vars<0){
        smart_throw(m_name+": Polynomials need to have a positive number of variables");
    }
    if(order<0){
        smart_throw(m_name+": Polynomials need to have a positive order");
    }

    if(a.size()!=0 && a.size() != (std::size_t) vars)
        smart_throw("Base polynomial: variables lower bound need to have the same size of the number of variables");
    if(a.size()!=0 && b.size() != (std::size_t) vars)
        smart_throw("Base polynomial: variables upper bound need to have the same size of the number of variables");
    if(a.size()>0){
        for(int i=0; i<vars; i++){
            if(a[i]>b[i])
                smart_throw("Base polynomial: variables bounds need to be a<b");
        }
    }

    int n = combination(vars,order);

    if(n>constants::MAX_POLYNOMIAL_ALGEBRA_SIZE){
        smart_throw(m_name+": The size of the algebra is too big. Reduce polynomial order rnumber of variables. You can incur in memory issues");
    }

    m_coeffs.resize(n);

    //save some info
    m_degree = order;
    m_nvar = vars;

    m_J.resize(vars+1);
    m_N.resize(vars+1);
    for(int i=0; i<=vars; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();

    if(m_a.size()==0)
        m_a = std::vector<T>(m_nvar,-1.0);
    if(m_b.size()==0)
        m_b = std::vector<T>(m_nvar,1.0);

}

template < class T >
base_polynomial<T>::base_polynomial(const int &vars, const int &order, const int &i, const T &a, const T &b): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_monomial_base(false), m_J(0), m_N(0){

    if(vars<0){
        smart_throw(m_name+": Polynomials need to have a positive number of variables");
    }
    if(order<=0){
        smart_throw(m_name+": Polynomials need to have a positive order");
    }
    if(i<0 || i>=vars){
        smart_throw(m_name+": First order Polynomial constructor need a variable index between [0,nvars-1]");
    }

    //allocate memory for coefficients vector

    int n = combination(vars,order);

    if(n>constants::MAX_POLYNOMIAL_ALGEBRA_SIZE){
        smart_throw(m_name+": The size of the algebra is too big. Reduce polynomial order rnumber of variables. You can incur in memory issues");
    }

    m_coeffs.resize(n);
    m_coeffs[i+1] = (b-a)/2.0;
    m_coeffs[0] = (b+a)/2.0;

    //save some info
    m_degree = order;
    m_nvar = vars;

    m_J.resize(vars+1);
    m_N.resize(vars+1);
    for(int i=0; i<=vars; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();

}


template < class T >
base_polynomial<T>::base_polynomial(const int &vars, const int &order, const T &value): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_monomial_base(false), m_J(0), m_N(0){

    if(vars<0){
        smart_throw(m_name+": Polynomials need to have a positive number of variables");
    }
    if(order<=0){
        smart_throw(m_name+": Polynomials need to have a positive order");
    }

    //allocate memory for coefficients vector

    int n = combination(vars,order);

    if(n>constants::MAX_POLYNOMIAL_ALGEBRA_SIZE){
        smart_throw(m_name+": The size of the algebra is too big. Reduce polynomial order rnumber of variables. You can incur in memory issues");
    }

    m_coeffs.resize(n);
    m_coeffs[0] = value;

    //save some info
    m_degree = order;
    m_nvar = vars;

    m_J.resize(vars+1);
    m_N.resize(vars+1);
    for(int i=0; i<=vars; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();


}

template <class T>
base_polynomial<T>::~base_polynomial(){

}


template <class T>
void base_polynomial<T>::interpolation(const std::vector<std::vector<T> > &x, const std::vector<T> &y, std::vector<std::vector<T> > &H){
    if(x.size()==0)
        smart_throw(m_name+": for polynomial interpolation non empty nodal values need to be provided");
    if(x.size()!=y.size())
        smart_throw(m_name+": for polynomial interpolation, the number of nodes and nodal values need to be the same");
    if(x[0].size()!= (std::size_t) m_nvar)
        smart_throw(m_name+": the number of variables is not the same as in the algebra");

    for(unsigned int i=0;i<H.size();i++)
        H[i].clear();

    int npoints = x.size();
    int ncoeffs = m_coeffs.size();

    if(npoints<ncoeffs)
        smart_throw(m_name+": the number of interpolation point need to be equal or greater than the size of the algebra");

    Eigen::VectorXd Y(npoints);
    Eigen::MatrixXd base_matrix (npoints,ncoeffs);
    Eigen::VectorXd coe(ncoeffs);

    //building matrix H
    for(int i=0; i<npoints; i++){
        base_matrix(i,0)=1.0;
        std::vector<T> basis = evaluate_basis(x[i]);
        for (int j=1;j<ncoeffs;j++){
            base_matrix(i,j)=basis[j];
        }
        Y[i] = y[i];
    }

    // solve by inversion, faster when we want a lot of representations
    if(npoints==ncoeffs){
        Eigen::MatrixXd base_inv (npoints,ncoeffs);
        base_inv=base_matrix.inverse();
        coe = base_inv*Y;

        for(int i=0;i<npoints;i++){
            std::vector<T> row(npoints);
            for(int j=0;j<npoints;j++)
                row[j]=base_inv(i,j);
            H.push_back(row);
        }
    }
    else //solve by Least Square
        coe = base_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y);

    std::vector<T> final_coeffs(ncoeffs);
    for(int i=0; i<ncoeffs; i++) final_coeffs[i] = coe[i];

    set_coeffs(final_coeffs);


}


template <class T>
void base_polynomial<T>::interpolation(const std::vector<std::vector<T> > &x, const std::vector<std::vector<T> >  &y, std::vector<std::vector<T> > &H, std::vector<std::vector<T> > &res_coeffs) const{
    if(x.size()==0)
        smart_throw(m_name+": for polynomial interpolation non empty nodal values need to be provided");
    if(x.size()!=y.size())
        smart_throw(m_name+": for polynomial interpolation, the number of nodes and nodal values need to be the same");
    if(x[0].size()!= (std::size_t) m_nvar)
        smart_throw(m_name+": the number of variables is not the same as in the algebra");

    res_coeffs.clear();

    for(std::size_t i=0; i<H.size(); i++)
        H[i].clear();

    std::size_t npoints = x.size();
    std::size_t ncoeffs = m_coeffs.size();

    if(npoints<ncoeffs)
        smart_throw(m_name+": the number of interpolation point need to be equal or greater than the size of the algebra");

    Eigen::VectorXd Y(npoints);
    Eigen::MatrixXd base_matrix (npoints,ncoeffs);
    Eigen::VectorXd coe(ncoeffs);
    std::vector<T> coeffs(ncoeffs);

    //building matrix H
    for(std::size_t i=0; i<npoints; i++){
        base_matrix(i,0)=1.0;
        std::vector<T> basis = evaluate_basis(x[i]);
        for(std::size_t j=1; j<ncoeffs; j++){
            base_matrix(i,j)=basis[j];
        }
    }

    // solve by inversion, faster when we want a lot of representations
    if(npoints==ncoeffs){
        Eigen::MatrixXd base_inv (npoints,ncoeffs);
        base_inv=base_matrix.inverse();
        for(std::size_t i=0; i<npoints; i++){
            std::vector<T> row(npoints);
            for(std::size_t j=0; j<npoints; j++)
                row[j]=base_inv(i,j);
            H.push_back(row);
        }

        for(std::size_t i=0; i<y[0].size(); i++){
            for(std::size_t j=0; j<npoints; j++){
                Y[j] = y[j][i];
            }
            coe = base_inv*Y;
            for(std::size_t j=0; j<ncoeffs; j++)
                coeffs[j] = coe[j];
            res_coeffs.push_back(coeffs);

            for(std::size_t i=0; i<npoints; i++){
                std::vector<T> row(npoints);
                for(std::size_t j=0; j<npoints; j++)
                    row[j]=base_inv(i,j);
                H.push_back(row);
            }
        }
    }
    else{ //solve by Least Square
        for(int i=0; i<m_nvar; i++){
            for(std::size_t j=0; j<npoints; j++){
                Y[j] = y[j][i];
            }
            coe = base_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y);
            for(std::size_t j=0; j<ncoeffs; j++)
                coeffs[j] = coe[j];
            res_coeffs.push_back(coeffs);
        }
    }

    return;

}

template <class T>
void base_polynomial<T>::solve(const std::vector<std::vector<T> > &H, const std::vector<std::vector<T> >  &y, std::vector<std::vector<T> > &res_coeffs) const{


    int nrows_H = H.size();
    if(nrows_H==0)
        smart_throw(m_name+": solving linear system, the number of rows of the matrix is zero");
    int ncolumns_H = H[0].size();
    if((std::size_t) ncolumns_H!=y.size())
        smart_throw(m_name+": cannot solve linear system, matrix-vector dimensions mismatch");

    std::size_t nvars = y[0].size();

    res_coeffs.clear();

    for(std::size_t i=0; i<nvars; i++)
        res_coeffs.push_back(std::vector<T>(nrows_H,0.0));

    for(std::size_t i=0; i<nvars; i++){
        for(int j=0; j<nrows_H; j++){
            for(int k=0; k<ncolumns_H; k++){
                res_coeffs[i][j]+= H[j][k]*y[k][i];
            }
        }
    }

}


/******************************/
/*Monomial multiplication     */
/******************************/
//initialisation of static terms for multiplication and so
template <class T>
std::vector<bool> base_polynomial<T>::m_Q(1,true);
template <class T>
std::vector<int> base_polynomial<T>::m_M(1,0);
template <class T>
int base_polynomial<T>::m_Mnvar = 0;
template <class T>
int base_polynomial<T>::m_Mdegree = 0;
template <class T>
void base_polynomial<T>::monomial_multiplication(const base_polynomial<T> &x1, const base_polynomial<T> &x2, base_polynomial<T> &res_poly){

    if(!x1.m_monomial_base || !x2.m_monomial_base){
        smart_throw("One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(x1.get_nvar()!=x2.get_nvar() || x1.get_nvar()!=res_poly.get_nvar()){
        smart_throw("Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(x1.get_degree()!=x2.get_degree() || x1.get_degree()!=res_poly.get_degree()){
        smart_throw("Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int degree = x1.get_degree();
    int nvar = x1.get_nvar();
    std::vector<std::vector<int> > J = x1.get_J();
    std::vector<std::vector<int> > N = x1.get_N();

    // use M instead of searching index for faster multiplication
    bool use_M = false;
    if(nvar == m_Mnvar && degree == m_Mdegree) use_M = true;

    std::vector<T> x1_coeffs = x1.get_coeffs();
    std::vector<T> x2_coeffs = x2.get_coeffs();
    std::vector<T> res_coeffs(x1_coeffs.size());
    int i_0, j_0, idx_0; //index offsets
    int idx;
    int coeff=0;


    for (int deg0=0; deg0<=degree; deg0++){ //loop over order of terms of poly0
        int max_i=J[nvar][deg0];
        if (deg0==0) i_0=0;
        else i_0=N[nvar][deg0-1];

        for (int deg1=0; deg1<=degree-deg0; deg1++){ //loop over other of terms of poly1
            int max_j=J[nvar][deg1];
            if (deg1==0) j_0=0;
            else j_0=N[nvar][deg1-1];

            int deg = deg0+deg1;//order of terms of res
            if (deg==0) idx_0=0;
            else idx_0=N[nvar][deg-1];

            for (int i=0;i<max_i;i++){
                for (int j=0;j<max_j;j++){
                    if(fabs(x1_coeffs[i_0+i])>ZERO && fabs(x2_coeffs[j_0+j])>ZERO){
                        if (use_M) idx = m_M[coeff];
                        else{
                        //find what term is the result contributing to
                            std::vector<int> row0=x1.get_row(i,deg0);
                            std::vector<int> row1=x1.get_row(j,deg1);
                            std::vector<int> row(nvar);
                            for (int k=0;k<nvar;k++) row[k]=row0[k]+row1[k];
                            idx = res_poly.get_idx(row);
                        }
                        //multiply and add contribution
                        res_coeffs[idx_0+idx]+= x1_coeffs[i_0+i]*x2_coeffs[j_0+j];
                    }
                    coeff++;
                }
            }

        }
    }

    res_poly.set_coeffs(res_coeffs);

}


/******************************/
/*J & N matrix manipulation   */
/******************************/
template < class T >
void base_polynomial<T>::initialize_J(){
    int i,j;

    //fill J
    for(j = 0; j <= m_degree; ++j)
        m_J[1][j] = 1;
    for(i = 2; i <= m_nvar; ++i) {
        m_J[i][0] = 1;
        for(j = 1; j <= m_degree; ++j)
            m_J[i][j] = m_J[i][j-1]+m_J[i-1][j];
    }

}

template < class T >
void base_polynomial<T>::initialize_N(){
    int i,j;
    //fill N
    for(i = 1; i <= m_nvar; ++i)
        m_N[i][0] = m_J[i][0];
    for(i = 1; i <= m_nvar; ++i) {
        for(j = 1; j <= m_degree; ++j)
            m_N[i][j] = m_N[i][j-1]+m_J[i][j];
    }

}

template < class T >
std::vector<int> base_polynomial<T>::get_row(const int &idx, const int &deg) const{

    int i, j, l;
    int idx_tmp = idx;
    std::vector<int> k(m_nvar);

    if(deg > 0)
        idx_tmp += m_N[m_nvar][deg-1];
    l = m_nvar;
    for(i = 0; i < m_nvar; ++i) {
        if(idx_tmp == 0) {
            for(j = i; j < m_nvar; k[j++] = 0);
            return k;
        }
        if(i == 0)
            j = deg;
        else {
            for(j = 0; idx_tmp >= m_N[l][j]; ++j);
            k[i-1] -= j;
        }
        k[i] = j;
        idx_tmp -= m_N[l--][j-1];
    }

    return k;
}

template < class T >
int base_polynomial<T>::get_idx(const std::vector<int> &k) const{
    int i, j, l;

    int idx = 0;
    for(i = m_nvar-1; (i>=0) && (k[i]==0); --i);

    if(i < 0)
        return idx;
    j = m_nvar-i;
    l = -1;
    while(i >= 0) {
        l += k[i--];
        idx += m_N[j++][l];
    }
    idx -= m_N[m_nvar][l];

    return idx;

}

//initialisation of M for faster multiplication 
//and Q for more accurate range estimation
template <class T>
void base_polynomial<T>::initialize_M(const int &nvar, const int &degree){

    std::vector<int> J=get_J()[nvar];
    std::vector<int> N=get_N()[nvar];
    std::vector<int> M;
    std::vector<bool> Q;
    for (int deg0=0; deg0<=degree; deg0++){ //loop over order of terms of poly0
        int max_i=J[deg0];
        for (int deg1=0; deg1<=degree-deg0; deg1++){ //loop over other of terms of poly1
            int max_j=J[deg1];
            //int deg = deg0+deg1; //order of terms of result
            for (int i=0;i<max_i;i++){
                std::vector<int> row0=get_row(i,deg0);

                if(deg1==0){
                    bool quadratic = 1;
                    for (int v=0;v<nvar;v++){
                        if (row0[v]%2!=0){
                            quadratic = 0;
                            break;
                        }
                    }
                    Q.push_back(quadratic);
                }

                for (int j=0;j<max_j;j++){
                    std::vector<int> row1=get_row(j,deg1);
                    //find what term is the result contributing to
                    std::vector<int> row(nvar);
                    for (int k=0;k<nvar;k++) row[k]=row0[k]+row1[k];
                    //store it to avoid computing it in every multiplication
                    M.push_back(get_idx(row));
                }
            }
        }
    }
    //modify static stuff
    m_Q=Q;
    m_M=M;
    m_Mdegree = degree;
    m_Mnvar = nvar;
}

template <class T>
void base_polynomial<T>::delete_M(){
    m_M.resize(1);
    m_Mnvar=0;
    m_Mdegree=0;
}

// template <class T>
// void base_polynomial<T>::evaluate_base1D_monomial(const int &index, const base_polynomial<T> &other, base_polynomial<T>  &out){

//     if(index==0){
//         std::vector<T> out_coeffs = out.get_coeffs();
//         out_coeffs[0] = 1.0;
//         for(int i=1; i<out_coeffs.size(); i++)
//             out_coeffs[i] = 0.0;
//     }
//     else if(index==1){
//         std::vector<T> other_coeffs = other.get_coeffs();
//         out.set_coeffs(other_coeffs);
//     }
//     else{
//         for(int i=2; i<index; i++){
//             std::vector<T> other_coeffs = other.get_coeffs();
//             out.set_coeffs(other_coeffs);
//             monomial_multiplication(out,other,out);
//         }
//     }
// }

template <class T>
std::vector<T> base_polynomial<T>::evaluate_basis_monomial(const std::vector<T> &x) const{

    if((std::size_t) m_nvar!=x.size()){
        smart_throw(m_name+": (evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }

    //evaluate the bases
    std::vector < std::vector <T> > base;
    base.resize(m_nvar);
    for (int i=0; i<m_nvar;i++){
        base[i].resize(m_degree+1);
        base[i][0]=1.0;
        for (int j=1; j<=m_degree; j++){
            base[i][j]=x[i]*base[i][j-1];
        }
    }

    //construct the full polynomial value
    std::vector<T> res(m_coeffs.size());
    int idx=0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            T prod = 1.0;
            //if (fabs(m_coeffs[idx])>ZERO){
                std::vector<int> row = this->get_row(i,deg);
                for(int j=0;j<m_nvar; j++){
                    prod*=base[j][row[j]];
                }
                res[idx] = prod;
            //}
            idx++;
        }
    }

    return res;
}

//private routine for 1d evaluation
template <class T>
T base_polynomial<T>::horner(T x, int i) const{
    if (i>m_degree) return 0;
    else return x*(m_coeffs[i]+horner(x,i+1));
}

template class base_polynomial<double>;
template class base_polynomial<float>;
template class base_polynomial<long double>;
