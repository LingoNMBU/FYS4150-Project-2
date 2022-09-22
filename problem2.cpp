#include "functions.hpp"
#include <armadillo>
#include <cmath>

// Needs to link with armadillo llpack and lblas
//g++ problem2.cpp -o problem2.exe -llapack -lblas -larmadillo
//ONLY g++ problem2.cpp -o problem2.exe -larmadillo WORKS

int main ()
{

    arma::mat m, eigvec, eigvec_ana;
    arma::vec eigval, eigval_ana;
    double a, b, c, lambda, h, xmin, xmax;
    int n;
    
    xmin = 0.;
    xmax = 1.;
    n = 6;
    h = 1./n;

    a = -1./pow(h,2.);
    b = 2./pow(h,2.);

    //std::cout << n;
    //std::cout << std::endl;
    //std::cout << h;
    //std::cout << std::endl;
    //std::cout << a;
    //std::cout << std::endl;
    //std::cout << b;
    //std::cout << std::endl;

    eigval = arma::vec(n);
    eigvec = arma::mat(n,n);
    eigval_ana = arma::vec(n).fill(0.);
    eigvec_ana = arma::mat(n,n).fill(0.);

    m = produce_tridiag(n, a, b, a);    

    eig_sym(eigval, eigvec, m);
    eigvec = arma::normalise(eigvec);

    
    analytical_sym_eig(n, a, b, a, eigval_ana, eigvec_ana);
    eigvec_ana = arma::normalise(eigvec_ana);
    

    std::cout << m;
    std::cout << std::endl;
    std::cout << eigvec;
    std::cout << std::endl;
    std::cout << eigvec_ana;
    std::cout << std::endl;
    std::cout << eigval;
    std::cout << std::endl;
    std::cout << eigval_ana;

    return 0;
}