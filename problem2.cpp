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
    double a, b, h, xmin, xmax;
    int N, nsteps;
    
    xmin = 0.;
    xmax = 1.;
    N = 6;
    nsteps = N-1;
    h = xmax-xmin/nsteps;

    a = -1./pow(h,2.);
    b = 2./pow(h,2.);

    eigval = arma::vec(N);
    eigvec = arma::mat(N,N);
    eigval_ana = arma::vec(N).fill(0.);
    eigvec_ana = arma::mat(N,N).fill(0.);

    //produce tridiagonal matrix A
    m = produce_tridiag(N, a, b, a);    

    //Solve using the eigsym function
    eig_sym(eigval, eigvec, m);
    eigvec = arma::normalise(eigvec);
    
    //Solve using the analytical solution
    analytical_sym_eig(N, a, b, eigval_ana, eigvec_ana);
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