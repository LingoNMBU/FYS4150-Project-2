
#include <armadillo>
#include <cmath>
#include "functions.hpp"



int main ()
{

    double h, xmax, xmin, eps, a, d;
    int nsteps, N, maxiter, iterations;
    arma::mat R, A, eigenvectors, eigvec, eigvec_ana;;
    arma::vec eigenvalues, eigval, eigval_ana;
    bool converged;

    
    //Initialization
    N = 6;
    nsteps = N - 1;
    eps = 10e-8;
    xmin = 0.0;
    xmax = 1.;
    h = (xmax-xmin)/nsteps;
    a = -1./pow(h,2.);
    d = 2./pow(h,2.);   
    converged = false;
    iterations = 0;
    maxiter = 1000000;

    //produce tridiagonal matrix A
    A = produce_tridiag(N, a, d, a);

    std::cout << __LINE__
              << std::endl
              << A.n_rows
              << std::endl;



    eigenvalues = arma::vec(N);
    eigenvectors = arma::mat(N,N);

    eigval = arma::vec(N);
    eigvec = arma::mat(N,N);
    eigval_ana = arma::vec(N).fill(0.);
    eigvec_ana = arma::mat(N,N).fill(0.);

        //std::cout << eigenvalues
          //        << std::endl;

    //solve using jacobi rotation
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    eigenvectors = arma::normalise(eigenvectors);

    //Solve using the eigsym function
    eig_sym(eigval, eigvec, A);
    eigvec = arma::normalise(eigvec);
    
    //Solve using the analytical solution
    analytical_sym_eig(N, a, d, eigval_ana, eigvec_ana);
    eigvec_ana = arma::normalise(eigvec_ana);

    std::cout << A;
    std::cout << std::endl;
    std::cout << eigvec;
    std::cout << std::endl;
    std::cout << eigvec_ana;
    std::cout << std::endl;
    std::cout << eigenvectors;
    std::cout << std::endl;
    std::cout << eigval;
    std::cout << std::endl;
    std::cout << eigval_ana;
    std::cout << std::endl;
    std::cout << eigenvalues;



    return 0;    
}