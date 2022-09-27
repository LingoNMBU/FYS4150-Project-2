
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
    N = 10;
    nsteps = N - 1;
    eps = 10e-8;
    xmin = 0.0;
    xmax = 1.;
    h = (xmax-xmin)/nsteps;
    a = -1./pow(h,2.);
    d = 2./pow(h,2.);   
    converged = false;
    iterations = 0;
    maxiter = N*N*N; //since iterations seem to scale about N^(2.12), this should be sufficient

    //produce tridiagonal matrix A
    A = produce_tridiag(N, a, d, a);

    //Initialize eigeinvectors and eigenvalues
    eigenvalues = arma::vec(N);
    eigenvectors = arma::mat(N,N);

    //solve using jacobi rotation
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    eigenvectors = arma::normalise(eigenvectors);

    std::cout << A;
    std::cout << std::endl;
    std::cout << eigenvectors;
    std::cout << std::endl;
    std::cout << eigenvalues;
    
    std::string filename1 = "Prob6_eigenvectors";
    eigenvectors.save(filename1, arma::csv_ascii);

    std::string filename2 = "Prob6_eigenvalues";
    eigenvalues.save(filename2, arma::csv_ascii);

    return 0;    
}