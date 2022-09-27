#include <armadillo>
#include <cmath>
#include "functions.hpp"

int main ()
{

    double h, xmax, xmin, eps, a, d;
    int nsteps, maxiter, iterations, NN, N, runs, run;
    arma::mat R, A, eigenvectors, eigvec, eigvec_ana;;
    arma::vec eigenvalues, eigval, eigval_ana, Ns, iterations_N;
    bool converged;

    
    //Initialization
    Ns = {10,20, 30, 40, 50 , 60, 70, 80, 90, 100};
    NN = Ns.n_rows;
    iterations_N = arma::vec(10);

    //timing variables
    runs = 10;
    run = 1;


    for (int i = 0; i < NN; i++){
        N = Ns(i);
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

        // Generate random N*N matrix
        arma::mat A2 = arma::mat(N, N).randn();  

        // Symmetrize the matrix by reflecting the upper triangle to lower triangle
        A2 = arma::symmatu(A);  

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

        iterations_N(i) = iterations;              

        std::cout << "Size: " << A2.n_rows << "\n";
        std::cout << "Number of iterations: " << iterations << "\n";
        std::cout << "Convergence: " << converged << "\n";
        //std::cout << "Max offdiag value: " << max_val << "\n";
        //std::cout << "Epsilon: " << eps << "\n";
        std::cout << "\n";
    }


    std::cout << "Number of iterations: " << "\n" <<iterations_N << "\n";

    std::string filename1 = "Prob5_Ns";
    Ns.save(filename1, arma::csv_ascii);

    std::string filename2 = "Prob5_iterations_Ns";
    iterations_N.save(filename2, arma::csv_ascii);

    


    return 0;    
}