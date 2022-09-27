
#include <armadillo>
#include <iostream>
#include <cmath>


arma::mat produce_tridiag(int n, double a, double b, double c)
{
    //Creates a tridiagonal n x n matrix with 
    //The superdiagonal contianing the value a
    //The diagonal containing the value b
    //The subdiagonal containg the value c
    //Order N operation

    arma::mat m = arma::mat(n,n).fill(0.);

    for (int i = 0; i < n; i++){

        m(i,i) = b;

        if (i < n-1){               

            m(i, i+1) = a;
            m(i+1, i) = c;
            
        }
    }

    return m;
}


//double ** produce_std_tridiag(int n, double a, double b, double c)
//{
//    //Creates a tridiagonal n x n matrix with 
//    //The superdiagonal contianing the value a
//    //The diagonal containing the value b
//    //The subdiagonal containg the value c
//    //Order N operation
//
//    double mat[n][n] {};
//
//    for (int i = 0; i < n; i++){
//
//        m[i,i] = b;
//
//        if (i < n-1){               
//
//            m[i, i+1] = a;
//            m[i+1, i] = c;
//            
//        }
//    }
//
//    return m;
//}

void analytical_sym_eig(int N, double a, double b, arma::vec &eigval, arma::mat &eigvec)
{
    //Solves the eigenvalue/eigenvector problem
    // of a a tridiag matrix of the form 
    //a, b, a

    double pi = 2.*acos(0.0);

    for (int i = 1; i < N+1; i++){
        //std::cout << i;

        eigval(i-1) = b + (2. * a * std::cos( (i * pi) / (N + 1.) ));
        //std::cout << eigval(i) 
        //          << std::endl;

       for (int j = 1; j < N+1; j++){
            //
            //std::cout << i
            //          << std::endl;
            //std::cout << j
            //          << std::endl;                                           

            eigvec(i-1,j-1) = std::sin( (j * i * pi) / (N + 1.));
            //std::cout << eigvec(i,j) 
            //          << std::endl;;
        }


    }
    //std::cout << eigval;
    //std::cout << std::endl;
    //std::cout << eigvec;
    //std::cout << std::endl;
}

double max_offdiag_element_sym(const arma::mat &m, int &k, int &l)
{
    //Find largest off-diagonal element in a symmetric matrix
    //Currently O(N^2), but its around 0.5n^2 checks
    double max_val = 0.;
    int N = m.n_rows;

    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){
            if (i < j){
                double absval = fabs(m(i,j));
                if (absval > max_val){
                    max_val = absval;
                    k = i;
                    l = j;
                }
            }        
        }
    }

    return max_val;
}


// "Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R" (from project description)
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
{
    double tau, t, s, c;
    int N = A.n_rows;

    //Compute tau
    

    if (A(k,l) != 0.0) {
        tau = ( A(l,l) - A(k,k) ) / ( 2 * A(k,l) );
        // Compute trigonometric parameters choosing the smallest rotation
        if (tau > 0){
            t =  1. / (  tau + sqrt( 1. + tau*tau ) );
        }
        else{
            t = -1. / ( -tau + sqrt( 1. + tau*tau ) );
        }

        c = 1./sqrt(1. + t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }

    //Transform A and R

    //Special l and k number elements
    double Akk_prev = A(k,k);
    double All_prev = A(l,l);

    A(k,k) = (Akk_prev * c * c) - (2. * A(k,l) * c * s) + (All_prev * s*s);
    A(l,l) = (All_prev * c * c) + (2. * A(k,l) * c * s) + (Akk_prev * s*s);
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    //Update A matrix
    for (int i = 0; i < N; i++){
        if (i != k && i !=l){

                double Aik_prev = A(i,k);
                double Ail_prev = A(i,l);

                A(i,k) = (Aik_prev * c) - (Ail_prev * s);
                A(k, i) = A(i,k);

                A(i,l) = (Ail_prev * c) + (Aik_prev * s);
                A(l,i) = A(i,l);            
        }
    }

    //Update Rotation matrix
    for (int i = 0; i < N; i++){

        double Rik_prev = R(i,k);
        double Ril_prev = R(i,l);

        R(i,k) = (Rik_prev * c) - (Ril_prev * s);

        R(i,l) = (Ril_prev * c) + (Rik_prev * s);                
    }
}


// "Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter"(from project description)
void jacobi_eigensolver(const arma::mat& A1, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged)
{
    
    int k,l, N;
    double max_val;
    arma::mat A = A1;
    N = A.n_rows;
    iterations = 0;

    arma::mat R = arma::mat(N, N, arma::fill::eye);
    converged = false;

    //First run through A for k and l    
    max_val = max_offdiag_element_sym(A, k, l); // Gives 0 start indexes of l and k

    while (converged == false && iterations < maxiter){

        jacobi_rotate(A,R,k,l);
        //rotate(A,R,k,l,N);
        //std::cout << std::endl
        //          << k
        //          << std::endl
        //          << l
        //          << std::endl
        //          << max_val
        //          << std::endl
        //          << A
        //          << std::endl
        //          << R
        //          << std::endl
        //          << std::endl;

        max_val = max_offdiag_element_sym(A, k, l); 
        iterations += 1;

        if (max_val < eps){
            converged = true;
        }

    }

    
    //std::cout << std::endl
    //          << max_val
    //          << std::endl
    //          << A
    //          << std::endl
    //          << R
    //          << std::endl
    //          << std::endl;


    for (int i = 0; i < N; i++){
        eigenvalues(i) = A(i,i);
    }
    eigenvectors = R;


}

