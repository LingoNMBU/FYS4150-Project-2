
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

void analytical_sym_eig(int n, double a, double b, double c, arma::vec &eigval, arma::mat &eigvec)
{
    //Solves the eigenvalue/eigenvector problem
    // of a a tridiag matrix of the form 
    //a, b, a

    double pi = 2.*acos(0.0);

    for (int i = 1; i < n+1; i++){
        //std::cout << i;

        eigval(i-1) = b + (2. * a * std::cos( (i * pi) / (n + 1.) ));
        //std::cout << eigval(i) 
        //          << std::endl;

       for (int j = 1; j < n+1; j++){
            //

            eigvec(i-1,j-1) = std::sin( (j * i * pi) / (n + 1.));
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
    double max_val = 0.;
    int n = m.n_rows;

    for (int i = 0; i < n-1; i++){
        for (int j = 1; j < n; j++){
            if (i < j){
                double absval = fabs(m(i,j));
                //std::cout << i << j ;
                //std::cout << std::endl;
                //std::cout << fabs(val) ;
                //std::cout << std::endl;
                //std::cout << std::endl;

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