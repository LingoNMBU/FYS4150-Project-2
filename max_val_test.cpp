#include "functions.hpp"

int main () 
{
arma::mat A;
int n, k, l;
double a, b, max_val, true_max;

n = 4;
a = 0.;
b = 1.;
true_max = 0.7;

A = produce_tridiag(n, a, b, a);

A(1,2) = -0.7;
A(0,3) = 0.5;
A(2,1) = -0.7;
A(3,0) = 0.5;

//std::cout << A;

k = 0;
l = 0;
max_val = max_offdiag_element_sym(A, k, l);

std::cout << A << std::endl
          << k << std::endl
          << l << std::endl
          << max_val<< std::endl;


if (max_val == true_max){
    return 0;
}
else{
    return 0;
}


}
