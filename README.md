# FYS4150-Project-2
Project 2 for the course FYS4150

The different programs in theproject are compiled and linked using the following commands

problem2.cpp:
g++ problem2.cpp -o problem2.exe  -larmadillo
comment: The program print the eigenvalues and eigenvectors of the two methods

max_val_test.cpp (problem3):
g++ max_val_test.cpp -o max_val_test.exe
comment: if the test fails the program will abort with an error

problem4.cpp:
g++ problem4.cpp -o problem4.exe -larmadillo
comment: The program print the eigenvalues and eigenvectors of the two trusted method and the jacobi rotation algorithm

problem5.cpp:
g++ problem5.cpp -o problem5.exe -larmadillo
comment: For checking the scaling of the dense matrix in problem 5, A has to be changed to A2 on line 59 in the eigensolver function

problem6.cpp:
g++ problem6.cpp -o problem6.exe -larmadillo
comment: For creating the data for N=100, N=10 has to be changed to N=100 on line 19

The C++ files create files that are then used in plotting in the python files problem5_plotting.py and problem6_plotting.py that create the figures used in the rapport

functions.hpp contains all the written functions used

