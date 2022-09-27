# !python
# -*- coding: utf-8 -*

__author__ = 'Erling Ween Eriksen'
__email__ = 'erlinge@nmbu.no'

import pyarma as arma
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math


Ns = pd.read_csv('Prob5_Ns', names = ['N'])
iterations = pd.read_csv('Prob5_iterations_Ns', names=['iterations'])
N2 = Ns**2
N212 = Ns**(2.12)
N22 = Ns**(2.2)

print(iterations)

plt.plot(Ns,iterations, label='measured', marker = 'o', linestyle='')
plt.plot(Ns,N2, label='N^2')
plt.plot(Ns,N212, label='N^2.12')
plt.plot(Ns,N22, label='N^2.2')
plt.title('Scaling of number of necessary rotations with eps = 10e-8')
plt.xlabel('Size of matrix [N x N]')
plt.ylabel('Necessary iterations')
plt.legend()
plt.savefig('Rotations_scaling.pdf')
plt.show()

# def transform_log(a,b):
#     a_array = np.array(a.values)
#     for i,a in enumerate(a_array):
#         a_array[i] = math.log(a,b)
#     return a_array
#
#
# plt.plot(Ns, transform_log(iterations,math.e), label='measured, transformed', marker = 'o', linestyle='')
# plt.title('Scaling of number of necessary rotations with eps = 10e-8')
# plt.xlabel('Size of matrix [N x N]')
# plt.xlabel('Necessary iterations')
# plt.legend()
# plt.show()