# !python
# -*- coding: utf-8 -*

__author__ = 'Erling Ween Eriksen'
__email__ = 'erlinge@nmbu.no'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

eigenvalues = pd.read_csv('Prob6_eigenvalues',
                          names=['eigenvalues'])
eigenvectors = pd.read_csv('Prob6_eigenvectors',
                           names=[f'eigenvector {i}' for i in range(1, len(eigenvalues)+1)])
eigenvectors.index = range(1, len(eigenvectors) + 1)
eigenvalues.index = range(1, len(eigenvalues) + 1)
print(eigenvectors)
N = len(eigenvalues)

min_inds = [eigenvalues['eigenvalues'].idxmin()]
eigenvalues.drop(eigenvalues['eigenvalues'].idxmin(), inplace=True)
min_inds.append(eigenvalues['eigenvalues'].idxmin())
eigenvalues.drop(eigenvalues['eigenvalues'].idxmin(), inplace=True)
min_inds.append(eigenvalues['eigenvalues'].idxmin())

# print(eigenvalues)

print(min_inds)
if N == 100:
    y1 = eigenvectors['eigenvector 53'].values
    y2 = eigenvectors['eigenvector 21'].values
    y3 = eigenvectors['eigenvector 69'].values
if N == 10:
    y1 = eigenvectors['eigenvector 5'].values
    y2 = eigenvectors['eigenvector 1'].values
    y3 = eigenvectors['eigenvector 9'].values

y1 = np.append(y1, 0)
y2 = np.append(y2, 0)
y3 = np.append(y3, 0)
y1 = np.insert(y1, 0, 0)
y2 = np.insert(y2, 0, 0)
y3 = np.insert(y3, 0, 0)

# y1[N + 1] = 0.0
# y2[N + 1] = 0.0
# y3[N + 1] = 0.0

x10 = np.linspace(0,1,12)
x100 = np.linspace(0,1,102)

if N == 10:
    plt.plot(x10, y1, label='eigenvector 5', marker='o', linestyle='-')
    plt.plot(x10, y2, label='eigenvector 1', marker='o', linestyle='-')
    plt.plot(x10, y3, label='eigenvector 9', marker='o', linestyle='-')
    plt.title('Buckling beam solutions for N=10')
    plt.xlabel('x-hat')
    plt.ylabel('u(x-hat)')
    plt.legend()
    plt.savefig('Buckling_n10.pdf')
if N == 100:
    plt.plot(x100, y1, label='eigenvector 53', marker='o', linestyle='-')
    plt.plot(x100, y2, label='eigenvector 21', marker='o', linestyle='-')
    plt.plot(x100, y3, label='eigenvector 69', marker='o', linestyle='-')
    plt.title('Buckling beam solutions for N=100')
    plt.xlabel('x-hat')
    plt.ylabel('u(x-hat)')
    plt.legend()
    plt.savefig('Buckling_n100.pdf')
plt.show()
