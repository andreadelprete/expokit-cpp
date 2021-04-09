# -*- coding: utf-8 -*-
"""
Created on Mon May 18 21:59:02 2020

Script to identify sparsity patterns of matrix powers

@author: student
"""


import numpy as np
from sympy import *
init_printing() #use_unicode=True)

C = zeros(6,6)
A, t = symbols('A, t')
C[0,0] = t*A
C[0,1] = t
C[0,3] = t
C[3,2] = t
C[4,0] = t
C[5,4] = t

print('C')
pprint(C)

print('C^2')
C2 = C @ C
pprint(C2)

C3 = C2 @ C
print('C^3')
pprint(C3)

C4 = C2 @ C2
print('C^4')
pprint(C4)

C5 = C3 @ C2
print('C^5')
pprint(C5)

b0, b1, b2, b3, b4, b5 = symbols('b0, b1, b2, b3, b4, b5')
D = b0*eye(6) + b1*C + b2*C2 + b3*C3 + b4*C4 + b5*C5
pprint(D)

D1, D2, D3, D4, D5 = symbols('D1, D2, D3, D4, D5')
D = zeros(6,6)
D[0,0] = D1
D[0,1] = D2
D[0,2] = D3
D[0,3] = D2
D[1,1] = b0
D[2,2] = b0
D[3,2] = b1
D[3,3] = b0
D[4,0] = D2
D[4,1] = D3
D[4,2] = D4
D[4,3] = D3
D[4,4] = b0
D[5,0] = D3
D[5,1] = D4
D[5,2] = D5
D[5,3] = D4
D[5,4] = b1
D[5,5] = b0
pprint(D)

Dinv = D**-1
print('D inverse')
pprint(D**-1)