# tableaus for gates

import numpy as np
import copy

def _swap(vec1,vec2):
    return vec2, vec1

def CNOT(tableau, control, target):
    n = (len(tableau[0,:]))//2
    za = tableau[:,control]
    xa = tableau[:,n+control]
    zb = tableau[:,target]
    xb = tableau[:,n+target]
    tableau[:,-1] ^= xa*zb*(xb^(za^1))
    tableau[:,n+target] = xb ^ xa
    tableau[:,control] = zb ^ za
    return tableau

def Hadamard(tableau, index):
    n = (len(tableau[0,:]))//2
    za = tableau[:,index]
    xa = tableau[:,n+index]
    tableau[:,-1] ^= (xa*za)
    temp = copy.copy(tableau[:,index])
    tableau[:,index] = tableau[:,n+index]
    tableau[:,n+index] = temp
    return tableau

def Phase(tableau,index):
    n = (len(tableau[0,:]))//2
    za = tableau[:,index]
    xa = tableau[:,n+index]
    tableau[:,-1] ^= xa*za
    tableau[:,index] = za^xa
    return tableau

def tab_row_to_str(rowvec):
    n = (len(rowvec))//2
    operator = '-' if rowvec[-1] else '+'
    for i in range(n):
        if rowvec[i] and rowvec[i+n]:
            operator += 'Y'
        elif rowvec[i]:
            operator += 'Z'
        elif rowvec[i+n]:
            operator += 'X'
        else:
            operator += 'I'
    return operator

def print_tab(tableau):
    n = len(tableau[0,:])//2
    stabilizers = [tab_row_to_str(tableau[i,:]) for i in range(n)]
    destabilizers = [tab_row_to_str(tableau[i+n,:]) for i in range(n)]
    print("Stabilizers:", stabilizers, "Destabilizers:", destabilizers) 