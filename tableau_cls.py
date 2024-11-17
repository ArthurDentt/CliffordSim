# Tableau Class file for efficient Clifford circuit simulation
# Author: Andrew Connelly

from numpy import shape
from copy import copy

class Tableau():
    """ Base Tableau Class which is used to simulate Clifford circuits """
    
    def __init__(self, data):
        """ Initialize a tableau from data, calculate number of qubits """
        self.data = data # base tableau data.
        self.num_qubits = (shape(data)[0])//2 # number of qubits being simulated. (2n+1)//2 = n 
    
    def __repr__(self):
        """ repr method for printing tableau stabs & destabs """
        stabilizers = [_row_to_str(self.data[i,:]) for i in range(self.num_qubits)]
        destabilizers = [_row_to_str(self.data[i+self.num_qubits,:]) for i in range(self.num_qubits)]
        return "Stabilizers: " + str(stabilizers) + "\nDestabilizers: " + str(destabilizers) 
    
    def H(self, index):
        """ Mutates the Tableau's data to perform a Hadamard gate on qubit in position `index` """
        
        # check for invalid qubit index, raise error if so.
        if index < 0 or index >= self.num_qubits:
            raise ValueError(rf"Invalid qubit index {index}")
        
        # adjust phase column (r = r^(x*z))
        self.data[:,-1] ^= (self.data[:,index]*self.data[:,self.num_qubits+index]) 

        # swap x_ind and z_ind columns 
        temp = copy(self.data[:,index]) # TODO: Cleaner swap?
        self.data[:,index] = self.data[:,self.num_qubits+index]
        self.data[:,self.num_qubits+index] = temp
    
    def S(self, index):
        """ Mutates the Tableau's data to perform an S gate on qubit in position `index` """
        # adjust phase column (r = r^(x*z))
        self.data[:,-1] ^= (self.data[:,self.num_qubits+index]*self.data[:,index])
        # adjust z column (z = z^x)
        self.data[:,index] = (self.data[:,index]^self.data[:,self.num_qubits+index])
    
    def CX(self, control, target):
        """ Mutates the Tableau's data to perform an CX gate acting from `control` to `target` """
        # adjust phase column (r = r ^ x_c * z_t * (x_t ^ z_c ^ 1))
        self.data[:,-1] ^= (self.data[:,self.num_qubits+control] * self.data[:,target] *
                            (self.data[:,self.num_qubits+target] ^ (self.data[:,control] ^ 1)))
        # adjust x_t column (x_t = x_t ^ x_c)
        self.data[:,self.num_qubits+target] ^= self.data[:,self.num_qubits+control]
        # adjust z_c column (z_c = z_c ^ z_t)
        self.data[:,control] = (self.data[:,target] ^ self.data[:,control])
    
    def CZ(self,control,target): # TODO: Analytical calc for less swaps
        """ Mutates the Tableau's data to perform a CZ gate acting from `control` to `target` """
        # naive CZ = H_t CX H_t
        self.H(target)
        self.CX(control,target)
        self.H(target)

# ------------------------------------- #
#           Helper Functions            #
# ------------------------------------- #

def _row_to_str(rowvec):
    """ row to string method for printing out a stab or destab from a tableau row """
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

def _swap(vec1,vec2):
    """ swap method to help with swapping data without temp variables (useful ??)"""
    return vec2, vec1