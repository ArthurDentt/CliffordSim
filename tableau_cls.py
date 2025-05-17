# Tableau Class file for efficient Clifford circuit simulation
# Author: Andrew Connelly

from numpy import shape, asarray

class Tableau():
    """ Base Tableau Class which is used to simulate Clifford circuits """
    
    def __init__(self, data):
        """ Initialize a tableau from data, calculate number of qubits """
        self._data = asarray(data,dtype=bool) # base tableau data.
        self.num_qubits = (shape(data)[0])//2 # number of qubits being simulated. (2n+1)//2 = n 
    
    def __repr__(self):
        """ Repr method for printing tableau stabs & destabs """
        stabilizers = [_row_to_str(self._data[i,:]) for i in range(self.num_qubits)]
        destabilizers = [_row_to_str(self._data[i+self.num_qubits,:]) for i in range(self.num_qubits)]
        return "Stabilizers: " + str(stabilizers) + "\nDestabilizers: " + str(destabilizers) 
    
    @property
    def data(self):
        """ Returns the base tableau data """
        return asarray(self._data,dtype=int) # cast back to integer datatype before returning

    def H(self, index):
        """ Mutates the Tableau's data to perform a Hadamard gate on qubit in position `index` """
        
        # check for invalid qubit index, raise error if so.
        self._check_index_invalid(index)
        
        # adjust phase column (r = r^(x*z))
        self._data[:,-1] ^= (self._data[:,index]&self._data[:,self.num_qubits+index]) 

        # swap x_ind and z_ind columns 
        self._data[:,[index,self.num_qubits+index]] = self._data[:,[self.num_qubits+index,index]]
    
    def S(self, index):
        """ Mutates the Tableau's data to perform an S gate on qubit in position `index` """
       
        # check for invalid qubit index, raise error if so.
        self._check_index_invalid(index)

        # adjust phase column (r = r^(x*z))
        self._data[:,-1] ^= (self._data[:,self.num_qubits+index]&self._data[:,index])

        # adjust z column (z = z^x)
        self._data[:,index] = (self._data[:,index]^self._data[:,self.num_qubits+index])
    
    def CX(self, control, target):
        """ Mutates the Tableau's data to perform an CX gate acting from `control` to `target` """
        
        # check for invalid qubit indices, raise error if so.
        self._check_valid_control_target(control,target)
        self._check_index_invalid(control)
        self._check_index_invalid(target)

        # adjust phase column (r = r ^ x_c * z_t * (x_t ^ z_c ^ 1))
        self._data[:,-1] ^= (self._data[:,self.num_qubits+control] & self._data[:,target] &
                            (self._data[:,self.num_qubits+target] ^ (self._data[:,control] ^ True)))
        # adjust x_t column (x_t = x_t ^ x_c)
        self._data[:,self.num_qubits+target] ^= self._data[:,self.num_qubits+control]
        # adjust z_c column (z_c = z_c ^ z_t)
        self._data[:,control] = (self._data[:,target] ^ self._data[:,control])
    
    def CZ(self,control,target): # TODO: Analytical calc for less swaps
        """ Mutates the Tableau's data to perform a CZ gate acting from `control` to `target` """
        
        # check for invalid qubit indices, raise error if so.
        self._check_valid_control_target(control,target)
        self._check_index_invalid(control)
        self._check_index_invalid(target)

        # naive CZ = H_t CX H_t
        self.H(target)
        self.CX(control,target)
        self.H(target)

    def _check_index_invalid(self, index):
        # check for invalid qubit index, raise error if so.
        if index < 0 or index >= self.num_qubits:
            raise ValueError(rf"Invalid qubit index {index}")
    
    def _check_valid_control_target(self, control, target):
        if control==target:
            raise ValueError("must have distinct control and target indices")


# ------------------------------------- #
#           Helper Functions            #
# ------------------------------------- #


def _row_to_str(rowvec):
    """ row to string method for printing out a stab or destab from a tableau row """
    n = (len(rowvec))//2
    operator = '-' if rowvec[-1] else '+'   # is negative bit present
    for i in range(n):
        if rowvec[i]:                       # if Z Stabilizer component is present
            if rowvec[i+n]:                 # if X Stabilizer component is also present
                operator += 'Y'
            else:
                operator += 'Z'
        elif rowvec[i+n]:                   # if only X Stabilizer component is present
            operator += 'X'
        else:                               # if none are present
            operator += 'I'
    return operator
