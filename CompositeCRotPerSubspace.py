''' 
Koen Groenland, December 2019

This code is intended to aid the implementation of subspace-dependent rotation gates using sequences 
of MS gates, as explained in the paper at [https://arxiv.org/abs/2001.05231].

In particular, define alphas as a list of rotations, where Rz(alphas[0]) occurs when control qubits 
have Hamming weight q=0, Rz(alphas[1]) when q=1, etc. The necessary intermediate gates are then obtained
through

    crot_per_subspace( alphas )

which returns the projectors P, Q (=1-P) in the formalism of Haah's algorithm (https://quantum-journal.org/papers/q-2019-10-07-190/)


The package relies on 'qspd' and 'mpmath', both of which have to be installed before use:

    pip install qspd
    
I found that errors may occur if an old version of Sympy is installed. To fix, use:
    
    pip install sympy==1.4


'''


# Import the required packages:
import qspd
from qspd.util import InputFunction, Parity

import numpy as np
from scipy.linalg import solve

def crot_per_subspace( alphas, return_mpmath = False  ):
    """ Find the intermediate gates E_P that implement C^(nqubits-1) Rz(alpha_q) (where q is the Hamming weigth 
    of the control qubits) in the composite gate formalism. 
    
        :param alphas:           (list of floats) the rotation angles alpha_q, whose length equals the number of qubits. 
        :param return_mpmath:   (boolean) should the function return angles as mpmath floats? 
                                (defaults to False, in which case numpy floats are returned). 
        :return:                (dictionary of floats) the intermediate gates E_P. The 0th entry is a unitary, entries 1...L 
                                are tuples (P,1-P) that define E_P(theta) = exp(i theta/2 ) P + exp( -i theta/2 ) (1-P).
    """     
    
    a_laurent, b_laurent = crot_per_subspace_series( alphas )
    
    # Encode these as InputFunctions to be used by the qspd module.
    big_a = InputFunction( a_laurent, Parity.EVEN )
    big_b = InputFunction( b_laurent, Parity.EVEN )
    
    # Retrieve the sequence of gates, encoded (P, 1-P) for gates E_P.
    gates = qspd.qspd( big_a, big_b, return_matrices = True )
    
    # Cast the sequqence of gates to numpy. 
    if not return_mpmath :
        gates = gate_sequence_to_numpy( gates )
    
    return(gates)

def crot_per_subspace_series( alphas ):    
    """ Find the laurent series for A(phi) and B(phi) that implements C^{N-1}Rz(alpha_q)
        
        The the functions A and B are chosen EVEN (i.e. are cosine series). 

        :param alphas:      (list of floats) the rotation angles alpha_q, whose length equals the number of qubits. 
                
        :return a_laurent:  (list of floats) the coefficients of the Laurent series A(theta)
        :return b_laurent:  (list of floats) the coefficients of the Laurent serues B(theta)
    """

    n_qubits = len(alphas)
    num_dof = 2 * n_qubits - 2
    ks = np.arange(0,num_dof)
    thetas = np.linspace(0, np.pi, n_qubits , endpoint=True)

    # Make the fitting matrix:
    fitting_matrix = [[ np.cos( k * theta  ) for k in ks ] for theta in thetas ]
    deriv_matrix = [[ np.sin( k * theta ) * -1 * k for k in ks ] for theta in thetas[1:-1] ]

    # Make the target values:
    big_g_fitting_a = [ np.cos( alpha /2 ) for alpha in alphas ]
    big_g_fitting_b = [ -np.sin( alpha /2 ) for alpha in alphas ]
    big_g_deriv = [ 0 ] * ( len( thetas ) - 2 )

    # Make the complete solving matrix:
    matrix_a = fitting_matrix + deriv_matrix
    matrix_b = fitting_matrix + deriv_matrix
    vector_a = big_g_fitting_a + big_g_deriv
    vector_b = big_g_fitting_b + big_g_deriv

    # Solve the system of equations
    a_list = solve( matrix_a, vector_a )
    b_list = solve( matrix_b, vector_b )
    
    # Encode as Laurent series, which is required by the qspd code. 
    a_laurent = cos_series_to_laurent( double_list_to_dict( a_list, ks ) )
    b_laurent = cos_series_to_laurent( double_list_to_dict( b_list, ks ) )
    
    return a_laurent, b_laurent
    
    
    
''' given an items and a keys list, create a dictionary ''' 
def double_list_to_dict( items, keys ): 
    return dict( zip( keys, items ) )

''' given a cosine series with coefficients 'series_dict', return the corresponding Laurent polynomial. '''
def cos_series_to_laurent( series_dict ): 
    laurent = {}
    for k, coeff in series_dict.items() :
        if k == 0 : 
            laurent[k] = coeff
        else :
            laurent[k] = coeff / 2.
            laurent[-k] = coeff / 2.
    return laurent

''' turn a gate sequence from qspd into numpy format '''
def gate_sequence_to_numpy( gate_sequence ):
    new_sequence = {}

    new_sequence[0] = mpmath_matrix_to_numpy( gate_sequence[0] )
    for j in range( 1, len(gate_sequence) ):
        new_sequence[j] = [
            mpmath_matrix_to_numpy( gate_sequence[j][0] ),
            mpmath_matrix_to_numpy( gate_sequence[j][1] ),
        ]
    return( new_sequence )

# Given an mpmath matrix object, cast this to a numpy object. 
# Thanks to the .tolist() functionality, no mpmath functions have to be imported. 
def mpmath_matrix_to_numpy( mat ) :
    return np.matrix( mat.tolist() , dtype = complex )


def evaluate_gate_sequence( gate_list, theta ) :
    ''' 
    Given a sequence of matices E_P (such as retrieved from the qspd algorithm), 
    construct the actual unitary matrix that it implements. 
    
    Note that we send theta -> theta/2 because the qspd algorithm assumes 4 pi periodicity 
    of theta (or equivalently, Laurent exponents k that change in steps of 2), whereas we 
    work with 2 pi periodic functions. 
    ''' 
    t = np.exp(+1j * theta / 2)
    circuit = gate_list[0]

    # Append each unitary sequentially:
    for x in range(1, len(gate_list)): 
        unitary =  t * gate_list[x][0] + 1/t * gate_list[x][1]
        circuit = circuit.dot( unitary )

    return(circuit)
