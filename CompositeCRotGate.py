''' 
Koen Groenland, December 2019

This code is intended to aid the implementation of controlled-rotation gates using sequences 
of MS gates, as explained in our paper [https://arxiv.org/abs/2001.05231]. In particular, use

    crot_angles( number_of_qubits,  rotation_angle_alpha )
    
to obtain a dictionary of N+1 angles, corresponding to the Z-rotation angles phi_j that 
interleave the MS gates. 

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



def crot_angles( n_qubits, alpha, return_mpmath = False ):
    """ Find the angles phi_j that implement C^(nqubits-1) Rz(alpha) in the composite gate formalism. 
    
        :param n_qubits:        (integer) the total number of qubits in the multiqubit gate.
        :param alpha:           (float) the rotation angle in Rz(alpha) = exp( i Z alpha / 2 )
        :param return_mpmath:   (boolean) should the function return angles as mpmath floats? 
                                (defaults to False, in which case numpy floats are returned). 
        :return:                (dictionary of floats) the rotation angles phi[j]. 
    """ 
    
    # Step 1: find the laurent series for A and B.
    a_laurent, b_laurent = crot_laurent_series( n_qubits, alpha )
    
    # Step 2: encode these as InputFunctions to be used by the qsdp module.
    big_a = InputFunction( a_laurent, Parity.EVEN )
    big_b = InputFunction( b_laurent, Parity.ODD )
    
    # Step 3: retrieve the corresponding sequence of angles phi_j
    phis = qspd.qspd( big_a, big_b, return_matrices = False )
    
    # Step 4: if needed, cast the angles phi_j to numpy.
    if not return_mpmath : 
        phis = mpmath_dict_to_numpy( phis )

    return( phis )


def crot_laurent_series( n_qubits, alpha ):
    """ Find the laurent series for A(theta) and B(theta) that implements 
        (Identity if theta != pi) and (Rz(alpha) if theta == pi).
        
        This implementation sets B(theta) = 0 and solves the coefficients of A(theta). 
        The function A is EVEN (i.e. is a cosine series), the function B is 0 and thus qualifies as ODD. 

        :param n_qubits:    (integer) the total number of qubits in the multiqubit gate, 
                            i.e. the number of values theta.
        :param alpha:       (float) the rotation angle in Rz(alpha) = exp( -i Z alpha / 2 )
                
        :return:            (list of floats) the coefficients of the *cosine series*, such that 
                            A(theta) = sum_{k}  a_list[k] * cos(k theta)

    """
    
    
    # Define the relevant angles theta, and the subsets needed to fix F^{phi}(theta) and its derivatives. 
    # By periodicity, it doesn't matter where the series start: we choose the special point -pi for simplicity.
    thetas = np.arange( -1 * np.pi, np.pi, 2*np.pi / n_qubits ) 
    thetas_fitting = thetas[0:n_qubits//2+1]
    thetas_deriv = thetas[1:(n_qubits+1)//2]
    
    # The degree of the resultant series depends on the number of free parameters (degrees of freedom DOF)
    num_dof = len(thetas_fitting) + len(thetas_deriv)  #big_n = num_dof - 1
    ks = np.arange(0,num_dof)
    
    # Define the function G(theta), which is the 'goal' function
    big_g = [1]*n_qubits  
    big_g[0] = np.cos(alpha/2) # choose the value at -pi (index 0) to be the special point. 
   
    # The relevant entries (throwing away entries that are unnecessary due to symmetries)
    big_g_fitting = big_g[0:len(thetas_fitting)]
    big_g_deriv = [0]*len(thetas_deriv)

    # Make matrix/vector combination that encodes the linear system we're solving.
    fitting_matrix_a = [[ np.cos( k * theta ) for k in ks ] for theta in thetas_fitting ]
    deriv_matrix_a = [[ np.sin( k * theta ) * -1 * k for k in ks] for theta in thetas_deriv ]

    matrix_a = fitting_matrix_a + deriv_matrix_a
    vector_a = np.real( big_g_fitting + big_g_deriv )

    # Solve the system of equations:
    a_list = solve( matrix_a, vector_a )
    
    # Create dictionaries for A and B:
    a_laurent = cos_series_to_laurent( double_list_to_dict( a_list, ks ) )
    b_laurent = double_list_to_dict( [0]*(num_dof), np.arange( 0, num_dof )) #create the all-0 Laurent poly. 
    
    return a_laurent, b_laurent



def mpmath_dict_to_numpy( mpmath_dict, cast_method = np.float ):
    ''' given a dictionary whose items are mpmath numbers, return a new object with the same items in numpy format
    
        :param mpmath_dict: the input dictionary, whose items are mpmath numbers.
        :param cast_method: the function used to turn the mpmath number into a numpy number. Defaults to np.float( ). 
        :return:            a dictionary of numpy numbers. 
    '''
    np_dict = {}
    for key, item in mpmath_dict.items():
        np_dict[key] = cast_method( item )
    return np_dict




# Check if a numpy matrix is unitary, up to tolerance. 
def matrix_is_unitary( mat, tolerance = 10**(-4) ):
    if np.linalg.norm( mat.dot( mat.H) - np.eye(2) ) < tolerance :
        return True
    return False 


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

''' evaluate the Laurent polynomial (at the point theta) that corresponds to the dictionary 'laurent'. '''
def evaluate_laurent_dict( laurent, theta ): 
    return sum( [ coeff * np.exp( 1j * k * theta ) for k, coeff in laurent.items() ])
    
    
    
    

def evaluate_angle_sequence( phi_dict, theta ):
    ''' Given a list of angles phi_j (as dictionary), evaluate the corresponding composte gate at angle theta. 
    
        :param phi_dict:  dictionary of angles, as received from the function 'crot_angles'. 
        :param theta:       (float) rotation angle, the argument of the composite gate F^{phi}(theta)
        
        :return:            approximately unitary matrix that represents the composite gate. 
        
    ''' 
    
    # The implementation first builds a list of all the 2x2 unitaries, 
    # and then multiplies them using np.linalg.multi_dot( ). 
    list_of_unitaries = [ rotz( phi_dict[0]  ) ] 
    x_rotation = rotx( theta )
    for j in range(1,len(phi_dict)):
        list_of_unitaries.append( rotz( phi_dict[j]  ) )
        list_of_unitaries.append( x_rotation )
        list_of_unitaries.append( rotz( -1 * phi_dict[j]  ) )
        
    return np.linalg.multi_dot( list_of_unitaries )

        
''' Rotation around the pauli-X ''' 
def rotx( angle ) :
    return np.cos(angle/2) * np.eye(2) - 1j * np.sin(angle/2) * np.eye(2)[[1,0]]

''' Rotation around the pauli-Z ''' 
def rotz( angle ):
    return np.cos(angle/2) * np.eye(2) - 1j * np.sin(angle/2) * (np.eye(2) * (1,-1))    
    
    
  
