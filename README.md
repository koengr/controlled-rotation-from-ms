# Controlled rotation quantum gates using trapped ions

This code complements the paper "Signal processing techniques for efficient compilation of controlled rotations in trapped ions
" (published as [Groenland et al. (2020) New J. Phys. *22* 063006](https://doi.org/10.1088/1367-2630/ab8830) and [on ArXiv](https://arxiv.org/abs/2001.05231)). 

The main command finds a list of the required *rotation angles* in a circuit that implements a controlled rotation gate (with N qubits and rotation angle alpha) using L=2N applications of the Molmer-Sorensen (MS) gate.

To obtain the L+1 rotation angles that sandwich the MS pulses, use
```python
angles = crot_angles( n_qubits, alpha )
```

# Installation
This code needs no installation. Simply move the file 'CompositeCRotGate.py' into your favorite folder, and from the same folder, start your Python 3 code with ``` from CompositeCRotGate import *```


Dependencies do need to be installed:
```
pip3 install qspd
pip3 install numpy
pip3 install sympy==1.4
```
