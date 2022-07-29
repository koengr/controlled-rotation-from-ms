# Controlled rotation quantum gates using trapped ions

This code complements the paper "Signal processing techniques for efficient compilation of controlled rotations in trapped ions
" (published as [Groenland et al. (2020) New J. Phys. *22* 063006](https://doi.org/10.1088/1367-2630/ab8830) and [on ArXiv](https://arxiv.org/abs/2001.05231)). 

The main command lists the required *rotation angles* to perform a controlled rotation gate (with ```n_qubits``` qubits and rotation angle ```angle```), as definted in the paper. They are returned in the form of a [Dictionary](https://www.w3schools.com/python/python_dictionaries.asp). Using these angles, one can construct the corresponding circuit that uses L=2N applications of the Molmer-Sorensen (MS) gate.

To obtain the L+1 rotation angles that sandwich the MS pulses, use
```python
angles = crot_angles( n_qubits, alpha )
# result example: {4: 4.8634, 3: 1.2011, 2: 1.4285, 1: 4.0494, 0: 0.511}
```

# Installation
This code needs no installation. Simply move the file 'CompositeCRotGate.py' into your favorite folder, and from the same folder, start your Python 3 code with ``` from CompositeCRotGate import *```


Dependencies do need to be installed:
```
pip3 install qspd
pip3 install numpy
pip3 install sympy==1.4
```
