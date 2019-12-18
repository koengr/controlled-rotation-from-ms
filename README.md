# Controlled rotation quantum gates using trapped ions

This code complements the paper "Sequences of Molmer-Sorensen gates can implement controlled rotations using signal processing techniques" (link coming soon). In particular, it finds the circuit that implements a controlled rotation gate (with N qubits and rotation angle alpha) using 2N applications of the Molmer-Sorensen (MS) gate.

To obtain the L+1 rotation angles that sandwich the MS pulses, use
```python
angles = crot_angles( n_qubits, alpha )
```

