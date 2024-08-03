# A simple and practical supersingularity test for elliptic curves over $\mathbb{F}_p$

A proof-of-concept implemtation using [SageMath] (https://www.sagemath.org) of version 10.3.
This implementation contains three algorithms for supersingularity tests:

1. Doliskani test in "Doliskani.sage"
2. Pairing-based supersingularity verification in "pairing.sage" (based on https://github.com/LinKaizhan/Pairingoptimizations)
3. Our new algorithm in "differential_addition_chain.sage".

Each algorithm count and output the number of field operations during the computation.
To execute the algorithm, one can use the command:
```
sage "algorithm_name.sage"
```
algorithm_name can be Doliskani, pairing or differential_addition_chain.