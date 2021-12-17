# QuadProg.jl

An allocation-free primal-dual interior point solver for convex QP's written in pure Julia. The Mehrotra predictor-corrector algorithm used is the same one detailed in [the cvxgen paper](https://stanford.edu/~boyd/papers/pdf/code_gen_impl.pdf). 

The linear systems are solved with [a fork of QDLDL](https://github.com/kevin-tracy/QDLDL.jl), with iterative refinement to ensure a successful factorization. 

```
code here 
```
