# BernoulliApproximation
BernoulliApproximation contains Julia code that uses Bernoulli matrix polynomials to approximate matrix exponentials. The algorithm is described in my paper "Implementing an algorithm using Bernoulli matrix polynomials to approximate matrix exponentials."
BernoulliApproximation.jl contains the function bernoulli(A, mk, s), which should be called by the user.
The other functions and the global variables in BernoulliApproximation.jl should not be directly used by the user.
The two scratch notebooks in this repository should not be run by the user.
Procedure for using the function bernoulli(A, mk, s):
1. Fork this repository.
2. Include the module BernoulliApproximation.jl in your Julia code.
3. Call the function bernoulli(A, mk, s) in your Julia code.
