"""
BernoulliApproximation.jl is a Julia module that contains the Julia functions P(mk, A) and bernoulli(A, mk, s).
bernoulli(A, mk, s) should be called by the user.
P(mk, A) should not be directly called by the user.
italicB and alpha also should not be directly used by the user.
"""
using OffsetArrays
using JLD2
"""
italicB: A precomputed array of Bernoulli numbers 0 through 50
alpha: A precomputed array of alpha values as specified in my paper
"""
italicB = load_object("./italicB.jld2")
alpha = load_object("./alpha.jld2")
"""
A helper function that computes P_m_k as described in my paper.
"""
function P(mk, A)
    q = Int64(sqrt(mk))
    Ap = OffsetArray{Matrix{Float64}}(undef, 0:q)
    Ap[0] = A^0
    for i in 1:1:q
        Ap[i] = Ap[i-1]*A
    end
    answer = alpha[mk,mk]*Ap[q]
    for i in 0:1:q-2
        for j in 1:1:q
            answer = answer + alpha[mk,mk-i*q-j]*Ap[q-j]
        end
        answer = answer * Ap[q]
    end
    for j in 1:1:q
        answer = answer + alpha[mk,q-j]*Ap[q-j]
    end
    return (ℯ-1)*answer
end
"""
Computes an approximation of e^A using Bernoulli polynomials and the parameters m_k and s specified in my paper.
Parameters:
A: A square matrix
mk: The parameter m_k described in my paper
s: The parameter s described in my paper
"""
function bernoulli(A, mk, s)
    B = P(mk, A/(2^s))
    for i in 1:1:s
        B = B*B
    end
    return B
end
