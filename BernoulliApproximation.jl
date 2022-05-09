using OffsetArrays
using JLD2
italicB = load_object("./italicB.jld2")
alpha = load_object("./alpha.jld2")
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
function bernoulli(A, mk, s)
    B = P(mk, A/(2^s))
    for i in 1:1:s
        B = B*B
    end
    return B
end
