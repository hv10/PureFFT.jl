module PureFFT
using Primes
import Base.show

export fft_cooley_tukey, dft
export plan_fft, fft, ifft

struct FFTPlan
    f::Int
    method::Function
    plan::Vector{FFTPlan}
end
function Base.show(io::IO, mime::MIME"text/plain", fp::FFTPlan)
    indent = get(io, :indent, 0)
    println(io, ' '^indent, "> $(fp.f) | $(fp.method)")
    for p in fp.plan
        show(IOContext(io, :indent => indent + 2), mime, p)
    end
end
Base.push!(fp::FFTPlan, args...) = push!(fp.plan, args...)
FFTPlan(f::Int, method::Function) = FFTPlan(f, method, Vector{FFTPlan}())

get_factor_arr(N::Int) =
    mapreduce(vcat, eachfactor(N)) do (f, c)
        repeat([f], c)
    end

plan_fft_min(N::Int; method::Symbol=:dit)::FFTPlan = begin
    @assert (method == :dit) || (method == :dif) "Method needs to be either :dit (decimation-in-time) or :dif (decimation-in-frequency)"
    # min rad planning
    # get prime factors of N
    N_fcts = sort(get_factor_arr(N), rev=true)
    plan = FFTPlan(N, fft_cooley_tukey)
    NN = N
    ptr = plan.plan
    for f in N_fcts
        NN = div(NN, f)
        if f == 1
            P = FFTPlan(f, dft_1)
        elseif f == 2
            P = FFTPlan(f, dft_2)
        elseif f == 4
            P = FFTPlan(f, dft_4)
        else
            P = FFTPlan(f, dft)
        end
        Q = FFTPlan(NN, fft_cooley_tukey)
        if method == :dit
            push!(ptr, Q, P)
        else
            push!(ptr, P, Q)
        end
        if NN == 1
            q_index = method == :dit ? 1 : 2
            ptr[q_index] = FFTPlan(NN, dft_1)
            break
        elseif NN == 2
            q_index = method == :dit ? 1 : 2
            ptr[q_index] = FFTPlan(NN, dft_2)
            break
        elseif NN == 4
            q_index = method == :dit ? 1 : 2
            ptr[q_index] = FFTPlan(NN, dft_4)
            break
        else
            ptr = Q.plan
        end
    end
    return plan
end

plan_fft_max(N::Int; method::Symbol=:dit) = begin
    @assert (method == :dit) || (method == :dif) "Method needs to be either :dit (decimation-in-time) or :dif (decimation-in-frequency)"
    throw(ErrorException("Not yet implemented"))
end

dft_1(x::AbstractArray, args...; kwargs...) = x
dft_2(x::AbstractArray, args...; kwargs...) = [x[1] + x[2], x[1] - x[2]]
dft_4(x::AbstractArray, args...; NN::Int=4, inverse::Bool=false, normalize::Bool=false) = begin
    @assert length(x) == 4 "DFT_4 is only valid for 4-element arrays"
    inv = inverse ? 1 : -1
    ret = [
        x[1] + x[2] + x[3] + x[4],
        x[1] + x[2] * twiddle(inv, 1, 4) + x[3] * twiddle(inv, 2, 4) + x[4] * twiddle(inv, 3, 4),
        x[1] + x[2] * twiddle(inv, 2, 4) + x[3] * twiddle(inv, 4, 4) + x[4] * twiddle(inv, 6, 4),
        x[1] + x[2] * twiddle(inv, 3, 4) + x[3] * twiddle(inv, 6, 4) + x[4] * twiddle(inv, 9, 4),
    ]
    if normalize
        ret /= one(eltype(x)) * 4
    end
    return ret
end

dft(x::AbstractArray, args...; NN::Int=length(x), inverse::Bool=false, normalize::Bool=inverse) = begin
    N = length(x)
    inv = inverse ? 1 : -1
    X = zeros(eltype(x), N)
    for k in 1:N
        for n in 1:N
            index = div((n - 1) * (k - 1) * NN, N) % NN
            twid = twiddle(inv, index, NN)
            X[k] += x[n] * twid
        end
        if inverse && normalize
            X[k] /= NN
        end
    end
    return X
end

twiddle(inv, ind, NN) = exp(-inv * 2 * π * 1im * ind / NN) # we could cache this

fft_cooley_tukey(x::AbstractArray, fft_plan::FFTPlan; NN::Int=length(x), inverse::Bool=false, normalize::Bool=false) = begin
    N = length(x)
    X = zeros(eltype(x), N)
    P = fft_plan.plan[1]
    Q = fft_plan.plan[2]
    # inner fft
    for q in 1:Q.f
        res = P.method(x[q:Q.f:end], P; NN, inverse=inverse, normalize=false)
        X[q:Q.f:end] = res
    end
    # mul twiddles
    for q in 1:Q.f
        for s in 1:P.f
            index = div((q - 1) * (s - 1) * NN, N)
            twid = twiddle(inverse ? 1 : -1, index, NN)
            X[Q.f*(s-1)+q] *= twid
        end
    end
    # outer fft
    for s in 1:P.f
        res = Q.method(X[Q.f*(s-1)+1:Q.f*s], Q; NN, inverse=inverse, normalize=false)
        X[Q.f*(s-1)+1:Q.f*s] = res
    end
    ret = zeros(eltype(X), N)
    for s in 1:P.f
        ret[s:P.f:end] = X[(s-1)*Q.f+1:s*Q.f]
    end
    if inverse && normalize
        ret ./= N
    end
    return ret
end

plan_fft(N::Int; method::Symbol=:dit, rad::Symbol=:min) = begin
    @assert (method == :dit) || (method == :dif) "Method needs to be either :dit (decimation-in-time) or :dif (decimation-in-frequency)"
    @assert (rad == :min) || (rad == :max) "Radix mode needs to be either :min or :max"
    if rad == :min
        return plan_fft_min(N; method=method)
    else
        return plan_fft_max(N; method=method)
    end
end

fft(x::AbstractArray; method::Symbol=:dit, rad::Symbol=:min) = begin
    plan = plan_fft(length(x); method, rad)
    fft(x, plan)
end

fft(x::AbstractArray, plan::FFTPlan) = fft_cooley_tukey(x, plan; inverse=false, normalize=false)

ifft(x::AbstractArray; method::Symbol=:dit, rad::Symbol=:min) = begin
    plan = plan_fft(length(x); method, rad)
    ifft(x, plan)
end
ifft(x::AbstractArray, plan::FFTPlan) = fft_cooley_tukey(x, plan; inverse=true, normalize=true)

#=
# Precompilation Statements for Common Types
=#
# precomputed functions
precompile(dft_1, (ComplexF16,))
precompile(dft_1, (ComplexF32,))
precompile(dft_1, (ComplexF64,))

precompile(dft_2, (Vector{ComplexF16},))
precompile(dft_2, (Vector{ComplexF32},))
precompile(dft_2, (Vector{ComplexF64},))

precompile(dft_4, (Vector{ComplexF16},))
precompile(dft_4, (Vector{ComplexF32},))
precompile(dft_4, (Vector{ComplexF64},))

# public interface
precompile(fft, (Vector{ComplexF16},))
precompile(fft, (Vector{ComplexF32},))
precompile(fft, (Vector{ComplexF64},))

precompile(ifft, (Vector{ComplexF16},))
precompile(ifft, (Vector{ComplexF32},))
precompile(ifft, (Vector{ComplexF64},))

precompile(plan_fft, (Int,))

end # module PureFFT
