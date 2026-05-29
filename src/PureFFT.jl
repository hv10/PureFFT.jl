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

#=
Helper Functions
=#
get_factor_arr(N::Int) =
    mapreduce(vcat, eachfactor(N)) do (f, c)
        repeat([f], c)
    end

twiddle(inv, ind, NN) = exp(-inv * 2 * π * 1im * ind / NN) # we could cache this

#= 
DFT Implementations
=#
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
    @inbounds for k in 1:N
        @inbounds for n in 1:N
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

#=
Codelet Dispatch
=#
const CODELET_MAP = Dict{Int,Function}(
    1 => dft_1,
    2 => dft_2,
    4 => dft_4,
)
# Sorted largest-first; updated by retune() / tuned_codelets.jl
const CODELET_SIZES = Ref{Vector{Int}}([4, 2, 1])

"""
    get_codelet(r) -> FFTPlan

Look up the best available codelet for size `r`.
Falls back to the generic `dft` for any size not in CODELET_MAP.
"""
get_codelet(r::Int)::FFTPlan = FFTPlan(r, get(CODELET_MAP, r, dft))

#=
Planning Strategies
=#
plan_fft_min(N::Int; method::Symbol=:dit)::FFTPlan = begin
    @assert (method == :dit) || (method == :dif) "Method must be :dit or :dif"

    N_fcts = sort(get_factor_arr(N), rev=true)   # largest primes first
    plan = FFTPlan(N, fft_cooley_tukey)
    NN = N
    ptr = plan.plan

    for f in N_fcts
        NN = div(NN, f)
        P = get_codelet(f)
        Q = NN <= 4 || haskey(CODELET_MAP, NN) ? get_codelet(NN) : FFTPlan(NN, fft_cooley_tukey)

        method == :dit ? push!(ptr, P, Q) : push!(ptr, Q, P) # recursion remainder in back, small radix first.

        (NN == 1 || haskey(CODELET_MAP, NN)) && break
        ptr = Q.plan
    end

    return plan
end

plan_fft_max(N::Int; method::Symbol=:dit)::FFTPlan = begin
    @assert (method == :dit) || (method == :dif) "Method must be :dit or :dif"

    plan = FFTPlan(N, fft_cooley_tukey)
    NN = N
    ptr = plan.plan

    while NN > 1
        divisors = [c for c in CODELET_SIZES[] if c > 1 && NN % c == 0]

        r = !isempty(divisors) ?
            maximum(divisors) :
            minimum(p for p in keys(factor(NN)) if p > 1)

        NN_next = div(NN, r)

        @assert NN_next < NN || NN == 1

        P = get_codelet(r)

        Q = (NN_next == 1 || haskey(CODELET_MAP, NN_next)) ?
            get_codelet(NN_next) :
            FFTPlan(NN_next, fft_cooley_tukey)

        method == :dit ? push!(ptr, Q, P) : push!(ptr, P, Q)

        (NN_next <= 1 || haskey(CODELET_MAP, NN_next)) && break

        ptr = Q.plan
        NN = NN_next
    end

    return plan
end


fft_cooley_tukey(x::AbstractArray, fft_plan::FFTPlan; NN::Int=length(x), inverse::Bool=false, normalize::Bool=false) = begin
    N = length(x)
    X = zeros(eltype(x), N)
    P = fft_plan.plan[1]
    Q = fft_plan.plan[2]
    # inner fft
    for q in 1:Q.f
        X[q:Q.f:end] = P.method(@view(x[q:Q.f:end]), P; NN, inverse=inverse, normalize=false)
    end
    # mul twiddles
    @inbounds for q in 1:Q.f, s in 1:P.f
        index = div((q - 1) * (s - 1) * NN, N)
        X[Q.f*(s-1)+q] *= twiddle(inverse ? 1 : -1, index, NN)
    end
    # outer fft
    @inbounds for s in 1:P.f
        X[Q.f*(s-1)+1:Q.f*s] = Q.method(@view(X[Q.f*(s-1)+1:Q.f*s]), Q; NN, inverse=inverse, normalize=false)
    end
    ret = vec(@views permutedims(reshape(X, Q.f, P.f)))
    inverse && normalize && (ret ./= N)
    return ret
end

plan_fft(N::Int; method::Symbol=:dit, rad::Symbol=:min) = begin
    @assert (method == :dit) || (method == :dif) "Method needs to be either :dit (decimation-in-time) or :dif (decimation-in-frequency)"
    @assert (rad == :min) || (rad == :max) "Radix mode needs to be either :min or :max"
    return rad == :min ? plan_fft_min(N; method=method) : plan_fft_max(N; method=method)
end

fft(x::AbstractArray; method::Symbol=:dit, rad::Symbol=:min) = fft(x, plan_fft(length(x); method, rad))
fft(x::AbstractArray, plan::FFTPlan) = fft_cooley_tukey(x, plan; inverse=false, normalize=false)

ifft(x::AbstractArray; method::Symbol=:dit, rad::Symbol=:min) = ifft(x, plan_fft(length(x); method, rad))
ifft(x::AbstractArray, plan::FFTPlan) = fft_cooley_tukey(x, plan; inverse=true, normalize=true)

#=
Machine-Specific Tuning
=#
include("tuning.jl")

# Load machine-tuned codelets if they exist; otherwise stay with hand-written fallbacks
const TUNED_FILE = joinpath(@__DIR__, "generated" , "tuned_codelets.jl")

if isfile(TUNED_FILE)
    include(TUNED_FILE)
    if TUNED_CPU != Sys.CPU_NAME
        @warn "PureFFT tuning cache was built for a different CPU. Run PureFFT.retune()."
    end
    merge!(CODELET_MAP, TUNED_CODELET_MAP)
    CODELET_SIZES[] = TUNED_CODELET_SIZES
    @info "PureFFT: loaded tuned codelets, sizes = $(TUNED_CODELET_SIZES)"
else
    @info "PureFFT: no tuned codelets found — run PureFFT.retune() to generate them"
end


#=
# Precompilation Statements for Common Types
=#
for T in [Vector{ComplexF64}, Vector{ComplexF32}, Vector{ComplexF16}]
    # codelets
    precompile(dft_1, (T,))
    precompile(dft_2, (T,))
    precompile(dft_4, (T,))
    # public interface
    precompile(fft, (T,))
    precompile(ifft, (T,))
end

precompile(plan_fft, (Int,))

end # module PureFFT
