dft_1(x::AbstractArray, args...; kwargs...) = x

dft_2(x::AbstractArray, args...; kwargs...) = [x[1] + x[2], x[1] - x[2]]

dft_4(x::AbstractArray, args...; NN::Int=4, inverse::Bool=false, normalize::Bool=false) = begin
    @assert length(x) == 4 "dft_4 is only valid for 4-element arrays"
    inv = inverse ? 1 : -1
    ret = [
        x[1] + x[2] + x[3] + x[4],
        x[1] + x[2] * twiddle(inv, 1, 4) + x[3] * twiddle(inv, 2, 4) + x[4] * twiddle(inv, 3, 4),
        x[1] + x[2] * twiddle(inv, 2, 4) + x[3] * twiddle(inv, 4, 4) + x[4] * twiddle(inv, 6, 4),
        x[1] + x[2] * twiddle(inv, 3, 4) + x[3] * twiddle(inv, 6, 4) + x[4] * twiddle(inv, 9, 4),
    ]
    normalize && (ret ./= 4)
    return ret
end
512, 256, 128, 64, 32, n_f

for T in [Vector{ComplexF64}, Vector{ComplexF32}, Vector{ComplexF16}]
    precompile(dft_1, (T,))
    precompile(dft_2, (T,))
    precompile(dft_4, (T,))
end
