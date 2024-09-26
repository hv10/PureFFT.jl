module PureFFT

function pad_to_power_of_two(a::AbstractArray{<:Number})
    n = length(a)
    next_pow = 2^ceil(Int, log2(n))  # Find the next power of 2
    if next_pow != n
        return vcat(a, zeros(eltype(a), next_pow - n))  # Pad with zeros
    else
        return a
    end
end

function _fft!(a::AbstractVector{<:Complex}, invert::Bool)
    n = length(a)
    if n == 1
        return
    end

    a0 = a[1:2:end]  # Get even indices
    a1 = a[2:2:end]  # Get odd indices

    _fft!(a0, invert)
    _fft!(a1, invert)

    ang = 2 * π / n * (invert ? -1 : 1)
    w = Complex(1.0, 0.0)
    wn = Complex(cos(ang), sin(ang))

    for i in 1:(n ÷ 2)
        a[i] = a0[i] + w * a1[i]
        a[i + n ÷ 2] = a0[i] - w * a1[i]
        if invert
            a[i] /= 2
            a[i + n ÷ 2] /= 2
        end
        w *= wn
    end
end

function fft(a::AbstractVector, invert::Bool = false)
    # Pad the array if necessary
    padded_a = pad_to_power_of_two(a)
	padded_a = padded_a .+ 0im
    
    # Perform the FFT in place
    _fft!(padded_a, invert)
    
    # # If performing inverse FFT, scale by the array size
    if invert
        padded_a ./= length(padded_a)
    end
    
    return padded_a[1:length(a)]
end


end # module PureFFT
