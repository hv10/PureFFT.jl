@testsnippet MakeData begin
    using Statistics
    make_data(len) =
        map(1:len) do t
            sinpi(t / 3.0) + 0im
        end
end

@testitem "Check twiddle Factors for small numbers" begin
    # check non-inverted
    @test PureFFT.twiddle(1, 0, 4) ≈ 1 + 0im
    @test PureFFT.twiddle(1, 1, 4) ≈ 0 - 1im
    @test PureFFT.twiddle(1, 2, 4) ≈ -1 + 0im
    @test PureFFT.twiddle(1, 3, 4) ≈ 0 + 1im
end

@testitem "FFT size 2" setup = [MakeData] begin
    a = make_data(2)
    plan = PureFFT.plan_fft_min(2; method=:dit)
    a_dft = PureFFT.dft(a)
    a_dft_r = PureFFT.dft(PureFFT.dft(a); inverse=true, normalize=true)
    # DFT^-1(DFT(a))=a
    @test a ≈ a_dft_r
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
end

@testitem "FFT size 3" setup = [MakeData] begin
    a = make_data(3)
    plan = PureFFT.plan_fft_min(3; method=:dit)
    a_dft = PureFFT.dft(a)
    a_dft_r = PureFFT.dft(PureFFT.dft(a); inverse=true, normalize=true)
    # DFT^-1(DFT(a))=a
    @test a ≈ a_dft_r
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
end

@testitem "FFT size 5" setup = [MakeData] begin
    a = make_data(5)
    plan = PureFFT.plan_fft_min(5; method=:dit)
    a_dft = PureFFT.dft(a)
    a_dft_r = PureFFT.dft(PureFFT.dft(a); inverse=true, normalize=true)
    # DFT^-1(DFT(a))=a
    @test a ≈ a_dft_r
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
end

@testitem "FFT size 7" setup = [MakeData] begin
    a = make_data(7)
    plan = PureFFT.plan_fft_min(7; method=:dit)
    a_dft = PureFFT.dft(a)
    a_dft_r = PureFFT.dft(PureFFT.dft(a); inverse=true, normalize=true)
    # DFT^-1(DFT(a))=a
    @test a ≈ a_dft_r
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
end

@testitem "FFT size 4" setup = [MakeData] begin
    a = make_data(4)
    plan = PureFFT.plan_fft_min(4; method=:dit)
    display(plan)
    a_dft = PureFFT.dft(a)
    a_dft_r = PureFFT.dft(PureFFT.dft(a); inverse=true, normalize=true)
    # DFT^-1(DFT(a))=a
    @test a ≈ a_dft_r
    a_ct_r = PureFFT.fft_cooley_tukey(PureFFT.dft(a), plan; inverse=true, normalize=true)
    @test a ≈ a_ct_r
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
end

@testitem "FFT size 6" setup = [MakeData] begin
    a = make_data(6)
    plan = PureFFT.plan_fft_min(6; method=:dit)
    display(plan)
    a_dft = PureFFT.dft(a)
    a_dft_r = PureFFT.dft(PureFFT.dft(a); inverse=true, normalize=true)
    # DFT^-1(DFT(a))=a
    @test a ≈ a_dft_r
    a_ct_r = PureFFT.fft_cooley_tukey(PureFFT.dft(a), plan; inverse=true, normalize=true)
    @test a ≈ a_ct_r
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
end

@testitem "DFT_1 matches DFT" setup = [MakeData] begin
    a = make_data(1)
    a_dft = PureFFT.dft(a)
    a_dft_1 = PureFFT.dft_1(a)
    @test a_dft ≈ a_dft_1
end

@testitem "DFT_2 matches DFT" setup = [MakeData] begin
    a = make_data(2)
    a_dft = PureFFT.dft(a)
    a_dft_2 = PureFFT.dft_2(a)
    @test a_dft ≈ a_dft_2
end

@testitem "DFT_4 matches DFT" setup = [MakeData] begin
    a = make_data(4)
    a_dft = PureFFT.dft(a)
    a_dft_4 = PureFFT.dft_4(a)
    @test a_dft ≈ a_dft_4
end


@testitem "Check Bigger Input" setup = [MakeData] begin
    a = make_data(360)
    plan = PureFFT.plan_fft(360)
    display(plan)
    a_dft = PureFFT.dft(a)
    a_ctfft = PureFFT.fft_cooley_tukey(a, plan)
    @test a_dft ≈ a_ctfft
    a_ctfft_r = PureFFT.fft_cooley_tukey(a_ctfft, plan; inverse=true, normalize=true)
    @test a ≈ a_ctfft_r
end