@testitem "fft function" setup = [MakeData] begin
    a = make_data(48)
    plan = PureFFT.plan_fft(48; method=:dit, rad=:min)
    a_fft = PureFFT.fft(a)
    a_fft_pplan = PureFFT.fft(a, plan)
    @test a_fft ≈ a_fft_pplan
end

@testitem "ifft function" setup = [MakeData] begin
    a = make_data(48)
    plan = PureFFT.plan_fft(48; method=:dit, rad=:min)
    a_fft = PureFFT.fft(a, plan)
    a_fft_r = PureFFT.ifft(a)
    a_fft_pplan_r = PureFFT.ifft(a, plan)
    @test a_fft_r ≈ a_fft_pplan_r
end