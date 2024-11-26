# PureFFT
> An implementation of FFT-Algorithm for General Factorizations by Cooley-Tukey, purely in Julia, licensed permissively.
> Contributions welcome.

## Why?
I am working on some projects which will be using FFT internally, but which will not be able to use GPL licensed code.
Most projects I am aware of use FFTW (dual licensed GPL + Commercial by MIT) or other GPL licensed implementations of FFT.

I wanted to learn how the Algorithm works, so I spend a few afternoons on understanding the literature about it.
Then I tried implementing it from scratch.

Long term, I would like to make this a drop in replacement for other FFT Packages, therefore we might need to commmit to matching an existing API.

## What?
This package provides:
- a $O(n^2)$ DFT implementation which works on arbitrary sized arrays
- a Type for holding a FFT-Plan (`FFTPlan`)
- a `plan_fft` function that plans the FFT for a known length
  - supporting Decimation-in-Frequency (`method=:dif`) and Decimation-in-Time (`method=:dit`)
  - supporting the Min-Radix planning style (`rad=:min`)
- a Cooley-Tukey FFT implementation for applying such plans to an array
- base-cases for arrays of sizes `[1,2,4]` to further improve execution speed

## TODO's
- [ ] Register as Julia Package
- [ ] Performance Benchmarking & Optimization
- [ ] Type Stability Checking
- [ ] Edge Case Hunting
- [ ] Max-Radix Plan not yet implemented
- [ ] Exhaustive Testing for Planning Method+Approach Combinations
- [ ] Write the one-stop shop `fft` function combining planning and execution.
