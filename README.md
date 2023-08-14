# A fast, accurate integer convolution

This is [an entry](https://codegolf.stackexchange.com/a/264008/78410) for Code Golf challenge
[Compute convolution quickly and accurately](https://codegolf.stackexchange.com/q/263748/78410).

The algorithm used is Number Theoretic Transform, which is a variation of FFT that does not suffer from
floating-point inaccuracies. The modulo used is `P = 71 * 2^57 + 1` which is slightly larger than `10^19`,
and I devised a specialized `a * b % P` function [here](https://github.com/Bubbler-4/cg263748/blob/5bb30c4e0f8d89125f4f2272b127f9dd94614f07/src/lib.rs#L8),
which gave ~4x speedup for the entire convolution compared to `(a as u128 * b as u128 % P as u128) as u64`.

The main entrypoint for the algorithm is [`solve_ntt`](https://github.com/Bubbler-4/cg263748/blob/5bb30c4e0f8d89125f4f2272b127f9dd94614f07/src/lib.rs#L291).

Additionally, three variants of the same algorithm is implemented as
[`solve_ntt_par`](https://github.com/Bubbler-4/cg263748/blob/5bb30c4e0f8d89125f4f2272b127f9dd94614f07/src/lib.rs#L313),
[`solve_ntt_par2`](https://github.com/Bubbler-4/cg263748/blob/5bb30c4e0f8d89125f4f2272b127f9dd94614f07/src/lib.rs#L335),
and [`solve_ntt_par3`](https://github.com/Bubbler-4/cg263748/blob/5bb30c4e0f8d89125f4f2272b127f9dd94614f07/src/lib.rs#L357),
which are some attempts at applying data-independent parallelism using `rayon`.
(Trying to parallelize individual butterflying turns out to be harmful, apparently because it causes much more cache misses.
Will it perform better if cache-line-sized chunks are used instead?)

Use `cargo bench` to run benchmarks locally.

I don't claim this to be the fastest code that supports up to `10^19` in the result of convolution.
Even the core algorithm itself (Cooley-Tukey) can be substituted for different ones such as Stockham for significant improvement.