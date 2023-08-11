use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cg263748::*;

fn bench_ntt(c: &mut Criterion) {
    let data = vec![10_000_000; 100000];
    c.bench_function("sequential ntt", |b| b.iter(||
        solve_ntt(black_box(&data))
    ));
}

fn bench_ntt_par(c: &mut Criterion) {
    let data = vec![10_000_000; 100000];
    c.bench_function("parallel ntt", |b| b.iter(||
        solve_ntt_par(black_box(&data))
    ));
}

fn bench_ntt_par2(c: &mut Criterion) {
    let data = vec![10_000_000; 100000];
    c.bench_function("parallel ntt2", |b| b.iter(||
        solve_ntt_par2(black_box(&data))
    ));
}

fn bench_ntt_par3(c: &mut Criterion) {
    let data = vec![10_000_000; 100000];
    c.bench_function("parallel ntt3", |b| b.iter(||
        solve_ntt_par3(black_box(&data))
    ));
}

criterion_group!(benches, bench_ntt, bench_ntt_par, bench_ntt_par2, bench_ntt_par3);
criterion_main!(benches);