use num_bigint::BigUint;
use dashu::integer::UBig;
use std::time::Instant;
use rayon::prelude::*;

pub const P: u64 = 10232178353385766913;
pub const P128: u128 = 10232178353385766913;
pub fn mulmod(a: u64, b: u64) -> u64 {
    let x = a as u128 * b as u128;
    let y = x >> 57;
    let y1 = (y >> 35) as u64;
    let y2 = y as u64 & ((1u64 << 35) - 1);
    let y_div_71 = y1 * 483939977 + (y1 + y2) / 71;
    let sub = y_div_71 as u128 * P as u128;
    let (res, borrow) = x.overflowing_sub(sub);
    if borrow { (res as u64).wrapping_add(P) } else { res as u64 }
}

pub fn mulmod128(a: u128, b: u128) -> u128 {
    a * b % P128
}

fn powmod(a: u64, p: u64) -> u64 {
    let mut cur = 1;
    let mut pow = a;
    let mut p = p;
    while p > 0 {
        if p % 2 > 0 {
            cur = mulmod(cur, pow);
        }
        p /= 2;
        pow = mulmod(pow, pow);
    }
    cur
}

fn powmod128(a: u128, p: u128) -> u128 {
    let mut cur = 1;
    let mut pow = a;
    let mut p = p;
    while p > 0 {
        if p % 2 > 0 {
            cur = mulmod128(cur, pow);
        }
        p /= 2;
        pow = mulmod128(pow, pow);
    }
    cur
}

fn ntt<const INV: bool>(input: &mut [u64], omega_table: &[u64], inv_p2: u64) {
    // length is a power of 2
    let len = input.len();
    let l = len.trailing_zeros() as usize;
    for i in 1..len {
        let j = i.reverse_bits() >> (64 - l);
        if i < j { input.swap(i, j); }
    }

    let mut root_pow = len / 2;
    for intvl_shift in 1..=l {
        let intvl = 1usize << intvl_shift;
        input.chunks_exact_mut(intvl).for_each(|chunk| {
            let mut root_idx = 0;
            let (left, right) = chunk.split_at_mut(intvl / 2);
            for (u, v) in left.into_iter().zip(right) {
                let u2 = *u;
                let v2 = mulmod(*v, omega_table[root_idx]);
                let (mut x, overflow) = u2.overflowing_sub(P - v2);
                if overflow { x = x.wrapping_add(P); }
                *u = x;
                let (mut y, overflow) = u2.overflowing_sub(v2);
                if overflow { y = y.wrapping_add(P); }
                *v = y;
                root_idx = root_idx + root_pow; // & len - 1;
            }
        });
        root_pow /= 2;
    }
    if INV {
        for x in input.into_iter() { *x = mulmod(*x, inv_p2); }
    }
}

fn ntt_par<const INV: bool>(input: &mut [u64], omega_table: &[u64], inv_p2: u64) {
    // length is a power of 2
    let len = input.len();
    let l = len.trailing_zeros() as usize;
    for i in 1..len {
        let j = i.reverse_bits() >> (64 - l);
        if i < j { input.swap(i, j); }
    }

    let mut root_pow = len / 2;
    for intvl_shift in 1..=l {
        let intvl = 1usize << intvl_shift;
        input.par_chunks_exact_mut(intvl).for_each(|chunk| {
            let mut root_idx = 0;
            let (left, right) = chunk.split_at_mut(intvl / 2);
            for (u, v) in left.into_iter().zip(right) {
                let u2 = *u;
                let v2 = mulmod(*v, omega_table[root_idx]);
                let (mut x, overflow) = u2.overflowing_sub(P - v2);
                if overflow { x = x.wrapping_add(P); }
                *u = x;
                let (mut y, overflow) = u2.overflowing_sub(v2);
                if overflow { y = y.wrapping_add(P); }
                *v = y;
                root_idx = root_idx + root_pow; // & len - 1;
            }
        });
        root_pow /= 2;
    }
    if INV {
        input.into_par_iter().for_each(|x| *x = mulmod(*x, inv_p2));
    }
}

fn ntt_par2<const INV: bool>(input: &mut [u64], omega_table: &[u64], inv_p2: u64) {
    // length is a power of 2
    let len = input.len();
    let l = len.trailing_zeros() as usize;
    for i in 1..len {
        let j = i.reverse_bits() >> (64 - l);
        if i < j { input.swap(i, j); }
    }

    let mut root_pow = len / 2;
    for intvl_shift in 1..=l {
        let intvl = 1usize << intvl_shift;
        input.par_chunks_exact_mut(intvl).flat_map(|chunk| {
            let (left, right) = chunk.split_at_mut(intvl / 2);
            let roots = (&omega_table).into_par_iter().step_by(root_pow);
            (left, right, roots)
        }).for_each(|(u, v, root)| {
            let u2 = *u;
            let v2 = mulmod(*v, *root);
            let (mut x, overflow) = u2.overflowing_sub(P - v2);
            if overflow { x = x.wrapping_add(P); }
            *u = x;
            let (mut y, overflow) = u2.overflowing_sub(v2);
            if overflow { y = y.wrapping_add(P); }
            *v = y;
        });
        root_pow /= 2;
    }
    if INV {
        input.into_par_iter().for_each(|x| *x = mulmod(*x, inv_p2));
    }
}

fn ntt_par3<const INV: bool>(input: &mut [u64], omega_table: &[u64], inv_p2: u64) {
    // length is a power of 2
    let len = input.len();
    let l = len.trailing_zeros() as usize;
    for i in 1..len {
        let j = i.reverse_bits() >> (64 - l);
        if i < j { input.swap(i, j); }
    }

    let mut root_pow = len / 2;
    for intvl_shift in 1..=l/2 {
        let intvl = 1usize << intvl_shift;
        input.par_chunks_exact_mut(intvl).for_each(|chunk| {
            let mut root_idx = 0;
            let (left, right) = chunk.split_at_mut(intvl / 2);
            for (u, v) in left.into_iter().zip(right) {
                let u2 = *u;
                let v2 = mulmod(*v, omega_table[root_idx]);
                let (mut x, overflow) = u2.overflowing_sub(P - v2);
                if overflow { x = x.wrapping_add(P); }
                *u = x;
                let (mut y, overflow) = u2.overflowing_sub(v2);
                if overflow { y = y.wrapping_add(P); }
                *v = y;
                root_idx = root_idx + root_pow; // & len - 1;
            }
        });
        root_pow /= 2;
    }
    for intvl_shift in l/2+1..=l {
        let intvl = 1usize << intvl_shift;
        input.chunks_exact_mut(intvl).for_each(|chunk| {
            let (left, right) = chunk.split_at_mut(intvl / 2);
            let roots = (&omega_table).into_par_iter().step_by(root_pow);
            (left, right, roots).into_par_iter().for_each(|(u, v, root)| {
                let u2 = *u;
                let v2 = mulmod(*v, *root);
                let (mut x, overflow) = u2.overflowing_sub(P - v2);
                if overflow { x = x.wrapping_add(P); }
                *u = x;
                let (mut y, overflow) = u2.overflowing_sub(v2);
                if overflow { y = y.wrapping_add(P); }
                *v = y;
            });
        });
        root_pow /= 2;
    }
    if INV {
        input.into_par_iter().for_each(|x| *x = mulmod(*x, inv_p2));
    }
}

fn ntt128(input: &mut [u128], invert: bool, omega_table: &[u128], _omega_table_prime: &[u128], inv_p2: u128) {
    // length is a power of 2
    let len = input.len();
    let l = len.trailing_zeros() as usize;
    let mut j = 0;
    for i in 1..len {
        let mut bit = len >> 1;
        while (j & bit) != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j { input.swap(i, j); }
    }

    let mut root_pow = len / 2;
    for intvl_shift in 1..=l {
        let intvl = 1usize << intvl_shift;
        for i in (0..len>>intvl_shift).map(|x| x << intvl_shift) {
            let mut root_idx = 0;
            for j in 0..intvl / 2 {
                let u = input[i+j];
                let v = input[i+j+intvl/2] * omega_table[root_idx] % P128;
                // input[i+j] = (u+v) % P128;
                let mut x = u+v;
                if x >= P128 { x -= P128; }
                input[i+j] = x;
                // input[i+j+intvl/2] = (u + P128 - v) % P128;
                let mut y = u + P128 - v;
                if y >= P128 { y -= P128; }
                input[i+j+intvl/2] = y;
                root_idx = root_idx + root_pow & len - 1;
            }
        }
        if invert { root_pow += len; }
        root_pow /= 2;
    }
    if invert {
        for x in input.into_iter() { *x = mulmod128(*x, inv_p2); }
    }
}

/// (omega_table, inv_p2_table)
fn ntt_precompute(len: usize) -> (Vec<u64>, Vec<u64>) {
    let l = len.trailing_zeros();

    let primitive_root = 3;
    let omega = powmod(primitive_root, P >> l);
    let mut omega_table = vec![1];
    // let mut omega_table_prime = vec![((1u128 << 64) / P as u128) as u64];
    for i in 1..len {
        // omega_table.push((omega_table[i-1] as u128 * omega as u128 % P as u128) as u64);
        omega_table.push(mulmod(omega_table[i-1], omega));
        // omega_table_prime.push((((omega_table[i] as u128) << 64) / P as u128) as u64);
    }
    let mut inv_p2 = 1;
    let mut inv_p2_table = vec![1];
    for _ in 0..l {
        if inv_p2 % 2 == 0 { inv_p2 /= 2; }
        else { inv_p2 = inv_p2 / 2 + P / 2 + 1; }
        inv_p2_table.push(inv_p2);
    }
    (omega_table, inv_p2_table)
}

/// (omega_table, omega_table_prime, inv_p2_table)
fn ntt128_precompute(len: usize) -> (Vec<u128>, Vec<u128>, Vec<u128>) {
    let l = len.trailing_zeros();
    let primitive_root = 3;
    let omega = powmod128(primitive_root, P128 >> l);
    let mut omega_table = vec![1];
    let mut omega_table_prime = vec![((1u128 << 64) / P128)];
    for i in 1..len {
        omega_table.push(omega_table[i-1] as u128 * omega as u128 % P128);
        omega_table_prime.push(((omega_table[i] as u128) << 64) / P128);
    }
    let mut inv_p2 = 1;
    let mut inv_p2_table = vec![1];
    for _ in 0..l {
        if inv_p2 % 2 == 0 { inv_p2 /= 2; }
        else { inv_p2 = (inv_p2 + P128) / 2; }
        inv_p2_table.push(inv_p2);
    }
    // println!("{:?} {:?}", omega_table, inv_p2_table);
    (omega_table, omega_table_prime, inv_p2_table)
}

pub fn solve_ntt(a: &[u64]) -> (Vec<u64>, u128) {
    let instant = Instant::now();

    let len_max = a.len().next_power_of_two() * 2;
    let (mut omega_table, inv_p2_table) = ntt_precompute(len_max);
    let inv_p2 = inv_p2_table[len_max.trailing_zeros() as usize];

    let mut a2 = a.to_vec();
    a2.resize(len_max, 0);
    ntt::<false>(&mut a2, &omega_table, inv_p2);
    for ax in a2.iter_mut() {
        *ax = mulmod(*ax, *ax);
    }
    omega_table[1..].reverse();
    ntt::<true>(&mut a2, &omega_table, inv_p2);
    a2.truncate((a.len() * 2).saturating_sub(1));

    let elapsed = instant.elapsed().as_micros();

    (a2, elapsed)
}

pub fn solve_ntt_par(a: &[u64]) -> (Vec<u64>, u128) {
    let instant = Instant::now();
    
    let len_max = a.len().next_power_of_two() * 2;
    let (mut omega_table, inv_p2_table) = ntt_precompute(len_max);
    let inv_p2 = inv_p2_table[len_max.trailing_zeros() as usize];

    let mut a2 = a.to_vec();
    a2.resize(len_max, 0);
    ntt_par::<false>(&mut a2, &omega_table, inv_p2);
    for ax in a2.iter_mut() {
        *ax = mulmod(*ax, *ax);
    }
    omega_table[1..].reverse();
    ntt_par::<true>(&mut a2, &omega_table, inv_p2);
    a2.truncate((a.len() * 2).saturating_sub(1));

    let elapsed = instant.elapsed().as_micros();

    (a2, elapsed)
}

pub fn solve_ntt_par2(a: &[u64]) -> (Vec<u64>, u128) {
    let instant = Instant::now();
    
    let len_max = a.len().next_power_of_two() * 2;
    let (mut omega_table, inv_p2_table) = ntt_precompute(len_max);
    let inv_p2 = inv_p2_table[len_max.trailing_zeros() as usize];

    let mut a2 = a.to_vec();
    a2.resize(len_max, 0);
    ntt_par2::<false>(&mut a2, &omega_table, inv_p2);
    for ax in a2.iter_mut() {
        *ax = mulmod(*ax, *ax);
    }
    omega_table[1..].reverse();
    ntt_par2::<true>(&mut a2, &omega_table, inv_p2);
    a2.truncate((a.len() * 2).saturating_sub(1));

    let elapsed = instant.elapsed().as_micros();

    (a2, elapsed)
}

pub fn solve_ntt_par3(a: &[u64]) -> (Vec<u64>, u128) {
    let instant = Instant::now();
    
    let len_max = a.len().next_power_of_two() * 2;
    let (mut omega_table, inv_p2_table) = ntt_precompute(len_max);
    let inv_p2 = inv_p2_table[len_max.trailing_zeros() as usize];

    let mut a2 = a.to_vec();
    a2.resize(len_max, 0);
    ntt_par3::<false>(&mut a2, &omega_table, inv_p2);
    for ax in a2.iter_mut() {
        *ax = mulmod(*ax, *ax);
    }
    omega_table[1..].reverse();
    ntt_par3::<true>(&mut a2, &omega_table, inv_p2);
    a2.truncate((a.len() * 2).saturating_sub(1));

    let elapsed = instant.elapsed().as_micros();

    (a2, elapsed)
}

pub fn solve_ntt128(a: &[u64]) -> (Vec<u64>, u128) {
    let len_max = a.len().next_power_of_two() * 2;
    let (omega_table, omega_table_prime, inv_p2_table) = ntt128_precompute(len_max);
    let inv_p2 = inv_p2_table[len_max.trailing_zeros() as usize];

    let instant = Instant::now();

    let mut a2 = a.into_iter().map(|&x| x as u128).collect::<Vec<_>>();
    a2.resize(len_max, 0);
    ntt128(&mut a2, false, &omega_table, &omega_table_prime, inv_p2);
    for ax in a2.iter_mut() {
        *ax = mulmod128(*ax, *ax);
    }
    ntt128(&mut a2, true, &omega_table, &omega_table_prime, inv_p2);
    a2.truncate((a.len() * 2).saturating_sub(1));
    let a2 = a2.into_iter().map(|x| x as u64).collect();

    let elapsed = instant.elapsed().as_micros();

    (a2, elapsed)
}

pub fn solve_num(input: &[u64]) -> (Vec<u64>, u128) {
    let mut num = Vec::with_capacity(input.len() * 2);
    for &x in input {
        num.push(x as u32);
        num.push(0);
    }
    let big = BigUint::from_slice(&num);
    let instant = Instant::now();
    let big2 = &big * &big;
    let elapsed = instant.elapsed().as_micros();
    (big2.to_u64_digits(), elapsed)
}

pub fn solve_dashu(input: &[u64]) -> (Vec<u64>, u128) {
    let big = UBig::from_words(input);
    let instant = Instant::now();
    let square = big.square();
    let elapsed = instant.elapsed().as_micros();
    (square.as_words().to_vec(), elapsed)
}