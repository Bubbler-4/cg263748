use cg263748::*;

fn main() {
    let input = (
        vec![10_000_000u64; 100_000],
        vec![1; 5]
    ).0;

    let (num_result, elapsed) = solve_num(&input);
    println!("elapsed: {}.{:06}", elapsed / 1000000, elapsed % 1000000);

    let (dashu_result, elapsed) = solve_dashu(&input);
    println!("elapsed: {}.{:06}", elapsed / 1000000, elapsed % 1000000);
    assert_eq!(num_result, dashu_result);

    let (ntt128_result, elapsed) = solve_ntt128(&input);
    println!("elapsed: {}.{:06}", elapsed / 1000000, elapsed % 1000000);
    assert_eq!(num_result, ntt128_result);

    println!("solve_ntt:");
    for _ in 0..20 {
        let (ntt_result, elapsed) = solve_ntt(&input);
        print!("{}.{:06} ", elapsed / 1000000, elapsed % 1000000);
        assert_eq!(num_result, ntt_result);
    }
    println!();

    println!("solve_ntt_par:");
    for _ in 0..20 {
        let (ntt_result, elapsed) = solve_ntt_par(&input);
        print!("{}.{:06} ", elapsed / 1000000, elapsed % 1000000);
        assert_eq!(num_result, ntt_result);
    }
    println!();

    println!("solve_ntt_par2:");
    for _ in 0..20 {
        let (ntt_result, elapsed) = solve_ntt_par2(&input);
        print!("{}.{:06} ", elapsed / 1000000, elapsed % 1000000);
        assert_eq!(num_result, ntt_result);
    }
    println!();

    println!("solve_ntt_par3:");
    for _ in 0..20 {
        let (ntt_result, elapsed) = solve_ntt_par3(&input);
        print!("{}.{:06} ", elapsed / 1000000, elapsed % 1000000);
        assert_eq!(num_result, ntt_result);
    }
    println!();
}

/*
P = 71 * 2^57 + 1
x <= P^2-2P+1
x % P = x - floor(x / P) * P
approximating floor(x / P) with floor(x / (P-1))
x / P < i, i + 1 <= x / (P-1)
x < Pi, x >= (P-1)(i+1) = Pi + P - i - 1
Pi + P - i - 1 <= x < Pi
i > P-1 -> contradiction because x / P <= P-2+1/P; i <= P-1

floor(x / (P-1)) = (x >> 57) / 71
y = x >> 57
y / 71 = (y1 * 2^35 + y2) / 71 = y1 * ((2^35-1)/71) + (y1+y2) / 71
 */

#[cfg(test)]
mod tests {
    use quickcheck::quickcheck;
    use crate::{P, mulmod, solve_ntt, solve_ntt128, solve_num};
    quickcheck! {
        fn prop(a: u64, b: u64) -> bool {
            mulmod(a % P, b % P) == (a as u128 * b as u128 % P as u128) as u64
        }

        fn prop2(a: Vec<u64>) -> bool {
            let a = a.into_iter().map(|x| x % 10000001).take(100000).collect::<Vec<_>>();
            let mut ntt = solve_ntt128(&a).0;
            let num = solve_num(&a).0;
            while ntt.last() == Some(&0) { ntt.pop(); }
            ntt == num
        }

        fn prop3(a: Vec<u64>) -> bool {
            let a = a.into_iter().map(|x| x % 10000001).take(100000).collect::<Vec<_>>();
            let mut ntt = solve_ntt(&a).0;
            let num = solve_num(&a).0;
            while ntt.last() == Some(&0) { ntt.pop(); }
            ntt == num
        }
    }

    #[test]
    fn fixed() {
        let arr = [0, 1, P-1, P/2, P/2+1, 1234567, P-1234568];
        for a in arr {
            for b in arr {
                assert_eq!(mulmod(a, b), (a as u128 * b as u128 % P as u128) as u64);
            }
        }
    }
}