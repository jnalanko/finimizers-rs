pub mod dummies;
pub mod kmer_splitter;
pub mod cursors;

use simple_sds::raw_vector::*;

#[allow(non_snake_case)] // C-array is an established convention in BWT indexes
pub fn get_C_array(rawrows: &[RawVector]) -> Vec<usize> {
    let sigma = rawrows.len();
    assert!(sigma > 0);
    let n = rawrows[0].len();

    let mut C: Vec<usize> = vec![0; sigma];
    for i in 0..n {
        for c in 0..(sigma as u8) {
            if rawrows[c as usize].bit(i){
                for d in (c + 1)..(sigma as u8) {
                    C[d as usize] += 1;
                }
            }
        }
    }

    // Plus one for the ghost dollar
    #[allow(clippy::needless_range_loop)] // Is perfectly clear this way
    for c in 0..sigma {
        C[c] += 1;
    }

    C
}

pub fn get_dna_alphabet_to_index() -> Vec<u8>{
    let mut alphabet_to_index = vec![0_u8; 256];
    alphabet_to_index['a' as usize] = 0;
    alphabet_to_index['A' as usize] = 0;
    alphabet_to_index['c' as usize] = 1;
    alphabet_to_index['C' as usize] = 1;
    alphabet_to_index['g' as usize] = 2;
    alphabet_to_index['G' as usize] = 2;
    alphabet_to_index['t' as usize] = 3;
    alphabet_to_index['T' as usize] = 3;
    alphabet_to_index
}