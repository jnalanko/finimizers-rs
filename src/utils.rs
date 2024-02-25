use rand::{self, RngCore, SeedableRng};
use rand::rngs::StdRng;

fn generate_random_dna_string(length: usize, seed: u64) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(seed);
    let alphabet: Vec<u8> = vec![b'A', b'C', b'G', b'T'];

    let mut s = Vec::<u8>::new();
    for _ in 0..length{
        s.push(alphabet[(rng.next_u64() % 4) as usize]);
    }
    s
}