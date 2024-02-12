use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::fmindex::BackwardSearchResult;
use clap::{Command, Arg};
use jseqio::reader::SeqRecordProducer;
use rand::seq::index::sample;
use rand::{self, Rng, RngCore, SeedableRng};
use rand::rngs::StdRng;

fn count(index: &FMIndex<&Vec<u8>, &Vec<usize>, &Occ>, pattern: &[u8]) -> usize {
    let res = index.backward_search(pattern.iter());
    let interval = match res {
        BackwardSearchResult::Complete(sai) => Some(sai.lower..sai.upper),
        BackwardSearchResult::Partial(sai, _l) => None,
        BackwardSearchResult::Absent => None,
    };

    if let Some(I) = interval{
        I.len()
    } else { 0 }
} 

fn generate_random_dna_string(length: usize, seed: u64) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(seed);
    let alphabet: Vec<u8> = vec![b'A', b'C', b'G', b'T'];

    let mut s = Vec::<u8>::new();
    for _ in 0..length{
        s.push(alphabet[(rng.next_u64() % 4) as usize]);
    }
    s
    
}

fn main() {

    let n = 1000_000_usize;
    let seed = 1234;
    let mut data = generate_random_dna_string(n, seed);

    data.push(b'$'); // The BWT library needs this?

    // Build the FM index
    eprintln!("Building FM index");
    let alphabet = dna::n_alphabet();
    let sa = suffix_array(data.as_slice());
    let bwt = bwt(data.as_slice(), &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let index: FMIndex<&Vec<u8>, &Vec<usize>, &Occ> = FMIndex::new(&bwt, &less, &occ);
    eprintln!("FM index built");

    let mut sampled_endpoints = Vec::<usize>::new();
    let mut lengths = Vec::<usize>::new();
    let k = 63;
    for i in 0..n-k+1{
        let kmer = &data[i..i+k];
        let mut finimizer: Option<&[u8]> = None;
        let mut f_start = 0;
        let mut f_end = kmer.len();
        for start in 0..kmer.len(){
            for end in start+1..kmer.len(){
                let x = &kmer[start..end];
                let freq = count(&index, x);
                if freq == 1 {
                    match finimizer{
                        None => {
                            finimizer = Some(x);
                            f_start = i + start;
                            f_end = i + end;
                        },
                        Some(cur) => {
                            if x.len() < cur.len() || (x.len() == cur.len() && x < cur){
                                finimizer = Some(x);
                                f_start = i + start;
                                f_end = i + end;
                            }
                        }
                    }
                    break; // No reason to check further endpoints because this one was unique already
                }
            }
        }
        match sampled_endpoints.last(){
            Some(prev) => {
                if *prev != f_end{
                    sampled_endpoints.push(f_end);
                    lengths.push(f_end - f_start);
                }
            },
            None => {
                sampled_endpoints.push(f_end);
                lengths.push(f_end - f_start);
            }
        }

        println!("{}, {}, {}", f_end, String::from_utf8(kmer.to_vec()).unwrap(), String::from_utf8(finimizer.unwrap().to_vec()).unwrap());
    }

    assert_eq!(sampled_endpoints.len(), lengths.len());
    println!("{}/{}", sampled_endpoints.len(), n);
    println!("Mean length: {}", lengths.iter().sum::<usize>() as f64 / lengths.len() as f64);


}
