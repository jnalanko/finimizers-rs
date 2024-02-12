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

fn get_dollar_concatenation(db: &jseqio::seq_db::SeqDB) -> Vec<u8>{
    let mut concat = Vec::<u8>::new();
    for i in 0..db.sequence_count(){
        let rec = db.get(i);
        concat.extend(rec.seq);
        concat.push(b'$');
    }
    concat
}

// Returns pair (finimizer endpoint vector, finimizer length vector)
fn get_finimizers(seq: &[u8], k: usize, index: &FMIndex<&Vec<u8>, &Vec<usize>, &Occ>) -> (Vec<usize>, Vec<usize>) {
    let mut sampled_endpoints = Vec::<usize>::new();
    let mut lengths = Vec::<usize>::new();
    let n = seq.len();
    for i in 0..n-k+1{
        let kmer = &seq[i..i+k];
        //println!("Processing {}", String::from_utf8(kmer.to_vec()).unwrap());
        let mut finimizer: Option<&[u8]> = None;
        let mut f_start = 0;
        let mut f_end = kmer.len();
        for start in 0..kmer.len(){
            for end in start+1..=kmer.len(){
                let x = &kmer[start..end];
                let freq = count(index, x);
                //eprintln!("{} {} {}", start, end, freq);
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

        //println!("{}, {}, {}", f_end, String::from_utf8(kmer.to_vec()).unwrap(), String::from_utf8(finimizer.unwrap().to_vec()).unwrap());
    }

    assert_eq!(sampled_endpoints.len(), lengths.len());

    (sampled_endpoints, lengths)

}

fn main() {

    // Read file path from argv
    let filepath = std::env::args().nth(1).unwrap();
    let k = std::env::args().nth(2).unwrap().parse::<usize>().unwrap();

    let reader = jseqio::reader::DynamicFastXReader::from_file(&filepath).unwrap();
    let db = reader.into_db().unwrap();
    let data = get_dollar_concatenation(&db);

    /*
    let n = 1000_000_usize;
    let seed = 1234;
    let mut data = generate_random_dna_string(n, seed);

    data.push(b'$'); // The BWT library needs this?
    */

    // Build the FM index
    eprintln!("Building FM index");
    let alphabet = dna::n_alphabet();
    let sa = suffix_array(data.as_slice());
    let bwt = bwt(data.as_slice(), &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let index: FMIndex<&Vec<u8>, &Vec<usize>, &Occ> = FMIndex::new(&bwt, &less, &occ);
    eprintln!("FM index built");

    let mut total_finimizer_count = 0_usize; // Number of endpoints that are at the end of a finimizer
    let mut total_seq_len = 0_usize;
    let mut total_finimizer_len = 0_usize;
    for i in 0..db.sequence_count(){
        let rec = db.get(i);
        eprintln!("Processing sequence {} of length {} (total processed: {}, density : {})", i, rec.seq.len(), total_seq_len, total_finimizer_count as f64 / total_seq_len as f64);
        let (ends, lengths) = get_finimizers(rec.seq, k, &index);
        total_finimizer_count += ends.len();
        total_seq_len += rec.seq.len();
        total_finimizer_len += lengths.iter().sum::<usize>();
    }

    println!("{}/{} = {}", total_finimizer_count, total_seq_len, total_finimizer_count as f64 /  total_seq_len as f64);
    println!("Mean length: {}", total_finimizer_len as f64 / total_finimizer_count as f64);


}
