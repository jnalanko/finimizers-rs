use clap::{Command, Arg};
use jseqio::reader::SeqRecordProducer;
use rand::{self, RngCore, SeedableRng};
use rand::rngs::StdRng;
use sbwt::sbwt::*;
use sbwt::subsetrank::*;
use bitvec::prelude::*;

fn sbwt_count(index: &Sbwt::<MatrixRank>, pattern: &[u8]) -> usize {
    match index.search(pattern){
        Some(interval) => interval.len(),
        None => 0,
    }
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

fn colex_smaller(a: &[u8], b: &[u8]) -> bool {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    for i in 0..n {
        if a[n-1-i] < b[n-1-i] { return true; }
        if a[n-1-i] > b[n-1-i] { return false; }
    }

    false
}

// Returns pair (finimizer endpoint vector, finimizer length vector
// Marks the lex positions of the minimizers in the given bit vector lex_marks
fn get_finimizers(seq: &[u8], k: usize, index: &Sbwt::<MatrixRank>, lex_marks: &mut BitVec) -> (Vec<usize>, Vec<usize>) {
    let mut sampled_endpoints = Vec::<usize>::new();
    let mut lengths = Vec::<usize>::new();

    // Bit set with k bits
    let n = seq.len();
    for i in 0..n-k+1{
        let kmer = &seq[i..i+k];
        //println!("Processing {}", String::from_utf8(kmer.to_vec()).unwrap());
        let mut finimizer: Option<&[u8]> = None;
        let mut f_start = 0;
        let mut f_end = kmer.len();
        let mut lex: Option<usize> = None;
        for start in 0..kmer.len(){
            for end in start+1..=kmer.len(){
                let x = &kmer[start..end];
                let range = index.search(x);
                let freq = match &range {
                    Some(interval) => interval.len(),
                    None => 0,
                };
                //eprintln!("{} {} {}", start, end, freq);
                if freq == 1 {
                    if finimizer.is_none() || finimizer.is_some_and(
                        |cur| x.len() < cur.len() || (x.len() == cur.len() && colex_smaller(x , cur))
                    ){
                        finimizer = Some(x);
                        lex = Some(range.unwrap().start);
                        f_start = i + start;
                        f_end = i + end;
                    }
                    break; // No reason to check further endpoints because this one was unique already
                }
            }
        }

        lex_marks.set(lex.unwrap(), true); // Should always exist because a full k-mer has freq 1
        let last = sampled_endpoints.last();
        if last.is_none() || last.is_some_and(|e| *e != f_end) {
            sampled_endpoints.push(f_end);
            lengths.push(f_end - f_start);
        }
    }

    assert_eq!(sampled_endpoints.len(), lengths.len());

    (sampled_endpoints, lengths)

}

#[allow(clippy::needless_range_loop, non_snake_case)]
fn get_streaming_finimizers(SS: &StreamingSupport<MatrixRank>, seq: &[u8], k : usize, lex_marks: &mut BitVec) -> (Vec<usize>, Vec<usize>) {
    assert!(seq.len() >= k);
    let mut sampled_endpoints = Vec::<usize>::new();
    let mut lengths = Vec::<usize>::new();
    let SFS = SS.shortest_freq_bound_suffixes(seq, 1);
    
    for start in 0..seq.len()-k+1 {
        // Figure out the finimizer

        let mut best = (usize::MAX, usize::MAX, -1_isize); // Length, colex, endpoint
        for end in start..start+k { // Inclusive end!
            if SFS[end].is_none() { continue } // No unique match ending here
            let (len, I) = SFS[end].as_ref().unwrap(); // Length, interval
            if end + 1 < start + len { continue } // Shortest unique match not fit in this k-mer window (end - len + 1 < start)
            if (*len, I.start, end as isize) < best {
                best = (*len, I.start, end as isize)
            }
        }
        assert!(best.2 >= 0); // Endpoint must be set by this point
        best.2 += 1; // Make the end exclusive

        // Report the finimizer

        lex_marks.set(best.1, true);

        let last = sampled_endpoints.last();
        if last.is_none() || last.is_some_and(|e| *e != best.2 as usize) {
            sampled_endpoints.push(best.2 as usize);
            lengths.push(SFS[best.2 as usize - 1].as_ref().unwrap().0); // -1: back to inclusive end for indexing SFS
        }
    }

    assert_eq!(sampled_endpoints.len(), lengths.len());

    (sampled_endpoints, lengths)

}

fn main() {

    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("finimizer-experiment")
        .arg_required_else_help(true)
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .help("Input fasta/fastq file")
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true))
        .arg(Arg::new("k")
            .short('k')
            .help("k-mer k")
            .value_parser(clap::value_parser!(usize))
            .required(true))
        .arg(Arg::new("threads")
            .help("Number of threads to use")
            .short('t')
            .value_parser(clap::value_parser!(usize))
            .default_value("4")
            .long("threads"))
        .arg(Arg::new("memory")
            .help("Memory budget in GB (not strictly enforced)")
            .short('m')
            .value_parser(clap::value_parser!(usize))
            .default_value("8")
            .long("mem-gb"));

    let matches = cli.get_matches();

    let filepath = matches.get_one::<std::path::PathBuf>("input").unwrap();
    let k = *matches.get_one::<usize>("k").unwrap();
    let nthreads = *matches.get_one::<usize>("threads").unwrap();
    let mem_gb= *matches.get_one::<usize>("memory").unwrap();

    let reader = jseqio::reader::DynamicFastXReader::from_file(&filepath).unwrap();

    // Choose the number of u64s in a k-mer based on the k
    let (sbwt, lcs) = match k {
        0..=32 => {
            Sbwt::<MatrixRank>::new::<1>(reader, k, mem_gb, nthreads, true, true)
        }
        33..=64 => {
            Sbwt::<MatrixRank>::new::<2>(reader, k, mem_gb, nthreads, true, true)
        }
        65..=96 => {
            Sbwt::<MatrixRank>::new::<3>(reader, k, mem_gb, nthreads, true, true)
        }
        97..=128 => {
            Sbwt::<MatrixRank>::new::<4>(reader, k, mem_gb, nthreads, true, true)
        }
        _ => {
            panic!("k > 128 not supported");
        }
    };
    let SS = StreamingSupport::new(&sbwt, lcs.unwrap());
   
    let mut total_finimizer_count = 0_usize; // Number of endpoints that are at the end of a finimizer
    let mut total_seq_len = 0_usize;
    let mut total_finimizer_len = 0_usize;
    let mut reader2 = jseqio::reader::DynamicFastXReader::from_file(&filepath).unwrap();
    let mut seq_id = 0_usize;
    let mut lex_marks = bitvec![0; sbwt.n_sets()];
    let mut total_kmers = 0_usize;
    // Print current time
    let start_time = std::time::Instant::now();
    println!("Starting finimizer search");
    while let Some(rec) = reader2.read_next().unwrap(){
        //println!("Processing sequence {} of length {} (total processed: {}, density : {})", seq_id, rec.seq.len(), total_seq_len, total_finimizer_count as f64 / total_seq_len as f64);
        //let (ends, lengths) = get_finimizers(rec.seq, k, &sbwt, &mut lex_marks);
        let (ends2, lengths2) = get_streaming_finimizers(&SS, rec.seq, k,  &mut lex_marks);
        total_finimizer_count += ends2.len();
        total_seq_len += rec.seq.len();
        total_finimizer_len += lengths2.iter().sum::<usize>();
        total_kmers += std::cmp::max(rec.seq.len() as isize - k as isize + 1, 0) as usize;
        seq_id += 1;
    }
    let now = std::time::Instant::now();
    println!("{} k-mers queried at {} us/k-mer", total_kmers, now.duration_since(start_time).as_micros() as f64 / total_kmers as f64);
    println!("{} lex marks (DBG density {})", lex_marks.count_ones(), lex_marks.count_ones() as f64 / sbwt.n_sets() as f64);
    println!("{}/{} = {}", total_finimizer_count, total_seq_len, total_finimizer_count as f64 /  total_seq_len as f64);
    println!("Mean length: {}", total_finimizer_len as f64 / total_finimizer_count as f64);

}
