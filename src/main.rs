mod brute_force_finimizers;

use clap::{Command, Arg};
use jseqio::reader::SeqRecordProducer;
use sbwt::sbwt::*;
use sbwt::subsetrank::*;
use bitvec::prelude::*;

// Returns a triple:
// 1) A Vec of distinct ranges i..j such that seq[i..j) is a finimizer
// 2) A Vec of lengths of finimizers that were seen for the first time in this sequence
// 2) A Vec of frequencies of finimizers that were seen for the first time in this sequence
#[allow(clippy::needless_range_loop, non_snake_case)]
fn get_streaming_finimizers(SS: &StreamingSupport<MatrixRank>, seq: &[u8], k: usize, t: usize, lex_marks: &mut BitVec, finimizer_strings_out: &mut Option<impl std::io::Write>) -> (Vec<std::ops::Range<usize>>, Vec<usize>, Vec<usize>) {
    assert!(seq.len() >= k);
    let mut finimizer_ranges = Vec::<std::ops::Range::<usize>>::new();
    let mut new_frequencies = Vec::<usize>::new();
    let mut new_lengths= Vec::<usize>::new();
    let SFS = SS.shortest_freq_bound_suffixes(seq, t);
    
    for start in 0..seq.len()-k+1 {
        // Figure out the finimizer

        let mut best = (usize::MAX, usize::MAX, usize::MAX, -1_isize); // Length, freq, colex, endpoint
        for end in start..start+k { // Inclusive end!
            if SFS[end].is_none() { continue } // No unique match ending here
            let (len, I) = SFS[end].as_ref().unwrap(); // Length, interval
            if end + 1 < start + len { continue } // Shortest unique match not fit in this k-mer window (end - len + 1 < start)
            if (*len, I.len(), I.start, end as isize) < best {
                best = (*len, I.len(), I.start, end as isize)
            }
        }
        assert!(best.3 >= 0); // Endpoint must be set by this point
        best.3 += 1; // Make the end exclusive

        // Write the finimizer string, if needed
        if let Some(writer) = finimizer_strings_out.as_mut() {
            if lex_marks.get(best.2).unwrap() == false { // First time seeing this
                let (len, end) = (best.0 as isize, best.3); // length, exclusive end
                writer.write_all(b">\n").unwrap();
                writer.write_all(&seq[(end-len) as usize .. end as usize]).unwrap();
                writer.write_all(b"\n").unwrap();
            }
        }
        
        // Report the finimizer
        if lex_marks.get(best.2).unwrap() == false { // First time seeing this
            new_frequencies.push(best.1); // Frequency
            new_lengths.push(best.0);
            lex_marks.set(best.2, true);
        }

        let last = finimizer_ranges.last();
        if last.is_none() || last.is_some_and(|I| I.end != best.3 as usize) {
            let (len, I) = SFS[best.3 as usize - 1].as_ref().unwrap(); // -1: back to inclusive end for indexing SFS
            finimizer_ranges.push((best.3 - *len as isize) as usize .. best.3 as usize);
        }
    }

    assert_eq!(new_lengths.len(), new_frequencies.len());
    (finimizer_ranges, new_lengths, new_frequencies)

}

#[allow(non_snake_case)]
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
        .arg(Arg::new("finimizer-out")
            .short('o')
            .long("finimizer-out")
            .help("Finimizer output fasta/fastq file")
            .value_parser(clap::value_parser!(std::path::PathBuf)))
        .arg(Arg::new("k")
            .short('k')
            .help("k-mer k")
            .value_parser(clap::value_parser!(usize))
            .required(true))
        .arg(Arg::new("t")
            .short('t')
            .help("Frequency bound t")
            .value_parser(clap::value_parser!(usize))
            .default_value("1"))
        .arg(Arg::new("threads")
            .help("Number of parallel threads to use")
            .short('p')
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
    let finimizer_outfile = matches.get_one::<std::path::PathBuf>("finimizer-out");
    let k = *matches.get_one::<usize>("k").unwrap();
    let t = *matches.get_one::<usize>("t").unwrap();
    let nthreads = *matches.get_one::<usize>("threads").unwrap();
    let mem_gb= *matches.get_one::<usize>("memory").unwrap();

    let reader = jseqio::reader::DynamicFastXReader::from_file(&filepath).unwrap();
    let mut finimizer_out = finimizer_outfile.map(|f| std::io::BufWriter::new(std::fs::File::create(f).unwrap()));

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
    let mut total_finimizer_len = 0_usize;

    let mut total_distinct_finimizer_count = 0_usize;
    let mut total_distinct_finimizer_len = 0_usize;
    let mut total_distinct_finimizer_freq= 0_usize;

    let mut reader2 = jseqio::reader::DynamicFastXReader::from_file(&filepath).unwrap();
    let mut lex_marks = bitvec![0; sbwt.n_sets()];
    let mut total_seq_len = 0_usize;
    let mut total_kmers = 0_usize;
    // Print current time
    let start_time = std::time::Instant::now();
    println!("Starting finimizer search");
    while let Some(rec) = reader2.read_next().unwrap(){
        total_seq_len += rec.seq.len();
        total_kmers += std::cmp::max(rec.seq.len() as isize - k as isize + 1, 0) as usize;

        let (text_ranges, new_lengths, new_freqs) = get_streaming_finimizers(&SS, rec.seq, k, t, &mut lex_marks, &mut finimizer_out);

        total_finimizer_count += text_ranges.len();
        total_finimizer_len += text_ranges.iter().fold(0_usize,|acc, R| acc + R.len());

        total_distinct_finimizer_count += new_lengths.len();
        total_distinct_finimizer_len += new_lengths.iter().sum::<usize>();
        total_distinct_finimizer_freq += new_freqs.iter().sum::<usize>();

    }
    let now = std::time::Instant::now();
    println!("{} k-mers queried at {} us/k-mer", total_kmers, now.duration_since(start_time).as_micros() as f64 / total_kmers as f64);
    println!("{} lex marks (DBG density {})", lex_marks.count_ones(), lex_marks.count_ones() as f64 / sbwt.n_sets() as f64);
    println!("{}/{} = {} streaming density", total_finimizer_count, total_seq_len, total_finimizer_count as f64 / total_seq_len as f64);
    println!("{} mean streaming length", total_finimizer_len as f64 / total_finimizer_count as f64);
    println!("{} mean distinct frequency", total_distinct_finimizer_freq as f64 / total_distinct_finimizer_count as f64);
    println!("{} mean distinct length", total_distinct_finimizer_len as f64 / total_distinct_finimizer_count as f64);


}
