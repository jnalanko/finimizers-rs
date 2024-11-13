use std::io::Write;
use sbwt::sbwt::*;
use sbwt::subsetrank::*;

use jseqio::reader::*;

fn build_command(matches: &clap::ArgMatches){

    let infile = matches.get_one::<std::path::PathBuf>("input").unwrap();
    let outfile = matches.get_one::<std::path::PathBuf>("output").unwrap();
    let k = *matches.get_one::<usize>("k").unwrap();
    let mem_gb = *matches.get_one::<usize>("mem-gb").unwrap();
    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let dedup_batches = matches.get_flag("dedup-batches");

    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();

    let reader = DynamicFastXReader::from_file(infile).unwrap();
    let mut out = std::io::BufWriter::new(std::fs::File::create(outfile).unwrap());
 
    log::info!("Building SBWT");
    let (sbwt, _) = match k {
        0..=32 => {
            Sbwt::<MatrixRank>::new::<1>(reader, k, mem_gb, n_threads, dedup_batches, false)
        }
        33..=64 => {
            Sbwt::<MatrixRank>::new::<2>(reader, k, mem_gb, n_threads, dedup_batches, false)
        }
        _ => {
            panic!("k > 64 not supported");
        }
    };
    log::info!("Serializing");
    let n_written = sbwt.serialize(&mut out).unwrap();
    log::info!("Wrote {} bytes ({:.2} bits / k-mer)", n_written, n_written as f64 * 8.0 / sbwt.n_kmers() as f64);

}

fn query_command(matches: &clap::ArgMatches){
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let outfile = matches.get_one::<std::path::PathBuf>("output").unwrap();
    let queryfile = matches.get_one::<std::path::PathBuf>("query").unwrap();

    let mut query_reader = DynamicFastXReader::from_file(queryfile).unwrap();
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let mut out = std::io::BufWriter::new(std::fs::File::create(outfile).unwrap());

    let sbwt = Sbwt::<MatrixRank>::load(&mut index_reader).unwrap();

    let start_time = std::time::Instant::now();
    let mut n_query_kmers = 0;
    let mut n_found = 0;
    while let Some(rec) = query_reader.read_next().unwrap(){
        let seq = rec.seq;
        for kmer in seq.windows(sbwt.k()){
            n_query_kmers += 1;
            if sbwt.search(kmer).is_some(){
                out.write_all(b"1").unwrap();
                n_found += 1;
            } else{
                out.write_all(b"0").unwrap();
            }
        }
    }
    let end_time = std::time::Instant::now();
    let elapsed = end_time - start_time;
    log::info!("Queried {} k-mers", n_query_kmers);
    log::info!("{:.2}% of queried k-mers found", n_found as f64 / n_query_kmers as f64 * 100.0);
    log::info!("Elapsed time: {:.2} seconds ({:.2} ns / k-mer)", elapsed.as_secs_f64(), elapsed.as_nanos() as f64 / n_query_kmers as f64);
    
}

fn main() {

    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = clap::Command::new("sbwt")
        .about("K-mer using using the SBWT.")
        .arg_required_else_help(true)
        .subcommand(clap::Command::new("build")
            .arg(clap::Arg::new("input")
                .help("Input fasta or fastq seuqence file")
                .short('i')
                .long("input")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output")
                .help("Output SBWT file")
                .short('o')
                .long("output")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("k")
                .help("k-mer length")
                .short('k')
                .required(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(clap::Arg::new("mem-gb")
                .help("An approximate memory budget in GB. The actual memory usage may be higher.")
                .short('m')
                .long("mem-gb")
                .required(false)
                .value_parser(clap::value_parser!(usize))
                .default_value("8")
            )
            .arg(clap::Arg::new("threads")
                .help("Number of threads to use")
                .long("threads")
                .short('t')
                .default_value("8")
                .global(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(clap::Arg::new("dedup-batches")
                .help("Slows down the construction, but saves disk if the data is very repetitive")
                .long("dedup-batches")
                .short('d')
                .action(clap::ArgAction::SetTrue)
            )
        )
        .subcommand(clap::Command::new("query")
            .arg(clap::Arg::new("index")
                .help("SBWT index file")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("query")
                .help("Query sequences in FASTA or FASTQ format, possibly gzipped")
                .short('q')
                .long("query")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output")
                .help("Output text file")
                .short('o')
                .long("output")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
        )
        ;

    let matches = cli.get_matches();
    match matches.subcommand(){
        Some(("build", sub_matches)) => build_command(sub_matches),
        Some(("query", sub_matches)) => query_command(sub_matches),
        _ => unreachable!(),
    }
    
}
