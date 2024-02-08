use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::fmindex::BackwardSearchResult;
use clap::{Command, Arg};
use jseqio::reader::SeqRecordProducer;

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

fn main() {
    let cli = Command::new("total-string-frequency")
        .arg_required_else_help(true)
        .arg(Arg::new("unitigs")
            .help("Fasta format")
            .short('u')
            .value_parser(clap::value_parser!(std::path::PathBuf)))
        .arg(Arg::new("patterns")
            .help("Fasta format")
            .short('p')
            .value_parser(clap::value_parser!(std::path::PathBuf)));

    let matches = cli.get_matches();
    let unitigs_file: &std::path::PathBuf = matches.get_one("unitigs").unwrap();
    let patterns_file: &std::path::PathBuf = matches.get_one("patterns").unwrap();

    // Read unitigs and concatenate them with dollar separators
    let mut unitigs_concat = Vec::<u8>::new();
    let mut unitig_reader = jseqio::reader::DynamicFastXReader::from_file(unitigs_file).unwrap();
    while let Some(rec) = unitig_reader.read_next().unwrap(){
        unitigs_concat.extend(rec.seq);
        unitigs_concat.push(b'$');
    }

    // Build the FM index
    let alphabet = dna::n_alphabet();
    let sa = suffix_array(unitigs_concat.as_slice());
    let bwt = bwt(unitigs_concat.as_slice(), &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let index: FMIndex<&Vec<u8>, &Vec<usize>, &Occ> = FMIndex::new(&bwt, &less, &occ);

    // Count occurrences of distinct patterns
    let mut searched_patterns = std::collections::HashSet::<Vec<u8>>::new();
    let mut pattern_reader = jseqio::reader::DynamicFastXReader::from_file(patterns_file).unwrap();
    let mut total_count = 0_usize;
    while let(Some(rec)) = pattern_reader.read_next().unwrap(){
        let pattern = rec.seq;
        if searched_patterns.contains(pattern){
            continue;
        } else{
            searched_patterns.insert(pattern.to_owned());
            total_count += count(&index, pattern);
        }
    }
    println!("Distinct patterns: {}", searched_patterns.len());
    println!("Total occurrences: {}", total_count);

}
