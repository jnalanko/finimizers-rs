use crate::construction::cursors;
use crate::subsetrank::*;
use crate::util;
use simple_sds::ops::Access;

use crate::construction::{dummies, kmer_splitter, cursors::DummyNodeMerger};

#[derive(Eq, PartialEq, Debug)]
#[allow(non_snake_case)] // C-array is an established convention in BWT indexes
pub struct Sbwt<SR: SubsetRank> {
    subset_rank: SR,
    alphabet: Vec<u8>,
    alphabet_to_index: Vec<u8>,
    n_kmers: usize,
    k: usize,
    C: Vec<usize>, // Cumulative character counts (includes one ghost dollar)
}

pub struct StreamingSupport<'a, SR: SubsetRank>{
    sbwt: &'a Sbwt<SR>,
    lcs: simple_sds::int_vector::IntVector,
}

impl <'a, SR: SubsetRank> StreamingSupport<'a, SR>{

    pub fn new(sbwt: &'a Sbwt<SR>, lcs: simple_sds::int_vector::IntVector) -> Self{
        Self{sbwt, lcs}
    }

    #[allow(non_snake_case)]
    pub fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        let mut new_start = I.start;
        let mut new_end = I.end;
        while new_start > 0 && self.lcs.get(new_start) >= target_len as u64 {
            new_start -= 1;
        }
        while new_end < self.sbwt.subset_rank.len() && self.lcs.get(new_end) >= target_len as u64 {
            new_end += 1;
        }
        new_start..new_end
    }

    #[allow(non_snake_case)]
    pub fn extend_right(&self, I: std::ops::Range<usize>, c: u8) -> std::ops::Range<usize> {
        let c_idx = self.sbwt.char_idx(c);
        self.sbwt.lf_step(I.start, c_idx)..self.sbwt.lf_step(I.end, c_idx)
    }

    // Returns pairs (match length, colex interval) (see paper "Finimizers: Variable-length bounded-frequency minimizers for k-mer sets")
    #[allow(non_snake_case)]
    pub fn shortest_freq_bound_suffixes(&self, query: &[u8], freq_bound: usize) -> Vec<Option<(usize, std::ops::Range<usize>)>> {
        let mut SFS = vec![None; query.len()];
        let mut I = 0..self.sbwt.subset_rank.len();
        let mut d = 0_usize; // String length of interval I
        for i in 0..query.len() {
            let mut Ic = self.extend_right(I.clone(), query[i]);
            while d > 0 && Ic.is_empty() {
                I = self.contract_left(I, d-1);
                d -= 1;
                Ic = self.extend_right(I.clone(), query[i]);
            }

            if !Ic.is_empty() {
                I = Ic;
                d = std::cmp::min(self.sbwt.k, d+1);
            }

            // Contract I to the sortest frequency-bounded suffix
            while I.len() <= freq_bound && d > 0 {
                let J = self.contract_left(I.clone(), d-1);
                if J.len() <= freq_bound {
                    I = J;
                    d -= 1;
                } else {
                    break
                }
            }

            if I.len() <= freq_bound{
                SFS[i] = Some((d, I.clone()));
            }
        }
        SFS
    }

    // Returns pairs (match length, colex interval)
    #[allow(non_snake_case)]
    pub fn matching_statistics(&self, query: &[u8]) -> Vec<(usize, std::ops::Range<usize>)> {
        let mut d = 0_usize;
        let mut MS = vec![(0_usize, 0..0); query.len()];
        let mut I = 0..self.sbwt.subset_rank.len();
        for i in 0..query.len() {
            while d > 0 && self.extend_right(I.clone(), query[i]).is_empty() {
                I = self.contract_left(I, d-1);
                d -= 1;
            }
            let Ic = self.extend_right(I.clone(), query[i]);
            if !Ic.is_empty() {
                I = Ic;
                d = std::cmp::min(self.sbwt.k, d+1);
            }
            MS[i] = (d, I.clone());
        }
        MS
    }
}

impl<SR: SubsetRank> Sbwt<SR> {
    pub fn n_kmers(&self) -> usize {
        self.n_kmers
    }

    pub fn n_sets(&self) -> usize {
        self.subset_rank.len()
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn char_idx(&self, c: u8) -> usize {
        self.alphabet_to_index[c as usize] as usize
    }

    #[allow(unused)]
    pub fn alphabet(&self) -> &[u8] {
        self.alphabet.as_slice()
    }

    #[allow(non_snake_case)]
    #[allow(unused)]
    pub fn C(&self) -> &[usize] {
        self.C.as_slice()
    }

    // Returns the number of bytes written
    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;

        n_written += self.subset_rank.serialize(out)?;

        // We're not using serde because we want full control over the bytes
        // in order to guarantee compatibility across languages
        n_written += util::write_bytes(out, &self.alphabet.len().to_le_bytes())?;
        n_written += util::write_bytes(out, &self.alphabet)?;

        n_written += util::write_bytes(out, &self.alphabet_to_index.len().to_le_bytes())?;
        n_written += util::write_bytes(out, &self.alphabet_to_index)?;

        n_written += util::write_bytes(out, &self.n_kmers.to_le_bytes())?;
        n_written += util::write_bytes(out, &self.k.to_le_bytes())?;
        n_written += util::write_bytes(out, &self.C.len().to_le_bytes())?; // TODO: check at build time that usize = u64, crash otherwise
        n_written += util::write_bytes(
            out,
            &self
                .C
                .iter()
                .flat_map(|x| x.to_le_bytes())
                .collect::<Vec<u8>>(),
        )?;

        Ok(n_written)
    }

    #[allow(non_snake_case)] // For C-array
    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let subset_rank = SR::load(input)?;

        let mut alphabet_len = [0_u8; 8];
        input.read_exact(&mut alphabet_len)?;
        let alphabet_len = u64::from_le_bytes(alphabet_len) as usize;

        let mut alphabet = vec![0_u8; alphabet_len];
        input.read_exact(&mut alphabet)?;

        let mut alphabet_to_index_len = [0_u8; 8];
        input.read_exact(&mut alphabet_to_index_len)?;
        let alphabet_to_index_len = u64::from_le_bytes(alphabet_to_index_len) as usize;

        let mut alphabet_to_index = vec![0_u8; alphabet_to_index_len];
        input.read_exact(&mut alphabet_to_index)?;

        let mut n_kmers = [0_u8; 8];
        input.read_exact(&mut n_kmers)?;
        let n_kmers = u64::from_le_bytes(n_kmers) as usize;

        let mut k = [0_u8; 8];
        input.read_exact(&mut k)?;
        let k = u64::from_le_bytes(k) as usize;

        let mut C_len = [0_u8; 8];
        input.read_exact(&mut C_len)?;
        let C_len = u64::from_le_bytes(C_len) as usize;

        let mut C_bytes = vec![0_u8; C_len * 8];
        input.read_exact(&mut C_bytes)?;

        let C = C_bytes
            .chunks_exact(8)
            .map(|chunk| {
                let mut arr = [0_u8; 8];
                arr.copy_from_slice(chunk);
                u64::from_le_bytes(arr) as usize
            })
            .collect::<Vec<usize>>();

        Ok(Self {
            subset_rank,
            alphabet,
            alphabet_to_index,
            n_kmers,
            k,
            C,
        })
    }


    pub fn lf_step(&self, i: usize, char_idx: usize) -> usize {
        self.C[char_idx] + self.subset_rank.rank(char_idx as u8, i)
    }

    // Returns the interval of the pattern, if found
    // The pattern is given in ascii characters
    pub fn search(&self, pattern: &[u8]) -> Option<std::ops::Range<usize>> {
        let mut left = 0_usize;
        let mut right = self.subset_rank.len();
        for chr in pattern.iter() {
            let c = self.alphabet_to_index[*chr as usize];
            left = self.lf_step(left, c as usize);
            right = self.lf_step(right, c as usize);
            if left >= right {
                return None;
            }
        }
        Some(left..right)
    }

    //pub fn streaming_search(pattern: &[u8], window_len: usize) {}

    // B is the number u64 words in a k-mer
    // Returns pair (SBWT, Optional LCS)
    pub fn new<const B: usize>(seqs: jseqio::reader::DynamicFastXReader, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, build_lcs: bool) -> (Self, Option<simple_sds::int_vector::IntVector>) {
        let sigma = 4_usize; // DNA alphabet size: todo: this is hard coded all over the place

        let tempdir = std::path::PathBuf::from("./temp"); // TODO: read from command line?
        std::fs::create_dir_all(&tempdir).unwrap();

        log::info!("Splitting k-mers into bins");
        let bin_files = kmer_splitter::split_to_bins::<B>(seqs, k, mem_gb, n_threads, dedup_batches);

        log::info!("Sorting and deduplicating bins");
        kmer_splitter::par_sort_and_dedup_bin_files::<B>(&bin_files, mem_gb, n_threads);

        let (kmers_filename, kmers_writer) = util::create_temp_file(&tempdir);
        kmer_splitter::concat_files(bin_files, kmers_writer);

        let n_kmers = std::fs::metadata(kmers_filename.path()).unwrap().len() as usize / crate::kmer::LongKmer::<B>::byte_size();

        log::info!("{} distinct k-mers found", n_kmers);

        let required_dummies = dummies::get_sorted_dummies::<B>(kmers_filename.path(), sigma, k);

        log::info!("{} dummy nodes needed", required_dummies.len());

        let n = n_kmers + required_dummies.len();

        // Write dummies to disk
        let (dummy_file, dummy_writer) = util::create_temp_file(&tempdir);
        dummies::write_to_disk(required_dummies, dummy_writer);
        
        log::info!("Constructing the sbwt subset sequence");

        let char_cursors = cursors::init_char_cursors::<B>(dummy_file.path(),kmers_filename.path(), k, sigma);

        let global_cursor = DummyNodeMerger::new(
            std::io::BufReader::new(std::fs::File::open(dummy_file.path()).unwrap()),
            std::io::BufReader::new(std::fs::File::open(kmers_filename.path()).unwrap()),
            k,
        );

        let (rawrows, lcs) = cursors::build_sbwt_bit_vectors(global_cursor, char_cursors, n, k, sigma, build_lcs);

        // Create the C array
        #[allow(non_snake_case)] // C-array is an established convention in BWT indexes
        let C: Vec<usize> = crate::construction::get_C_array(&rawrows);

        let alphabet_to_index = crate::construction::get_dna_alphabet_to_index();

        log::info!("Building the subset rank structure");
        (Self {
            subset_rank: SR::new_from_bit_vectors(rawrows.into_iter().map(simple_sds::bit_vector::BitVector::from).collect()),
            alphabet: vec![b'A', b'C', b'G', b'T'],
            alphabet_to_index,
            n_kmers,
            k,
            C,
        }, lcs)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use simple_sds::ops::Access;

    use super::*;

    #[allow(non_snake_case)]
    fn ACGT_to_0123(c: u8) -> u8 {
        match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!("Invalid DNA character {} ", c as char),
        }
    }

    fn encode(s: &str) -> Vec<u8> {
        s.as_bytes()
            .iter()
            .map(|c| ACGT_to_0123(*c))
            .collect::<Vec<u8>>()
    }

    fn to_fasta(seqs: Vec<Vec<u8>>) -> Vec<u8>{
        let mut ans = Vec::<u8>::new();
        for s in seqs.iter(){
            ans.push(b'>');
            ans.push(b'\n');
            ans.extend_from_slice(s);
            ans.push(b'\n');
        }
        ans
    }

    #[test]
    #[allow(non_snake_case)]
    fn LCS_paper_example() {
        if std::env::var("RUST_LOG").is_err(){
            std::env::set_var("RUST_LOG", "info");
        }
    
        env_logger::init();

        let seqs = vec![b"AGGTAAA".to_vec(), b"ACAGGTAGGAAAGGAAAGT".to_vec()];

        let fasta = to_fasta(seqs);
        let reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(fasta)).unwrap();

        let (sbwt, lcs) = Sbwt::<MatrixRank>::new::<2>(reader, 4, 3, 2, false, true);
        let lcs = lcs.unwrap();

        assert_eq!(sbwt.subset_rank.len(), 18);

        assert_eq!(sbwt.subset_rank.access(0), encode("A")); //   $$$$
        assert_eq!(sbwt.subset_rank.access(1), encode("C")); //   $$$A
        assert_eq!(sbwt.subset_rank.access(2), encode("G")); //   GAAA
        assert_eq!(sbwt.subset_rank.access(3), encode("")); //    TAAA
        assert_eq!(sbwt.subset_rank.access(4), encode("A")); //   GGAA
        assert_eq!(sbwt.subset_rank.access(5), encode("A")); //   GTAA
        assert_eq!(sbwt.subset_rank.access(6), encode("G")); //   $ACA
        assert_eq!(sbwt.subset_rank.access(7), encode("A")); //   AGGA
        assert_eq!(sbwt.subset_rank.access(8), encode("AG")); //  GGTA
        assert_eq!(sbwt.subset_rank.access(9), encode("A")); //   $$AC
        assert_eq!(sbwt.subset_rank.access(10), encode("GT")); // AAAG
        assert_eq!(sbwt.subset_rank.access(11), encode("G")); //  ACAG
        assert_eq!(sbwt.subset_rank.access(12), encode("G")); //  GTAG
        assert_eq!(sbwt.subset_rank.access(13), encode("AT")); // AAGG
        assert_eq!(sbwt.subset_rank.access(14), encode("")); //   CAGG
        assert_eq!(sbwt.subset_rank.access(15), encode("")); //   TAGG
        assert_eq!(sbwt.subset_rank.access(16), encode("")); //   AAGT
        assert_eq!(sbwt.subset_rank.access(17), encode("A")); //  AGGT

        assert_eq!(sbwt.search(b""), Some(0..18));
        assert_eq!(sbwt.search(b"AGG"), Some(13..16));
        assert_eq!(sbwt.search(b"AAGT"), Some(16..17));

        let true_lcs = [0,0,1,3,2,2,1,1,1,0,0,2,2,1,3,3,0,2];
        for (i, x) in lcs.iter().enumerate(){
            println!("LCS {}", x);
            assert_eq!(true_lcs[i], x);
        }

        // Test streaming support
        let SS = StreamingSupport{sbwt: &sbwt, lcs};
        let mut I = 0..sbwt.subset_rank.len();
        I = SS.extend_right(I, b'A') ;
        I = SS.extend_right(I, b'G') ;
        I = SS.extend_right(I, b'G') ;
        assert_eq!(I, 13..16);

        let I3 = SS.contract_left(I.clone(), 3);
        assert_eq!(I3, 13..16);

        let I2 = SS.contract_left(I.clone(), 2);
        assert_eq!(I2, 13..16);

        let I1 = SS.contract_left(I.clone(), 1);
        assert_eq!(I1, 10..16);

        let I0 = SS.contract_left(I.clone(), 0);
        assert_eq!(I0, 0..18);
        dbg!(I3, I2, I1, I0);
    }

    #[test]
    #[allow(non_snake_case)]
    fn finimizer_paper_example(){

        // Using the example from the finimizer paper
        let seqs = vec![b"GTAAGTCT".to_vec(), b"AGGAAA".to_vec(), b"ACAGG".to_vec(), b"GTAGG".to_vec(), b"AGGTA".to_vec()];
        let k = 4;
        let fasta = to_fasta(seqs);
        let reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(fasta)).unwrap();

        let (sbwt , lcs) = Sbwt::<MatrixRank>::new::<2>(reader, k, 3, 2, true, true);
        let SS = StreamingSupport{sbwt: &sbwt, lcs: lcs.unwrap()};
        let MS = SS.matching_statistics(b"AAGTAA");
        eprintln!("{:?}", MS);
        let true_MS_lengths: [usize; 6] = [1,2,3,4,3,4];
        let true_MS_freqs: [usize; 6] = [7,3,1,1,1,1];
        let true_MS_colex_starts: [usize; 6] = [2,3,11,17,8,5];
        assert_eq!(MS.len(), true_MS_lengths.len());
        for i in 0..MS.len(){
            assert_eq!(MS[i].0, true_MS_lengths[i]);
            assert_eq!(MS[i].1.len(), true_MS_freqs[i]);
            assert_eq!(MS[i].1.start + 1, true_MS_colex_starts[i]); // Paper has 1-indexing
        }

        let mut SFS = SS.shortest_freq_bound_suffixes(b"AAGTAA", 1);

        let true_SFS = [None, None, Some((3,11..12)), Some((3, 17..18)), Some((2, 8..9)), Some((3, 5..6))];

        // Put our SFS in 1-based indexing
        for item in SFS.iter_mut(){
            if let Some(X) = item{
                *item = Some((X.0, X.1.start + 1..X.1.end + 1));
            }

        }
        eprintln!("{:?}", SFS);

        assert_eq!(SFS, true_SFS);
    }

    #[test]
    fn serialize_and_load() {
        let seqs = vec![b"AGGTAAA".to_vec(), b"ACAGGTAGGAAAGGAAAGT".to_vec()];
        let fasta = to_fasta(seqs);
        let reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(fasta)).unwrap();

        let (sbwt , _) = Sbwt::<MatrixRank>::new::<2>(reader, 4, 3, 2, true, false);

        let mut buf = Vec::<u8>::new();
        sbwt.serialize(&mut buf).unwrap();
        let sbwt2 = Sbwt::<MatrixRank>::load(&mut buf.as_slice()).unwrap();

        assert_eq!(sbwt, sbwt2);
    }

    #[test]
    #[allow(non_snake_case)]
    fn non_ACGT(){
        let seqs = vec![b"AGGTAAA".to_vec(), b"ACAGGTAGGANAAGGAAAGT".to_vec()];           
        //.............................................................^...................

        let fasta = to_fasta(seqs);
        let reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(fasta)).unwrap();
        let (sbwt, _) = Sbwt::<MatrixRank>::new::<2>(reader, 4, 3, 2, true, false);

        let mut buf = Vec::<u8>::new();
        sbwt.serialize(&mut buf).unwrap();
        let sbwt2 = Sbwt::<MatrixRank>::load(&mut buf.as_slice()).unwrap();

        assert_eq!(sbwt, sbwt2);
    }
}
