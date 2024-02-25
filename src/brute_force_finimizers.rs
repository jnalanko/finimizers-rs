use sbwt::sbwt::*;
use sbwt::subsetrank::*;
use bitvec::prelude::*;

#[allow(dead_code)]
fn sbwt_count(index: &Sbwt::<MatrixRank>, pattern: &[u8]) -> usize {
    match index.search(pattern){
        Some(interval) => interval.len(),
        None => 0,
    }
} 

#[allow(dead_code)]
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
#[allow(dead_code)]
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