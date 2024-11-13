use std::{io::BufWriter, fs::File};

// Returns the number of bytes written
pub fn write_bytes<W: std::io::Write>(out: &mut W, bytes: &[u8]) -> std::io::Result<usize>{
    out.write_all(bytes)?;
    Ok(bytes.len() + 8)
}

// TODO TODO TODO THIS HAS A RARE RACE CONDITION: FIX IT!!!
pub fn create_temp_file(dir: &std::path::Path) -> (tempfile::NamedTempFile<()>, BufWriter<File>) {
    let f = tempfile::Builder::new().rand_bytes(10).prefix("sbwt-temp-").make_in(dir, |path| {
        if path.is_file(){
            return Err(std::io::ErrorKind::AlreadyExists.into());
        }
        
        Ok(())
    }).unwrap();

    let writer = std::io::BufWriter::new(std::fs::File::create(f.path()).unwrap());
    (f, writer)
}

// Searcher the range [0..n)
// Return the index of the answer, or n if does not exist
pub fn binary_search_leftmost_that_fulfills_pred<T, Access: Fn(usize) -> T, Pred: Fn(T) -> bool>(access: Access, pred: Pred, n: usize) -> usize {
    let mut ans = n;
    let mut step = n;
    while step > 0 {
        while ans as isize - step as isize >= 0 && pred(access(ans-step)) {
            ans -= step;
        }
        step /= 2;
    }
    ans
}