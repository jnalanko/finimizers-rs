use simple_sds::bit_vector::*;
use simple_sds::raw_vector::*;
use simple_sds::ops::*;
use simple_sds::serialize::*;

#[allow(clippy::len_without_is_empty)]
pub trait SubsetRank{

    // Todo: make char type into a generic unsigned integer.
    // Issues with that: it can't seem to figure out how to use such
    // generic integers as array indexes. Lol.

    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self;
    fn new_from_bit_vectors(vecs: Vec<simple_sds::bit_vector::BitVector>) -> Self;

    fn len(&self) -> usize; // Number of sets in the sequence
    fn rank(&self, c: u8, i: usize) -> usize;
    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>);

    fn access(&self, i: usize) -> Vec<u8>{
        let mut v = Vec::new();
        self.append_set_to_buf(i, &mut v);
        v
    }

    // Returns the number of bytes written
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>;

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> where Self: Sized;
}

#[derive(Eq, PartialEq, Debug)]
pub struct MatrixRank{
    rows: Vec<BitVector>
}

impl SubsetRank for MatrixRank{
    
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>{
        let n_rows = self.rows.len() as u64;
        let mut n_written = 0_usize;
        out.write_all(&n_rows.to_le_bytes())?;
        for row in self.rows.iter(){
            row.serialize(out)?;
            n_written += row.size_in_bytes();
        }
        Ok(n_written)
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>{
        let n_rows = u64::load(input)? as usize;

        let mut rows = Vec::<BitVector>::new();
        for _ in 0..n_rows{
            rows.push(BitVector::load(input)?);
        }
        Ok(Self{rows})
    }

    fn new_from_bit_vectors(mut rows: Vec<simple_sds::bit_vector::BitVector>) -> Self{
        
        for row in rows.iter_mut(){
            row.enable_rank();
        }

        Self{rows}
    }

    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self{
        let n = subset_seq.len();
        let mut rawrows = Vec::<RawVector>::new();
        for _ in 0..sigma{
            rawrows.push(RawVector::with_len(n, false));
        }
        for (i, set) in subset_seq.iter().enumerate(){
            for c in set.iter(){
                rawrows[*c as usize].set_bit(i, true)
            }
        }

        let rows: Vec<BitVector> = rawrows.into_iter().map(BitVector::from).collect();
        Self::new_from_bit_vectors(rows)
    }

    fn len(&self) -> usize{
        if self.rows.is_empty() { 0 }
        else { self.rows[0].len() }
    }

    fn rank(&self, c: u8, i: usize) -> usize{
        self.rows[c as usize].rank(i)
    }

    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>){
        for c in 0..self.rows.len(){
            if self.rows[c].get(i){
                buf.push(c as u8);
            }
        }
    }
}

impl std::fmt::Display for MatrixRank{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.rows.iter(){
            for i in 0..row.len(){
                if row.get(i){
                    write!(f, "1")?;
                } else{
                    write!(f, "0")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serialize_and_load(){
        let sets: Vec<Vec<u8>> = vec![vec![1,2,3], vec![0,2], vec![0,1,3,4], vec![], vec![0,1,2]];
        let sr = MatrixRank::new(sets, 5);
        let mut buf = Vec::<u8>::new();
        sr.serialize(&mut buf).unwrap();
        let sr2 = MatrixRank::load(&mut buf.as_slice()).unwrap();
        assert_eq!(sr, sr2);
    }
}