use crate::kmer::*;
use crate::kmer_chunk::*;
use crate::util;
use tempfile::NamedTempFile;
use std::io::BufWriter;
use std::thread;
use jseqio::reader::*;

pub fn split_to_bins<const B: usize>(mut seqs: jseqio::reader::DynamicFastXReader, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool) -> Vec<NamedTempFile<()>>{

    // Suppose we have a memory budget of m bytes and t threads.
    // Suppose each k-mer takes s bytes and there are 64 bins.
    // Let b be the number of k-mers in each splitter thread bin buffer.
    // A splitter thread uses 64bs bytes
    // In total the splitter threads use 64bst threads.
    // So, we need:
    // 64bbt = m
    // b = m / (64bt)

    let bin_prefix_len = 3; // If you update this you must update all the logic below
    let n_bins = (4_usize).pow(bin_prefix_len); // 64
    let producer_buf_size = 1_000_000_usize; // TODO: respect this
    let encoder_bin_buf_size = mem_gb * (1_usize << 30) / ((n_bins * LongKmer::<B>::byte_size()) * n_threads);

    log::info!("Splitting k-mers into {} bins", n_bins);
    log::info!("Bin buffer size: {}", encoder_bin_buf_size);

    use crossbeam::crossbeam_channel::bounded;
    let (parser_out, encoder_in) = bounded(4);
    let (encoder_out, writer_in) = bounded(4);

    let temp_dir = std::path::PathBuf::from("temp"); // TODO: read from config

    // Create producer
    let producer_handle = thread::spawn(move || {
        let mut buf = Vec::<Box<[u8]>>::new();
        let mut current_total_buffer_size = 0_usize;
        while let Some(rec) = seqs.read_next().unwrap(){
            current_total_buffer_size += rec.seq.len();
            let mut seq = rec.seq.to_owned();
            seq.reverse(); // Reverse to get colex sorting
            buf.push(seq.into_boxed_slice());
            if current_total_buffer_size > producer_buf_size{
                let mut sendbuf = Vec::<Box<[u8]>>::new();
                std::mem::swap(&mut sendbuf, &mut buf);
                parser_out.send(sendbuf).unwrap();
                current_total_buffer_size = 0;
            }
        }
        parser_out.send(buf).unwrap();
        drop(parser_out);
    });

    // Create encoder-splitters
    let mut encoder_handles = Vec::<thread::JoinHandle::<()>>::new();
    for _ in 0..n_threads{
        let receiver_clone = encoder_in.clone();
        let sender_clone = encoder_out.clone();
        encoder_handles.push(std::thread::spawn(move || {
            let mut bin_buffers = vec![Vec::<LongKmer::<B>>::new(); n_bins];
            for buf in bin_buffers.iter_mut(){
                buf.reserve_exact(encoder_bin_buf_size);
            }
            while let Ok(batch) = receiver_clone.recv(){
                for seq in batch{
                    for kmer in seq.windows(k){
                        match LongKmer::<B>::from_ascii(kmer) {
                            Ok(kmer) => {
                                let bin_id = kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize; // Interpret nucleotides in base-4
                                bin_buffers[bin_id].push(kmer);
                                if bin_buffers[bin_id].len() == encoder_bin_buf_size{
                                    if dedup_batches{
                                        bin_buffers[bin_id].sort_unstable();
                                        bin_buffers[bin_id].dedup();
                                    }
                                    sender_clone.send(bin_buffers[bin_id].clone()).unwrap();
                                    bin_buffers[bin_id].clear();
                                }
                            }
                            Err(KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
                            Err(KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
                        }        
                    }
                }
            }

            // Send remaining buffers
            for mut b in bin_buffers{
                if dedup_batches{
                    b.sort_unstable();
                    b.dedup();
                }
                sender_clone.send(b).unwrap();
            }
        }));
    }

    // Create writer
    let mut bin_filenames = 
        Vec::<NamedTempFile<()>>::new();
    let mut bin_writers = 
        Vec::<std::io::BufWriter::<std::fs::File>>::new();

    for _ in 0..n_bins{
        let (f, w) = util::create_temp_file(&temp_dir);
        bin_filenames.push(f);
        bin_writers.push(w);
    }


    let writer_handle = thread::spawn( move || {
        while let Ok(batch) = writer_in.recv(){
            if !batch.is_empty() {
                let bin_id = batch[0].get_from_left(0) as usize * 16 + batch[0].get_from_left(1) as usize * 4 + batch[0].get_from_left(2) as usize; // Intepret nucleotides in base-4
                let bin_file = &mut bin_writers[bin_id];
                for kmer in batch{
                    kmer.serialize(bin_file).unwrap(); // Todo: write all at once
                }
            }
        }
    });

    producer_handle.join().unwrap();
    drop(encoder_in); // Close the channel
    for h in encoder_handles{
        h.join().unwrap();
    }
    drop(encoder_out); // Close the channel
    writer_handle.join().unwrap();

    bin_filenames
}

// https://stackoverflow.com/questions/59205184/how-can-i-change-the-number-of-threads-rayon-uses
/*pub fn create_pool(n_threads: usize) -> rayon::ThreadPool {
    match rayon::ThreadPoolBuilder::new()
       .num_threads(n_threads)
       .build()
    {
       Err(e) => panic!("Failed to create thread pool: {}", e)
       Ok(pool) => pool,
    }
 }*/

// Overwrite the files with sorted and deduplicates files
pub fn par_sort_and_dedup_bin_files<const B: usize>(bin_files: &[NamedTempFile<()>], mem_gb: usize, n_threads: usize){
    let mut files_and_sizes: Vec<(std::path::PathBuf, usize)> = 
        bin_files.iter().map(|f| (f.path().to_owned(), f.path().metadata().unwrap().len() as usize)).collect();
    files_and_sizes.sort_by_key(|(_, size)| *size);

    let max_mem = mem_gb * (1_usize << 30);

    log::info!("Sorting k-mer bins");

    use crossbeam::unbounded;

    // A work queue
    let (queue_sender, queue_recvr) = unbounded::<(std::path::PathBuf, usize)>();

    // A queue to notify the producer that a bin has been processed.
    // The usize in the channel is the size of the bin that was processed.
    let (producer_notify, producer_recv_notify) = unbounded::<usize>();

    // Wrap in mutex to share between threads
    let mut total_size_in_processing = 0_usize;

    // Start the producer
    let producer_handle = thread::spawn(move || {
        while !files_and_sizes.is_empty() {
            // Push as much work to the queue as possible
            while !files_and_sizes.is_empty(){    
                let (f,s) = files_and_sizes.last().unwrap();
                let (f,s) = (f.clone(), *s);
                if total_size_in_processing == 0 || total_size_in_processing + s <= max_mem {
                    queue_sender.send((f, s)).unwrap();
                    files_and_sizes.pop();
                    total_size_in_processing += s;
                } else {
                    break;
                }
            }

            let s_done = producer_recv_notify.recv().unwrap(); // Wait for a notification
            total_size_in_processing -= s_done;
        }

        // All files have been pushed to the channel
        drop(queue_sender); // Close the channel
    });

    let mut consumer_handles = Vec::<thread::JoinHandle<()>>::new();

    // Spawn consumers
    for _ in 0..n_threads{
        let recv_clone = queue_recvr.clone();
        let producer_notify = producer_notify.clone();

        consumer_handles.push(std::thread::spawn( move || {
            while let Ok((f, s)) = recv_clone.recv(){
                log::info!("Sorting bin {} of size {}", f.to_str().unwrap(), s);
                let mut reader = std::io::BufReader::new(std::fs::File::open(&f).unwrap());
                let chunk = KmerChunk::<B>::load(&mut reader).unwrap();
                drop(reader); // Flush
        
                let mut chunk = chunk.sort();
                chunk.dedup();
                let chunk_out = std::io::BufWriter::new(std::fs::File::create(&f).unwrap()); // Overwrite
                chunk.serialize(chunk_out).unwrap();

                // Notify the producer that s bytes have been processed and
                // new work can possibly be pushed to the queue.
                let _ = producer_notify.send(s); // This may fail if the producer has already exited. That is ok.
            }
        }));
    }

    producer_handle.join().unwrap();
    for h in consumer_handles{
        h.join().unwrap();
    }

}

// The original files are deleted
pub fn concat_files(infiles: Vec<NamedTempFile<()>>, mut out_writer: BufWriter<std::fs::File>){
    for fp in infiles {
        let mut reader = std::io::BufReader::new(std::fs::File::open(fp).unwrap());
        std::io::copy(&mut reader, &mut out_writer).unwrap();
        // fp is dropped here, which deletes the file
    }
}
