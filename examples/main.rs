use bio::io::fastq;
use minima::SupermerIter;
use nthash::NtHashIterator;
use std::env;

// This is a really dull example. We just generate all the supermers
// for each read.
//
// We should do something a bit more useful.
fn main() {
    let args: Vec<String> = env::args().collect();
    let file = args.get(1).expect("Please provide a file");
    let reader = fastq::Reader::from_file(file).expect("A readable fq");
    let k = 13;
    let s = 40;
    for res in reader.records() {
        let rec = res.unwrap();
        let h = NtHashIterator::new(rec.seq(), k).unwrap();
        for _ in SupermerIter::new(h, k, s) {
	}
    }
}
