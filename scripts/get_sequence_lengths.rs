use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: get_sequence_lengths <filename>");
        std::process::exit(1);
    }

    let file = File::open(&args[1]).expect("Could not open file");

    let mut sequence = String::new();
    let mut lengths = Vec::new();

    for line in BufReader::new(file).lines() {
        let line = line.expect("Could not read line");
        if line.starts_with('>') {
            if !sequence.is_empty() {
                lengths.push(sequence.len());
                sequence.clear();
            }
        } else {
            sequence.push_str(&line);
        }
    }

    if !sequence.is_empty() {
        lengths.push(sequence.len());
    }

    for length in lengths {
        println!("{}", length);
    }
}