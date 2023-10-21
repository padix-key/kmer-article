use std::convert::TryInto;
use std::time::Instant;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use json;
use rand::prelude::*;
use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha8Rng;


#[derive(Debug, Clone)]
struct Alphabet {
    symbols: Vec<u8>,
    symbol_to_code: [u8; 256],
    illegal_code: u8
}

impl Alphabet {
    fn from(symbols: Vec<u8>) -> Alphabet{
        // The last symbol code of the alphabet + 1 is always illegal
        // Since this code cannot occur from symbol encoding,
        // it can be later used to check for illegal symbols
        let illegal_code: u8 = symbols.len().try_into().unwrap();
        // An array based map that maps from symbol to code
        // Since the maximum value of a char is 256 the size of the map is known at compile time
        // Initially fill the map with the illegal symbol
        // Consequently, the map will later return the illegal symbol,
        // when indexed with a character that is not part of the alphabet
        let mut symbol_to_code: [u8; 256] = [illegal_code; 256];
        // Then fill in entries for the symbols of the alphabet
        for (i, symbol) in symbols.iter().enumerate() {
            symbol_to_code[*symbol as usize] = i.try_into().unwrap();
        }
        Alphabet {symbols, symbol_to_code, illegal_code}
    }

    fn encode(&self, sequence: &[u8]) -> Result<Vec<u8>, String>{
        let mut code: Vec<u8> = Vec::with_capacity(sequence.len());
        for symbol in sequence.iter() {
            unsafe {
                // This is actually safe: 'symbol_to_code' comprises all 256 u8 values
                let c = self.symbol_to_code.get_unchecked(*symbol as usize);
                if *c == self.illegal_code {
                    return Err(format!("Symbol {} is not in the alphabet", *symbol as char));
                }
                code.push(*c)
            }
        }
        Ok(code)
    }

    fn size(&self) -> usize {
        return self.symbols.len()
    }
}


#[derive(Debug, Clone)]
struct KmerAlphabet {
    alphabet: Alphabet,
    k: usize,
    radix_multipliers: Vec<usize>
}

impl KmerAlphabet {
    fn new(alphabet: Alphabet, k: usize) -> KmerAlphabet {
        let mut radix_multipliers: Vec<usize> = Vec::new();
        for i in 0..k {
            radix_multipliers.insert(0, alphabet.size().pow(i as u32));
        }
        KmerAlphabet {alphabet, k, radix_multipliers}
    }

    fn decompose_naive(&self, sequence: &[u8]) -> Result<Vec<usize>, String> {
        let seq_code = self.alphabet.encode(sequence)?;
        let mut kmers: Vec<usize> = Vec::with_capacity(sequence.len() - self.k + 1);

        for i in 0..kmers.capacity() {
            let mut c: usize = 0;
            for j in 0..self.k {
                c += seq_code[i+j] as usize * self.radix_multipliers[j];
            }
            kmers.push(c);
        }

        Ok(kmers)
    }

    fn decompose_fast(&self, sequence: &[u8]) -> Result<Vec<usize>, String> {
        let seq_code = self.alphabet.encode(sequence)?;
        let mut kmers: Vec<usize> = Vec::with_capacity(sequence.len() - self.k + 1);
        let end_radix_multiplier: usize = self.alphabet.size().pow((self.k - 1) as u32);

        // Compute initial kmer code using the naive way
        let mut c: usize = 0;
        for j in 0..self.k {
            c += seq_code[j] as usize * self.radix_multipliers[j];
        }
        kmers.push(c);

        // Compute subsequent kmer codes using the information from the previous one
        let mut prev_kmer: usize = kmers[0];
        for i in 1..kmers.capacity() {
            let kmer: usize =
                (
                    // Remove first symbol
                    (prev_kmer - seq_code[i - 1] as usize * end_radix_multiplier)
                    // Shift k-mer to left
                    * self.alphabet.size()
                )
                // Add new symbol
                + seq_code[i + self.k - 1] as usize;
            kmers.push(kmer);
            prev_kmer = kmer;
        }

        Ok(kmers)
    }

    fn decompose_bitshift(&self, sequence: &[u8]) -> Result<Vec<usize>, String> {
        if self.alphabet.size() != 4 {
            return Err(String::from("Bit shift is only implemented for alphabets of size 4"))
        }

        let seq_code = self.alphabet.encode(sequence)?;
        let mut kmers: Vec<usize> = Vec::with_capacity(sequence.len() - self.k + 1);
        let end_radix_multiplier: usize = self.alphabet.size().pow((self.k - 1) as u32);

        // Compute initial kmer code using the naive way
        let mut c: usize = 0;
        for j in 0..self.k {
            c += seq_code[j] as usize * self.radix_multipliers[j];
        }
        kmers.push(c);

        // Compute subsequent kmer codes using the information from the previous one
        let mut prev_kmer: usize = kmers[0];
        for i in 1..kmers.capacity() {
            let kmer: usize =
                (
                    // Remove first symbol
                    (prev_kmer - seq_code[i - 1] as usize * end_radix_multiplier)
                    // Shift k-mer to left
                    << 2
                )
                // Add new symbol
                + seq_code[i + self.k - 1] as usize;
            kmers.push(kmer);
            prev_kmer = kmer;
        }

        Ok(kmers)
    }
}


const REP: usize = 100_000;
const SEQ_LEN: usize = 1000;
const MAX_K: usize = 32;


fn benchmark_decomposition<F>(decompose_func: F, sequence: &[u8], repetitions: usize) -> usize
where
    F: Fn(&[u8]) -> Result<Vec<usize>, String>
{
    let now = Instant::now();
    for _ in 0..repetitions {
        decompose_func(&sequence).unwrap();
    }
    now.elapsed().as_nanos() as usize / repetitions
}


fn main() {
    let alphabet = Alphabet::from(vec![b'A', b'C', b'G', b'T']);
    let mut rng = ChaCha8Rng::seed_from_u64(0);
    let sequence: Vec<u8> = (0..SEQ_LEN)
        .map(|_| {
            let idx = rng.gen_range(0..alphabet.symbols.len());
            alphabet.symbols[idx] as u8
        })
        .collect();

    // Perform benchmark
    let mut k_list: Vec<usize> = Vec::new();
    let mut naive_time_ns_list: Vec<usize> = Vec::new();
    let mut fast_time_ns_list: Vec<usize> = Vec::new();
    let mut bitshift_time_ns_list: Vec<usize> = Vec::new();
    for k in 1..=MAX_K {
        let kmer_alphabet = KmerAlphabet::new(alphabet.clone(), k);

        k_list.push(k);
        naive_time_ns_list.push(
            benchmark_decomposition(|seq| kmer_alphabet.decompose_naive(seq), &sequence, REP)
        );
        fast_time_ns_list.push(
            benchmark_decomposition(|seq| kmer_alphabet.decompose_fast(seq), &sequence, REP)
        );
        bitshift_time_ns_list.push(
            benchmark_decomposition(|seq| kmer_alphabet.decompose_bitshift(seq), &sequence, REP)
        );
    }

    // Prepare JSON
    let mut json_content = HashMap::new();
    json_content.insert(String::from("k"), k_list);
    json_content.insert(String::from("naive"), naive_time_ns_list);
    json_content.insert(String::from("fast"), fast_time_ns_list);
    json_content.insert(String::from("bitshift"), bitshift_time_ns_list);
    // Write JSON
    let mut file = File::create("benchmark.json").unwrap();
    file.write_all(&json::stringify_pretty(json_content, 4).as_bytes()).unwrap()
}