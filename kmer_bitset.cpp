#include "kmer.hpp"
/*
Implementation of k-mers using bitsets

Warning: Lots of bit hacks are used
*/




// For each length 0 <= l <= MAX_KMER_LENGTH, store kmer_bitset with exactly 2l 1s
kmer_bitset contiguous_kmer_array[MAX_KMER_LENGTH+1]; 

// Function to initialise contiguous_kmer_array
// TO DO: Initialise this at compile time
void initialise_contiguous_kmer_array(){
    if (LOGGING) std::clog << INFO_LOG << "Initialising contiguous_kmer_array" << std::endl;
    for (int i = 0; i <= MAX_KMER_LENGTH; ++i){
        contiguous_kmer_array[i] = kmer_bitset(KMER_BITSET_SIZE);
        if (LOGGING) std::clog << "Initialising contiguous_kmer_array " << i << " with size " << contiguous_kmer_array[i].size() << std::endl;

    }
    if (LOGGING) std::clog << INFO_LOG << "Initialising prefixes" << std::endl;
    for (int i = 1; i <= MAX_KMER_LENGTH; ++i){
        contiguous_kmer_array[i] |= contiguous_kmer_array[i-1];
        for (int j = 0; j < NUCLEOTIDE_BIT_SIZE; ++j){
            contiguous_kmer_array[i][j + (i-1) * NUCLEOTIDE_BIT_SIZE] = 1;
        }
    }
    if (LOGGING) std::clog << INFO_LOG << "contiguous_kmer_array initialised" << std::endl;
}

// Function to obtain contiguous kmers
kmer_bitset contiguous_kmer(const int kmer_length){
    if (kmer_length > MAX_KMER_LENGTH) throw std::runtime_error("Given k-mer length exceeds maximum k-mer length");
    return contiguous_kmer_array[kmer_length];
}

// For each power of two from 1 to KMER_BITSET_SIZE / 2, 
// compute a bitset consisting of alternating 1s and 0s 
// with run length equal to that power of two.
// This is for reversing the kmer in LOG_KMER_BITSET_SIZE operations
kmer_bitset reversing_kmer_array[LOG_KMER_BITSET_SIZE]; 
kmer_bitset invert_reversing_kmer_array[LOG_KMER_BITSET_SIZE]; 

// Function to initialise reversing_kmer_array
// TO DO: Initialise this at compile time
void initialise_reversing_kmer_array(){
    if (LOGGING) std::clog << INFO_LOG << "Initialising reversing_kmer_array" << std::endl;

    for (int i = 0; i < LOG_KMER_BITSET_SIZE; ++i){
        reversing_kmer_array[i] = kmer_bitset(KMER_BITSET_SIZE);
        invert_reversing_kmer_array[i] = kmer_bitset(KMER_BITSET_SIZE);
    }

    int array_idx = 0;
    for (int gap_size = NUCLEOTIDE_BIT_SIZE; gap_size < KMER_BITSET_SIZE; gap_size *= 2, ++array_idx){
        if (LOGGING) std::clog << INFO_LOG << "Initialising reversing_kmer_array " << array_idx << " with size " << reversing_kmer_array[array_idx].size() << std::endl;  
        if (LOGGING) std::clog << std::flush;
        for (uint32_t gap_idx = 0; gap_idx < KMER_BITSET_SIZE / gap_size; ++gap_idx){
            for (int pos_idx = gap_idx * gap_size; pos_idx < (gap_idx + 1) * gap_size; ++pos_idx){
                reversing_kmer_array[array_idx][pos_idx] = (gap_idx & 0x1); // Alternate between blocks of 0s and 1s
                invert_reversing_kmer_array[array_idx][pos_idx] = ~reversing_kmer_array[array_idx][pos_idx] ; // Alternate between blocks of 0s and 1s
            }
        }
        if (DEBUG) std::cout << "Reversing array " << array_idx << " set to " << reversing_kmer_array[array_idx] << std::endl;
        if (DEBUG) std::cout << "Invert reversing array " << array_idx << " set to " << invert_reversing_kmer_array[array_idx] << std::endl;
    }
}

// Function to reverse a kmer_bitset (used for reverse complementation)
kmer_bitset reverse_kmer_bitset(kmer_bitset kbs){
    int array_idx = 0;
    kmer_bitset cur = kbs;
    if (DEBUG) std::cout << "kbs    = " << kbs << std::endl;

    if (DEBUG) std::cout << "cur    = " << cur << std::endl;
    // log KMER_BITSET_SIZE iterations of reversing in blocks of powers of 2
    
    for (int gap_size = NUCLEOTIDE_BIT_SIZE; gap_size < KMER_BITSET_SIZE; gap_size *= 2, ++array_idx){
        if (DEBUG) std::cout << "f half = " << (cur & reversing_kmer_array[array_idx]) << std::endl;
        if (DEBUG) std::cout << "s half = " << (cur & invert_reversing_kmer_array[array_idx]) << std::endl;

        cur = (
            ((cur & reversing_kmer_array[array_idx]) >> gap_size) 
            | ((cur & (invert_reversing_kmer_array[array_idx])) << gap_size)
        );
        if (DEBUG) std::cout << "cur    = " << cur << std::endl;
    }
    return cur;
}
