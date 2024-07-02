#include <string>
#include <vector>
#include <cstdint>
#include <iostream>
#include <fstream>

typedef std::vector<uint8_t> acgt_string;

std::vector<std::string> strings_from_fasta(const char fasta_filename[]);
void add_nucleotide_strings(std::vector<acgt_string> &return_strings, const std::string &raw_string);
std::vector<acgt_string> cut_nucleotide_strings(const std::vector<std::string> &raw_strings);
std::vector<acgt_string> nucleotide_strings_from_fasta_file(const char fasta_filename[]);