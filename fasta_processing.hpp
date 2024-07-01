#include <string>
#include <vector>
#include <cstdint>
#include <iostream>
#include <fstream>

std::vector<std::string> strings_from_fasta(const char fasta_filename[]);
void add_nucleotide_strings(std::vector<std::vector<uint8_t>> &return_strings, const std::string &raw_string);
std::vector<std::vector<uint8_t>> cut_nucleotide_strings(const std::vector<std::string> &raw_strings);