/**
 * @file fasta_processing.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-07-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */
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
