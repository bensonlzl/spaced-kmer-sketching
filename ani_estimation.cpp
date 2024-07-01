#include "ani_estimator.hpp"

double containment(int intersection, int set_size){
    return ((double) intersection) / ((double) set_size);
}


double binomial_estimator(double containment, int kmer_num_ones){
    return std::pow(containment,(((double) 1.0)/((double) kmer_num_ones)));
}