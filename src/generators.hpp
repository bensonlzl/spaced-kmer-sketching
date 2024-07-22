/**
 * @file generators.hpp
 * @author Benson Lin (bensonlinzl@gmail.com)
 * @brief 
 * @date 2024-07-15
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "stl_includes.hpp"

/**
 * @brief 
 * Helper function to compute adjacent pairwise combinations in a list (vector)
 * 
 * @tparam T 
 * @param v list of elements (vector)
 * @return a pair of vectors that whose corresponding pairs of elements represent every pair of elements in v
 */
template<typename T>
std::pair<std::vector<T>, std::vector<T>> generate_pairwise_from_vector(
    const std::vector<T> &v
){
    std::vector<T> v1, v2;

    int sz = v.size();

    for (size_t i = 0; i < v.size(); ++i){
        v1.push_back(v[i]);
        v2.push_back(v[(i+1)%sz]);
    }

    return make_pair(v1,v2);
}

/**
 * @brief 
 * Helper function to compute all pairwise combinations in a list (vector)
 * 
 * @tparam T 
 * @param v list of elements (vector)
 * @return a pair of vectors that whose corresponding pairs of elements represent every pair of elements in v
 */
template<typename T>
std::pair<std::vector<T>, std::vector<T>> generate_all_pairs_from_vector(
    const std::vector<T> &v
){
    std::vector<T> v1, v2;

    for (size_t i = 0; i < v.size(); ++i){
        for (size_t j = 0; j < v.size(); ++j){
            v1.push_back(v[i]);
            v2.push_back(v[j]);
        }
    }

    return make_pair(v1,v2);
}
