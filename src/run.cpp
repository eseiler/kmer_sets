// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "run.hpp"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

static constexpr uint64_t kmer_size{32u};

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

std::vector<uint64_t> get_kmers(std::filesystem::path const & path)
{
    std::vector<uint64_t> result{};
    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>> fin{path};
    auto adaptor = seqan3::views::minimiser_hash(seqan3::ungapped{kmer_size}, seqan3::window_size{kmer_size}, seqan3::seed{0u});
    // auto adaptor = seqan3::views::minimiser_hash(seqan3::ungapped{kmer_size},
    //                                              seqan3::window_size{kmer_size},
    //                                              seqan3::seed{0x8F'3F'73'B5'CF'1C'9A'DEULL});

    for (auto && [seq] : fin)
    {
        result.append_range(adaptor(seq));
    }

    std::ranges::sort(result);
    auto const [first, last] = std::ranges::unique(result);
    result.erase(first, last);
    result.shrink_to_fit();
    return result;
}

void convert_kmer_to_dna4(uint64_t kmer, std::vector<seqan3::dna4> & vector)
{
    assert(vector.size() == kmer_size);
    for (size_t i = 0; i < kmer_size; ++i)
    {
        vector[kmer_size - 1 - i].assign_rank(kmer & 0b11);
        kmer >>= 2;
    }
}

void convert_dna_to_string(std::vector<seqan3::dna4> const & vector, std::string & str)
{
    assert(vector.size() == str.size());
    std::ranges::transform(vector,
                           str.begin(),
                           [](seqan3::dna4 letter)
                           {
                               return seqan3::to_char(letter);
                           });
}

void run(configuration const & config)
{
    std::vector<uint64_t> const sample_A = get_kmers(config.sample_A);
    std::vector<uint64_t> const sample_B = get_kmers(config.sample_B);
    std::vector<uint64_t> A_not_in_B{};

    std::ranges::set_difference(sample_A, sample_B, std::back_inserter(A_not_in_B));

    seqan3::sequence_file_output file_output{config.output};
    seqan3::dna4_vector sequence(kmer_size);
    std::string id(kmer_size, char{});
    for (auto && elem : A_not_in_B)
    {
        convert_kmer_to_dna4(elem, sequence);
        convert_dna_to_string(sequence, id);
        file_output.emplace_back(sequence, id);
    }
}
