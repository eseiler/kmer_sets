// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>

#include "run.hpp"

void setup_parser(sharg::parser & parser, configuration & config)
{
    parser.info.author = "Enrico Seiler";
    parser.info.version = "1.0.0";

    parser.add_option(config.sample_A,
                      sharg::config{.short_id = '\0',
                                    .long_id = "sample_A",
                                    .description = "Sample A.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fastq"}}});

    parser.add_option(config.sample_B,
                      sharg::config{.short_id = '\0',
                                    .long_id = "sample_B",
                                    .description = "Sample B.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fastq"}}});

    parser.add_option(config.output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "Output.",
                                    .required = true,
                                    .validator = sharg::output_file_validator{{"fasta"}}});
}

int main(int argc, char ** argv)
{
    configuration config{};

    sharg::parser parser{"k-mer_set_difference", {argv, argv + argc}, sharg::update_notifications::off};

    try
    {
        setup_parser(parser, config);
        parser.parse();
        run(config);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR]: " << ext.what() << '\n';
        return -1;
    }

    return 0;
}
