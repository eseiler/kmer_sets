// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>

struct configuration
{
    std::filesystem::path sample_A{};
    std::filesystem::path sample_B{};
    std::filesystem::path output{};
};
