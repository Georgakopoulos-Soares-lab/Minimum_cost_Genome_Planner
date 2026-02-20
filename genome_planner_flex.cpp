#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cmath>
#include <map>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>

namespace fs = std::filesystem;
using fm_index_t = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<256>>, 512, 1024>;

// --- FUNCTION PROTOTYPES ---
std::map<std::string, std::string> read_fasta_and_clean(const std::string& path);
std::vector<fm_index_t> load_single_fm_index(const std::string& index_path);
bool query_kmer(const std::string& kmer, const std::vector<fm_index_t>& indexes);
double cost_synth(int length, double cost_per_base);

struct PlannerStats {
    double cost = 0.0;
    std::uint64_t reuse_moves = 0;
    std::uint64_t synth_moves = 0;
    std::uint64_t joins = 0;
    std::uint64_t segments = 0;
    std::uint64_t reuse_bases = 0;
    std::uint64_t synth_bases = 0;
    std::uint64_t length = 0;
};

PlannerStats solve_dp_for_chromosome(const std::string& chrom_seq, int max_kmer_len, const std::vector<fm_index_t>& indexes, double cost_pcr, double cost_join, double cost_synth_linear, double cost_synth_quad);

// --- MAIN PROGRAM ---
int main(int argc, char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: " << argv[0]
                  << " <W> <target.fasta> <pcr> <join> <synth_linear> [synth_quad] <source_index.fm>\n\n"
                  << "Optimal (DP) minimum-cost genome construction planner.\n"
                  << "Partitions the target genome into blocks of length <= W, choosing reuse\n"
                  << "(PCR) or synthesis for each block to minimise total cost.\n\n"
                  << "Arguments:\n"
                  << "  W                Max block length (bp). Reflects experimental PCR/synthesis limits.\n"
                  << "  target.fasta     FASTA file of the genome to construct (target).\n"
                  << "  pcr              Fixed cost per reused (PCR-amplified) block, regardless of length.\n"
                  << "                   A block is reusable if it occurs as an exact substring of the source.\n"
                  << "  join             Fixed cost per junction between adjacent blocks.\n"
                  << "                   Not charged for the first block (no preceding junction).\n"
                  << "  synth_linear     Per-base synthesis cost coefficient (linear term c_s).\n"
                  << "                   Synthesis cost = c_s * L  (for a block of length L).\n"
                  << "  synth_quad       [optional] Quadratic synthesis cost coefficient (c_s2).\n"
                  << "                   Synthesis cost = c_s * L + c_s2 * L^2.\n"
                  << "                   Omit (or set to 0) for purely linear synthesis cost.\n"
                  << "  source_index.fm  FM-index file built over the source genome (via create_index).\n\n"
                  << "Output (CSV, one row per chromosome/record plus a TOTAL row):\n"
                  << "  filename, chromosome, length_bp, total_cost\n\n"
                  << "Examples:\n"
                  << "  # Linear synthesis cost (virus experiments):\n"
                  << "  ./genome_planner_flex 500 target.fasta 5 1.5 0.2 source.fm\n\n"
                  << "  # Nonlinear synthesis cost (bacterial sweep, W=1000):\n"
                  << "  ./genome_planner_flex 1000 target.fasta 5 1.5 0.2 1e-4 source.fm\n\n"
                  << "  # Nonlinear synthesis cost (eukaryotic sweep, W=800):\n"
                  << "  ./genome_planner_flex 800 target.fasta 5 1.5 0.2 1e-4 source.fm\n"
                  << std::endl;
        return 0;
    }
    if (argc != 7 && argc != 8) {
        std::cerr << "Usage: " << argv[0]
                  << " <W> <target.fasta> <pcr> <join> <synth_linear> [synth_quad] <source_index.fm>"
                  << "  (use --help for details)" << std::endl;
        return 1;
    }
    int W = std::stoi(argv[1]);
    std::string fasta_path = argv[2];
    double cost_pcr_arg = std::stod(argv[3]);
    double cost_join_arg = std::stod(argv[4]);
    double cost_synth_linear_arg = std::stod(argv[5]);
    double cost_synth_quad_arg = 0.0;
    std::string index_path_arg;
    if (argc == 7) {
        index_path_arg = argv[6];
    } else {
        cost_synth_quad_arg = std::stod(argv[6]);
        index_path_arg = argv[7];
    }

    std::vector<fm_index_t> indexes = load_single_fm_index(index_path_arg);
    if (indexes.empty()) { return 1; }
    
    std::map<std::string, std::string> target_chromosomes = read_fasta_and_clean(fasta_path);

    PlannerStats total;
    total.cost = 0.0;

    for (const auto& pair : target_chromosomes) {
        std::string chrom_header = pair.first; // Make a mutable copy
        const std::string& chrom_seq = pair.second;
        if (chrom_seq.empty()) continue;

        // --- NEW: Clean the header for safe CSV output ---
        for (char &c : chrom_header) {
            if (c == ' ' || c == ',') {
                c = '_';
            }
        }

        PlannerStats stats = solve_dp_for_chromosome(chrom_seq, W, indexes, cost_pcr_arg, cost_join_arg, cost_synth_linear_arg, cost_synth_quad_arg);
        
        // Output in CSV format
        std::cout << fs::path(fasta_path).filename().string() << ","
                  << chrom_header << ","
                  << chrom_seq.length() << ","
                  << stats.cost << std::endl;

        total.cost += stats.cost;
        total.reuse_moves += stats.reuse_moves;
        total.synth_moves += stats.synth_moves;
        total.joins += stats.joins;
        total.segments += stats.segments;
        total.reuse_bases += stats.reuse_bases;
        total.synth_bases += stats.synth_bases;
        total.length += static_cast<std::uint64_t>(chrom_seq.length());
    }

    // Stats line intended for scripts to parse (not CSV of same schema).
    std::cout << "STATS_TOTAL,"
              << total.reuse_moves << ","
              << total.synth_moves << ","
              << total.joins << ","
              << total.segments << ","
              << total.reuse_bases << ","
              << total.synth_bases
              << std::endl;

    // Final TOTAL line keeps the historical 4-column CSV schema and is last.
    std::cout << fs::path(fasta_path).filename().string() << ",TOTAL,"
              << total.length << "," << total.cost << std::endl;

    return 0;
}

// --- FUNCTION IMPLEMENTATIONS ---

std::map<std::string, std::string> read_fasta_and_clean(const std::string& path) {
    std::map<std::string, std::string> sequences;
    std::ifstream fasta_file(path);
    if (!fasta_file.is_open()) { exit(1); }
    std::string line, header, current_sequence;
    while (std::getline(fasta_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) sequences[header] = current_sequence;
            header = line.substr(1);
            current_sequence.clear();
        } else {
            for (char c : line) {
                char uc = std::toupper(static_cast<unsigned char>(c));
                if (uc == 'A' || uc == 'T' || uc == 'C' || uc == 'G') current_sequence += uc;
            }
        }
    }
    if (!header.empty()) sequences[header] = current_sequence;
    return sequences;
}

std::vector<fm_index_t> load_single_fm_index(const std::string& index_path) {
    std::vector<fm_index_t> indexes;
    fm_index_t index;
    if (sdsl::load_from_file(index, index_path)) {
        indexes.push_back(std::move(index));
    } else {
        std::cerr << "ERROR: Could not load index file: " << index_path << std::endl;
    }
    return indexes;
}

bool query_kmer(const std::string& kmer, const std::vector<fm_index_t>& indexes) {
    for (const auto& index : indexes) {
        if (sdsl::count(index, kmer.begin(), kmer.end()) > 0) return true;
    }
    return false;
}

// Fast check for all substring lengths ending at position i (exclusive):
// Uses incremental backward_search over the FM-index (CSA), so querying
// lengths 1..W costs O(W) rank operations total (not O(W^2)).
static inline void update_min_cost_for_i(
    const std::string& chrom_seq,
    long long i,
    int W,
    const fm_index_t& index,
    const std::vector<double>& DP,
    std::vector<long long>& pred,
    std::vector<uint16_t>& chosen_len,
    std::vector<uint8_t>& chosen_is_reuse,
    double cost_pcr,
    double cost_join,
    double cost_synth_linear,
    double cost_synth_quad,
    double& min_cost_for_i
) {
    const int max_w = static_cast<int>(std::min<long long>(W, i));
    if (max_w <= 0) return;

    // Interval for empty pattern is the full suffix array range.
    fm_index_t::size_type l = 0;
    fm_index_t::size_type r = index.size() - 1;

    for (int w = 1; w <= max_w; ++w) {
        const long long j = i - w;
        const char c = chrom_seq[static_cast<size_t>(j)];

        fm_index_t::size_type l2 = 0, r2 = 0;
        const auto occ = sdsl::backward_search(index, l, r, static_cast<fm_index_t::char_type>(c), l2, r2);
        l = l2;
        r = r2;

        const bool is_reuse = (occ > 0);
        const double acquisition_cost = is_reuse
                            ? cost_pcr
                            : (cost_synth_linear * static_cast<double>(w) + cost_synth_quad * static_cast<double>(w) * static_cast<double>(w));
        const double path_cost = DP[static_cast<size_t>(j)] + acquisition_cost + (j > 0 ? cost_join : 0.0);
        if (path_cost < min_cost_for_i) {
            min_cost_for_i = path_cost;
            pred[static_cast<size_t>(i)] = j;
            chosen_len[static_cast<size_t>(i)] = static_cast<uint16_t>(w);
            chosen_is_reuse[static_cast<size_t>(i)] = static_cast<uint8_t>(is_reuse ? 1 : 0);
        }

        // Once no match exists for length w, longer strings (w+1,...) also cannot match.
        if (occ == 0) {
            for (int w2 = w + 1; w2 <= max_w; ++w2) {
                const long long j2 = i - w2;
                const double ww = static_cast<double>(w2);
                const double path_cost2 = DP[static_cast<size_t>(j2)] + (cost_synth_linear * ww + cost_synth_quad * ww * ww) + (j2 > 0 ? cost_join : 0.0);
                if (path_cost2 < min_cost_for_i) {
                    min_cost_for_i = path_cost2;
                    pred[static_cast<size_t>(i)] = j2;
                    chosen_len[static_cast<size_t>(i)] = static_cast<uint16_t>(w2);
                    chosen_is_reuse[static_cast<size_t>(i)] = 0;
                }
            }
            return;
        }
    }
}

double cost_synth(int length, double cost_per_base) {
    return static_cast<double>(length) * cost_per_base;
}

PlannerStats solve_dp_for_chromosome(const std::string& chrom_seq, int W, const std::vector<fm_index_t>& indexes, double cost_pcr, double cost_join, double cost_synth_linear, double cost_synth_quad) {
    PlannerStats stats;
    const long long N = static_cast<long long>(chrom_seq.length());
    stats.length = static_cast<std::uint64_t>(N);
    if (N == 0) return stats;
    std::vector<double> DP(N + 1, 1e18);
    std::vector<long long> pred(N + 1, -1);
    std::vector<uint16_t> chosen_len(N + 1, 0);
    std::vector<uint8_t> chosen_is_reuse(N + 1, 0);
    DP[0] = 0.0;
    // Common case in this codebase: one source FM-index per run.
    // If multiple indexes are provided, we fall back to the original logic
    // (still without substr allocations) to preserve semantics.
    const bool single_index = (indexes.size() == 1);

    for (long long i = 1; i <= N; ++i) {
        double min_cost_for_i = 1e18;

        if (single_index) {
            update_min_cost_for_i(chrom_seq, i, W, indexes[0], DP, pred, chosen_len, chosen_is_reuse, cost_pcr, cost_join, cost_synth_linear, cost_synth_quad, min_cost_for_i);
        } else {
            const int max_w = static_cast<int>(std::min<long long>(W, i));
            for (int w = 1; w <= max_w; ++w) {
                const long long j = i - w;
                bool reusable = false;
                for (const auto& index : indexes) {
                    if (sdsl::count(index,
                                    chrom_seq.begin() + j,
                                    chrom_seq.begin() + i) > 0) {
                        reusable = true;
                        break;
                    }
                }
                const double ww = static_cast<double>(w);
                const double acquisition_cost = reusable ? cost_pcr : (cost_synth_linear * ww + cost_synth_quad * ww * ww);
                const double path_cost = DP[static_cast<size_t>(j)] + acquisition_cost + (j > 0 ? cost_join : 0.0);
                if (path_cost < min_cost_for_i) {
                    min_cost_for_i = path_cost;
                    pred[static_cast<size_t>(i)] = j;
                    chosen_len[static_cast<size_t>(i)] = static_cast<uint16_t>(w);
                    chosen_is_reuse[static_cast<size_t>(i)] = static_cast<uint8_t>(reusable ? 1 : 0);
                }
            }
        }

        DP[static_cast<size_t>(i)] = min_cost_for_i;
    }

    stats.cost = DP[static_cast<size_t>(N)];

    // Backtrack to count moves.
    std::uint64_t segments = 0;
    std::uint64_t reuse_moves = 0;
    std::uint64_t synth_moves = 0;
    std::uint64_t reuse_bases = 0;
    std::uint64_t synth_bases = 0;

    long long cur = N;
    while (cur > 0) {
        const uint16_t len = chosen_len[static_cast<size_t>(cur)];
        const uint8_t is_reuse = chosen_is_reuse[static_cast<size_t>(cur)];
        const long long p = pred[static_cast<size_t>(cur)];

        if (p < 0 || len == 0) {
            // Should not happen, but avoid infinite loops.
            break;
        }

        segments++;
        if (is_reuse) {
            reuse_moves++;
            reuse_bases += len;
        } else {
            synth_moves++;
            synth_bases += len;
        }
        cur = p;
    }

    stats.segments = segments;
    stats.joins = (segments > 0) ? (segments - 1) : 0;
    stats.reuse_moves = reuse_moves;
    stats.synth_moves = synth_moves;
    stats.reuse_bases = reuse_bases;
    stats.synth_bases = synth_bases;
    return stats;
}