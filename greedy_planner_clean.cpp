#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <omp.h>

namespace fs = std::filesystem;
using fm_index_t = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<256>>, 512, 1024>;

// Function Prototypes
std::map<std::string, std::string> read_fasta_and_clean(const std::string& path);
std::vector<fm_index_t> load_single_fm_index(const std::string& index_path);
bool query_kmer(const std::string& kmer, const std::vector<fm_index_t>& indexes);
double cost_synth(int length, double cost_per_base);
double solve_greedy_for_chromosome(const std::string& chrom_seq, int max_kmer_len, const std::vector<fm_index_t>& indexes, double cost_pcr, double cost_join, double cost_synth_per_base);

struct GreedyStats {
    double cost = 0.0;
    std::uint64_t reuse_moves = 0;
    std::uint64_t synth_moves = 0;
    std::uint64_t joins = 0;
    std::uint64_t segments = 0;
    std::uint64_t reuse_bases = 0;
    std::uint64_t synth_bases = 0;
    std::uint64_t length = 0;
};

GreedyStats solve_greedy_for_chromosome_stats(const std::string& chrom_seq, int W, const std::vector<fm_index_t>& indexes, double cost_pcr, double cost_join, double cost_synth_linear, double cost_synth_quad);

// Main Program
int main(int argc, char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: " << argv[0]
                  << " <W> <target.fasta> <pcr> <join> <synth_linear> [synth_quad] <source_index.fm>\n\n"
                  << "Replication-First greedy genome construction planner.\n"
                  << "At each position, greedily selects the longest reusable block (up to W bp);\n"
                  << "falls back to synthesis if no reusable block is found.\n\n"
                  << "Arguments:\n"
                  << "  W                Max block length (bp).\n"
                  << "  target.fasta     FASTA file of the genome to construct (target).\n"
                  << "  pcr              Fixed cost per reused (PCR-amplified) block.\n"
                  << "  join             Fixed cost per junction (not charged for the first block).\n"
                  << "  synth_linear     Per-base synthesis cost (linear term c_s). Cost = c_s * L.\n"
                  << "  synth_quad       [optional] Quadratic term c_s2. Cost = c_s*L + c_s2*L^2.\n"
                  << "                   Omit for purely linear synthesis cost.\n"
                  << "  source_index.fm  FM-index over the source genome (built with create_index).\n\n"
                  << "Output (CSV): filename, chromosome, length_bp, total_cost\n\n"
                  << "Examples:\n"
                  << "  # Linear synthesis cost:\n"
                  << "  ./greedy_planner_clean 500 target.fasta 5 1.5 0.2 source.fm\n\n"
                  << "  # Nonlinear synthesis cost (bacterial/eukaryotic experiments):\n"
                  << "  ./greedy_planner_clean 1000 target.fasta 5 1.5 0.2 1e-4 source.fm\n"
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

    GreedyStats total;
    total.cost = 0.0;
    for (const auto& pair : target_chromosomes) {
        std::string chrom_header = pair.first;
        const std::string& chrom_seq = pair.second;
        if (chrom_seq.empty()) continue;
        for (char &c : chrom_header) { if (c == ' ' || c == ',') { c = '_'; } }

        GreedyStats stats = solve_greedy_for_chromosome_stats(chrom_seq, W, indexes, cost_pcr_arg, cost_join_arg, cost_synth_linear_arg, cost_synth_quad_arg);
        
        std::cout << fs::path(fasta_path).filename().string() << "," << chrom_header << "," << chrom_seq.length() << "," << stats.cost << std::endl;

        total.cost += stats.cost;
        total.reuse_moves += stats.reuse_moves;
        total.synth_moves += stats.synth_moves;
        total.joins += stats.joins;
        total.segments += stats.segments;
        total.reuse_bases += stats.reuse_bases;
        total.synth_bases += stats.synth_bases;
        total.length += static_cast<std::uint64_t>(chrom_seq.length());
    }

    std::cout << "STATS_TOTAL,"
              << total.reuse_moves << ","
              << total.synth_moves << ","
              << total.joins << ","
              << total.segments << ","
              << total.reuse_bases << ","
              << total.synth_bases
              << std::endl;

    std::cout << fs::path(fasta_path).filename().string() << ",TOTAL," << total.length << "," << total.cost << std::endl;
    return 0;
}

// Function Implementations
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
double cost_synth(int length, double cost_per_base) {
    return static_cast<double>(length) * cost_per_base;
}

static inline double cost_synth_nonlinear(int length, double linear_per_base, double quad_coeff) {
    const double x = static_cast<double>(length);
    return (linear_per_base * x) + (quad_coeff * x * x);
}
double solve_greedy_for_chromosome(const std::string& chrom_seq, int W, const std::vector<fm_index_t>& indexes, double cost_pcr, double cost_join, double cost_synth_per_base) {
    long long N = chrom_seq.length();
    if (N == 0) return 0.0;
    double total_cost = 0.0;
    long long i = 0;
    while (i < N) {
        int best_w = 0;
        for (int w = std::min((long long)W, N - i); w >= 1; --w) {
            std::string block = chrom_seq.substr(i, w);
            if (query_kmer(block, indexes)) {
                best_w = w;
                break;
            }
        }
        if (best_w > 0) {
            total_cost += cost_pcr;
            if (i > 0) { total_cost += cost_join; }
            i += best_w;
        } else {
            int synth_len = 1;
            total_cost += cost_synth_nonlinear(synth_len, cost_synth_per_base, 0.0);
            if (i > 0) { total_cost += cost_join; }
            i += synth_len;
        }
    }
    return total_cost;
}

GreedyStats solve_greedy_for_chromosome_stats(const std::string& chrom_seq, int W, const std::vector<fm_index_t>& indexes, double cost_pcr, double cost_join, double cost_synth_linear, double cost_synth_quad) {
    GreedyStats stats;
    const long long N = static_cast<long long>(chrom_seq.length());
    stats.length = static_cast<std::uint64_t>(N);
    if (N == 0) return stats;

    double total_cost = 0.0;
    std::uint64_t segments = 0;
    std::uint64_t reuse_moves = 0;
    std::uint64_t synth_moves = 0;
    std::uint64_t reuse_bases = 0;
    std::uint64_t synth_bases = 0;

    long long i = 0;
    while (i < N) {
        int best_w = 0;
        for (int w = static_cast<int>(std::min<long long>(W, N - i)); w >= 1; --w) {
            std::string block = chrom_seq.substr(static_cast<size_t>(i), static_cast<size_t>(w));
            if (query_kmer(block, indexes)) {
                best_w = w;
                break;
            }
        }

        segments++;
        if (i > 0) { total_cost += cost_join; }

        if (best_w > 0) {
            total_cost += cost_pcr;
            reuse_moves++;
            reuse_bases += static_cast<std::uint64_t>(best_w);
            i += best_w;
        } else {
            const int synth_len = 1;
            // Nonlinear term has no effect for synth_len=1, but keep the model consistent.
            total_cost += cost_synth_nonlinear(synth_len, cost_synth_linear, cost_synth_quad);
            synth_moves++;
            synth_bases += static_cast<std::uint64_t>(synth_len);
            i += synth_len;
        }
    }

    stats.cost = total_cost;
    stats.segments = segments;
    stats.joins = (segments > 0) ? (segments - 1) : 0;
    stats.reuse_moves = reuse_moves;
    stats.synth_moves = synth_moves;
    stats.reuse_bases = reuse_bases;
    stats.synth_bases = synth_bases;
    return stats;
}