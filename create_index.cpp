#include <iostream>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <sdsl/csa_wt.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/util.hpp> // Required for register_tmp_file

using namespace sdsl;
using fm_index_t = csa_wt<wt_huff<bit_vector_il<256>>, 512, 1024>;

int main(int argc, char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: " << argv[0] << " <input.fasta> <output.fm>\n\n"
                  << "Build an FM-index (SDSL csa_wt) over the nucleotide sequence(s)\n"
                  << "contained in a FASTA file and serialise it to a binary .fm file.\n\n"
                  << "Arguments:\n"
                  << "  input.fasta   Path to the source genome FASTA (single or multi-record).\n"
                  << "                Non-ACGT characters are stripped before indexing.\n"
                  << "  output.fm     Destination path for the serialised FM-index.\n\n"
                  << "Environment variables:\n"
                  << "  SDSL_CACHE_DIR   Directory for SDSL temporary construction files\n"
                  << "                   (defaults to SLURM_TMPDIR, then '.' if unset).\n\n"
                  << "Example:\n"
                  << "  ./create_index source.fasta source.fm\n"
                  << std::endl;
        return 0;
    }
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.fasta> <output.fm>  (use --help for details)" << std::endl;
        return 1;
    }
    std::string input_file = argv[1];
    std::string output_file = argv[2];

    const char* cache_dir_env = std::getenv("SDSL_CACHE_DIR");
    if (cache_dir_env == nullptr || std::string(cache_dir_env).empty()) {
        cache_dir_env = std::getenv("SLURM_TMPDIR");
    }
    std::string cache_dir = (cache_dir_env == nullptr || std::string(cache_dir_env).empty())
                                ? std::string(".")
                                : std::string(cache_dir_env);
    try {
        std::filesystem::create_directories(cache_dir);
    } catch (...) {
        cache_dir = ".";
    }

    cache_config config(false, cache_dir, util::basename(output_file));
    fm_index_t index;
    construct(index, input_file, config, 1); 
    
    if (store_to_file(index, output_file)) {
        std::cout << "✅ Successfully created index '" << output_file << "' from '" << input_file << "'" << std::endl;
        
        // --- CORRECTED FUNCTION NAME ---
        util::delete_all_files(config.file_map); // Changed from delete_files
        
        return 0;
    } else {
        std::cerr << "Error: Could not write index to " << output_file << std::endl;
        return 1;
    }
}