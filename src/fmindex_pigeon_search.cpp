#include <sstream>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

//this should also be okay
std::vector <std::vector <seqan3::dna5>> splice(int numSlice, std::vector <seqan3::dna5> query){
    int elementNum = query.size() / numSlice;
    std::vector <std::vector <seqan3::dna5>> result;
    std::vector <seqan3::dna5> temp;
    temp.push_back(query[0]);
    for (int i=1; i< query.size(); i++){
        if(i % elementNum == 0){
            result.push_back(temp);
            temp.clear();
        }
        temp.push_back(query[i]);
    }
    if (!temp.empty()) result.push_back(temp);
    return result;
}

//reference begin position returned will be the first position
//i think this should be good
// --------------------------copy pasted into mismatch
// int hamVerify(auto refPart, auto query){
//     int count = 0;
//     for (int i = 0; i < query.size(); i++){
//         if (refPart[i] != query[i]) count++;
//     }
//     return count;
// }


void mismatch(std::vector<std::vector<seqan3::dna5>> const& ref, std::vector<seqan3::dna5> const& query, auto& index, int k){
    std::vector <std::vector <seqan3::dna5>> qParts = splice(k+1, query);
    int totalLen = query.size();
    int shiftSize = totalLen / (k+1);
    for(int i = 0; i < qParts.size(); i++){
        auto results = seqan3::search(qParts[i], index);
        for(auto& res : results){
            auto shift = res.reference_begin_position() - i * shiftSize;
            if (shift >= 0){
                int count = 0;
                for (int j = 0; j < query.size(); j++){
                    if (ref[res.reference_id()][j + shift] != query[j]) count++;
                }
                if (count <= k){
                    std:: cout<< "found query at " << shift;
                }
            }
            
        }
    }
}

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read query n ref into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream) {
        reference.push_back(record.sequence());
    }

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    seqan3::configuration const cfg0 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};

    //!TODO here adjust the number of searches
    queries.resize(100); // will reduce the amount of searches
    int k=2;
    //only using 1 reference, for multiple queries
    for (auto& q : queries){
        mismatch(reference, q, index, k);
    }
    

    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text
    // Pseudo code (might be wrong):
    // for query in queries:
    //      parts[3] = cut_query(3, query);
    //      for p in {0, 1, 2}:
    //          for (pos in search(index, part[p]):
    //              if (verify(ref, query, pos +- ....something)):
    //                  std::cout << "found something\n"

    return 0;
}
