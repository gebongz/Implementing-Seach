#include <sstream>
#include <chrono>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

//this should also be okay
std::vector <std::vector <seqan3::dna5>> splice(int numSlice, std::vector <seqan3::dna5> query){
    int qSize = query.size() ;
    int elementNum = qSize / numSlice;
    std::vector <std::vector <seqan3::dna5>> result;
    std::vector <seqan3::dna5> temp;

    temp.push_back(query[0]);
    for (int i=1; i< qSize; i++){
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


// void mismatch(std::vector<std::vector<seqan3::dna5>> const& ref, std::vector<seqan3::dna5> const& query, auto& index, int k){
//     std::vector <std::vector <seqan3::dna5>> qParts = splice(k+1, query);
//     int totalLen = query.size();
//     int shiftSize = totalLen / (k+1);
//     int numQParts = qParts.size();
//     for(int i = 0; i < numQParts; i++){
//         auto results = seqan3::search(qParts[i], index);
//         //seqan3::debug_stream << "search finished\n";
//         for(auto& res : results){
//             auto beginPos = res.reference_begin_position();
//             auto shift = beginPos - i * shiftSize;
//             int count = 0;
//             int j=0;
//             //match left
//             while (j < (i * shiftSize) && count <= k){
//                 if (ref[res.reference_id()][j + shift] != query[j]) count++;
//                 j++;
//             }
//             if (count<=k){
//                 //match right
//                 j += qParts[i].size()-1;
//                 while (j < totalLen && count <= k ){
//                     if (ref[res.reference_id()][j + shift] != query[j]) count++;
//                     j++;
//                 }
//             }
//             if (count <= k){
//                 seqan3::debug_stream << "found query at " << shift <<" with "<< count <<" errors\n";
//             }
//         }
//     }
// }

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

    auto queryNum = int{100};
    parser.add_option(queryNum, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto k = int{0};
    parser.add_option(k, '\0', "errors", "number of allowed hamming distance errors");

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
    queries.resize(queryNum); // will reduce the amount of searches
    
    const std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    for (auto& q : queries){
        //seqan3::debug_stream << "query ID: "<<q <<"\n" ;
        //mismatch(reference, q, index, k);
        std::vector <std::vector <seqan3::dna5>> qParts = splice(k+1, q);
        int totalLen = q.size();
        int shiftSize = totalLen / (k+1);
        int numQParts = k+1;
        for(int i = 0; i < numQParts; i++){
            auto results = seqan3::search(qParts[i], index);
            //seqan3::debug_stream << "search finished\n";
            for(auto& res : results){
                auto beginPos = res.reference_begin_position();
                auto shift = beginPos - i * shiftSize;
                int count = 0;
                int j=0;
                //match left
                while (j < (i * shiftSize) && count <= k){
                    if (reference[res.reference_id()][j + shift] != q[j]) count++;
                    j++;
                }
                if (count<=k){
                    //match right
                    j += qParts[i].size()-1;
                    while (j < totalLen && count <= k ){
                        if (reference[res.reference_id()][j + shift] != q[j]) count++;
                        j++;
                    }
                }
                if (count <= k){
                    seqan3::debug_stream << "found query at " << shift <<" with "<< count <<" errors\n";
                }
            }
        }
    }
    const auto end = std::chrono::steady_clock::now();

    seqan3::debug_stream << "The Program took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms.\n";

    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text
    // Pseudo code (might be wrong):
    // for query in queries:
    //      parts[3] = cut_query(3, query);
    //      for p in {0, 1, 2}:
    //          for (pos in search(index, part[p]):
    //              if (verify(ref, query, pos +- ....something)):
    //                  std::cout << "found something\n"
    seqan3::debug_stream << "function end";
    return 0;
}
