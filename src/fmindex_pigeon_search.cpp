#include <sstream>
#include <cmath>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

//need to be checked
std::vector <std::vector <seqan3::dna5>> splice(int numSlice, std::vector <seqan3::dna5> query){
    float aa = query.size();
    int elementNum = (int) std:: ceil( (float) aa / numSlice); //check this again
    std::vector <std::vector <seqan3::dna5>> result;
    std::vector <seqan3::dna5> temp;
    temp.push(query[0]);
    for (int i=1, i< query.size(), i++){
        if(i % elementNum ==0){
            result.push(temp);
            temp.clear();
        }
        temp.push(query[i]);
    }
    return result;
}

//reference begin position returned will be the first position
//i think this should be good
std :: vector <seqan3::search_result> indelVerify(auto res1, auto res2, int k, int partLen){
    std :: vector <seqan3::search_result> res;
    for (auto && result: res1){
        for (int i =0, i < res2.size(), i++){
            r2 = res2[i].reference_begin_position();
            if (result.reference_begin_position()+ partLen <= r2 && result.reference_begin_position()+ partLen + k >= r2){
                res.push(result);
            }
        }
    }
    return res;
}

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

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
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
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
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    //!TODO here adjust the number of searches
    queries.resize(100); // will reduce the amount of searches


    //--------------------------------------check splice funct first ------------------------------------

    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text
    // Pseudo code (might be wrong):
    std :: vector<auto> results;
    int k=2;
    for (int i = 0; i< queries.size(); i++){
        //for each query
        std::vector <std::vector <seqan3::dna5>> div = splice (k+1, queries[i]);
        for (int j = 0 ; j < k+1 ; j++){
            results.push(seqan3::search(queries, index, cfg0));
        }
        std::vector <seqan3::search_result> res = indelVerify(results[0], results[1], 1, div[0].size());
        if (res.empty()) std::cout<<"lol no matches";
        else{
            std::vector <seqan3::search_result> res2 = indelVerify(res, results[2], 1, div[1].size());
            if (res2.empty()){
                std::cout<<"lol no matches";
            }
            else
            {
                for (int k = 0, k <res2.size(), k++) seqan3::debug_stream << "\tFound query " << res[i].query_id() << " at " << res[i].reference_begin_position() << "\n";
            }
        }
    }
    
    // for query in queries:
    //      parts[3] = cut_query(3, query);
    //      for p in {0, 1, 2}:
    //          for (pos in search(index, part[p]):
    //              if (verify(ref, query, pos +- ....something)):
    //                  std::cout << "found something\n"

    return 0;
}
