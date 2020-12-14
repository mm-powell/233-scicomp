#include "method.h"
#include <map>
#include <math.h>
#include <numeric>
#include <algorithm>

using namespace std;

method::method()
{

}






std::vector<double> method::polarity(std::vector<std::string> pos, std::vector<std::string> neg, std::vector<std::string> tweets) {
    double pos_count = compare(pos, tweets);
    double neg_count = compare(neg, tweets);
    std::vector<double> pol_scores;
    pol_scores.push_back(pos_count/ (pos_count+neg_count) );
    pol_scores.push_back(neg_count/ (pos_count+neg_count) );
    return pol_scores;
}







void method::remove_punct(std::vector<std::string> tweets) {
    for (int j = 0, len = tweets.size(); j < len; j++) {
        for (int i = 0, len = tweets[j].size(); i < len; i++)
        {
            if (ispunct(tweets[j][i]))
            {
                tweets[j].erase(i--, 1);
                len = tweets[j].size();
            }
        }
    }
}
//https://stackoverflow.com/questions/19138983/c-remove-punctuation-from-string





void method::remove_stop(std::vector<std::string> stop_words, std::vector<std::string> tweets) {
    for (int k = 0; k < stop_words.size(); ++k) {
        tweets.erase(std::remove(tweets.begin(), tweets.end(), stop_words[k]), tweets.end());
    }

}





std::vector<std::string> method::tokenize(std::vector<std::string> tweets) {
    std::vector<std::string> token_vec;
    for (int i = 0; i < tweets.size(); ++i) {
        char* str;
        std::string str_obj(tweets[i]);
        str = &str_obj[0];
        //char str[] = char(tweets[i]);
        char *token = strtok(str, " ");
        while (token != NULL) {
            token_vec.push_back(token);
            token = strtok(NULL, " ");
        }
    }
    return token_vec;
}
//https://www.geeksforgeeks.org/tokenizing-a-string-cpp/






int method::compare(std::vector<std::string> emo_words, std::vector<std::string> tweets) {
    int k = 0;
    for (int i = 0; i < emo_words.size(); ++i) {
        for (int j = 0; j < tweets.size(); ++j) {
            if (emo_words[i] == tweets[j]) {
                k = k+1;
            }
        }
    }
    return k;
}






std::vector<std::string> method::compare_with_vec(std::vector<std::string> emo_words, std::vector<std::string> tweets) {
    std::vector<std::string> word_vec;
    for (int i = 0; i < emo_words.size(); ++i) {
        for (int j = 0; j < tweets.size(); ++j) {
            if (emo_words[i] == tweets[j]) {
                word_vec.push_back(emo_words[i]);
            }
        }
    }
    return word_vec;
}






double method::jsd(std::vector<std::string> candidate1, std::vector<std::string> candidate2) {


    // keep NONUNIQUE lists, with repeats
    std::vector<std::string> candidate1_orig = candidate1;
    std::vector<std::string> candidate2_orig = candidate2;


    // create vector of ALL tokens, with repeats
    std::vector<std::string> total_vec;
    total_vec.insert( total_vec.end(), candidate1.begin(), candidate1.end() );
    total_vec.insert( total_vec.end(), candidate2.begin(), candidate2.end() );



    // compute pi parameters
    double pi_cand1 = (double (candidate1.size()) / double (total_vec.size()) );
    double pi_cand2 = (double (candidate2.size()) / double (total_vec.size()) );


    // create NONUINQUE vector for counts
    std::vector<std::string> total_vec_c;
    total_vec_c.insert( total_vec_c.end(), candidate1.begin(), candidate1.end() );
    total_vec_c.insert( total_vec_c.end(), candidate2.begin(), candidate2.end() );
    sort( total_vec_c.begin(), total_vec_c.end() );



    // sort
    sort( candidate1.begin(), candidate1.end() );
    candidate1.erase( unique( candidate1.begin(), candidate1.end() ), candidate1.end() );

    sort( candidate2.begin(), candidate2.end() );
    candidate2.erase( unique( candidate2.begin(), candidate2.end() ), candidate2.end() );




    // keep only unique tokens
    sort( total_vec.begin(), total_vec.end() );
    total_vec.erase( unique( total_vec.begin(), total_vec.end() ), total_vec.end() );



    // counting occurences
    std::map<std::string, int> total_map;
    for (int i = 0; i < total_vec.size(); ++i) {
        total_map.insert ( std::pair<std::string, int> ( total_vec[i], count(total_vec_c.begin(), total_vec_c.end(), total_vec[i]) ) );
    }

    std::map<std::string, int> cand1_map;
    for (int i = 0; i < candidate1.size(); ++i) {
        cand1_map.insert ( std::pair<std::string, int> ( candidate1[i], count(candidate1_orig.begin(), candidate1_orig.end(), candidate1[i]) ) );
    }

    std::map<std::string, int> cand2_map;
    for (int i = 0; i < candidate2.size(); ++i) {
        cand2_map.insert ( std::pair<std::string, int> ( candidate2[i], count(candidate2_orig.begin(), candidate2_orig.end(), candidate2[i]) ) );
    }


    // compute probabilities

    std::map<std::string, double> total_probs;
    for (int j = 0; j < total_vec.size(); ++j) {
        total_probs.insert ( std::pair<std::string, double> (total_vec[j], double(total_map[total_vec[j]]) / double(total_vec_c.size()) ) );
    }

    std::map<std::string, double> cand1_probs;
    for (int j = 0; j < candidate1.size(); ++j) {
        cand1_probs.insert ( std::pair<std::string, double> (candidate1[j], double(cand1_map[candidate1[j]]) / double(candidate1_orig.size()) ) );
    }

    std::map<std::string, double> cand2_probs;
    for (int j = 0; j < candidate2.size(); ++j) {
        cand2_probs.insert ( std::pair<std::string, double> (candidate2[j], double(cand2_map[candidate2[j]]) / double(candidate2_orig.size()) ) );
    }



    // creating mixed distribution

    std::vector<double> cand1_mix;
    for (int k = 0; k < cand1_probs.size(); ++k) {
        cand1_mix.push_back( double((double)pi_cand1 * (double)( (double)cand1_probs[candidate1[k]] * (log2( (double)cand1_probs[candidate1[k]] / (double)total_probs[candidate1[k]])) ) ) );
    }


    std::vector<double> cand2_mix;
    for (int k = 0; k < cand2_probs.size(); ++k) {
        cand2_mix.push_back( (double)pi_cand2 * (double)( (double)cand2_probs[candidate2[k]] * (log2( (double)cand2_probs[candidate2[k]] / (double)total_probs[candidate2[k]])) ) );
    }


    return (double)std::accumulate(cand1_mix.begin(), cand1_mix.end(), 0.) + (double)std::accumulate(cand2_mix.begin(), cand2_mix.end(), 0.);




}



