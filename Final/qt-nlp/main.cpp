#include <iostream>
#include <vector>
#include <fstream>
#include <method.h>
#include <string>
#include <istream>
#include <iterator>
#include <map>
#include <utility>
#include <math.h>
#include <numeric>


using namespace std;

int main ()
{

    // IMPORTING STOP WORDS
    string myText_stop;
    std::vector<std::string> stop_vec;
    ifstream stop_words_vec("/Users/mpowell2/Desktop/SciComp/final/data/stop.txt");
    while (getline (stop_words_vec, myText_stop)) {
      stop_vec.push_back(myText_stop);
    }
    stop_words_vec.close();






    //////////// IMPORTING EMOTION WORDS
    ///
    // NEGATIVE
    string myText_neg;
    std::vector<std::string> neg_vec;
    ifstream neg_words("/Users/mpowell2/Desktop/SciComp/final/test-1/neg_words.txt");
    while (getline (neg_words, myText_neg)) {
      neg_vec.push_back(myText_neg);
    }
    neg_words.close();


    // POSITIVE
    string myText_pos;
    std::vector<std::string> pos_vec;
    ifstream pos_words("/Users/mpowell2/Desktop/SciComp/final/test-1/pos_words.txt");
    while (getline (pos_words, myText_pos)) {
      pos_vec.push_back(myText_pos);
    }
    pos_words.close();


    // ANGER
    string myText_ang;
    std::vector<std::string> ang_vec;
    ifstream ang_words("/Users/mpowell2/Desktop/SciComp/final/test-1/ang_words.txt");
    while (getline (ang_words, myText_ang)) {
      ang_vec.push_back(myText_ang);
    }
    ang_words.close();


    // ANTICIPATION
    string myText_ant;
    std::vector<std::string> ant_vec;
    ifstream ant_words("/Users/mpowell2/Desktop/SciComp/final/test-1/ant_words.txt");
    while (getline (ant_words, myText_ant)) {
      ant_vec.push_back(myText_ant);
    }
    ant_words.close();


    // DISGUST
    string myText_dis;
    std::vector<std::string> dis_vec;
    ifstream dis_words("/Users/mpowell2/Desktop/SciComp/final/test-1/dis_words.txt");
    while (getline (dis_words, myText_dis)) {
      dis_vec.push_back(myText_dis);
    }
    dis_words.close();


    // FEAR
    string myText_fea;
    std::vector<std::string> fea_vec;
    ifstream fea_words("/Users/mpowell2/Desktop/SciComp/final/test-1/fea_words.txt");
    while (getline (fea_words, myText_fea)) {
      fea_vec.push_back(myText_fea);
    }
    fea_words.close();


    // JOY
    string myText_joy;
    std::vector<std::string> joy_vec;
    ifstream joy_words("/Users/mpowell2/Desktop/SciComp/final/test-1/joy_words.txt");
    while (getline (joy_words, myText_joy)) {
      joy_vec.push_back(myText_joy);
    }
    joy_words.close();


    // SADNESS
    string myText_sad;
    std::vector<std::string> sad_vec;
    ifstream sad_words("/Users/mpowell2/Desktop/SciComp/final/test-1/sad_words.txt");
    while (getline (sad_words, myText_sad)) {
      sad_vec.push_back(myText_sad);
    }
    sad_words.close();


    // SURPRISE
    string myText_sur;
    std::vector<std::string> sur_vec;
    ifstream sur_words("/Users/mpowell2/Desktop/SciComp/final/test-1/sur_words.txt");
    while (getline (sur_words, myText_sur)) {
      sur_vec.push_back(myText_sur);
    }
    sur_words.close();


    // TRUST
    string myText_tru;
    std::vector<std::string> tru_vec;
    ifstream tru_words("/Users/mpowell2/Desktop/SciComp/final/test-1/tru_words.txt");
    while (getline (tru_words, myText_tru)) {
      tru_vec.push_back(myText_tru);
    }
    tru_words.close();



    //////////// IMPORTING TWEETS

    // JOE BIDEN
    string myText_jb;
    std::vector<std::string> jb_vec;
    ifstream jb_tweets("/Users/mpowell2/Desktop/SciComp/final/test-1/dem_tweets/jb_tweets.txt");
    while (getline (jb_tweets, myText_jb)) {
      jb_vec.push_back(myText_jb);
    }
    jb_tweets.close();

    // BERNIE SANDERS
    string myText_bs;
    std::vector<std::string> bs_vec;
    ifstream bs_tweets("/Users/mpowell2/Desktop/SciComp/final/test-1/dem_tweets/bs_tweets.txt");
    while (getline (bs_tweets, myText_bs)) {
      bs_vec.push_back(myText_bs);
    }
    bs_tweets.close();

    // ELIZABETH WARREN
    string myText_ew;
    std::vector<std::string> ew_vec;
    ifstream ew_tweets("/Users/mpowell2/Desktop/SciComp/final/test-1/dem_tweets/ew_tweets.txt");
    while (getline (ew_tweets, myText_ew)) {
        ew_vec.push_back(myText_ew);
    }
    ew_tweets.close();

    // KAMALA HARRIS
    string myText_kh;
    std::vector<std::string> kh_vec;
    ifstream kh_tweets("/Users/mpowell2/Desktop/SciComp/final/test-1/dem_tweets/kh_tweets.txt");
    while (getline (kh_tweets, myText_kh)) {
      kh_vec.push_back(myText_kh);
    }
    kh_tweets.close();













    // CLEAN DATA
    // Remove punctuation, remove stop words, and tokensize

    method met;

    met.remove_punct(jb_vec);
    met.remove_stop(stop_vec, jb_vec);
    std::vector<std::string> jb_tweets_clean = met.tokenize(jb_vec);

    met.remove_punct(bs_vec);
    met.remove_stop(stop_vec, bs_vec);
    std::vector<std::string> bs_tweets_clean = met.tokenize(bs_vec);

    met.remove_punct(ew_vec);
    met.remove_stop(stop_vec, ew_vec);
    std::vector<std::string> ew_tweets_clean = met.tokenize(ew_vec);

    met.remove_punct(kh_vec);
    met.remove_stop(stop_vec, kh_vec);
    std::vector<std::string> kh_tweets_clean = met.tokenize(kh_vec);






    // Polarity

    std::vector<double> pol_jb = met.polarity(pos_vec, neg_vec, jb_tweets_clean);





    // Sentiment Analysis

    std::vector<std::string> bs_emo_trust_vec = met.compare_with_vec(tru_vec, bs_tweets_clean);
    double bs_emo_trust = met.compare(tru_vec, bs_tweets_clean);





    // Jensen-Shannon Divergence

    double div_KE = met.jsd(kh_tweets_clean, ew_tweets_clean);
    double div_KJ = met.jsd(kh_tweets_clean, jb_tweets_clean);
    double div_KB = met.jsd(kh_tweets_clean, bs_tweets_clean);
    double div_BE = met.jsd(bs_tweets_clean, ew_tweets_clean);
    double div_BJ = met.jsd(bs_tweets_clean, jb_tweets_clean);
    double div_JE = met.jsd(jb_tweets_clean, ew_tweets_clean);


    std::cout << "Kamala/Elizabeth = " << div_KE << std::endl;
    std::cout << "Kamala/Joe = " << div_KJ << std::endl;
    std::cout << "Kamala/Bernie = " << div_KB << std::endl;
    std::cout << "Bernie/Elizabeth = " << div_BE << std::endl;
    std::cout << "Bernie/Joe = " << div_BJ << std::endl;
    std::cout << "Joe/Elizabeth = " << div_JE << std::endl;












}

