#ifndef METHOD_H
#define METHOD_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <istream>

class method
{
public:
    method();
    std::vector<double> polarity(std::vector<std::string> pos, std::vector<std::string> neg, std::vector<std::string> tweets);
    void remove_punct(std::vector<std::string> tweets);
    void remove_stop(std::vector<std::string> stop_words, std::vector<std::string> tweets);
    std::vector<std::string> tokenize(std::vector<std::string> tweets);
    int compare(std::vector<std::string> emo_words, std::vector<std::string> tweets);
    std::vector<std::string> compare_with_vec(std::vector<std::string> emo_words, std::vector<std::string> tweets);
    double jsd(std::vector<std::string> candidate1, std::vector<std::string> candidate2);
};

#endif // METHOD_H
