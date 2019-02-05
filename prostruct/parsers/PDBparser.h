//
// Created by gil on 05/04/18.
//

#ifndef PROSTRUCT_PDBPARSER_H
#define PROSTRUCT_PDBPARSER_H

#include "prostruct/struct/chain.h"

struct AASequenceOrder
{
    bool operator()(std::string const& left, std::string const& right) const
    {
        std::string::size_type p_r  = right.find('-');
        std::string::size_type pp_r  = right.find('-', p_r+1);
        auto right_str = right.substr(p_r + 1, pp_r - p_r - 1);
        int right_int = std::stoi(right_str);

        std::string::size_type p_l = left.find('-');
        std::string::size_type pp_l  = left.find('-', p_l+1);
        auto left_str = left.substr(p_l + 1, pp_l - p_l - 1);
        int left_int = std::stoi(left_str);

        if (pp_r + 1 == right.size() && pp_l + 1 == left.size()) return right_int > left_int;// compare non-insertions -> ALA2 > ALA1
        else if (pp_r + 1 < right.size() && pp_l + 1 == left.size() && right_int == left_int) return true; // compare non-insertion with insertion -> ALA1A > ALA1
        else if (pp_l + 1 < left.size() && pp_r + 1 == right.size() && right_int == left_int) return false;// compare insertion with non-insertion -> ALA1 > ALA1A
        else if (pp_l + 1 < left.size() && pp_r + 1 < right.size() && right_int == left_int) return right.substr(pp_r+1, right.npos) > left.substr(pp_l+1, left.npos); // compare insertions -> ALA2A > ALA2A
        else return right_int > left_int;
    }
};

void createMap(std::string &, std::map<std::string, std::map<std::string, std::vector<std::shared_ptr<Atom>>, AASequenceOrder>> &,
               std::vector<std::string> &);


#endif //PROSTRUCT_PDBPARSER_H
