//
// Created by Gil Hoben on 2019-02-06.
//

#include "prostruct/parsers/PDBparser.h"

#include <cstdio>
#include <fstream>
#include <sstream>

inline void remove_whitespace(std::string &str) {
  str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
}

template <typename T>
void createMap(
    const std::string &fname,
    std::map<std::string, std::map<std::string, atomVector<T>, AASequenceOrder>> &  chainResMap,
    std::vector<std::string> &chainOrder) {

  std::ifstream file(fname);

  std::string line;

  if (!file.is_open())
    throw "File does not exist!";

  while (std::getline(file, line)) {

    auto record = line.substr(0, 5);

    remove_whitespace(record);

    if (record == "ATOM") {

      auto serialStr = line.substr(6, 5);
      //            int serial = std::stoi(serialStr);
      auto name = line.substr(12, 4);
      //            auto altLocation = line.at(16);
      auto residue = line.substr(17, 3);
      auto chainID = line.substr(21, 1);
      auto resSeqStr = line.substr(22, 4);
      //            int resSeq = std::stoi(resSeqStr);
      auto insCode = line.substr(26, 1);
      T x = std::stod(line.substr(30, 8));
      T y = std::stod(line.substr(38, 8));
      T z = std::stod(line.substr(46, 8));
      //            double occupancy = std::stod(line.substr(54, 6));
      //            double tempFactor = std::stod(line.substr(60, 6));
      std::string element = line.substr(76, 2);
      //            auto charge = line.substr(78, 2);

      remove_whitespace(name);
      remove_whitespace(serialStr);
      remove_whitespace(element);
      remove_whitespace(resSeqStr);
      remove_whitespace(insCode);

      if (std::find(chainOrder.begin(), chainOrder.end(), chainID) ==
          chainOrder.end()) {
        chainOrder.push_back(chainID);
      }

      auto residueID = residue + "-" + resSeqStr + "-" + insCode;
      auto atom = std::make_shared<Atom<T>>(element, name, x, y, z);

      chainResMap[chainID][residueID].emplace_back(atom);
    }
  }
}

template void createMap(const std::string &, std::map<std::string, std::map<std::string, atomVector<float>, AASequenceOrder>> &,
                        std::vector<std::string> &);

template void createMap(const std::string &, std::map<std::string, std::map<std::string, atomVector<double>, AASequenceOrder>> &,
                        std::vector<std::string> &);
