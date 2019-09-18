#ifndef NUCLEARUTILITIES
#define NUCLEARUTILITIES

#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace nme {

namespace GeneralUtilities {

inline std::vector<std::vector<std::string> > GetCSVData(std::string filename, std::string delimeter) {
  std::ifstream file(filename);

  std::vector<std::vector<std::string> > dataList;

  std::string line = "";
  // Iterate through each line and split the content using delimeter
  while (getline(file, line)) {
    std::vector<std::string> vec;
    boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
    dataList.push_back(vec);
  }
  // Close the File
  file.close();

  return dataList;
}
}
}
#endif
