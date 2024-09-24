#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <fstream>

namespace utils {
  std::string parse_input(std::string const& input_path) {
    std::cout << "Parsing input...\n";
    std::ifstream input(input_path);
    std::string parsed = "";
    std::string in;
    while(std::getline(input, in)) {
      if(in[0] == '>') continue;
      parsed = parsed + in;
    }
    input.close();
    //std::cout << parsed << '\n';
    std::cout << "done.\n";
    return parsed;
  }
}

#endif //UTILS_HPP
