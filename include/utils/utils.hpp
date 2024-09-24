#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <fstream>

namespace utils {
  std::string parse_input(std::string const& input_path) {
    std::cout << "Parsing input...\n";
    uint32_t limit = 100000;
    std::ifstream input(input_path);
    std::string parsed = "";
    std::string in;
    uint32_t c = 0;
    while(std::getline(input, in) and c<limit) {
      if(in[0] == '>') continue;
      parsed = parsed + in;
      c++;
    }
    input.close();
    //std::cout << parsed << '\n';
    std::cout << "done.\n";
    return parsed;
  }
}

#endif //UTILS_HPP
