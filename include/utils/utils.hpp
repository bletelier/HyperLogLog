#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <hyperloglog/hyperloglog.hpp>

namespace utils {
  void parse_input(sketch::hyperloglog &sketch, std::string const& input_path) {
    std::ifstream input(input_path);
    std::string in;
    while(std::getline(input, in)) {
      if(in[0] == '>') continue;
      uint32_t in_size = (uint32_t) in.size();
      for(uint32_t i = 0; i < in_size; ++i) {
        char c = in[i];
        if(c == 'A' or c == 'C' or c == 'G' or c == 'T')
          sketch.read_stream(c);
      }
    }
    input.close();
  }

  double estimate_jaccard(sketch::hyperloglog &A, sketch::hyperloglog &B) {
    sketch::hyperloglog AuB(A,B);
    double est_card_A = A.estimate_cardinality();
    double est_card_B = B.estimate_cardinality();
    double est_card_AuB = AuB.estimate_cardinality();
    return (est_card_A + est_card_B - est_card_AuB)/(est_card_AuB); 
  }
  
  double real_jaccard(sketch::hyperloglog &A, sketch::hyperloglog &B) {
    sketch::hyperloglog AuB(A,B);
    uint64_t real_card_A = A.real_cardinality();
    uint64_t real_card_B = B.real_cardinality();
    uint64_t real_card_AuB = AuB.real_cardinality();
    return (1.0*(real_card_A + real_card_B - real_card_AuB))/(1.0 * real_card_AuB);
  }

  std::pair<double,double> calculate_errors(std::vector<sketch::hyperloglog> &A_i) {
    uint8_t N = A_i.size();
    std::cout << "Calculating errors for " << +((N*(N-1))/2) << " combinations of genomes\n";
    double erm = 0.0;
    double eam = 0.0;
    double prom_jh = 0.0, prom_jr = 0.0;
    uint16_t tot = 0;
    for(uint8_t i = 0; i < N; ++i) {
      for(uint8_t j = i+1; j < N; ++j) {
        double jaccard_hat = std::max(0.0, estimate_jaccard(A_i[i], A_i[j]));
        double jaccard_real = real_jaccard(A_i[i], A_i[j]);
        double diff = 1.0 * std::abs(jaccard_hat - jaccard_real);
        erm += jaccard_real == 0 ? diff : diff/jaccard_real;
        eam += diff;
        tot++;
        std::cout << "\tCombination " << tot << " done (gen" << +(i+1) << ", gen" << +(j+1)<< ")\n";
      } 
    }
    std::cout << "done\n";
    return {(1.0/tot)*erm, (1.0/tot)*eam};
  } 

}

#endif //UTILS_HPP
