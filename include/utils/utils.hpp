#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <hyperloglog/hyperloglog.hpp>

namespace utils {
  void parse_input(sketch::hyperloglog &sketch, std::string const& input_path) {
    std::cout << "Reading stream " << input_path << '\n';
    std::ifstream input(input_path);
    std::string in;
    int limit = 200000;
    int c = 0;
    while(std::getline(input, in) and c < limit) {
      if(in[0] == '>') continue;
      c++;
      uint32_t in_size = (uint32_t) in.size();
      for(uint32_t i = 0; i < in_size; ++i) {
        char c = in[i];
        if(c == 'A' or c == 'C' or c == 'G' or c == 'T')
          sketch.read_stream(c);
      }
    }
    input.close();
    std::cout << "done.\n";
  }

  double estimate_jaccard(sketch::hyperloglog &A, sketch::hyperloglog &B) {
    sketch::hyperloglog AuB(A,B);
    double est_card_A = A.estimate_cardinality();
    double est_card_B = B.estimate_cardinality();
    double est_card_AuB = AuB.estimate_cardinality();
    std::cout << "Esti: " << (uint64_t) est_card_A << ' ' << (uint64_t) est_card_B << ' ' << (uint64_t) est_card_AuB << '\n';
    return (est_card_A + est_card_B - est_card_AuB)/(est_card_AuB); 
  }
  
  double real_jaccard(sketch::hyperloglog &A, sketch::hyperloglog &B) {
    sketch::hyperloglog AuB(A,B);
    uint64_t real_card_A = A.real_cardinality();
    uint64_t real_card_B = B.real_cardinality();
    uint64_t real_card_AuB = AuB.real_cardinality();
    std::cout << "Real: " << real_card_A << ' ' << real_card_B << ' ' << real_card_AuB << '\n';
    return (1.0*(real_card_A + real_card_B - real_card_AuB))/(1.0 * real_card_AuB);
  }

  std::pair<double,double> calculate_errors(std::vector<sketch::hyperloglog> &A_i) {
    uint8_t N = A_i.size();
    std::cout << "Calculating Errors with N = " << +N << '\n';
    double erm = 0.0;
    double eam = 0.0;
    uint16_t tot = 0;
    for(uint8_t i = 0; i < N; ++i) {
      std::cout << "i = " << +(i+1) << '\n';
      for(uint8_t j = i+1; j < N; ++j) {
        double jaccard_hat = std::max(0.0, estimate_jaccard(A_i[i], A_i[j]));
        double jaccard_real = real_jaccard(A_i[i], A_i[j]);
        std::cout << "hat: " << jaccard_hat << " | real: " << jaccard_real << '\n';
        double diff = 1.0 * std::abs(jaccard_hat - jaccard_real);
        std::cout << diff << '\n';
        erm += jaccard_real == 0 ? diff : diff/jaccard_real;
        eam += diff;
        tot++;
      } 
    }
    std::cout << "done\n";
    return {(1.0/tot)*erm, (1.0/tot)*eam};
  } 

}

#endif //UTILS_HPP
