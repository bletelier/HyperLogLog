#include<iostream>
#include<hyperloglog/hyperloglog.hpp>
#include<utils/utils.hpp>

uint8_t const N = 5;

int main (int argc, char* argv[]) {
  uint16_t K;
  uint32_t W;
  uint8_t p = 14;
  uint32_t seed = 15;

  if(argc < 2) {
    std::cout << "Error... execute as: ./build/main <K (kmer size)> <W (window size, default 0)>";
    return 1;
  }
  if(argc == 2) W = 0;
  else W = atoi(argv[2]);
  K = atoi(argv[1]);
  std::cout << "Executing with p = " << +p << ", K = " << K << " and W = " << W << '\n'; 

  std::vector<sketch::hyperloglog> A_i;

  for(uint8_t i = 0; i < N; ++i) {
    sketch::hyperloglog A(p, K, W, seed);
    utils::parse_input(A,"input/gen"+std::to_string(i+1)+".fna");
    A_i.push_back(A);
  }
 
  std::pair<double,double> errors = utils::calculate_errors(A_i);
  std::cout << "ERM = " << errors.first << '\n';
  std::cout << "EAM = " << errors.second << '\n';
  return 0;
}
