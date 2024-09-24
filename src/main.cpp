#include<iostream>
#include<hyperloglog/hyperloglog.hpp>
#include<utils/utils.hpp>

const uint16_t K = 20;
int main () {
  hll::hyperloglog<K> a (14);
  //std::string gen2 = "ACGTGACCCG";
  //std::string gen4 = "ATCGGATCCCCCTA";

  hll::hyperloglog<K> aa (14);
  std::string gen = utils::parse_input("input/gen1.fna");
  for(char c: gen)
    a.read_stream(c);
  std::cout << "cardinality = " << a.estimate_cardinality() << '\n';
  std::cout << "real cardinality = " << a.real_cardinality() << '\n';

  std::string gen3 = utils::parse_input("input/gen2.fna");
  for(char c: gen3)
    aa.read_stream(c);
  std::cout << "cardinality = " << aa.estimate_cardinality() << '\n';
  std::cout << "real cardinality = " << aa.real_cardinality() << '\n';

  hll::hyperloglog<K> auaa(a, aa);
  std::cout << "cardinality = " << auaa.estimate_cardinality() << '\n';
  std::cout << "real cardinality = " << auaa.real_cardinality() << '\n';
  /*hll::hyperloglog<25> auaa(a, aa);
  std::cout << "cardinality = " << auaa.estimate_cardinality() << '\n';
  hll::hyperloglog<25, 100> b (14);
  for(char c: gen)
    b.read_stream(c);
  std::cout << "cardinality = " << b.estimate_cardinality() << '\n';*/
  return 0;
}
