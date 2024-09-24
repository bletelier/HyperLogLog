#include<iostream>
#include<hyperloglog/hyperloglog.hpp>
#include<utils/utils.hpp>

int main () {
  hll::hyperloglog<25> a (16);
  std::string aa = "ACGTGACCG";

  std::string gen = utils::parse_input("input/gen1.fna");
  for(char c: gen)
    a.read_stream(c);
  std::cout << "cardinality = " << a.estimate_cardinality() << '\n';
  hll::hyperloglog<25, 100> b (16);
  for(char c: gen)
    b.read_stream(c);
  std::cout << "cardinality = " << b.estimate_cardinality() << '\n';
  return 0;
}
