#include<iostream>
#include<hyperloglog/hyperloglog.hpp>

int main () {
  hll::hyperloglog<3> a;
  std::string aa = "ACGTGACCG";

  for(char c: aa)
    a.read_stream(c);

  std::cout << '\n';
  hll::hyperloglog<3, 6> b;
  for(char c: aa)
    b.read_stream(c);
  return 0;
}
