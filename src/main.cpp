#include<iostream>
#include<hyperloglog/hyperloglog.hpp>

int main () {
  hll::hyperloglog a(76);
  std::cout << a.get_val() << '\n';

  hll::hyperloglog b;
  std::cout << b.get_val() << '\n';
  return 0;
}
