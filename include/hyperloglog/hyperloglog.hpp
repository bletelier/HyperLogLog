#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <string>

namespace hll {
  template<uint16_t K, uint32_t W = 0>
  class hyperloglog {
    public:
      typedef uint32_t t_value;
    protected:
      bool is_first;
      t_value kmer_size;
      t_value window_size;
      std::string kmer;
      std::string window;
    public:
      hyperloglog() {
        is_first = false;
        kmer_size = 0;
        window_size = 0;
        kmer.resize(K, 'Z');
        window.resize(W, 'Z');
      };

      void read_stream(const char &c) {
        if(W == 0) return read_stream_kmers(c);
        return read_stream_minimizers(c);
      }         
    private:
      void compute_hyperloglog() {
        if(kmer_size < K or window_size < W) return;
        std::cout << kmer << '\n';
      }

      void read_stream_kmers(const char &c) {
        if(kmer_size >= K) 
          kmer = kmer.substr(1) + c;  
        else
          kmer[kmer_size++] = c;
        compute_hyperloglog();
      }

      void read_stream_minimizers(const char &c) {
        if(window_size >= W)
          window = window.substr(1) + c;
        else
          window[window_size++] = c;
        get_minimizer();
        compute_hyperloglog();
      }

      void get_minimizer() {
        if(window_size < W) return;
        std::string aux1(K, 'Z');
        std::string aux2 = aux1;
        kmer_size = 0;
        for(char const& c: window) {
          if(kmer_size >= K) {
            aux2 = aux1.compare(aux2) <= 0 ? aux1 : aux2; 
            aux1 = aux1.substr(1) + c;
          }
          else
            aux1[kmer_size++] = c;

        }
        aux2 = aux1.compare(aux2) <= 0 ? aux1 : aux2; 
        kmer = aux2;
      }
  };
}

#endif //HYPERLOGLOG_HPP
