#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <string>
#include <MurmurHash1.h>
#include <cassert>
#include <bitset>
#include <cmath>

namespace hll {
  template<uint16_t K, uint32_t W = 0>
  class hyperloglog {
    public:
      typedef uint32_t t_32;
      typedef uint8_t t_8;
    protected:
      bool is_first;
      t_32 kmer_size;
      t_32 window_size;
      std::string kmer;
      std::string window;

      t_32 seed;
      t_8 p; 
      t_8 b;

      double alpha_m;
      t_32 m;
      t_8 *M;

    public:
      hyperloglog (t_8 _p, t_32 _seed = 5) {
        is_first = false;
        kmer_size = 0;
        window_size = 0;
        kmer.resize(K, 'Z');
        window.resize(W, 'Z');

        if(_p > 16) p = 16; 
        if(_p < 4) p = 4; 
        p = _p;
        m = (1 << p) - 1; //Tamaño M y es <=> (2^(p) - 1)
        b = 32 - p;
        seed = _seed;

        alpha_m = 0.673;
        if(m == 32) alpha_m = 0.697;
        else if(m == 64) alpha_m = 0.709;
        else if(m >= 128) alpha_m = 0.7213/(1 + (1.079/m));

        M = (t_8*) calloc(m, sizeof(t_8));
      }

      void read_stream (const char &c) {
        if(W == 0) return read_stream_kmers(c);
        return read_stream_minimizers(c);
      }         

      double estimate_cardinality() {
        double Z = 0;
        t_32 zero_registers = 0;
        for(t_32 i = 0; i < m; ++i) {
          if(M[i] == 0) zero_registers++;
          t_32 val = (1 << M[i]);
          Z += 1.0/val;
        }
        double E = alpha_m*m*m*(1.0/Z); 
        double E_star = E;

        long pow2_32 = 1LL<<32;

        if(E <= (5/2)*m) {
          if(zero_registers != 0)
            E_star = m*log2((m/zero_registers));
        }
        else if(E > pow2_32) {
          E_star = -(pow2_32 * log2(1 - (E/pow2_32)));
        }        
        return E_star;
      }

    private:
      void compute_hyperloglog () {
        if(kmer_size < K or window_size < W) return;
        t_32 hash_kmer = MurmurHash1(kmer.c_str(), K, seed);
        t_32 index = hash_kmer >> (32 - p);
        t_32 value_part = hash_kmer << p;
        t_8 value = __builtin_clz(value_part) + 1;
        //std::cout << +M[index] << ' ' << +value << '\n';
        M[index] = std::max(M[index] , value);
      }

      void read_stream_kmers (const char &c) {
        if(kmer_size >= K) 
          kmer = kmer.substr(1) + c;  
        else
          kmer[kmer_size++] = c;
        compute_hyperloglog();
      }

      void read_stream_minimizers (const char &c) {
        if(window_size >= W)
          window = window.substr(1) + c;
        else
          window[window_size++] = c;
        get_minimizer();
        compute_hyperloglog();
      }

      void get_minimizer () {
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
