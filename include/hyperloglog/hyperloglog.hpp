#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <string>
#include <MurmurHash2.h>
#include <cassert>
#include <set>
#include <vector>
#include <bitset>
#include <cmath>

namespace sketch {
  class hyperloglog {
    public:
      typedef uint64_t t_64;
      typedef uint32_t t_32;
      typedef uint16_t t_16;
      typedef uint8_t t_8;
    protected:
      t_16 K;
      t_32 W;

      t_32 kmer_size;
      t_32 window_size;
      std::string kmer;
      std::string window;

      std::set<std::string> real;

      t_32 seed;
      t_8 p; 
      t_8 b;

      double alpha_m;
      t_32 m;
      std::vector<t_8> M;

    public:
      //constructor para la union de dos hyperloglog
      hyperloglog (hyperloglog const& A, hyperloglog const& B) { 
        if(A.K != B.K) {
          std::cout << "error, distinct K\n";
          return;
        }
        if(A.W != B.W) {
          std::cout << "error, distinct W\n";
          return;
        }
        if(A.p != B.p) {
          std::cout << "error, distinct p\n";
          return;
        }
        if(A.seed != B.seed) {
          std::cout << "error, distinct seed\n";
          return;
        }
        real.clear();
        p = A.p;
        b = A.b;
        seed = A.seed;
        m = A.m;
        alpha_m = A.alpha_m;
        K = A.K;
        W = A.W;
        real.insert(A.real.begin(), A.real.end());
        real.insert(B.real.begin(), B.real.end());
        M.resize(m, 0);
        for(t_32 i = 0; i < m; ++i) {
          M[i] = std::max(A.M[i], B.M[i]);
          assert(M[i] >= A.M[i]);
          assert(M[i] >= B.M[i]);
        }
      }

      //constructor para 1 hyperloglog
      hyperloglog (t_8 _p, t_16 _K, t_32 _W = 0, t_32 _seed = 5) {
        K = _K;
        W = _W;
        kmer_size = 0;
        window_size = 0;
        kmer.resize(K, 'Z');
        window.resize(W, 'Z');

        real.clear();

        if(_p > 16) p = 16; 
        if(_p < 4) p = 4; 
        p = _p;
        m = 1 << p; //Tamaño M y es <=> (2^(p))
        b = 32 - p;
        seed = _seed;

        alpha_m = 0.673;
        if(m == 32) alpha_m = 0.697;
        else if(m == 64) alpha_m = 0.709;
        else if(m >= 128) alpha_m = 0.7213/(1 + (1.079/m));
        M.resize(m, 0);
      }

      t_64 real_size() {
        return (t_64) real.size() * (t_64) sizeof(*real.begin());
      }

      t_64 estimated_size() {
        return m*(t_64) sizeof(t_8);
      }

      void read_stream (const char &c) {
        if(W == 0) return read_stream_kmers(c);
        return read_stream_minimizers(c);
      }         

      t_64 real_cardinality() {
        return (uint64_t) real.size();
      }

      double estimate_cardinality() {
        double Z = 0.0;
        t_32 zero_registers = 0;
        for(t_32 i = 0; i < m; ++i) {
          if(M[i] == 0) zero_registers++;
          double val = (1 << M[i]);
          Z += 1.0/val;
        }
        double E = alpha_m*m*m*(1.0/Z); 
        double E_star = E;

        t_64 pow2_32 = 1LL<<32;

        if(E <= (5.0/2.0)*m) {
          if(zero_registers != 0) {
            double log = log2(((1.0*m)/zero_registers));
            E_star = m*log;
          }
        }
        else if(E > (1.0/30.0)*pow2_32) {
          double log = log2((1.0 - E)/(1.0*pow2_32));
          E_star = -pow2_32 * log;
        }        
        return E_star;
      }

    private: 
      //11010100100100010010001100110101
      //p = 14, b = 32 - 14 = 18
      //index = 11010100100100 => M[13604] |  value = 010010001100110101 ==> M[13604] = 2
      //11010100100100000000000000000000
      //p = 14, b = 18
      //index = 11010100100100 => M[13604] |  value = 000000000000000000 ==>  M[13604] = 19
      void compute_hyperloglog () {
        if(kmer_size < K or window_size < W) return;
        real.insert(kmer);
        t_32 hash_kmer = MurmurHash2(kmer.c_str(), K, seed);
        t_32 index = hash_kmer >> b;  
        t_32 value_part = hash_kmer << p; 
        t_8 value = __builtin_clz(value_part) + 1;
        if (value > b) value = b+1;
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
