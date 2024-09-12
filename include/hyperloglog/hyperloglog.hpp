#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

namespace hll {
  class hyperloglog {
    public:
      typedef uint32_t t_value; 
    protected:
      t_value val = 0;
    public:
      hyperloglog() {}; //default
      hyperloglog(t_value v) {
        val = v;
      };

      t_value get_val() {
        return val;
      }
                        
    private:


  };
}

#endif //HYPERLOGLOG_HPP
