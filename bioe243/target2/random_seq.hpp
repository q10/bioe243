#ifndef __RANDOM_SEQ_HPP
#define __RANDOM_SEQ_HPP

#include "math.hpp"
#include "random.hpp"
#include "random_seq.hpp"
#include "bisect.hpp"
#include "builtin.hpp"
#include "time.hpp"

using namespace __shedskin__;
namespace __random_seq__ {

extern str *const_0;

using __bisect__::bisect_left;



extern list<__ss_int> *orig_seq, *seen_seq;
extern __ss_int __4, __5, i, num;
extern str *__name__;

__ss_int get_seq_int(list<__ss_int> *seq);
__ss_bool check_in(__ss_int i, __ss_int key);
list<__ss_int> *fastest_create_random_sequence();

PyMODINIT_FUNC initrandom_seq(void);

PyMODINIT_FUNC addrandom_seq(void);

} // module namespace
namespace __shedskin__ { /* XXX */

}
#endif
