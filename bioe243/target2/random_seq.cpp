#include "random_seq.hpp"

namespace __random_seq__ {

str *const_0;

list<__ss_int> *orig_seq, *seen_seq;
__ss_int __4, __5, i, num;
str *__name__;


__ss_int get_seq_int(list<__ss_int> *seq) {
    __ss_int __0, __1, i, seq_int;

    seq_int = 0;

    FAST_FOR(i,0,30,1,0,1)
        if (orig_seq->__getfast__(i)) {
            seq_int = (seq_int+(1<<i));
        }
    END_FOR

    return seq_int;
}

__ss_bool check_in(__ss_int i, __ss_int key) {
    __ss_bool __2, __3;

    if (((i!=len(seen_seq)) and (seen_seq->__getfast__(i)==key))) {
        return True;
    }
    return False;
}

list<__ss_int> *fastest_create_random_sequence() {
    
    __random__::shuffle(orig_seq);
    return orig_seq;
}

void __init() {
    const_0 = new str("__main__");

    __name__ = new str("random_seq");

    orig_seq = (((new list<__ss_int>(1, 0)))->__mul__(15))->__add__(((new list<__ss_int>(1, 1)))->__mul__(15));
    seen_seq = (new list<__ss_int>());
    if (__eq(__name__, const_0)) {
        num = 5000;

        FAST_FOR(i,0,num,1,4,5)
            fastest_create_random_sequence();
        END_FOR

    }
}

} // module namespace

/* extension module glue */

extern "C" {
#include <Python.h>
#include <structmember.h>

PyObject *__ss_mod_random_seq;

namespace __random_seq__ { /* XXX */
PyObject *Global_random_seq_get_seq_int(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        list<__ss_int> *arg_0 = __ss_arg<list<__ss_int> *>("seq", 0, 0, NULL, args, kwargs);

        return __to_py(__random_seq__::get_seq_int(arg_0));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->msg)?(e->msg->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_random_seq_check_in(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        __ss_int arg_0 = __ss_arg<__ss_int >("i", 0, 0, NULL, args, kwargs);
        __ss_int arg_1 = __ss_arg<__ss_int >("key", 1, 0, NULL, args, kwargs);

        return __to_py(__random_seq__::check_in(arg_0, arg_1));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->msg)?(e->msg->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_random_seq_fastest_create_random_sequence(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {

        return __to_py(__random_seq__::fastest_create_random_sequence());

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->msg)?(e->msg->unit.c_str()):""));
        return 0;
    }
}

static PyNumberMethods Global_random_seq_as_number = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
};

static PyMethodDef Global_random_seqMethods[] = {
    {(char *)"__newobj__", (PyCFunction)__ss__newobj__, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"get_seq_int", (PyCFunction)Global_random_seq_get_seq_int, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"check_in", (PyCFunction)Global_random_seq_check_in, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"fastest_create_random_sequence", (PyCFunction)Global_random_seq_fastest_create_random_sequence, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {NULL}
};

PyMODINIT_FUNC initrandom_seq(void) {
    __shedskin__::__init();
    __math__::__init();
    __time__::__init();
    __random__::__init();
    __bisect__::__init();
    __random_seq__::__init();

    __ss_mod_random_seq = Py_InitModule((char *)"random_seq", Global_random_seqMethods);
    if(!__ss_mod_random_seq)
        return;


    addrandom_seq();

}

PyMODINIT_FUNC addrandom_seq(void) {
    PyModule_AddObject(__ss_mod_random_seq, (char *)"num", __to_py(__random_seq__::num));
    PyModule_AddObject(__ss_mod_random_seq, (char *)"i", __to_py(__random_seq__::i));
    PyModule_AddObject(__ss_mod_random_seq, (char *)"orig_seq", __to_py(__random_seq__::orig_seq));
    PyModule_AddObject(__ss_mod_random_seq, (char *)"seen_seq", __to_py(__random_seq__::seen_seq));

}

} // namespace __random_seq__

} // extern "C"
int main(int, char **) {
    __shedskin__::__init();
    __math__::__init();
    __time__::__init();
    __random__::__init();
    __bisect__::__init();
    __shedskin__::__start(__random_seq__::__init);
}
