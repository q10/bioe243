CC=g++
#CCFLAGS=-O2 -msse2 -fomit-frame-pointer -pipe -Wno-deprecated $(CPPFLAGS)
CCFLAGS=-O3 -march=native -Wno-deprecated $(CPPFLAGS) -I. -I/usr/lib/python2.5/site-packages/shedskin/lib -g -fPIC -D__SS_BIND -I/usr/include/python2.5 -I/usr/include/python2.5
LFLAGS=-lgc -lpcre $(LDFLAGS) -shared -Xlinker -export-dynamic -lpthread -ldl  -lutil -lm -lpython2.5 -L/usr/lib/python2.5/config

all:	random_seq.so

CPPFILES=random_seq.cpp /usr/lib/python2.5/site-packages/shedskin/lib/bisect.cpp /usr/lib/python2.5/site-packages/shedskin/lib/random.cpp /usr/lib/python2.5/site-packages/shedskin/lib/builtin.cpp /usr/lib/python2.5/site-packages/shedskin/lib/time.cpp /usr/lib/python2.5/site-packages/shedskin/lib/math.cpp /usr/lib/python2.5/site-packages/shedskin/lib/re.cpp
HPPFILES=random_seq.hpp /usr/lib/python2.5/site-packages/shedskin/lib/bisect.hpp /usr/lib/python2.5/site-packages/shedskin/lib/random.hpp /usr/lib/python2.5/site-packages/shedskin/lib/builtin.hpp /usr/lib/python2.5/site-packages/shedskin/lib/time.hpp /usr/lib/python2.5/site-packages/shedskin/lib/math.hpp /usr/lib/python2.5/site-packages/shedskin/lib/re.hpp

random_seq.so:	$(CPPFILES) $(HPPFILES)
	$(CC)  $(CCFLAGS) $(CPPFILES) $(LFLAGS) -o random_seq.so

random_seq.so_prof:	$(CPPFILES) $(HPPFILES)
	$(CC) -pg -ggdb $(CCFLAGS) $(CPPFILES) $(LFLAGS) -o random_seq.so_prof

random_seq.so_debug:	$(CPPFILES) $(HPPFILES)
	$(CC) -g -ggdb $(CCFLAGS) $(CPPFILES) $(LFLAGS) -o random_seq.so_debug

clean:
	rm -f random_seq.so random_seq.so_prof random_seq.so_debug

.PHONY: all run clean

