CXX=g++
CPPFLAGS=-I edlib/include -std=c++0x -lpthread
edlib-benchmarks: edlib-benchmarks.cpp edlib/src/edlib.cpp
	$(CXX) $^ -o $@ $(CPPFLAGS)
edlib-benchmarks-dbg: edlib-benchmarks.cpp edlib/src/edlib.cpp
	$(CXX) $^ -o $@ $(CPPFLAGS) -ggdb
clean:
	rm -f edlib-benchmarks edlib-benchmarks-dbg
