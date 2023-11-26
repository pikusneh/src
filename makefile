CXX = g++
CXXFLAGS = $(shell root-config --cflags) -O2 -Wall -fPIC
LIBS = $(shell root-config --libs)

main.o: main.cpp
	$(CXX) -c main.cpp $(CXXFLAGS)

# Link
GenRecoAnalyzer: main.o GenRecoAnalyzer.o
	$(CXX) -o GenRecoAnalyzer main.o GenRecoAnalyzer.o $(LIBS)

# CXX = g++
# CXXFLAGS = $(shell root-config --cflags) -O2 -Wall -fPIC
# LIBS = $(shell root-config --libs)

# GenRecoAnalyzer: main.o GenRecoAnalyzer.o
# 	$(CXX) -o GenRecoAnalyzer main.o GenRecoAnalyzer.o $(LIBS)

# main.o: main.cpp
# 	$(CXX) -c main.cpp $(CXXFLAGS)

# GenRecoAnalyzer.o: GenRecoAnalyzer.cpp GenRecoAnalyzer.h
# 	$(CXX) -c GenRecoAnalyzer.cpp $(CXXFLAGS)

# clean:
# 	rm -f *.o GenRecoAnalyzer
