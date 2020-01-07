run : bench.cpp
	g++ -I ~/Documents/external_package/EIGENDIR -O2 -mfma bench.cpp

.PHONY : clean
clean :
	-rm a.out
