Msuite2: bin/preprocessor.pe bin/preprocessor.se bin/T2C.pe.m3 bin/T2C.pe.m4 bin/T2C.se.m3 bin/T2C.se.m4 bin/rmdup.w.pe bin/rmdup.c.pe bin/rmdup.w.se bin/rmdup.c.se bin/meth.caller.CpG bin/meth.caller.CpH bin/pair.CpG bin/pair.CpH bin/profile.DNAm.around.TSS util/bed2wig util/extract.meth.in.region
	@echo Build Msuite2 done.

cc=g++
## note that the g++ MUST support c++11 standard (i.e., version 4.8 or higher)
options=-std=c++11 -O3
multithread=-fopenmp
gzsupport=-lz

bin/preprocessor.pe: src/preprocessor.pe.cpp src/common.h src/util.h
	$(cc) $(options) $(multithread) $(gzsupport) -o bin/preprocessor.pe src/preprocessor.pe.cpp src/util.cpp

bin/preprocessor.se: src/preprocessor.se.cpp src/common.h src/util.h
	$(cc) $(options) $(multithread) $(gzsupport) -o bin/preprocessor.se src/preprocessor.se.cpp src/util.cpp

bin/T2C.pe.m3: src/T2C.pe.mode3.cpp src/common.h src/util.h src/util.cpp
	$(cc) $(options) $(multithread) -o bin/T2C.pe.m3 src/T2C.pe.mode3.cpp src/util.cpp

bin/T2C.pe.m4: src/T2C.pe.mode4.cpp src/common.h src/util.h src/util.cpp
	$(cc) $(options) $(multithread) -o bin/T2C.pe.m4 src/T2C.pe.mode4.cpp src/util.cpp

bin/T2C.se.m3: src/T2C.se.mode3.cpp src/common.h src/util.h src/util.cpp
	$(cc) $(options) $(multithread) -o bin/T2C.se.m3 src/T2C.se.mode3.cpp src/util.cpp

bin/T2C.se.m4: src/T2C.se.mode4.cpp src/common.h src/util.h src/util.cpp
	$(cc) $(options) $(multithread) -o bin/T2C.se.m4 src/T2C.se.mode4.cpp src/util.cpp

bin/rmdup.w.pe: src/rmdup.w.pe.cpp src/util.h src/util.cpp
	$(cc) $(options) -o bin/rmdup.w.pe src/rmdup.w.pe.cpp src/util.cpp

bin/rmdup.w.se: src/rmdup.w.se.cpp src/util.h src/util.cpp
	$(cc) $(options) -o bin/rmdup.w.se src/rmdup.w.se.cpp src/util.cpp

bin/rmdup.c.pe: src/rmdup.c.pe.cpp src/util.h src/util.cpp
	$(cc) $(options) -o bin/rmdup.c.pe src/rmdup.c.pe.cpp src/util.cpp

bin/rmdup.c.se: src/rmdup.c.se.cpp src/util.h src/util.cpp
	$(cc) $(options) -o bin/rmdup.c.se src/rmdup.c.se.cpp src/util.cpp

bin/meth.caller.CpG: src/meth.caller.CpG.cpp src/common.h src/util.h
	$(cc) $(options) -o bin/meth.caller.CpG src/meth.caller.CpG.cpp src/util.cpp

bin/meth.caller.CpH: src/meth.caller.CpH.cpp src/common.h src/util.h
	$(cc) $(options) -o bin/meth.caller.CpH src/meth.caller.CpH.cpp src/util.cpp

bin/pair.CpG: src/pair.CpG.cpp src/common.h src/util.h src/util.cpp
	$(cc) $(options) -o bin/pair.CpG src/pair.CpG.cpp src/util.cpp

bin/pair.CpH: src/pair.CpH.cpp src/common.h src/util.h src/util.cpp
	$(cc) $(options) -o bin/pair.CpH src/pair.CpH.cpp src/util.cpp

bin/profile.DNAm.around.TSS: src/profile.DNAm.around.TSS.cpp
	$(cc) $(options) -o bin/profile.DNAm.around.TSS src/profile.DNAm.around.TSS.cpp

util/bed2wig: util/bed2wig.cpp
	$(cc) $(options) -o util/bed2wig util/bed2wig.cpp

util/extract.meth.in.region: util/extract.meth.in.region.cpp
	$(cc) $(options) -o util/extract.meth.in.region util/extract.meth.in.region.cpp

clean:
	rm -f bin/preprocessor.pe bin/preprocessor.se bin/T2C.pe.m3 bin/T2C.pe.m4 bin/T2C.se.m3 bin/T2C.se.m4 bin/rmdup.w.pe bin/rmdup.c.pe bin/rmdup.w.se bin/rmdup.c.se bin/meth.caller.CpG bin/meth.caller.CpH bin/pair.CpG bin/pair.CpH bin/profile.DNAm.around.TSS util/bed2wig util/extract.meth.in.region

