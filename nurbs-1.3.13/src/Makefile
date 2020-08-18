OCTFILES=basisfun.oct \
basisfunder.oct \
bspderiv.oct \
bspeval.oct \
curvederivcpts.oct \
nrb_srf_basisfun__.oct \
nrb_srf_basisfun_der__.oct \
nrbsurfderiveval.oct \
surfderivcpts.oct \
surfderiveval.oct \
tbasisfun.oct 

MKOCTFILE ?= mkoctfile

all: $(OCTFILES)

low_level_functions.o: low_level_functions.cc
	$(MKOCTFILE) -c $<

%.oct:  %.cc low_level_functions.o
	$(MKOCTFILE) $< low_level_functions.o

clean:
	-rm -f *.o core octave-core *.oct *~

