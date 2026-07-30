CC      = cc
CFLAGS  = -O3 -Wall -Wextra
OMP_INC = -I/opt/homebrew/opt/libomp/include
OMP_LIB = -L/opt/homebrew/opt/libomp/lib -lomp
OMPFLAGS = -Xpreprocessor -fopenmp $(OMP_INC) $(OMP_LIB)

all: aak_parse stage2_filter libwildkernel.dylib

aak_parse: aak_parse.c
	$(CC) $(CFLAGS) -o aak_parse aak_parse.c

stage2_filter: stage2_filter.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -o stage2_filter stage2_filter.c

libwildkernel.dylib: wild_kernel.c
	$(CC) $(CFLAGS) -shared -fPIC -o libwildkernel.dylib wild_kernel.c

clean:
	rm -f aak_parse stage2_filter libwildkernel.dylib

.PHONY: all clean
