CC      = cc
CFLAGS  = -O3 -Wall -Wextra

all: aak_parse libwildkernel.dylib

aak_parse: aak_parse.c
	$(CC) $(CFLAGS) -o aak_parse aak_parse.c

libwildkernel.dylib: wild_kernel.c
	$(CC) $(CFLAGS) -shared -fPIC -o libwildkernel.dylib wild_kernel.c

clean:
	rm -f aak_parse libwildkernel.dylib

.PHONY: all clean
