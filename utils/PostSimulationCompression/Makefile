CC = gcc
CFLAGS = -O3 -march=native -ffast-math -funroll-loops -std=c11
LDFLAGS = -lm

all: encode_trj decode_trj

encode_trj: encode_trj.c
	$(CC) $(CFLAGS) -o encode_trj encode_trj.c

decode_trj: decode_trj.c
	$(CC) $(CFLAGS) -o decode_trj decode_trj.c $(LDFLAGS)

clean:
	rm -f encode_trj decode_trj