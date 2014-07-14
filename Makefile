CC:=gcc
CFLAGS:= -Wall -g
.PHONY:all clean test

all: genprog

genprog : genprog.c genprog.h
	$(CC) $(CFLAGS) -o $@ $< -lm

test : genprog test.tsv
	 ./genprog --random-seed 1 test.tsv

test.tsv: 
	tr "\0" "\n" < /dev/zero | head -n 3000 | awk '{printf("%f\n",rand());}' |\
	paste - - - | awk '{printf("%s\t%s\t%s\t%f\n",$$1,$$2,$$3,2.0 *($$1 - $$3 /($$1 + $$3)));}' > $@

clean:
	rm -f genprog test.tsv

