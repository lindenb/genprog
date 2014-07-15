CC:=gcc
CFLAGS:= -Wall -O3
.PHONY:all clean test

all: genprog

genprog : genprog.c genprog.h
	$(CC) $(CFLAGS) -o $@ $< -lm

test : genprog test.tsv
	 ./genprog --min-bases 1 --max-bases 30  test.tsv

test.tsv: 
	tr "\0" "\n" < /dev/zero | head -n 3000 | awk '{printf("%f\n",rand());}' |\
	paste - - - | awk '{printf("%s\t%s\t%s\t%f\n",$$1,$$2,$$3,2.0 *(($$1 - $$2) /($$1 + $$3)));}' > $@

clean:
	rm -f genprog test.tsv

