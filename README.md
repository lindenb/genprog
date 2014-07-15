genprog
=======

Genetic programming implementation in C 

Author: Pierre Lindenbaum PhD @yokofakun

> In artificial intelligence, genetic programming (GP) is an evolutionary algorithm-based methodology inspired by biological evolution to find computer programs that perform a user-defined task.
> Essentially GP is a set of instructions and a fitness function to measure how well a computer has performed a task.

See https://en.wikipedia.org/wiki/Genetic_programming

## Compilation

```
$ make
```

## Usage

```
$ genprog [options] (stdin|file)
```

## Input

* Input is a table delimited file. It contains only numeric values.
* It must contains at least 2 columns.
* It must contains at least 2 rows.
* All columns but last one are the observed values.
* The last column is the expected result.

### Example:
```tsv
$ head test.tsv 
0.066221	0.984740	0.526733	-3.098112
0.804358	0.838329	0.799516	-0.042361
0.457574	0.442436	0.014716	0.064105
0.324609	0.703579	0.053475	-2.004687
0.752822	0.681952	0.559829	0.107980
0.040027	0.728014	0.724929	-1.798762
0.881736	0.343190	0.991670	0.574938
0.005937	0.787656	0.133316	-11.227320
0.644971	0.032298	0.838944	0.825752
0.133133	0.559989	0.730044	-0.989035
```

## Options

TODO

## Example


```bash
tr "\0" "\n" < /dev/zero | head -n 3000 | awk '{printf("%f\n",rand());}' |\
	paste - - - | awk '{printf("%s\t%s\t%s\t%f\n",$1,$2,$3,2.0 *(($1 - $2) /($1 + $3)));}' > test.tsv
./genprog --min-bases 3 --max-bases 20 --min-genomes 3 --max-genomes 50 test.tsv
```
output
```
FITNESS=2.078584E+03	GENERATION=0	SECONDS=0	(Add(-0.599715,${1}));
FITNESS=1.927753E+03	GENERATION=1	SECONDS=0	(Minus(${1},Minus(Invert(0.544634),${1})));
FITNESS=1.759247E+03	GENERATION=2	SECONDS=1	(Minus(${1},${2}));
FITNESS=1.473560E+03	GENERATION=7	SECONDS=1	(Minus(Add(Minus(Negate(0.274741),-0.661806),Minus(${1},${2})),${2}));
FITNESS=1.149314E+03	GENERATION=21	SECONDS=4	(Mul(Invert(Minus(Negate(0.274741),-0.661806)),Minus(${1},${2})));
FITNESS=8.915954E+02	GENERATION=47	SECONDS=8	(Mul(Invert(Minus(Negate(${1}),${3})),${2}));
FITNESS=8.241330E+02	GENERATION=54	SECONDS=9	(Mul(Invert(Minus(Negate(${1}),${3})),Div(${2},0.899698)));
FITNESS=7.760294E+02	GENERATION=56	SECONDS=10	(Mul(Invert(Minus(Negate(${1}),${3})),Invert(Mul(Invert(${2}),0.682955))));
FITNESS=7.641078E+02	GENERATION=64	SECONDS=11	(Mul(Invert(Minus(Negate(${1}),${3})),Invert(Mul(Invert(${2}),0.734345))));
FITNESS=6.547658E+02	GENERATION=149	SECONDS=28	(Mul(Invert(Minus(Negate(${1}),${3})),Minus(${2},${1})));
FITNESS=2.959147E+02	GENERATION=153	SECONDS=28	(Mul(Invert(Minus(Negate(${1}),${3})),Minus(${2},Minus(${1},${2}))));
FITNESS=8.254130E-11	GENERATION=220	SECONDS=42	(Mul(Invert(Minus(Negate(${1}),${3})),Minus(Minus(${2},Minus(${1},${2})),${1})));
#min fitness reached
```
