/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef GENPROG_H
#define GENPROG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include <getopt.h>

typedef double floating_t;
typedef int boolean_t;
enum nodeType {CONSTANT,OPERATOR,COLUMN};

#define MIN(a,b) (a<b?a:b)

/** SpreadSheet */
typedef struct spreadsheet_t
	{
	/** number of columns */
	size_t columns;
	/** number of cells */
	size_t size;
	/** all cells */
	floating_t* data;
	/** normalized data for last column */
	floating_t* normalized;
	}SpreadSheet,*SpreadSheetPtr;


/**
 * Operator
 */
typedef struct operator_t
	{
	char name[100];
	size_t num_children;
	floating_t (*eval)(floating_t* );
	int weight;
	size_t index;
	} Operator,*OperatorPtr;


/**
 * OperatorList
 */
typedef struct operator_list_t
	{
	Operator** operators;
	size_t size;
	int total_weight;
	/* below, shortcuts to operator to create silent mutations e.g x/1.0 ; x*1.0   x+0.0; y-0.0 */
	OperatorPtr plus;
	OperatorPtr minus;
	OperatorPtr div;
	OperatorPtr mul;
	} OperatorList,*OperatorListPtr;



/** Configuration */
typedef struct config_t
	{
	SpreadSheetPtr spreadsheet;
	OperatorListPtr operators;
	long max_generations;
	long curr_generations;
	int min_genomes_per_generation;
	int max_genomes_per_generation;
	int min_base_per_genome;
	int max_base_per_genome;
	unsigned int seedp;
	int weight_constant;
	int weight_op;
	int weight_column;
	float max_fraction_of_errors;
	float probability_mutation;
	boolean_t remove_introns;
	boolean_t best_will_survive;
	boolean_t enable_self_self;
	boolean_t normalize_data;
	long massive_extinction_every;
	boolean_t sort_on_genome_size;
	boolean_t remove_clone;
	} Config,*ConfigPtr;
	
#define RANDOM_FLOAT(cfg) ((double)rand_r(&(cfg->seedp))/(double)RAND_MAX)
#define RANDOM_SIZE_T(cfg,LEN)  (rand_r(&(cfg->seedp))% (LEN) )

/** Node */
typedef struct node_t
	{
	union
		{
		floating_t constant;
		size_t column;
		size_t operator;
		} core;
	enum nodeType type;
	} Node,*NodePtr;	
	





/**
 * A Genome
 */
typedef struct genome_t
	{
	/** associated config */
	ConfigPtr config;
	/** all nodes */
	NodePtr nodes;
	/* number of nodes */
	size_t node_count;
	/** fitness */
	floating_t fitness;
	/** valid ? */
	boolean_t bad_flag;
	/** associated generation */
	long  generation;
	}Genome,*GenomePtr;

/*
 * Generation
 */
typedef struct generation_t
	{
	/* generation count */
	long generation;
	/* genomes */ 
	Genome** genomes;
	/* number of genomes */
	size_t genome_count;
	} Generation,*GenerationPtr;

SpreadSheetPtr SpreadSheetRead(FILE* in);
size_t SpreadSheetColumns(const SpreadSheetPtr ptr);
size_t SpreadSheetRows(const SpreadSheetPtr ptr);
floating_t SpreadSheetAt(const SpreadSheetPtr ptr,size_t y,size_t x);


void GenomeFree(GenomePtr ptr);
GenomePtr GenomeNew(ConfigPtr cfg);
NodePtr GenomeAt(GenomePtr g,size_t i);
size_t GenomeSize(const GenomePtr g);
GenomePtr GenomeNew1(ConfigPtr cfg);
GenomePtr GenomeClone(const GenomePtr src);
int GenomeCompare(const GenomePtr g1,const GenomePtr g2);
void GenomePrint(const GenomePtr g,FILE* out);
void GenomeMute(GenomePtr cfg);
boolean_t GenomeEquals(const GenomePtr g1,const GenomePtr g2);

OperatorListPtr OperatorsListNew();
OperatorPtr OperatorListAt(OperatorListPtr list,size_t index);
size_t OperatorListSize(OperatorListPtr list);
#endif

