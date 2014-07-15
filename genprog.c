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

#include <unistd.h>
#include <errno.h>
#include "genprog.h"


#define DEBUG(msg) do { fprintf(stderr,"[%d]DEBUG ",__LINE__); fputs(msg,stderr);fputc('\n',stderr); } while(0)
#define THROW_ERROR(a) do {fputs(a,stderr);fputc('\n',stderr);assert(0);} while(0)



static const size_t NPOS=(size_t)-1UL;


floating_t SpreadSheetAt(const SpreadSheetPtr ptr,size_t y,size_t x)
	{
	if( y >= SpreadSheetRows(ptr) )
		{
		fprintf(stderr,"Y=%d >= rows=(%d)\n",y, SpreadSheetRows(ptr) );
		THROW_ERROR("y>rows");
		}
	if( x >= SpreadSheetColumns(ptr) )
		{
		fprintf(stderr,"X=%d >= cols=(%d)\n",x, SpreadSheetColumns(ptr) );
		THROW_ERROR("x>=col");
		}
	return ptr->data[ y*SpreadSheetColumns(ptr) + x];
	}




static void NodeInitConstant(ConfigPtr cfg,NodePtr op)
	{
	op->type = CONSTANT;
	if(RANDOM_FLOAT(cfg)<0.5)
		{
		op->core.constant = RANDOM_FLOAT(cfg);
		}
	else
		{
		op->core.constant = (double)rand_r(&(cfg->seedp));
		}
	if(RANDOM_FLOAT(cfg)<0.5) op->core.constant *=-1;
	}
	
static void NodeInitColumn(ConfigPtr cfg,NodePtr op)
	{
	op->type = OPERATOR;
	op->core.operator= RANDOM_SIZE_T(cfg, OperatorListSize(cfg->operators) );
	}
	
static void NodeInitOperator(ConfigPtr cfg,NodePtr op)
	{
	op->type = COLUMN;
	op->core.column = RANDOM_SIZE_T(cfg, SpreadSheetColumns(cfg->spreadsheet)-1);

	}

static void NodeMut(ConfigPtr cfg,NodePtr op)
	{
	switch(op->type)
		{
		case CONSTANT: NodeInitConstant(cfg,op);break;
		case OPERATOR : NodeInitOperator(cfg,op);break;
		case COLUMN: NodeInitColumn(cfg,op);break;
		default: THROW_ERROR("BOUM");
		}
	}

static void NodeInit(ConfigPtr cfg,NodePtr op)
	{
	float r= RANDOM_FLOAT(cfg);

	//int weight = rand_r(&(cfg->seedp)) % cfg->weight_constant + cfg->weight_column + cfg->weight_op ;
	
	if( r < 0.333f)
		{
		NodeInitConstant(cfg,op);
		}
	else if(  r < 0.666f)
		{
		NodeInitOperator(cfg,op);
		}
	else
		{
		NodeInitColumn(cfg,op);
		}

	}

GenomePtr GenomeNew1(ConfigPtr cfg)
	{
	GenomePtr g=(GenomePtr)calloc(1,sizeof(Genome));
	if(g==NULL) THROW_ERROR("boum");

	g->config = cfg;
	g->node_count = 0;
	g->nodes = NULL;
	g->bad_flag=0;
	g->fitness = NAN;
	g->generation = -1L;
	g->creation = time(NULL);
	return g;
	}

GenomePtr GenomeClone(const GenomePtr src)
	{
	GenomePtr g=(GenomePtr)calloc(1,sizeof(Genome));
	if(g==NULL) THROW_ERROR("boum");
	g->config = src->config;
	g->node_count = src->node_count;
	g->nodes = (NodePtr)calloc(src->node_count,sizeof(Node));
	if(g->nodes==NULL) THROW_ERROR("boum");
	memcpy((void*)g->nodes,(void*)src->nodes,sizeof(Node)*(g->node_count));
	g->bad_flag = src->bad_flag;
	g->fitness = src->fitness;
	g->generation = src ->generation;
	g->creation = src->creation; 
	return g;
	}	
	


size_t GenomeSize(const GenomePtr g)
	{
	return g->node_count;
	}

NodePtr GenomeAt(GenomePtr g,size_t i)
	{
	if( i>=GenomeSize(g) ) THROW_ERROR("boum");
	return &(g->nodes[i]);
	}


GenomePtr GenomeNew(ConfigPtr cfg)
	{
	int i;
	GenomePtr g = GenomeNew1(cfg);

	g->node_count = cfg->min_base_per_genome +
		RANDOM_SIZE_T(cfg, abs(cfg->max_base_per_genome - cfg->min_base_per_genome)) ;
		
	g->nodes = (NodePtr)calloc(g->node_count,sizeof(Node));
	if(g->nodes==NULL) THROW_ERROR("boum");
	for(i=0;i <  GenomeSize(g) ; ++i)
		{
		NodeInit(cfg,&g->nodes[i]);
		}
	
	return g;
	}

void GenomeFree(GenomePtr ptr)
	{
	if(ptr==NULL) return;
	free( ptr-> nodes );
	free(ptr);
	}

int GenomeCompare(const GenomePtr g1,const GenomePtr g2)
	{
	if(g1 -> bad_flag)
		{
		if(g2->bad_flag)
			{
			if( g1->config->sort_on_genome_size)
				{
				if(g1->node_count < g2->node_count) return -1;
				if(g1->node_count > g2->node_count) return 1;
				}
			return 0;
			}
		return 1;
		}
	if(g2->bad_flag ) return -1;
	if( g1->fitness < g2->fitness) return -1;
	if( g1->fitness > g2->fitness) return 1;
	
	if( g1->config->sort_on_genome_size)
		{
		if(g1->node_count < g2->node_count) return -1;
		if(g1->node_count > g2->node_count) return 1;
		}
	return 0;
	}
static int _GenomeCompare(const void* a,const void* b)
	{
	const GenomePtr g1=(const GenomePtr)(*((const Genome**)a));
	const GenomePtr g2=(const GenomePtr)(*((const Genome**)b));
	return GenomeCompare(g1,g2);
	}



static size_t GenerationCount(const GenerationPtr g)
	{
	return g->genome_count;
	}

static GenomePtr GenerationAt(GenerationPtr g,size_t i)
	{
	if( i>= GenerationCount(g)) THROW_ERROR("boum");
	return g->genomes[i];
	}

static GenerationPtr GenerationNew1(ConfigPtr cfg)
	{
	GenerationPtr g=(GenerationPtr)calloc(1,sizeof(Generation));
	if(g==NULL) THROW_ERROR("boum");
	g->generation = cfg->curr_generations;
	g->genome_count = 0;
	g->genomes = NULL;
	return g;
	}

static GenerationPtr GenerationNew(ConfigPtr cfg)
	{
	int i;
	GenerationPtr g = GenerationNew1(cfg); 
	g->genome_count =
		cfg->min_genomes_per_generation +
		RANDOM_SIZE_T(cfg, abs( cfg->max_genomes_per_generation - cfg->min_genomes_per_generation ))
		;
	g->genomes=(Genome**)calloc(g->genome_count,sizeof(Genome*));
	if(g->genomes==NULL) THROW_ERROR("boum");
	
	for(i=0;i< GenerationCount(g) ;++i)
		{
		g->genomes[i] = GenomeNew(cfg);
		}
	
	return g;
	}

void GenerationAdd(GenerationPtr g,GenomePtr genome)
	{
	g->genomes=(Genome**)realloc(
		g->genomes,
		(g->genome_count+1)*sizeof(GenomePtr)
		);
	if(g->genomes==NULL) THROW_ERROR("boum");
	g->genomes[g->genome_count] = genome;
	g->genome_count++;
	}


static void GenerationFree(GenerationPtr g)
	{
	size_t i;
	if(g==NULL) return;
	for(i=0;i< GenerationCount(g) ;++i)
		{
		GenomeFree( g->genomes[i] );
		}
	free(g->genomes);
	free(g);
	}



#define MAX_ARITY 2

static floating_t _plus(floating_t* array) { return array[0] + array[1];}
static floating_t _minus(floating_t* array) { return array[0] - array[1];}
static floating_t _mul(floating_t* array) { return array[0] * array[1];}
static floating_t _div(floating_t* array) { return array[1]==0?NAN:array[0]/array[1];}
static floating_t _sqrt(floating_t* array) { return array[0]<=0?NAN:sqrt(array[0]);}
static floating_t _negate(floating_t* array) { return  - array[0];}
static floating_t _invert(floating_t* array) { return array[0]==0?NAN:1.0/array[0];}

#define NEW_OPERATOR (OperatorPtr)(OperatorPtr)calloc(1,sizeof(Operator));\
	if(op==NULL) THROW_ERROR("BOUM");\
	ALL_OPERATORS->operators = (Operator**)realloc(ALL_OPERATORS->operators,sizeof(OperatorPtr)*(ALL_OPERATORS->size+1));\
	if(ALL_OPERATORS->operators==NULL) THROW_ERROR("BOUM");\
	ALL_OPERATORS->operators[ALL_OPERATORS->size]=op;\
	op->index = ALL_OPERATORS->size;\
	ALL_OPERATORS->size++
	

OperatorListPtr OperatorsListNew()
	{

	OperatorPtr op=NULL;
	OperatorListPtr ALL_OPERATORS=(OperatorListPtr)calloc(1,sizeof(OperatorList));
	if(ALL_OPERATORS==NULL) THROW_ERROR("");
	
	/** PLUS **/
	op =NEW_OPERATOR;
	strcpy(op->name,"Add");
	op->num_children = 2UL;
	op->eval = _plus;
	ALL_OPERATORS->plus=op;

	
	/** MINUS **/
	op = NEW_OPERATOR;
	strcpy(op->name,"Minus");
	op->num_children = 2UL;
	op->eval = _minus;
	ALL_OPERATORS->minus=op;

	
	/** MUL **/
	op = NEW_OPERATOR;
	strcpy(op->name,"Mul");
	op->num_children = 2UL;
	op->eval = _mul;
	ALL_OPERATORS->mul=op;

	
	/** DIV **/
	op = NEW_OPERATOR;
	strcpy(op->name,"Div");
	op->num_children = 2UL;
	op->eval = _div;
	ALL_OPERATORS->div=op;
	
	
	/** Negate **/
	op = NEW_OPERATOR;
	strcpy(op->name,"Negate");
	op->num_children = 1UL;
	op->eval = _negate;

	/** Invert **/
	op = NEW_OPERATOR;
	strcpy(op->name,"Invert");
	op->num_children = 1UL;
	op->eval = _invert;

	
	/** SQRT **/
	/*
	op =NEW_OPERATOR;
	strcpy(op->name,"Sqrt");
	op->num_children = 1UL;
	op->eval = _sqrt;
	*/

	return ALL_OPERATORS;
	} 

OperatorPtr OperatorListAt(OperatorListPtr list,size_t index)
	{
	if(index>=OperatorListSize(list))
		{
		fprintf(stderr,"%d >= %d\n",	index, OperatorListSize(list));
		THROW_ERROR("BOUM");
		}
	return list->operators[index];
	}

size_t OperatorListSize(OperatorListPtr list)
	{
	return list->size;
	}


static void GenomeEval1(
		const GenomePtr genome,
		const size_t rowIndex,
	 	size_t *nodeIndex,
		floating_t *value,
		boolean_t  *error
		)
	{
	NodePtr node=NULL;
	if( *error ) return;
	if( *nodeIndex >= GenomeSize(genome))
		{
		*error = 1;
		return;
		}
	
	node = GenomeAt(genome, *nodeIndex );
	switch(node->type)
		{
		case CONSTANT:
			{
			*value = node->core.constant;
			break;
			}
		case COLUMN:
			{
			*value = SpreadSheetAt(
				genome->config->spreadsheet,
				rowIndex,
				node->core.column
				);
			/* fprintf(stderr,"Value[%d][%d] = %f \n", rowIndex,node->core.column,*value ); */
			break;
			}
		case OPERATOR:
			{
			size_t i;
			floating_t values[MAX_ARITY];
			OperatorPtr op = OperatorListAt(genome->config->operators,node->core.operator);
			
			//fprintf(stderr,"Operator = %s \n",op->name);
			
			for(i=0;i< op->num_children;++i)
				{
				*nodeIndex=*nodeIndex + 1;
				GenomeEval1(
					genome,
					rowIndex,
					nodeIndex,
					&values[i],
					error
					);
				if(isnan(values[i]))
					{
					*error=1;
					break;
					}
				//fprintf(stderr,"%s[%d]  = %f\n",op->name,i,values[i]);
				}
			if( i != op->num_children )
				{
				//fprintf(stderr,"ERRORR nm children\n");
				*error=1;
				break;
				}
			*value = op->eval(values);
			//GenomePrint(genome,stderr);
			//fprintf(stderr,"eval==%f %s %f %f\n",*value,op->name,values[0],values[1]);
			
			if(isnan(*value))
				{
				*error=1;
				break;
				}
			break;
			}
		default:break;
		}
	}	


void _GenomePrint(const GenomePtr g,size_t *nodeIndex,FILE* out)
	{
	NodePtr node;
	if( *nodeIndex >= GenomeSize(g))
		{
		fprintf(out,"<EOF>");
		return;
		}
	
	node = GenomeAt(g, *nodeIndex );
	switch(node->type)
		{
		case CONSTANT:
			{
			fprintf(out,"%f", node->core.constant);
			break;
			}
		case COLUMN:
			{
			fprintf(out,"${%d}", node->core.column + 1);
			break;
			}
		case OPERATOR:
			{
			size_t i;
			OperatorPtr op = OperatorListAt(g->config->operators, node->core.operator);
			fprintf(out,"%s(", op->name);
			
			for(i=0;i< op->num_children;++i)
				{
				if(i>0) fprintf(out,",");
				*nodeIndex=*nodeIndex + 1;
				_GenomePrint(g,nodeIndex,out);
				}
			fprintf(out,")");
			break;
			}
		default:break;
		}

	
	}

void GenomePrint(const GenomePtr g,FILE* out)
	{
	size_t node_index=0;
	fprintf(stderr,"FITNESS=%E\tGENERATION=%ld\tSECONDS=%d\t(",g->fitness,g->generation,(int)difftime(g->creation,g->config->startup));
	_GenomePrint(g,&node_index,out);
	fputs(");\n",out);
	}

boolean_t GenomeEquals(const GenomePtr g1,const GenomePtr g2)
	{
	if(g1==g2) return 1;
	if(g1->node_count != g2->node_count ) return 0;
	return memcmp(
		(const void*)g1->nodes,
		(const void*)g2->nodes,
		g2->node_count*sizeof(Node)
		)==0;
	}


static void GenomeEval(GenomePtr g)
	{
	
	floating_t min_value= DBL_MAX;
	floating_t max_value=-DBL_MAX;
	size_t rowIndex;
	size_t num_errors=0UL;
	size_t max_errors=(size_t)(g->config->max_fraction_of_errors)*SpreadSheetRows(g->config->spreadsheet);
	floating_t  *norms=(floating_t*)calloc(SpreadSheetRows(g->config->spreadsheet),sizeof(floating_t));
	if(norms==NULL) THROW_ERROR("BOUM");
	
	if( SpreadSheetRows(g->config->spreadsheet)==0) THROW_ERROR("BOUM");
	
	for(rowIndex=0;
		rowIndex< SpreadSheetRows(g->config->spreadsheet);
		++rowIndex)
		{
		size_t nodeIndex=0;
		floating_t  value=0.0;
		boolean_t  error=0;
		
		
		GenomeEval1(g,rowIndex,&nodeIndex,&value,&error);
		if(error)
			{
			num_errors++;
			norms[rowIndex]=NAN;
			if( num_errors>max_errors )
				{
				//fprintf(stderr,"too many errors\n");
				break;
				}
			continue;
			}
		if( g->config->remove_introns )
			{
			g->node_count=1+nodeIndex;
			}
		norms[rowIndex]=value;
		if(value < min_value) min_value = value;
		if(value > max_value) max_value = value;
		}
	
	
	if(max_value==min_value || num_errors>max_errors)
		{
		g->bad_flag=1;
		g->fitness=NAN;
		}
	else
		{
		g->fitness=0.0;
		for(rowIndex=0;
			rowIndex< SpreadSheetRows(g->config->spreadsheet);
			++rowIndex)
			{
			double diff=0.0;
			if( isnan(norms[rowIndex]) ) continue;
			if(g->config->normalize_data)
				{
				norms[rowIndex]=(norms[rowIndex]-min_value)/(max_value-min_value);
				diff = fabs( norms[rowIndex] - g->config->spreadsheet->normalized[rowIndex] );
				}
			else
				{
				diff = fabs( norms[rowIndex] - SpreadSheetAt(g->config->spreadsheet,rowIndex,SpreadSheetColumns(g->config->spreadsheet)-1) );
				}
			g->fitness += pow(diff,2);
			}
		//fprintf(stderr,"fitness =%f\n",g->fitness);
		}
	free(norms);
	
	}


size_t SpreadSheetColumns(const SpreadSheetPtr ptr)
	{
	return ptr->columns;
	}

size_t SpreadSheetRows(const SpreadSheetPtr ptr)
	{
	return ptr->size/SpreadSheetColumns(ptr);
	}


/**
 * read a SpreadSheet from a FILE
 */
SpreadSheetPtr SpreadSheetRead(FILE* in)
	{
	size_t i=0UL;
	floating_t min_value,max_value;
	size_t curr_cols=0,curr_lines=0UL;
	SpreadSheetPtr p=(SpreadSheetPtr)calloc(1,sizeof(SpreadSheet));
	if(p==NULL) THROW_ERROR("Out of memory");
	
	while(!feof(in))
		{
		int c=fgetc(in);
		if(c==EOF) break;
		ungetc(c,in);
		for(;;)
			{
			p->data=(floating_t*)realloc(p->data,(p->size+1)*sizeof(floating_t));
			if(p->data==NULL) THROW_ERROR("OUT OF MEMORY");	
			if(fscanf(in,"%lf",&(p->data)[p->size])!=1)
				{
				THROW_ERROR("BAD FLOAT");	
				}
			p->size++;
			curr_cols++;
			
			c=fgetc(in);
			if(c=='\n' || c==EOF)
				{
				if(curr_cols<2)
					{
					fprintf(stderr,"Not enough column line %d\n",
						SpreadSheetRows(p));
					exit(EXIT_FAILURE);
					}
				else if(curr_lines == 0UL )
					{
					p->columns = curr_cols;
					}
				else if( curr_cols != p->columns)
					{
					fprintf(stderr,"Inconsistent column number %d\n",1+SpreadSheetRows(p));
					exit(EXIT_FAILURE);
					}
				curr_lines++;
				curr_cols=0UL;
				break;
				}
			if(c!='\t')
				{
				fprintf(stderr,"Expected a tab\n");
				exit(EXIT_FAILURE);
				}
			ungetc(c,in);
			}
		}
	
	if( SpreadSheetRows(p)==0) THROW_ERROR("Empty rows");
	
	min_value = DBL_MAX;
	max_value =-DBL_MAX;
	p->normalized=calloc(SpreadSheetRows(p),sizeof(floating_t));
	for(i=0;i< SpreadSheetRows(p);++i)
		{
		floating_t v= SpreadSheetAt(p,i,SpreadSheetColumns(p)-1);
		if(v < min_value) min_value=v;
		if(v > max_value) max_value=v;
		}
	
	if(max_value==min_value)
		{
		fprintf(stderr,"Min==Max\n");
		exit(EXIT_FAILURE);
		}
	
	for(i=0;i< SpreadSheetRows(p);++i)
		{
		floating_t v=SpreadSheetAt(p,i,SpreadSheetColumns(p)-1);
		p->normalized[i] = (v-min_value)/(max_value-min_value);
		}
	
	
	
	return p;
	}

void GenomeMute(GenomePtr g)
	{
	size_t i;
	ConfigPtr cfg=g->config;
	
	while(  GenomeSize(g) >0 &&
		RANDOM_FLOAT(cfg) < cfg->probability_mutation)
		{
		i= RANDOM_SIZE_T(cfg,g->node_count);
		
		float rnd= RANDOM_FLOAT(cfg);
		if( rnd < 0.05 )
			{
			//insert
			size_t j,insert_size= 1 + RANDOM_SIZE_T( cfg , 5); 
			g->nodes = (NodePtr)realloc(g->nodes, (g->node_count + insert_size)*sizeof(Node));
			if(g->nodes==NULL) THROW_ERROR("boum");
			memmove(
				(void*)&g->nodes[i+insert_size],//dest
				(void*)&g->nodes[i],//src
				(GenomeSize(g) - i)*sizeof(Node)
				);
			for(j=0;j< insert_size;++j)
				{
				NodeInit(cfg,&(g->nodes[i+j]));
				}
			g->node_count+=insert_size;
			}
		else if(rnd< 0.1)
			{
			//delete
			size_t del_size= 1 + RANDOM_SIZE_T( cfg , 5); 
			if( i+del_size > GenomeSize(g))
				{
				del_size = GenomeSize(g) - i;
				}
				
			memmove(
				(void*)&g->nodes[i],//dest
				(void*)&g->nodes[i+del_size],//src
				(GenomeSize(g) - (i+del_size))*sizeof(Node)
				);
				
			g->node_count -= del_size;
				
			}
		else if(rnd < 0.2)
			{
			//silent mutation
			size_t insert_size=2; 
			g->nodes = (NodePtr)realloc(g->nodes, (g->node_count + insert_size)*sizeof(Node));
			if(g->nodes==NULL) THROW_ERROR("boum");
			
			
			memmove(
				(void*)&g->nodes[i+1+insert_size],//dest
				(void*)&g->nodes[i+1],//src
				(GenomeSize(g) - (i+1))*sizeof(Node)
				);
			memcpy(
				(void*)&g->nodes[i+1],//dest
				(void*)&g->nodes[i],//src
				sizeof(Node)
				);
			memcpy(
				(void*)&g->nodes[i+2],//dest
				(void*)&g->nodes[i],//src
				sizeof(Node)
				);
			g->node_count+=insert_size;
			
			switch(  RANDOM_SIZE_T( cfg , 4) )
				{
				case 0: {
					g->nodes[i].type=OPERATOR;
					g->nodes[i].core.operator= cfg->operators->plus->index;
					
					g->nodes[i+2].type=CONSTANT;
					g->nodes[i+2].core.constant= 0.0;
					break;
					};
				case 1: {
					g->nodes[i].type=OPERATOR;
					g->nodes[i].core.operator= cfg->operators->minus->index;
					
					g->nodes[i+2].type=CONSTANT;
					g->nodes[i+2].core.constant= 0.0;
					break;
					};
				case 2: {
					g->nodes[i].type=OPERATOR;
					g->nodes[i].core.operator= cfg->operators->mul->index;
					
					g->nodes[i+2].type=CONSTANT;
					g->nodes[i+2].core.constant= 1.0;
					break;
					};
				default:
					{
					
					g->nodes[i].type=OPERATOR;
					g->nodes[i].core.operator= cfg->operators->mul->index;
					
					g->nodes[i+2].type=CONSTANT;
					g->nodes[i+2].core.constant= 1.0;
					
					break;
					}
				}
			}
		else if( rnd < 0.3)
			{
			NodeInit(cfg,&(g->nodes[i]));
			}
		else
			{
			NodeMut(cfg,&(g->nodes[i]));
			}
		}
	}

static GenomePtr GenomeXCross(
	ConfigPtr cfg,
	GenomePtr p1,
	GenomePtr p2
	)
	{
	GenomePtr g = GenomeNew1(cfg);
	size_t xcross1;
	size_t xcross2;
	
	if(RANDOM_FLOAT(cfg)<0.1)
		{
		xcross1 = RANDOM_SIZE_T(cfg,GenomeSize(p1));
		xcross2 = RANDOM_SIZE_T(cfg,GenomeSize(p2));
		}
	else
		{
		xcross1 = RANDOM_SIZE_T(cfg, MIN( GenomeSize(p1), GenomeSize(p2) ));
		xcross2 = xcross1;
		}
	
	
	g->node_count=  xcross1 + (GenomeSize(p2)-xcross2);
	
	
	g->nodes = (NodePtr)calloc(g->node_count,sizeof(Node));
	if(g->nodes==NULL) THROW_ERROR("boum");
	memcpy((void*)g->nodes,(void*)p1->nodes,xcross1*sizeof(Node));
	memcpy((void*)&g->nodes[xcross1],(void*)&(p2->nodes[xcross2]),(GenomeSize(p2)-xcross2)*sizeof(Node));
	
	GenomeMute(g);
	return g;
	}

static void doWork(ConfigPtr config)
	{
	GenerationPtr gen=NULL;
	config->curr_generations=0L;
	GenomePtr best=NULL;
	/* create initial family */
	gen = GenerationNew(config);
	
	while( config->max_generations==-1L || config->curr_generations < config->max_generations )
		{
		size_t i,j;
		GenerationPtr gen1= GenerationNew1(config);
		GenerationPtr tmp = NULL;
		
		while( GenerationCount(gen) < config->max_genomes_per_generation )
			{
			GenerationAdd(gen,GenomeNew(config));
			}
		
		if(best!=NULL && config->best_will_survive)
			{
			GenomePtr copy=GenomeClone(best);
			GenomeMute(copy);
			GenerationAdd(gen,copy);
			}
		
		for(i=0 ; i < GenerationCount(gen) ; ++i)
			{
			GenomePtr gi=GenerationAt(gen,i);
			for(j=0 ; j < GenerationCount(gen) ; ++j)
				{
				GenomePtr newgen=NULL;
				GenomePtr gj=GenerationAt(gen,j);
				if( i == j && !config->enable_self_self) continue;
				
				newgen =  GenomeXCross(config,gi,gj);
				
				if(newgen==NULL ||
					GenomeSize(newgen)==0 ||
					GenomeSize(newgen) < config->min_base_per_genome ||
					GenomeSize(newgen) > config->max_base_per_genome ||
					(config->remove_clone && (GenomeEquals(newgen,gi) || GenomeEquals(newgen,gj)) )
					)
					{
					GenomeFree(newgen);
					continue;
					}
				GenomeEval(newgen);
				if( newgen->bad_flag )
					{
					GenomeFree(newgen);
					continue;
					}
				newgen->generation = config->curr_generations;
				
				gen1->genomes = (Genome**)realloc(gen1->genomes,(gen1->genome_count+1)*sizeof(Genome*));
				if(gen1->genomes==NULL) THROW_ERROR("boum");
				gen1->genomes[ gen1->genome_count ] = newgen ;
				gen1->genome_count++;
				}
			}

		
		/* sort on fitness */
		qsort(	(void*)gen1->genomes,
			gen1->genome_count,
			sizeof(GenomePtr),
			_GenomeCompare
			);	
		
		/* remove duplicate children */
		i=0;
		while(i+1< GenerationCount(gen1))
			{
			if( GenomeEquals(GenerationAt(gen1,i),GenerationAt(gen1,i+1) ))
				{
				GenomeFree( GenerationAt(gen1,i) );
				memmove(
					(void*)&gen1->genomes[i],
					(const void*)&gen1->genomes[i+1],
					sizeof(GenomePtr)*(GenerationCount(gen1)-(i+1))
					);
				gen1->genome_count--;
				}
			else
				{
				++i;
				}
			}
				
		while( GenerationCount(gen1) > config->min_genomes_per_generation ||
			(config->massive_extinction_every!=-1L && config->curr_generations % config->massive_extinction_every==0 &&  GenerationCount(gen1) >0 ))
			{
			GenomeFree( gen1->genomes[ GenerationCount(gen1) -1 ] );
			gen1->genome_count--;
			}
		


		
		
		if(GenerationCount(gen1)>0 && !(GenerationAt(gen1,0)->bad_flag) )
			{
			//GenomePrint(GenerationAt(gen1,0),stderr);
			if(best==NULL || best->fitness > GenerationAt(gen1,0)->fitness)
				{
				GenomeFree(best);
				best=GenomeClone( GenerationAt(gen1,0) );

				GenomePrint(best,stdout);
				if(best->fitness < config->min_fitness)
					{
					fprintf(stderr,"min fitness reached\n");
					break;
					}
				}
			}
		else
			{
			fprintf(stderr,"too many errors\n");
			}
		tmp = gen;
		gen = gen1;
		GenerationFree(tmp);
		config->curr_generations++;
		}
	GenerationFree(gen);
	}

int main(int argc,char** argv)
	{
	Config config;
	memset((void*)&config,0,sizeof(Config));
	config.max_generations = -1L;
	config.min_genomes_per_generation=5;
	config.max_genomes_per_generation=50;
	config.min_base_per_genome=3;
	config.max_base_per_genome=30;
	
	config.probability_mutation=0.01;
	config.max_fraction_of_errors=0.1f;
	config.weight_constant=10;
	config.weight_op=10;
	config.weight_column=10;
	
	config.remove_introns=0;
	config.best_will_survive=0;
	config.enable_self_self=0;
	config.normalize_data=0;
	config.massive_extinction_every= -1L;
	config.sort_on_genome_size=0;
	config.remove_clone=0;
	
	config.min_fitness = 1E-6;
	config.startup=time(NULL);
	config.operators = OperatorsListNew();
	config.seedp=(int)time(NULL);
	srand(time(NULL));
	
	for(;;)
		{
		struct option long_options[] =
		     {
		       /* These options set a flag. */
		       //{"verbose", no_argument,       &verbose_flag, 1},
		       //{"brief",   no_argument,       &verbose_flag, 0},
		       /* These options don't set a flag.
		          We distinguish them by their indices. */
		       //{"add",     no_argument,       0, 'a'},
		       //{"append",  no_argument,       0, 'b'},
		       //{"delete",  required_argument, 0, 'd'},
		       //{"create",  required_argument, 0, 'c'},
		       {"enable-self-self",  no_argument , &config.enable_self_self , 1},
		       {"enable-best-survives",  no_argument , &config.best_will_survive , 1},
		       {"enable-remove-introns",  no_argument , &config.remove_introns , 1},
		       {"enable-remove-clone",  no_argument , &config.remove_clone , 1},
		       {"genome-size-matters",  no_argument , &config.sort_on_genome_size , 1},
		       {"normalize-data",  no_argument , &config.normalize_data , 1},
		       {"generations",    required_argument, 0, 'g'},
		       {"random-seed",    required_argument, 0, 's'},
		       {"min-bases",    required_argument, 0, 'b'},
		       {"max-bases",    required_argument, 0, 'B'},
		       {"min-genomes",    required_argument, 0, 'n'},
		       {"max-genomes",    required_argument, 0, 'N'},
		     
		       {0, 0, 0, 0}
		     };
		 /* getopt_long stores the option index here. */
		int option_index = 0;
	     	int c = getopt_long (argc, argv, "g:s:b:B:n:N:",
		                    long_options, &option_index);
		if(c==-1) break;
		switch(c)
			{
			case 's':
				{
				config.seedp=atoi(optarg);
				break;
				};
			case 'g':
				{
				config.max_generations=atol(optarg);
				break;
				};
			case 'b':
				{
				config.min_base_per_genome=strtoul(optarg,NULL,10);
				break;
				};
			case 'B':
				{
				config.max_base_per_genome=strtoul(optarg,NULL,10);
				break;
				};
			case 'n':
				{
				config.min_genomes_per_generation=atoi(optarg);
				break;
				};
			case 'N':
				{
				config.max_genomes_per_generation=atoi(optarg);
				break;
				};
			case 0: break;
			case '?': break;
			default: exit(EXIT_FAILURE); break;
			}
		}
	
	if( config.max_base_per_genome < config.min_base_per_genome)
		{
		fprintf(stderr," config.max_base_per_genome < config.min_base_per_genome\n");
		return EXIT_FAILURE;
		}
	
	if( config.min_genomes_per_generation<1)
		{
		fprintf(stderr," bad  config.min_genomes_per_generation\n");
		return EXIT_FAILURE;
		}
	
	if( config.max_genomes_per_generation < config.min_genomes_per_generation)
		{
		fprintf(stderr," config.max_genomes_per_generation < config.min_genomes_per_generation\n");
		return EXIT_FAILURE;
		}
	
	if(optind==argc)
		{
		config.spreadsheet=SpreadSheetRead(stdin);
		}
	else
		{
		FILE* in=fopen(argv[optind],"r");
		if(in==NULL)
			{
			fprintf(stderr,"Cannot open %s %s\n",
				argv[optind],
				strerror(errno)
				);
			return EXIT_FAILURE;
			}
		config.spreadsheet=SpreadSheetRead(in);
		fclose(in);
		}
	doWork(&config);
	return EXIT_SUCCESS;
	}
