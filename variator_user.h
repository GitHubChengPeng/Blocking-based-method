/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
 
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
 
  ========================================================================
  PISALIB 

  Pisa basic functionality that have to be implemented by the user  
   
  Header file.
  
  file: variator_user.h
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  last change: $date$
  
  ========================================================================
*/

#ifndef VARIATOR_USER_H
#define VARIATOR_USER_H

#define PISA_WIN		/**** replace with PISA_WIN if compiling for Windows */
						/**** replace with PISA_UNIX if compiling for UNIX */

/* maximal length of filenames */
#define FILE_NAME_LENGTH 128  /**** change the value if you like */

/* maximal length of entries in local cfg file */
#define CFG_NAME_LENGTH 128   /**** change the value if you like */




//----------BA算法-------------

extern	double SM;			//	Safety margin
extern	double A01;			//	The proportion of trans which change 0's into ?.  0 < A01 < 1

//-----------------------------




/*---| declaration of global variables (defined in variator_user.c) |-----*/

extern char *log_file; /* file to log to */

extern char paramfile[]; /* file with local parameters */

extern int length; /* length of the binary string */


/*-----------------------------------------------------------------------*/

struct individual_t
{
     /**********| added for PPDM |**************/
     int *bit_string; /* the binary decision variables */
     int length;      /* length of the bit_string */
     double *objective_value; /* objective values */
	 //int *objective_value;
     
     /**********| addition for PPDM end |*******/
	 // 对于PPDM, objective_value应是整数;但考虑到最终的fitness func. 还是要将到 reference point 的距离合并过来，
	 // 且double型数据也可以进行计数工作（alpha,beta,gamma）， 更重要的是Slector里的目标向量是double型的
	 // 因此这里仍然采用double型
};

/*-------------------------| individual |-------------------------------*/




typedef struct individual_t individual; 
/* 'individual_t' has to be defined in variator_user.h */




/*-------------------| functions for individual struct |----------------*/

void free_individual(individual *ind);
/* Frees the memory for one indiviual.

   post: memory for ind is freed
*/


double get_objective_value(int identity, int i); 
/* Gets objective value of an individual.

   pre: 0 <= i <= dim - 1 (dim is the number of objectives)

   post: Return value == the objective value number 'i' in individual '*ind'.
         If no individual with ID 'identity' return value == -1. 
*/   

/*-------------------------| statemachine |-----------------------------*/


int state0(); 
/* Do what needs to be done in state 0.

   pre: The global variable 'paramfile' contains the name of the
        parameter file specified on the commandline.
        The global variable 'alpha' contains the number of indiviuals
        you need to generate for the initial population.
                
   post: Optionally read parameter specific for the module.
         Optionally do some initialization.
         Initial population created.
         Information about initial population written to the ini file
         using write_ini().
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/


int state2();
/* Do what needs to be done in state 2.

   pre: The global variable 'mu' contains the number of indiviuals
        you need to read using 'read_sel()'.
        The global variable 'lambda' contains the number of individuals
        you need to create by variation of the individuals specified the
        'sel' file.
        
   post: Optionally call read_arc() in order to delete old uncessary
         individuals from the global population.
         read_sel() called
         'lambda' children generated from the 'mu' parents
         Children added to the global population using add_individual().
         Information about children written to the 'var' file using
         write_var().
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/


int state4();
/* Do what needs to be done in state 4.

   pre: State 4 means the variator has to terminate.

   post: Free all memory.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/


int state7();
/* Do what needs to be done in state 7.

   pre: State 7 means that the selector has just terminated.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/


int state8();
/* Do what needs to be done in state 8.

   pre: State 8 means that the variator needs to reset and get ready to
        start again in state 0.

   post: Get ready to start again in state 0, this includes:
         Free all memory.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/


int state11();
/* Do what needs to be done in state 11.

   pre: State 11 means that the selector has just reset and is ready
        to start again in state 1.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/

int is_finished();
/* Tests if ending criterion of your algorithm applies.

   post: return value == 1 if optimization should stop
         return value == 0 if optimization should continue
*/


int read_local_parameters();
/* read local parameters from file */

individual *new_individual();

individual *copy_individual(individual* ind);

int variate(int *parents, int *offspring);


/* flips one bit in ind */
int one_bit_mutation(individual *ind);

/* flips each bit with probability bit_turn_prob/ind.lenght */
int indep_bit_mutation(individual *ind);

/* xover is stored in ind1 and ind2 */
int one_point_crossover(individual *ind1, individual *ind2);

/* xover is stored in ind1 and ind2 */
int uniform_crossover(individual *ind1, individual *ind2);

int irand(int range);
/* Generate a random integer. */

double drand(double range);
/* Generate a random double. */



int eval(individual *p_ind);
/* Determines the objective value. PISA always minimizes. */

void cal_fitness(individual *x);	// 最基本的“模式匹配” 速度慢

void cal_fitness_two(individual *x);  // 对于1-item频繁 和 2-item 频繁 采用cal_fitness_two()评价个体， 而不是cal_fitness( )

void cal_fitness_two_rules(individual *x);

void cal_fitness_by_trie_FIM( individual *x );    //【通用】, 不管itemset包含的item有多少个   【精彩， 利用树和递归】

void cal_fitness_by_trie_AR( individual *x ) ;    //【通用】, 不管rules包含的item有多少个   【精彩， 利用树和递归】




void read_sensitive_item();		//	打开敏感项文件，读取到vector<item>类型的 s_set中

void filter_sen_trans();		//  过滤出包含敏感项的 transactions ; 同时计算每个敏感项频次

void cal_length();				//  计算length(染色体编码长度),根据敏感项频次和sup_l  

void filter_prelarge_trans();   //  过滤出包含prelarge项的 transactions ; 同时计算每个prelarge项频次

void cal_overlap();				//  计算每个频繁项的overlap degree，写入文件

void write_stat();				//  向文件写入统计信息

void init_random_shuffle();		// 借用STL中random_shuffle生成不重复的随机数序列，初始化。
								// 在函数new_individual(),crossover(), mutation()中用到不重复的随机序列

// 生成的子代染色体中的gene取值 不重复
// 新的交叉算子！！！！！！
// 新的交叉算子！！！！！！
int shuffle_crossover(individual *ind1, individual *ind2);


// 生成的子代染色体中的gene取值 不重复
// 新的变异算子！！！！！！
// 新的变异算子！！！！！！
int shuffle_mutation(individual *ind);




void write_output_file();

void write_output_obj();

void write_ini_obj();

/**********| addition for PPDM end |*******/

#endif /* VARIATOR_USER.H */
