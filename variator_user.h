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




//----------BA�㷨-------------

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
	 // ����PPDM, objective_valueӦ������;�����ǵ����յ�fitness func. ����Ҫ���� reference point �ľ���ϲ�������
	 // ��double������Ҳ���Խ��м���������alpha,beta,gamma���� ����Ҫ����Slector���Ŀ��������double�͵�
	 // ���������Ȼ����double��
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

void cal_fitness(individual *x);	// ������ġ�ģʽƥ�䡱 �ٶ���

void cal_fitness_two(individual *x);  // ����1-itemƵ�� �� 2-item Ƶ�� ����cal_fitness_two()���۸��壬 ������cal_fitness( )

void cal_fitness_two_rules(individual *x);

void cal_fitness_by_trie_FIM( individual *x );    //��ͨ�á�, ����itemset������item�ж��ٸ�   �����ʣ� �������͵ݹ顿

void cal_fitness_by_trie_AR( individual *x ) ;    //��ͨ�á�, ����rules������item�ж��ٸ�   �����ʣ� �������͵ݹ顿




void read_sensitive_item();		//	���������ļ�����ȡ��vector<item>���͵� s_set��

void filter_sen_trans();		//  ���˳������������ transactions ; ͬʱ����ÿ��������Ƶ��

void cal_length();				//  ����length(Ⱦɫ����볤��),����������Ƶ�κ�sup_l  

void filter_prelarge_trans();   //  ���˳�����prelarge��� transactions ; ͬʱ����ÿ��prelarge��Ƶ��

void cal_overlap();				//  ����ÿ��Ƶ�����overlap degree��д���ļ�

void write_stat();				//  ���ļ�д��ͳ����Ϣ

void init_random_shuffle();		// ����STL��random_shuffle���ɲ��ظ�����������У���ʼ����
								// �ں���new_individual(),crossover(), mutation()���õ����ظ����������

// ���ɵ��Ӵ�Ⱦɫ���е�geneȡֵ ���ظ�
// �µĽ������ӣ�����������
// �µĽ������ӣ�����������
int shuffle_crossover(individual *ind1, individual *ind2);


// ���ɵ��Ӵ�Ⱦɫ���е�geneȡֵ ���ظ�
// �µı������ӣ�����������
// �µı������ӣ�����������
int shuffle_mutation(individual *ind);




void write_output_file();

void write_output_obj();

void write_ini_obj();

/**********| addition for PPDM end |*******/

#endif /* VARIATOR_USER.H */
