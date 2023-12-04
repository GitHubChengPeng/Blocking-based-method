/*========================================================================
  PPDM based on PISA
 
  ========================================================================

  HIT-- Harbin Institue of Technology
 
  ========================================================================
  PISALIB 
  
  Pisa basic functions for use in variator_user.c
  
  C file.
  
  file: variator.c
  author: Peng CHENG, chengpeng_mailbox@163.com

  revision by: 
  last change: $date$
  
  ========================================================================
*/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <math.h>
//#include <cmath>

#include <string.h>

#include "variator.h"
#include "variator_user.h"


/*==========================Add for PPDM_dt=========================*/


#include "old_apriori.h"


/*--------| global variable definitions |-------*/

/*==== declared in variator_user.h used in other files as well ====*/

char *log_file = "PPDM_dt_error.log"; /**** Changed for PPDM_dt.*/

char paramfile[FILE_NAME_LENGTH]; /* file with local parameters */


/*========================== END for PPDM =========================*/




/*==== only used in this file ====*/

double *weights;
double *profits;
double *profitSums;
double *capacities;
int *selectOrder;


/* local parameters from paramfile*/
int seed;   /* seed for random number generator */

//int db_size; /* ���� apriori.cpp�ļ��ж���*/


/*=================================Add for PPDM ===================================*/

	int length; /* length of the binary string */  // MAX_del_trans�ǵ��ߣ� ������sup_l,��������ɾ��������length(�����볤��)
	int maxgen; /* maximum number of generations (stop criterion) */


	//------------------------------------------------------------------------------[void cal_length()����]
	int sum1_prelarge =0,  sum2_sensitive = 0;
	int length1;	// length1��ֵ�����ˡ����α��롱�ķֶ�λ�ã� ǰ��γ���Ϊlength1��ǰ�����prelarge_trans_array[]������
	int length2;	// ���γ���Ϊlength2, ���α����� sen_trans_array[] ������
	//------------------------------------------------------------------------------


	//----------------------------------------------
	//			�����漰sensitive item sets��trans
	//----------------------------------------------

		
	int sen_trans_array_size=0;				// sen_trans_array[]���鳤��
 											// ��sen_trans_setת������������˰�������������ݿ��¼�� ���ݿ��е����кţ����кŷ�Χ[0,database.size()-1];
	int ** sen_trans_array;					// ��filter_sen_trans()�������ڴ�ռ�,����ÿ��sen item�漰��trans���


	int * sen_frequency_array=NULL;			// ������ÿ��������ĳ���Ƶ�Σ�Ҳ����filter_sen_trans()�������ڴ�ռ�

	//----------------------------------------------
	//			�����漰prelarge item sets��trans
	//----------------------------------------------
	set<int> prelarge_trans_set;				// pre_large���漰��trans���ϣ� ʹ��set<int>���Ϳ��Է�ֹ �����ظ���Ԫ��

	int prelarge_trans_array_size=0;			// prelarge_trans_array[]���鳤��
	int * prelarge_trans_array=NULL;			// ��prelarge_trans_setת������������˰�������������ݿ��¼�� ���ݿ��е����кţ����кŷ�Χ[0,database.size()-1];
												// ��filter_prelarge_trans()�������ڴ�ռ�

	int * prelarge_frequency_array=NULL;		// ������ÿ��prelarge��ĳ���Ƶ�Σ�Ҳ����filter_prelarge_trans()�������ڴ�ռ�

	
	//----------------------------------------------
	//			����overlap for each large item set
	//----------------------------------------------

	double sel_overlap_percent=0.1;		// ȷ���ڸ߳�sup_l�İٷ�֮���ٱ�����trans��������overlap,��Ϊ��ƵƵ�����Ե�Ƶ��overlap����Ӱ��ϴ�
										// ����ģ����ƵƵ����֮���overlap����ҳ���sup_l����ĸ�ƵƵ�������ʺ����������Ϊ��������ɾ����trans̫��

	int *encoding_len;					// ��̬�������飬 ���ڴ���ÿ��sensitive item��Ӧ��Ⱦɫ���еı��볤�ȡ���ο�PPT�С���α��롱���֡�


	vector<int> v1;					// ���ڴ�Ŵ���������е�0~prelarge_trans_array_size-1 ���У� �����void init_random_shuffle()

	vector<vector<int>> v2;			// ÿһ�����ڴ�Ŵ���������е�0~sen_trans_vec[i].size()-1���У� �����void init_random_shuffle()

/*================================================================================*/


/*********************************************/
char outfile[FILE_NAME_LENGTH]; /* output file for last population */

int gen;

int mutation_type;		/*	0 = no mutation
							1 = one bit mutation
							2 = independent bit mutation */

int recombination_type; /*	0 = no recombination
							1 = one point crossover
							2 = uniform crossover */

double mutat_prob;		/* probability that individual is mutated */

double recom_prob;		/* probability that 2 individual are recombined */

double bit_turn_prob;	/* probability, that bit is turned when mutation occurs only used for independent bit mutation */



/*-------------------------| individual |-------------------------------*/

void free_individual(individual *ind) 
/* Frees the memory for one indiviual.

   post: memory for ind is freed
*/
{
     /**********| added for KNAPSACK |**************/
     
     if (ind == NULL)
          return;
     
     free(ind->bit_string);
     free(ind->objective_value);
     
     /**********| addition for KNAPSACK end |*******/
     
     free(ind);
}

double get_objective_value(int identity, int i)
/* Gets objective value of an individual.

   pre: 0 <= i <= dimension - 1 (dimension is the number of objectives)

   post: Return value == the objective value number 'i' in individual '*ind'.
         If no individual with ID 'identity' return value == -1. 
*/   
{
     /**********| added for KNAPSACK |**************/
     individual *temp_ind;
     /**********| addition for KNAPSACK end |*******/
     
     double objective_value = -1.0;

     assert(0 <= i && i < dimension); /* asserting the pre-condition */
     
     /**********| added for KNAPSACK |**************/
    
     if (i < 0 || i > (dimension - 1))
	  return(-1);
     
     temp_ind = get_individual(identity);
     if (temp_ind == NULL)
		return(-1);
     
     objective_value = temp_ind->objective_value[i];
     
     /**********| addition for KNAPSACK end |*******/
  
     return (objective_value);
}

/*-------------------------| statemachine functions |-------------------*/

int state0() 
{
     int i;
     
     int result; /* stores return values of called functions */
     int *initial_population; /* storing the IDs of the individuals */

	 gen = 1;			// genΪȫ�ֱ����� state0()Ϊ��һ��

     initial_population = (int *) malloc(alpha * sizeof(int)); 
     if (initial_population == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }

     /*=============| added for PPDM |======================*/

     result = read_local_parameters();

	 apriori_bodon();			// ����apriori_bodon()����õ����С���������д���ļ�(��������bodon_ruleset)�� �ٴ���ѡ�����й���
								// ��ȡ��1-item���͡�2-item(����)��Ƶ�������ݽṹ����ȡ ������main_trie_copy��

	 cout<<"\nPlease define the sensitive rules in the file.  When finished press [Enter] to continue! \n ";
	 getchar();

	 read_sensitive_item();		// ���������ļ�����ȡ��vector<item>���͵� s_set_pair��

	 apriori_large_sets();		// ���þ�apriori�㷨 [����ѵ���apriori_bodon,��apriori�㷨ֻ����sup_u��sup_l]
								// ��һ����cal_support()����sup_u��sup_l������sup_l-->��֪��MAX_del_trans-->��֪��"ÿ���������support"-->����filter_sen_trans();
								// ע������� ������filter_sen_trans()

	 //if(AR == FIM_or_AR && 2 == SIZE_threshold)
	 //	apriori_gen_rules();	// �Լ�д�Ĳ�����������!!!��ֻ�ܲ���2-item���� apriori_bodon()�Ѿ������˹������򲢱���

								
	 //filter_sen_trans();		// �ҳ���������������ݿ��¼��ţ����õ����������Ƶ��
								// ��Ϊ��old_apriori.cpp�ļ� cal_support()-->cal_MAX_delTrans()�����е���

	 //filter_prelarge_trans();	// �ҳ�����prelarge������ݿ��¼���漰prelarge��
								// ɾ��item�Ĳ��� �Ͳ���Ҫ �˺�����

	 //cal_length();				// ����length(Ⱦɫ����볤��),����������Ƶ�κ�sup_l ���㣻read_local_parameters()������lengthֵ���滻 
								// ���漰prelarge��

	 write_stat();				//  ���ļ�д��ͳ����Ϣ���漰prelarge��

	 //cal_overlap();			// ����ÿ��Ƶ�����overlap degree��д���ļ�
								// (����һ�����ݼ� ����һ�ξ͹��ˣ�)

	 //init_random_shuffle();		// ����STL��random_shuffle���ɲ��ظ�����������У���ʼ�����ں���new_individual(),crossover(), mutation()���õ����ظ����������
								// ���漰prelarge��

	 /*====================================================*/
     
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't read local parameters");
          return (1);
     }
     
     /* create first alpha individuals */
     for(i = 0; i < alpha; i++)
     {
          initial_population[i] = add_individual(new_individual());  // new_individual()Ϊ�¸��������ڴ棬�����ʼ��
												//���漰prelarge��	 // ���¸����ַ ������Ⱥ global_population.individual_array[]
																	 // add_individual��������ֵΪ ��������Ⱥ�е�λ�ã������֤�ţ�
          if(initial_population[i] == -1)
               return(1);
     } 

     result = write_ini(initial_population);	//����һ��Ⱥ������֤�ţ�����������Ⱥ�еĲ���λ�ã��� ��άĿ��ȡֵ д��ini_file�ļ�

     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write ini");
          free(initial_population);
          return (1);
     }

     free(initial_population);

	 cout<<"\n>>> State 0() finished! \n";
     return (0);
}


int state2()
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
{
     int *parent_identities, *offspring_identities; /* array of identities */
     int result; /* stores return values of called functions */

     parent_identities = (int *) malloc(mu * sizeof(int)); 
     if (parent_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }

     offspring_identities = (int *) malloc(lambda * sizeof(int)); 
     if (offspring_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }
     
     result = read_sel(parent_identities);
     if (result != 0) /* if some file reading error occurs, return 2 */
          return (2);

     result = read_arc(); 
     if (result != 0) /* if some file reading error occurs, return 2 */
          return (2);

	 if(1==gen)
		 	write_ini_obj();   // add for PPDM�� ֻд���һ��Ⱥ������š�Ŀ������

     /*=================| added for PPDM |=========================*/

     gen++;		// ����������1
	 cout<<"\n\n"<<"gen = "<<gen<<endl;


     result = variate(parent_identities, offspring_identities);  //���漰prelarge��
     if (result != 0)
          return (1);


	 printf("\n\n\n Current generation: %d", gen);

     /*================| addition for PPDM end |===================*/


     result = write_var(offspring_identities);
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write var");
          free(offspring_identities);
          free(parent_identities);
          return (1);
     }

     free(offspring_identities);
     free(parent_identities);
     return (0);
}
 

int state4() 
/* Do what needs to be done in state 4.

   pre: State 4 means the variator has to terminate.

   post: Free all memory.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     
     int result;
     result = read_arc();
   
     if (0 == result) /* arc file correctly read
                         this means it was not read before,
                         e.g., in a reset. */
     {
        write_output_file();  // д������š�Ŀ����������������

		write_output_obj();   // add for PPDM�� ֻд������š�Ŀ������
     }


	 /*-----------�ͷš���α��롱�漰�Ķ�̬����Ķ����ڴ�--------------*/
	 //sen_frequency_array=(int *)new int[s_set.size()];	

	 for(unsigned int i=0;i<s_set.size();i++)
	 {
			delete [] sen_trans_array[i];
	 }	
	 delete [] sen_trans_array;

	 delete [] sen_frequency_array;

	 delete [] prelarge_trans_array;

	 delete [] encoding_len;


     
     return (0);
}


int state7()
/* Do what needs to be done in state 7.

   pre: State 7 means that the selector has just terminated.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return(0);  
}


int state8()
/* Do what needs to be done in state 8.

   pre: State 8 means that the variator needs to reset and get ready to
        start again in state 0.

   post: Get ready to start again in state 0. 
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     /**********| added for KNAPSACK |**************/

   int result;
   
   gen = 1;
     
   result = read_arc();   //ɾ��������Ⱥ�еķǴ浵����

   if (0 == result) /* arc file correctly read
                       this means it was not read before */
   {
      write_output_file();
   }
   
   /*all of the following are allocated again in read_local_parameters() */
   /*
   free(weights);
   weights = NULL;
   
   free(profits);
   profits = NULL;
   
   free(profitSums);
   profitSums = NULL;
   
   free(capacities);
   capacities = NULL;
   
   free(selectOrder);
   selectOrder = NULL;
   */
   
     /**********| addition for KNAPSACK end |*******/
     
     return (0);
}


int state11()
/* Do what needs to be done in state 11.

   pre: State 11 means that the selector has just reset and is ready
        to start again in state 1.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return (0);  
}


int is_finished()
/* Tests if ending criterion of your algorithm applies.

   post: return value == 1 if optimization should stop
         return value == 0 if optimization should continue
*/
{
     /**********| added for KNAPSACK |**************/
     return (gen >= maxgen);
     /**********| addition for KNAPSACK end |*******/
}


/**********| added for KNAPSACK |**************/

/* temporary array used for sorting items */
static double* profitWeightRatios;


/* function used for sorting items in decreasing order */
// ����profitWeightRatios�������ݶ��±��������

int cmpItems(const void*  itemPtr1, const void*  itemPtr2)
{
    if (profitWeightRatios[*((int*) itemPtr1)] > profitWeightRatios[*((int*) itemPtr2)])  
		return -1;
    if (profitWeightRatios[*((int*) itemPtr2)] > profitWeightRatios[*((int*) itemPtr1)])  
		return 1;
    return 0;
}

int RandomInt(int  min, int  max)
/* generates a random integer */
{
  return (int) (min + (rand() % (max - min + 1)));
} /* RandomInt */



//ɾ���ܶ�knapsack��ش���
int read_local_parameters()
{
     FILE *fp;
     //int result;
     char str[CFG_NAME_LENGTH];

	 cout<<"\nA new [Blocking Algorithm] with fewer undesirable side effects and more desirable side effects.\n";

	 cout<<"\nBEGIN to read local parameters!\n";
     /* reading parameter file with parameters for selection */
     fp = fopen(paramfile, "r"); 
     assert(fp != NULL);


     fscanf(fp, "%s", str);
     assert(strcmp(str, "seed") == 0);
     fscanf(fp, "%d", &seed);
     

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "SM") == 0);
     fscanf(fp, "%lf", &SM);				//����safety margin
	 cout<<str<<" == "<<SM<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "A01") == 0);
     fscanf(fp, "%lf", &A01);				//����The proportion of 0's and 1's which are blocked
	 cout<<str<<" == "<<A01<<endl;

	 /************* add for PPDM ***************************/

	 // ��param_file�ж�������Դ���������ļ�����suppor��ֵ

     fscanf(fp, "%s", str);
     assert(strcmp(str, "datafrom") == 0);
	 fscanf(fp, "%s", datafrom);
	 cout<<endl<<str<<" == "<<datafrom<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "sensfile") == 0);
	 fscanf(fp, "%s", sensfile);
	 cout<<str<<" == "<<sensfile<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "FIM_or_AR") == 0);
     fscanf(fp, "%u", &FIM_or_AR);
	 cout<<str<<" == "<<FIM_or_AR<<endl;	

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "SIZE_threshold") == 0);
     fscanf(fp, "%u", &SIZE_threshold);
	 cout<<str<<" == "<<SIZE_threshold<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "MST") == 0);
     fscanf(fp, "%lf", &MST);
	 cout<<str<<" == "<<MST<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "MCT") == 0);
     fscanf(fp, "%lf", &MCT);
	 cout<<str<<" == "<<MCT<<endl;


	 fscanf(fp, "%s", str);
     assert(strcmp(str, "MST_discount") == 0);
     fscanf(fp, "%lf", &MST_discount);
	 cout<<str<<" == "<<MST_discount<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "MCT_discount") == 0);
     fscanf(fp, "%lf", &MCT_discount);
	 cout<<str<<" == "<<MCT_discount<<endl;


		//fscanf(fp, "%s", str);
		//   assert(strcmp(str, "MAX_del_trans_perc") == 0);
		//   fscanf(fp, "%lf", &MAX_del_trans_perc);
		//cout<<str<<" == "<<MAX_del_trans_perc<<endl;


	 /*****************************************************/

     
     srand(seed); /* seeding random number generator */

     fclose(fp);

	 cout<<"\nEND to read local parameters!\n";

     return (0);
}


//****************************************************************************
//       ���ܣ� ����ÿ��Ƶ�����overlap degree�� д���ļ� 
//		 �ȼ����һ��Ƶ�����漰��������transactions��Ҳ��������Ƶ�������Ƶ�Σ�
//		 coverage(large_item_i) = ��Ƶ��/large_item_i�漰����transaction����
//****************************************************************************

void cal_overlap()
{
	//int count = 0;
	//int sum_item = 0;

	////double sel_percent=0.1;	//����Ϊȫ�ֱ�����ȷ���ڸ߳�sup_l�İٷ�֮���ٱ�����trans��������overlap,��Ϊ��ƵƵ������overlapӰ���
	//							//����ģ����ƵƵ����֮���overlap��𡣲���ѡ���ҳ���sup_l����ĸ�ƵƵ�������ʺ������������ɾ����trans̫��

	////map<item,double> item_overlap;
	////typedef map<item, double>::const_iterator itemsetCI_double;
	////typedef map<item, double>::value_type newdata_double;

	//cout<<"\n>>>Begin to calculate overlap degree for each large item set\n\n";

	//ofstream outfile;
	//char filename[80];
	//sprintf(filename, "%s_%.3f_overlap.txt",datafrom, MST);

	//outfile.open(filename, ios::out);



	//for (unsigned int i=0; i<large_sets.size(); i++)								// ͳ��Ƶ��������������reallarge��prelarge
	//{
	//	for (itemsetCI tempCI=large_sets[i].begin(); tempCI!=large_sets[i].end(); tempCI++)	 
	//	{
	//			if(tempCI->second< sup_u + (int)(db_size* sel_overlap_percent ) )
	//			{				  //ԭ����sup_l
	//				item goalitem = tempCI->first;
	//				sum_item = 0;     // ����ͳ�� overlap degree;

	//				vector<int> Titem;		//������¼ ֧�� ��ǰƵ���� �����ݿ��¼ ���

	//				for (unsigned int j=0; j<database.size();j++)					// �ҳ�֧�ֵ�ǰƵ�����trans,������Ŵ���Titem
	//				{
	//					transaction goaltrans = database[j];				
	//					itemCI itemCI = goalitem.begin();							//  itemΪset<int>, �����Ԫ���Զ���������
	//					transactionCI transCI = goaltrans.begin();	
	//				
	//					while (itemCI!=goalitem.end() && transCI!=goaltrans.end())  //  �����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���
	//																				//  ��set����set<int>���Զ�����Ĺ���,vector<int>û�У�����initial()�����˹�����
	//					{
	//						if (*transCI == *itemCI)
	//						{
	//							count++;
	//							itemCI++;
	//							transCI++;
	//						}
	//						else
	//						{
	//							transCI++;
	//						}
	//					}
	//					if (count == goalitem.size())								// if������ʾƥ��
	//					{
	//						Titem.push_back(j);				 
	//					}
	//					count = 0;
	//				}

	//				cout<<"�ҳ���Titem   Titem.size()="<<Titem.size()<<endl;

	//				for(unsigned int j=0;j<Titem.size();j++)						// ΪTitem��ÿ��trans, Ѱ��ƥ���large item, ÿ�ҵ�1����������sum_item��1
	//				{
	//					transaction goaltrans2 = database[ Titem[j] ];		
	//				
	//					for (unsigned int k=0; k<large_sets.size(); k++)						
	//					{
	//						for (itemsetCI tempCII=large_sets[k].begin(); tempCII!=large_sets[k].end(); tempCII++)
	//						{
	//								if(tempCII->second< sup_u + (int)(db_size* sel_overlap_percent ) )
	//								{
	//									int count2=0;
	//									if(tempCII->first == tempCI->first)   continue; 
	//									else
	//									{
	//										transactionCI transCII = goaltrans2.begin();
	//										item goalitem2 = tempCII->first;
	//										itemCI itemCII = goalitem2.begin();	

	//										while (itemCII!=goalitem2.end() && transCII!=goaltrans2.end())  //  �����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���
	//																			//  ��set����set<int>���Զ�����Ĺ���,vector<int>û�У�����initial()�����˹�����
	//										{
	//											if (*transCII == *itemCII)
	//											{
	//												count2++;
	//												itemCII++;
	//												transCII++;
	//											}
	//											else{transCII++;}
	//										}
	//										if (count2 == goalitem2.size())								// if������ʾƥ��
	//										{
	//											sum_item = sum_item + 1; 		 
	//										}
	//										count2 = 0;
	//									}
	//								}
	//					
	//						}
	//					}
	//				}

	//				//item_overlap.insert( newdata_double(tempCI->first, (double)(sum_item)/(double)(tempCI->second)  ) );
	//				//����ǰƵ���������overlapֵ������item_overlapӳ��

	//				{
	//							item goalitem3 = tempCI->first;		
	//							itemCI itemC3 = goalitem3.begin();	
	//							outfile<<"<";
	//							while( itemC3!=goalitem3.end() )
	//							{
	//								outfile<<" "<<*itemC3<<" ";
	//								itemC3++;
	//							}
	//							outfile<<"> :  ";
	//				}

	//			
	//				outfile<<" "<<tempCI->second<<" \t"<< (double)(sum_item)/(double)(tempCI->second) <<endl;

	//				Titem.clear();  //���Titem

	//				cout<<"...";
	//			}

	//	}
	//}


	//outfile.close();

	//cout<<"\n>>>End of calculating overlap degree. Writing into file! \n\n";
}



//*******************************************************************
//       ���ļ�д��ͳ����Ϣ�� ������ 
//		 database_size, MST, sup_u, sup_l,datafile, sen_file, 
//*******************************************************************

void write_stat()
{
	ofstream outfile;
	char filename[80];

	cout<<"\n>>> Begin to write the statistics information into file!\n";

	sprintf(filename,"%s_MST_%.3f_MCT_%.2f_MLT_%.2f_statistic.dat",datafrom, MST,MCT,MLT);
	outfile.open(filename, ios::out);

	outfile<<"��α���--ɾ��items\n\n";

	outfile<<"datafile:                 "<<datafrom<<endl;
	outfile<<"database size:            "<<db_size<<endl;
	outfile<<"Frequent Itemset(1) or Association Rule(2):"<<FIM_or_AR<<endl;;   //��FIM��ARѡ�񿪹ء������������
	outfile<<"Minimum support threshold (MST):      "<<MST<<endl;
	outfile<<"Minimum confidence threshold (MCT):   "<<MCT<<endl;				// �����������
	outfile<<"Minimum lift threshold (MLT):         "<<MLT<<endl;				// �����������

	outfile<<"sup_u =               "<<sup_u<<endl;
	outfile<<"sup_l =               "<<sup_l<<endl<<endl;
	outfile<<"MAX_del_trans =       "<<MAX_del_trans<<endl;
	outfile<<"MAX_del_trans_perc =  "<<MAX_del_trans_perc<<endl;
	outfile<<"ALPHA              =  "<<ALPHA<<endl;
	outfile<<"Avg_Dist_MST       =  "<<Avg_Dist_MST<<endl;

	outfile<<endl;
	outfile<<"Max generation               =   "<<maxgen<<endl; 
	outfile<<"Population size /alpha       =   "<<alpha<<endl; 
	outfile<<"Parent size / mu             =   "<<mu<<endl; 
	outfile<<"Children size /lambda        =   "<<lambda<<endl; 
	outfile<<"Dimension (objective space)  =   "<<dimension<<endl; 

	outfile<<endl;
	outfile<<"sel_overlap_percent =	 "<<sel_overlap_percent<<endl;
	char filename2[80];
	sprintf(filename2, "%s_%.3f_overlap.txt",datafrom, MST);
	outfile<<"overlap degree for each large item = "<<filename2<<endl; 

	char tempStr[128];
	if(0==mutation_type) strcpy(tempStr, "no mutation");
	else if(1== mutation_type) strcpy(tempStr,"one bit mutation");
	else if(2== mutation_type) strcpy(tempStr,"independent bit mutation"); 	 
	else if(3== mutation_type) strcpy(tempStr,"shuffle_random mutation"); 
	else strcpy(tempStr,"INVALID mutation type");

	outfile<<"mutation_type	= \t"<<tempStr<<endl; 

	if(0==recombination_type) strcpy(tempStr, "no recombination");
	else if(1== recombination_type) strcpy(tempStr,"one point crossover");
	else if(2== recombination_type) strcpy(tempStr,"uniform crossover");
	else if(3== recombination_type) strcpy(tempStr,"shuffle_random crossover");
	else strcpy(tempStr,"INVALID crossover type");
	//����ƵĽ������ӣ� ϴ�Ʒ�ʽ���� shuffle random��, ��֤�Ӵ��������ظ�

	outfile<<"recombination_type    =  "<<tempStr<<endl; 
	outfile<<"mutatation_prob       =  "<<mutat_prob<<endl;
	outfile<<"recombination_prob    =  "<<recom_prob<<endl; 
	outfile<<"bit_turn_prob         =  "<<bit_turn_prob<<endl; 

	outfile<<"Number of strong rules  =   "<<num_strong_rules<<endl; 
	outfile<<"Number of strong rules calculated by trie tree =   "<<num_strong_rules_by_trie<<endl; 

	outfile<<"\n/*-----------------------------------------------------------------------------*/\n\n";

	outfile<<"sensitive file: \t"<<sensfile<<endl;
	outfile<<"Count of sensitive item sets:"<<s_set.size()<<endl;
	outfile<<"sensitive item sets LIST:"<<endl;
	int j=0;
	for(vector<pair<item,item>>::const_iterator sen_CI=s_set_pair.begin();sen_CI!=s_set_pair.end();sen_CI++)				
	{
		item condition, consequence;

		condition = sen_CI->first;
		consequence = sen_CI->second;

		set<int>::const_iterator CI1,CI2;

		CI1=condition.begin();
		for(;CI1!=condition.end();CI1++)
		{
			outfile<< *CI1 <<" ";
		}

		if(AR==FIM_or_AR)
		{
			outfile<<"-> ";

			CI2=consequence.begin();
			for(;CI2!=consequence.end();CI2++)
			{
				outfile<< *CI2 <<" ";
			}
		}

		outfile<<"  : ("<<sen_frequency_array[j]<<")"<<endl;
		j++;
	}
	if(j!=s_set_pair.size()) {cout<<" ERROR in write_stat()!!! "; exit(-1);}
	outfile<<endl<<endl;

	outfile<<"\n/*-----------------------------------------------------------------------------*/\n\n";

    /*  //ɾ��item�������trans������ ��û��prelarge itemsets
	outfile<<"\nprelarge item sets LIST: \n";
	for (unsigned int i=0; i<prelarge_sets.size(); i++)												//����n+1Ƶ���
	{

		outfile<<"\nprelarge"<<i+1<<"-item sets LIST (num="<< prelarge_sets[i].size()<<") :\n";
		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	//����n+1Ƶ�����ÿһ��Ƶ����
		{
				item goalitem = tempCI->first;
				itemCI itemCI = goalitem.begin();	
				outfile<<"<";
				while(itemCI!=goalitem.end())
				{
					outfile<<" "<<*itemCI<<" ";
					itemCI++;
				}
				outfile<<"> : "<<(tempCI->second)<<endl;
		}
	}
	*/
	outfile<<"\n/*-----------------------------------------------------------------------------*/\n\n";
	
	outfile<<"The distance of sensitive or prelarge itemsets to absolute MST: \n";
	outfile<<"sum1_prelarge = "<<sum1_prelarge<<"     sum2_sensitive = "<<sum2_sensitive<<endl;
	outfile<<"The length of chromosome: "<<length<<"   length1 = "<<length1<<"  length2 = "<<length2<<endl;
	outfile<<"The actual sup_l = (db_size-length)*MST = "<<(db_size-length)*MST<<endl<<endl;

	int count=0;
	for (unsigned int i=0; i<large_sets.size(); i++)									//ͳ��Ƶ��������������reallarge��prelarge
	{
		for (itemsetCI tempCI=large_sets[i].begin(); tempCI!=large_sets[i].end(); tempCI++)	 
		{
			count++;
		}
	}
	outfile<<"\nThe total number of large item sets (including reallarge and prelarge): "<<count<<endl;
	outfile<<"\nDetailed list of large item sets is storded in another file!"<<endl;

	outfile.close();

	cout<<"\n>>> END of writing statistical information into file!\n";

}



/* Performs variation. */
int variate(int *selected, int *result_ids)
{
     int result, i, k;

     result = 1;

     /* copying all individuals from selected */
     for(i = 0; i < mu; i++)
     {
          result_ids[i] = add_individual(copy_individual(get_individual(selected[i])));  //����ֵΪ�¸��������Ⱥ��λ�ã����֤�ţ�

          if(result_ids[i] == -1)
          {
               log_to_file(log_file, __FILE__, __LINE__,
                           "copying + adding failed");
               return (1);
          }
     }
 
     /* if odd number of individuals, last one is
        left as is */
     if((((double)mu/2) - (int)(mu/2)) != 0) k = mu - 1; 
     else k = mu;

     /* do recombination */
     for(i = 0; i < k; i+= 2)
     {  
          result = 1;
          if (drand(1) <= recom_prob)
          {
               result = 1;
               if (recombination_type == 1)
               {
                    result = one_point_crossover(  get_individual(result_ids[i]),get_individual(result_ids[i + 1])  );
               }
               else if (recombination_type == 2)
               {
                    result = uniform_crossover(  get_individual(result_ids[i]), get_individual(result_ids[i + 1])  );
               }
			   else if (recombination_type == 3)   //����ƵĽ������ӣ� ϴ�Ʒ�ʽ���� shuffle random��, ��֤�Ӵ��������ظ�
			   {
					result = shuffle_crossover(  get_individual(result_ids[i]), get_individual(result_ids[i + 1])  );
			   }
               else if (recombination_type == 0)
                    result = 0;

               if (result != 0)
                    log_to_file(log_file, __FILE__, __LINE__, "recombination failed!");
          }
     }
     
     /* do mutation */
     for(i = 0; i < mu; i++)
     {
          result = 1;
          if (drand(1) <= mutat_prob) /* only mutate with mut.probability */
          { 
               if(mutation_type == 1)
               {
                    result = one_bit_mutation(get_individual(result_ids[i]));  //one_bit_mutation()���������ޣ� �Ƽ� indep_bit_mutation()
               }
               else if(mutation_type == 2)
               {
                    result = indep_bit_mutation(get_individual(result_ids[i]));
               }
			   else if(mutation_type == 3)
               {
                    result = shuffle_mutation(get_individual(result_ids[i]));
               }
               else if(mutation_type == 0)
               {
                    result = 0;
               }
    
               if(result_ids[0] == -1)
                    log_to_file(log_file, __FILE__, __LINE__,
                                "mutation failed!");
          }
     }
     
     /* do evaluation */
     for(i = 0; i < mu; i++)
     {
		int result;
		result = eval(get_individual(result_ids[i]));
     }
     
     return (0);
}



// ���PPDM����, mutation������Ҫ����ǿ�������޸ģ�

// �������ɵ�Ⱦɫ���еı�Ų��ظ��� (ͨ��<set>���Ϳ���ʵ��)
/* flip one bit at random position */
int one_bit_mutation(individual *ind)
{
     int position1, position2;

     if(ind == NULL)
          return (1);
     /* flip bit at position */

	position1 = irand(length1);
	ind->bit_string[position1] = prelarge_trans_array[ irand(prelarge_trans_array_size)] ;  

	 //ind->bit_string[length1+position2] = sen_trans_array[ irand(sen_trans_array_size)] ;  // ������֤��Ⱦɫ���������ı�Ų��ظ���

	 int offset=0;
	 for(int j=0;j<num_sen_item;j++)
	 {
			position2 = irand(encoding_len[j]);
			ind->bit_string[ length1+offset+position2 ] = sen_trans_array[j][ irand( sen_frequency_array[j]) ] ;	// ������֤��Ⱦɫ���������ı�Ų��ظ���
																										// ������֤��Ⱦɫ���������ı�Ų��ظ�	
			offset=offset+encoding_len[j];
	 }
											

     return (0);
}


// �������ɵ�Ⱦɫ���еı�Ų��ظ��� (ͨ��<set>���Ϳ���ʵ��)
/* flip all bits with certain probability */
int indep_bit_mutation(individual *ind)
{

     double probability;
     /* absolute probability */
     probability = bit_turn_prob;  //param_file��Ϊ0.01,�Ƿ�̫С�ˣ� ʵ��ȷ�����ȸĳ�0.1��
     if(ind == NULL)
          return (1);

	 //ǰ��Σ�for prelarge_item ����
	 random_shuffle(v1.begin(), v1.end());

     for(int i = 0; i < length1; i++)
     {
          if(drand(1) < probability)
          {
			  ind->bit_string[i] = prelarge_trans_array[ v1[i] ] ;			// ������֤��Ⱦɫ���������ı�Ų��ظ���random_shuffle()����													
          }
     }

	 int offset=0;
	 for(int j=0;j<num_sen_item;j++)
	 {
		for(int i = 0; i < encoding_len[j]; i++)
		{
			if(drand(1) < probability)
			{
				ind->bit_string[ length1+offset+i ] = sen_trans_array[j][ irand( sen_frequency_array[j]) ] ;	// ������֤��Ⱦɫ���������ı�Ų��ظ���																									// ������֤��Ⱦɫ���������ı�Ų��ظ�
			}
		}
		offset=offset+encoding_len[j];
	 }
     
     return (0);
}


// �������ɵ�Ⱦɫ���еı�Ų��ظ��� (ͨ��<set>���Ϳ���ʵ��)
/* do a one point crossover on ind1 and 2, the individual are
   overwritten! */
int one_point_crossover(individual *ind1, individual *ind2)
{
     int position1, i;
	 int *position2;
     int *bit_string_ind2;
     bit_string_ind2 = (int *) malloc(sizeof(int) * ind2->length);
  
     for(i = 0; i < ind2->length; i++)
     {
          bit_string_ind2[i] = ind2->bit_string[i];
     }

     position1 = irand(length1);
	 for(i=0;i< position1; i++)
	 {
          ind2->bit_string[i] = ind1->bit_string[i];		// �������ɵ�Ⱦɫ���еı�Ų��ظ���
          ind1->bit_string[i] = bit_string_ind2[i];   		// �������ɵ�Ⱦɫ���еı�Ų��ظ���		  
	 }

	 position2 = new int(num_sen_item);

	 for(int j=0;j<num_sen_item;j++)
		position2[j] = irand(encoding_len[j]);

	 int offset=0;
     for(int j=0;j<num_sen_item;j++) 
	 {
		for(i = length1+ offset ; i < length1+ offset+ position2[j]; i++) 
		{
			ind2->bit_string[i] = ind1->bit_string[i];		// �������ɵ�Ⱦɫ���еı�Ų��ظ���
			ind1->bit_string[i] = bit_string_ind2[i];   		// �������ɵ�Ⱦɫ���еı�Ų��ظ���
		}  

		offset = offset+encoding_len[j];
	 }

     free(bit_string_ind2);
	 free(position2);

     return(0);
}


// �������ɵ�Ⱦɫ���еı�Ų��ظ���
/* do a uniform crossover on ind1 and 2, the individual are
   overwritten! */
int uniform_crossover(individual *ind1, individual *ind2)
{

     int choose, i;
     int *bit_string_ind2;
     bit_string_ind2 = (int *) malloc(sizeof(int) * ind2->length);
  
     for(i = 0; i < length; i++)
     {
          bit_string_ind2[i] = ind2->bit_string[i];
     }

     for(i = 0; i < length; i++) {							// ����"���α���"�� uniform_crossover�ƺ�����Ҫ�޸�
          choose = irand(2);
          if(choose == 1) /* switch around bits */
          { 
               ind2->bit_string[i] = ind1->bit_string[i];			// �������ɵ�Ⱦɫ���еı�Ų��ظ��� (ͨ��<set>���Ϳ���ʵ��)
               ind1->bit_string[i] = bit_string_ind2[i];			// �������ɵ�Ⱦɫ���еı�Ų��ظ��� (ͨ��<set>���Ϳ���ʵ��)
          } /* else leave bit as is */   
     }  

     free(bit_string_ind2);

     return (0);
}



// ���ɵ��Ӵ�Ⱦɫ���е�geneȡֵ ���ظ�
// �µĽ������ӣ�����������
// �µĽ������ӣ�����������
int shuffle_crossover(individual *ind1, individual *ind2)
{
     int  i;

	 set<int> pre_combine_gene_set;
	 vector<int> pre_combine_gene_vec;

	 cout<<"\nshuffle_crossover\n";

	 //��ǰ��Ρ�Pre_large����ϲ�,��setת��vector
	 for(i=0;i< length1; i++)
	 {
		 pre_combine_gene_set.insert( ind1->bit_string[i] );
		 pre_combine_gene_set.insert( ind2->bit_string[i] );
	 }
	 set<int>::const_iterator setCI= pre_combine_gene_set.begin();
	 while(setCI != pre_combine_gene_set.end())
	 {
		 pre_combine_gene_vec.push_back(*setCI);
		 setCI++;
	 }
	 random_shuffle(pre_combine_gene_vec.begin(), pre_combine_gene_vec.end());
	 for(i=0;i< length1; i++)
	 {
		 ind1->bit_string[i] =  pre_combine_gene_vec[i];
	 }
	 random_shuffle(pre_combine_gene_vec.begin(), pre_combine_gene_vec.end());
	 for(i=0;i< length1; i++)
	 {
		 ind2->bit_string[i] =  pre_combine_gene_vec[i];
	 }
	 pre_combine_gene_set.clear();
	 pre_combine_gene_vec.clear();



	 //�����Ρ��е�ÿ��section�е�gene�ϲ�
	 int offset=0;
     for(int j=0;j<num_sen_item;j++) 
	 {
		set<int> sen_combine_gene_set;
		vector<int> sen_combine_gene_vec;

		for(i = length1+ offset ; i < length1+ offset + encoding_len[j]; i++) 
		{
			sen_combine_gene_set.insert(ind1->bit_string[i]);			//���������塰���Ρ���j������Ƭ���е�gene�ϲ�
			sen_combine_gene_set.insert(ind2->bit_string[i]);
		}  

		set<int>::const_iterator setCI2= sen_combine_gene_set.begin();	//��set����vector
		while(setCI2 != sen_combine_gene_set.end())
		{
			sen_combine_gene_vec.push_back(*setCI2);
			setCI2++;
		}

		random_shuffle(sen_combine_gene_vec.begin(), sen_combine_gene_vec.end());	// ���¶Ժϲ���gene����ϴ�� ����֤���ظ��� �ؼ����ڣ���

		for(i = length1+ offset ; i < length1+ offset + encoding_len[j]; i++) 
		{
			ind1->bit_string[i] = sen_combine_gene_vec[ i-(length1+ offset) ];
		}

		random_shuffle(sen_combine_gene_vec.begin(), sen_combine_gene_vec.end());	// �ٴζԺϲ���gene����ϴ�� ����֤���ظ��� �ؼ����ڣ���

		for(i = length1+ offset ; i < length1+ offset + encoding_len[j]; i++) 
		{
			ind2->bit_string[i] = sen_combine_gene_vec[ i-(length1+ offset) ];
		}

		sen_combine_gene_set.clear();
		sen_combine_gene_vec.clear();

		offset = offset+encoding_len[j];
	 }



     return (0);
}



// ���ɵ��Ӵ�Ⱦɫ���е�geneȡֵ ���ظ�
// �µı������ӣ�����������
// �µı������ӣ�����������
int shuffle_mutation(individual *ind)
{
	cout<<"shuffle_mutation"<<endl;
	//ǰ��α���
	random_shuffle(v1.begin(), v1.end());

	for(int i=0;i<length1;i++)
	{
		ind->bit_string[i] = prelarge_trans_array[ v1[i] ];
	}

	//���α���
		
	for(unsigned int i=0;i<s_set.size(); i++)
	{
		random_shuffle( v2[i].begin(),v2[i].end() );
	}

	int offset=0;
	for ( int i=0; i < num_sen_item; i++)
	{
		 for( int j=0; j < encoding_len[i]; j++ )
		 {
				ind->bit_string[ length1 + offset + j] = sen_trans_array[i][ v2[i][j] ] ;   
													
		 }
		 offset = offset + encoding_len[i];	
	}
	cout<<endl;

	if(offset!=length2){cout<<"\n ERROR! in    int shuffle_mutation(individual *ind):  offset!=length2 "; getchar(); exit(-1);}



	return 0;
}




// �����������������д��Ľ���
// �����������������д��Ľ���
// �����������������д��Ľ���
/* Generate a random integer. */
//�������������Ҫ
int irand(int range)
{
     int j;

	 //Paramfile�ļ����Ѿ��������ӣ�

	 //srand((unsigned)time(NULL));		//�Ե�ǰʱ����Ϊ�������

     j=(int) ((double) range * (double) rand() / (RAND_MAX + 1.0)); // RAND_MAX��32767

     return (j);
}


/* Generate a random double. */
double drand(double range)
{
     double j;
     j=(range * (double) rand() / (RAND_MAX + 1.0));
     return (j);
}




//*******************************************************************
//       ����length(Ⱦɫ����볤��),����������Ƶ�κ�sup_l 
//		 length = ( SIGMA ( frequency_i - sup_l) )*2
//*******************************************************************

void cal_length()
{
	//int sum1_prelarge =0,  sum2_sensitive = 0;
	//int length1=0,length2=0;							// ����ȫ�ֱ�������Ϊlength1��ֵ�����ˡ����α��롱�ķֶ�λ��

	cout<<">>>Begin to calculate length!\n";

	encoding_len = new int[ s_set.size() ];				// ���ڴ�š���α��롱��ÿ�α���ĳ��ȣ�һ�α����Ӧһ��sensitive item.

	//�ȼ��㣺 ÿ��sensitive item��Ӧ�ı��볤��
	for(unsigned int i=0;i<s_set.size();i++)
	{
		if(sen_frequency_array[i]>sup_u)
		{
			encoding_len[i] = sen_frequency_array[i]-sup_l + 1;
			sum2_sensitive = sum2_sensitive + encoding_len[i] ;
		}
		else
		{
			cout<<"\nERROR in void cal_length(): ������"<<i<<"����Ƶ����"; getchar(); exit(-1);
		}
	}

	length2 = sum2_sensitive;

	if( MAX_del_trans < (unsigned)length2)  {
		cout<<"MAX_del_transȡֵ ��С �� ��������볤�ȳ�����\n";
		cout<<"MAX_del_trans = "<<MAX_del_trans <<"    length2 = "<<length2<<endl;

		write_stat();

		getchar();

		length2 = MAX_del_trans ;  

		exit(-1);
	}


	//ǰ��α��볤�� //ɾ��item���� �Ͳ���Ҫ�˲����ˣ�
	/*
	for (unsigned int i=0; i<prelarge_sets.size(); i++)												//����n+1Ƶ���
	{
		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	//����n+1Ƶ�����ÿһ��Ƶ����
		{
			if(tempCI->second>=sup_l && tempCI->second<sup_u)
				sum1_prelarge = sum1_prelarge + (tempCI->second - sup_l) + 1;
			else{
				cout<< "prelarge��Ƶ�� Ӧλ�� [sup_l,sup_u)���䡣 ������"; getchar(); exit(-1);
			}
		}
	}
	*/


	/*
	if (sum1*1.5 < (double)prelarge_trans_array_size)			length1 = (int)(sum1*1.5);
	else if (sum1*1.3 < (double)prelarge_trans_array_size)  	length1 = (int)(sum1*1.3);
	else if (sum1*1.1 < (double)prelarge_trans_array_size)		length1 = (int)(sum1*1.1);
	else {
		cout<<"\nMST ���ù�С\n"; 
		getchar(); exit(-1);
	}
	*/

	/*
	if (sum1_prelarge  <= (double)prelarge_trans_array_size && sum1_prelarge<= MAX_del_trans - length2 )			
		length1 = (int)(sum1_prelarge);
	else {

		cout<<"\n pre_large�εı��볤�� Խ�磬 ��������ֵ��" ;

		int temp_min =  ( (MAX_del_trans - length2)< prelarge_trans_array_size) ? (MAX_del_trans - length2): prelarge_trans_array_size;

		cout<<temp_min<<endl;

		length1 = temp_min;

		//getchar(); 
		
		//exit(-1);
	}
	*/
	sum1_prelarge = 0;
	length1 = 0;

	length = length1 + length2;   // Ⱦɫ������ܳ��ȣ� 3������ȫ�ֱ���

	cout<<"\nprelarge���ֵsup_l Ƶ���ܼ� sum1="<<sum1_prelarge<<"  length1="<<length1<<endl;
	cout<<"�������ֵsup_l Ƶ���ܼ�     sum2="<<sum2_sensitive<< "  length2="<<length2<<endl;
	cout<<"�����ܳ���  length = length1+length2 ="<<length<<endl;
	cout<<"\nactual_sup_l = "<<(db_size*MST)<<"\n";
	cout<<"original sup_l = "<<sup_l<<endl;

	sup_l_actual = (int)(db_size*MST);


	getchar();

	cout<<"\n>>>END of calculating length!\n";
}




//*************************************************************************************
//       ���ܣ����˳�����pre_large_set��item�ļ�¼; ����ÿ��pre_large��Ƶ��
//
//		 ���: ���� ������pre_large���¼�ļ�¼�š���prelarge_trans_array[]��prelarge_trans_set,
//*************************************************************************************

/*
       ֻ����� pre_large 1-itemset, ���������pre_large 2-item, 
	   ��Ϊ����pre_large 2-item ��transactions һ��Ҳ��������2-item��1-item.

*/
void filter_prelarge_trans()	//ʹ���ˡ��������ݿ�
{
	//apriori �㷨 �õ��� �����: large_sets, reallarge_sets, prelarge_sets

	int count =0 ;

	cout<<"\n\nBEGIN:  void filter_prelarge_trans()\n";

	cout<<"\nprelarge_sets.size()="<<prelarge_sets.size()<<"\n";

	for (unsigned int i=0; i<prelarge_sets.size(); i++)												//����n+1Ƶ���	
	{
		if( prelarge_trans_set.size() >= db_size/5 ) break;


		cout<<"\n--------prelarge_sets["<<i<<"].size = "<<prelarge_sets[i].size()<<"\n";

		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	//����n+1Ƶ�����ÿһ��Ƶ����
		{
				if( prelarge_trans_set.size() >= db_size/4 ) break;

				item goalitem = tempCI->first;

				for(unsigned int j=0;j<thin_DB.size();j++)// thin_DB �������ݿ�
				{
					transaction goaltrans = thin_DB[j];   // thin_DB �������ݿ�

					itemCI itemCI = goalitem.begin();										//itemΪset<int>, �����Ԫ���Զ���������
					transactionCI transCI = goaltrans.begin();	

					while (itemCI!=goalitem.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
					{
						if (*transCI == (copy_new_code[*itemCI]-1) ){  // ��ע�⡿��copy_new_code[]�����롱
							count++;
							itemCI++;
							transCI++;
						}
						else { transCI++;}
					}
					if (count == goalitem.size() )
					{
						//cout<<"ƥ����һ�� pre_large �� "<< j  ;
						prelarge_trans_set.insert(j);										//ʹ��set���Ϳ��Է�ֹ �����ظ���¼��
						//cout<<"   ��ǰprelarge trans ����Ϊ: "<< prelarge_trans_set.size()<<endl;
					}
					count = 0;
				}

		}
	}

	//��prelarge��¼ �ɼ���prelarge_trans_setת�浽����prelarge_trans_array[]
	prelarge_trans_array_size = prelarge_trans_set.size();
	prelarge_trans_array = (int *)new int[prelarge_trans_array_size];

	int k=0;
	set<int>::const_iterator  prelarge_trans_setCI=prelarge_trans_set.begin();
	while(prelarge_trans_setCI!=prelarge_trans_set.end())
	{
		prelarge_trans_array[k] = *prelarge_trans_setCI;
		k++;
		prelarge_trans_setCI++;
	}
	if(k!=prelarge_trans_array_size)		//�������Ƿ����
	{
		cout<<endl<<"�漰prelarge�� �� transaction�������󣡣���"; getchar(); exit(-1);
	}

	//getchar();

	// ��ʾ����prelarge�Ƶ�Σ� �漰������� ���ݿ��¼��
	cout<<"\nprelarge������:  "<<  prelarge_trans_set.size()<<endl;
	for (unsigned int i=0; i<prelarge_sets.size(); i++)										
	{
		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	
		{
			cout<< "prelarge�<";
			item goalitem = (tempCI->first) ;
			itemCI itemCI = goalitem.begin();
			while(itemCI!=goalitem.end())
			{
				cout<<" "<<*itemCI<<" ";
				itemCI++;
			}
			cout<<">    frquencyƵ�Σ�"<< (tempCI->second)<<endl;
			goalitem.clear();
		}
	}
		
	/*
	cout<<"\n��ʾ����prelarge������ݿ��¼�ţ�"<<endl<<"������Ϊ����"<< prelarge_trans_set.size() <<endl;
	for(unsigned int j=0;j<prelarge_trans_set.size();j++)
	{
		if(j%20==0) cout<<endl;
		cout<<" "<<prelarge_trans_array[j];
	}
	*/
	cout<<"\n\nEND:  void filter_prelarge_trans()\n\n";

	//getchar();

	//getchar();

}



//**************************************************************************************
//       ���ܣ����˳������������ transactions������ÿ��������Ƶ��  
//
//		 ���:  ���� �������������¼�ļ�¼�š� sen_trans_array[]��sen_trans_set, 
//				�Լ������˸����������Ƶ�ε�����sen_frequency_array[]
//**************************************************************************************


void filter_sen_trans()
{													
	int count=0;

	vector<set<int>> sen_trans_vec;			// �������漰��trans���ϣ� ʹ��set<int>���Ϳ��Է�ֹ �����ظ���Ԫ��
											// ÿ���������Ӧһ��sensitive transaction����

	sen_frequency_array=(int *)new int[s_set_pair.size()];								//sen_frequency_array[]����ÿ��������ĳ���Ƶ��

	cout<<"\n=======BEGIN to filter out transactions containing sensitive rules or itemsets=======\n";

	//ģʽƥ�䣬 �ҳ����м�¼����¼�Ŵ��뼯��sen_trans_set 
	//����ÿ��������ĳ���Ƶ�Σ�����sen_frequency_array[]
	for(unsigned int j=0;j<s_set_pair.size();j++)				
	{																				//�������м���(vector)�е�ÿһ�У�ÿһ����һ��item(��set<int>����)	
		int sen_frequency_count=0;
		set<int>  sen_trans_set;

		item goalitem;
		pair<item,item>  goalitem_pair = s_set_pair[j];
		set<int>::const_iterator it;
		for( it=goalitem_pair.first.begin(); it!=goalitem_pair.first.end(); it++)
		{
			goalitem.insert(*it);
		}
		for( it=goalitem_pair.second.begin(); it!=goalitem_pair.second.end();it++)
		{
			goalitem.insert(*it);
		}

		//s_set.push_back(goalitem);      //����condition ��consequence�ϲ����� �������s_set
										// s_set�����ǣ�  vector<set<int>> s_set;	
		item goalitem_encoding;
		basket_recode3(goalitem, goalitem_encoding);		//��goalitem��Ϊʹ���ڲ����ʾ��trie�����������ݿⶼʹ����Ļ�룩

		//for(unsigned int i=0;i<database.size(); i++)				
		for(unsigned int i=0;i<thin_DB.size(); i++)				//ʹ�á��������ݿ�(thin_DB��ÿ�о�����������)
		{
			//transaction goaltrans = database[i];
			transaction goaltrans = thin_DB[i];					//ʹ�á��������ݿ�			

			itemCI itemCI = goalitem_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
			transactionCI transCI = goaltrans.begin();	

			while (itemCI!=goalitem_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )	//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
				//if ( *transCI  == *itemCI )//ԭ���ݿ�
				if( *transCI  == *itemCI )   //������ʽһ���� ����������ʹ�õ����ڲ���
				{
					count++;
					itemCI++;
					transCI++;
				}
				else { transCI++;}
			}
			if (count == goalitem_encoding.size() && goalitem.size()!=0)
			{
				//cout<<"������ƥ����һ��Trans:  "<< i  ;
				sen_trans_set.insert(i);											//ʹ��set���Ϳ��Է�ֹ �����ظ���¼��
				
				//cout<<"   ��ǰsensitive trans ����Ϊ: "<< sen_trans_set.size()<<endl;
				sen_frequency_count++;
			}
			count = 0;
		}

		sen_trans_vec.push_back(sen_trans_set);			// ��֧�� ��j�������� �� sensitive transactions ID���� ����sen_trans_vec

		sen_frequency_array[j]=sen_trans_set.size();	// ��֧�� ��j�������� �� transaction���� ����ȫ������ sen_frequency_array[ ]

		cout<<"\n\n>>>��"<<j<<"���������Ƶ��Ϊ��"<<sen_frequency_count<< "  ��ǰsen_trans_set.size()="<<sen_trans_set.size()<<"\n\n";

		sen_trans_set.clear();
	}

	if(s_set_pair.size()!=s_set.size()) {cout<<"\n ������������������� ��С��һ�£�\n"; getchar();exit(-1);}

	for(unsigned int i=0;i<s_set_pair.size();i++)		// �������
	{
		if(sen_trans_vec[i].size()!=sen_frequency_array[i])		{cout<<"\nERROR 2 in void filter_sen_trans()"; getchar();exit(-1);}
		if(sen_frequency_array[i]<sup_u)	
		{
			cout<<"sen_frequency_array[i]="<<sen_frequency_array[i]<<"  sup_u="<<sup_u;
			cout<<"\nERROR 3 in void filter_sen_trans():������"<<i<<"����Ƶ����. ����sensitive data file�Ƿ��������"; 
			getchar();exit(-1);
		}
	}

	//�����м�¼ �ɼ���sen_trans_setת�浽����sen_trans_array[]
	int k=0;

	//int **sen_trans_array;			//��ӵ�����ֲ��������պ����ˣ� ��ϸ���������κ�������
	sen_trans_array_size = s_set_pair.size();
	sen_trans_array = (int **)new int*[sen_trans_array_size];

	for(unsigned int i=0;i<s_set_pair.size();i++)
	{
		sen_trans_array[i]=(int *)new int[ sen_frequency_array[i] ];
	}
	for(unsigned int i=0;i<s_set_pair.size();i++)
	{
		unsigned int j=0;

		itemCI sen_trans_setCI=sen_trans_vec[i].begin();
		while(sen_trans_setCI!=sen_trans_vec[i].end())
		{
			sen_trans_array[i][j]=*sen_trans_setCI;			//  �õ�һ����Ҫ���ݽṹ�� ��ά��̬���飨ÿ�в��ȳ���sen_trans_array[i][j]
			sen_trans_setCI++;								//  ÿ�ж�Ӧһ�����й��� ÿ�������Ƕ�Ӧ���й����sensitive transactions ID
			j++;
		}
		if( j!= sen_frequency_array[i] )  {cout<<"\nsen_trans_array[][]ÿ�еļ�����������\n "; getchar(); exit(-1);}
	}



	// ��ʾ���������Ƶ�Σ� �漰������� ���ݿ��¼��

	cout<<"\n����������:  "<<  s_set_pair.size()<<endl;
	for(unsigned int j=0;j<s_set_pair.size();j++)
	{
		itemCI tmp_itemCI = s_set_pair[j].first.begin();
		for (; tmp_itemCI != s_set_pair[j].first.end(); tmp_itemCI++ )
		{
			cout<< *tmp_itemCI <<" ";
		}

		if(!s_set_pair[j].second.empty())
		{
			cout<<" -> ";
			tmp_itemCI = s_set_pair[j].second.begin();
			for(;tmp_itemCI!=s_set_pair[j].second.end(); tmp_itemCI++ )
			{
				cout<< *tmp_itemCI <<"  ";
			}
		}
		cout<<  "  Ƶ��Ϊ�� "<< sen_frequency_array[j]<<endl;
	}


	int total_sen_trans_num=0;
	for(unsigned int i=0;i<s_set_pair.size();i++)
	{
		cout<<"\n\n��"<<i<<"���������漰��trans������"<<sen_trans_vec[i].size();
		cout<<"\n��"<<i<<"���������漰��trans���£�\n";

		for(unsigned int j=0;j<sen_trans_vec[i].size();j++)
		{		
			if(j%20==0) cout<<endl;
			cout<<" "<<sen_trans_array[i][j];
			total_sen_trans_num++;
		}
	}
	cout<<"\n\nSensitive trans����Ϊ��"<<total_sen_trans_num<<endl;

	cout<<"\n=======END of sensitive_trans filtering!=================\n";

	//getchar();
}


//****************************************************************
//			���������ļ�����ȡ��vector<item>���͵� s_set��
//****************************************************************/
void read_sensitive_item()
{
	ifstream sensitive;
	string line;
	string token;
	int data = 0;
	set<int> condition, consequence;

	cout<<"\n>>>Begin to read sensitive item sets from file!\n";

	//sensitive.open("sens.txt", ios::in);
	sensitive.open (sensfile, ios::in);

	if (!sensitive.is_open())
		cout << "�޷���dataset\n";	 

	while (!sensitive.eof())
	{
		getline(sensitive, line);					//��ȡһ�У� ע���ļ���ʹ�õ��ǡ���Ʒ��š���ʾ��Ʒ
		stringstream ss(line);	

		int consequence_flag = 0;

		while (ss >> token)
		{
			if(token.compare("->") == 0 )
			{
				consequence_flag = 1;
				cout<<"\n������ ���й��� A->B  \n";
				continue;
			}
			if(0==consequence_flag)
			{
				data = atoi(token.c_str());	 
				condition.insert(data);
			}
			else
			{
				data = atoi(token.c_str());	 
				consequence.insert(data);
			}
		}

		//vector<pair<item,item>> s_set_pair;	

		s_set_pair.push_back( pair<item,item>(condition, consequence) );

		condition.clear();
		consequence.clear();
	}

	/*=========================================
	[0](  (234) -> (258) )
	[1]( ��245,247��-> (277) )
	[2]( ��281��-> (190)  )
	===========================================*/

	//��һ���ҵ���������itemsets����������itemsets.�磺�û�ָ��AB���У���ABC,ABD������
	if( FIM_or_AR == FIM )		
	{
		set<set<int>> all_include_condition;

		vector<pair<item,item>>::const_iterator sen_CI;
		for( sen_CI = s_set_pair.begin(); sen_CI != s_set_pair.end(); sen_CI++ )
		{
			condition.clear();
			consequence.clear();
			basket_recode3(sen_CI->first, condition); //����
			if(sen_CI->first.size()!=condition.size())  {cout<<"read_sensitive_item()�����������Ƶ����"; getchar(); exit(-1);}
				
			p_apriori->get_pointer_to_Apriori_Trie()->find_all_include_itemset(condition, all_include_condition, sup_u);

			//���ڵõ���all_include_condition��ʹ���ڲ�����ģ� ��Ҫ���룡

			cout<<"all_include_condition.size() = "<<all_include_condition.size()<<endl;
		}

		for(set<set<int>>::const_iterator it=all_include_condition.begin(); it!=all_include_condition.end(); it++ )
		{   
					set<int> temp_itemset = *it;
					set<int> decode_itemset;
					basket_decode3( temp_itemset, decode_itemset );

					s_set_pair.push_back( pair<item,item>( decode_itemset, consequence ));
					//decode_itemset.clear();
		}
	}

	//if(1==gen)			// �ڵ�һ����֤�������Ƿ����ɹ�
	{
		cout<<"\n��֤sensative items sets read from file! \n";
		
		vector<pair<item,item>>::const_iterator sen_CI;

		for( sen_CI = s_set_pair.begin(); sen_CI != s_set_pair.end(); sen_CI++ )
		{
			condition.clear();
			consequence.clear();

			condition = sen_CI->first;
			consequence = sen_CI->second;

			cout<<"( ";
			set<int>::const_iterator it1,it2;
			for(it1 = condition.begin(); it1!=condition.end();it1++)
				cout<<*it1<<" ";
			cout<<") ";

			if(consequence.empty() && FIM_or_AR == AR)
			{
				cout<<"����������ò��� �������� �� FIM_or_AR == AR\n";
				getchar();
				exit(-1);
			}

			if(!consequence.empty())
			{
				cout<<"-> ( ";
				for(it2 = consequence.begin(); it2!=consequence.end(); it2++)
					cout<<*it2<<" ";
				cout<<")";
			}
			cout<<endl;

		}
	}

	//��s_set_pair�е�condition��consequence�ϲ��� ����s_set
	for(unsigned int j=0;j<s_set_pair.size();j++)				
	{															
		item goalitem;
		pair<item,item>  goalitem_pair = s_set_pair[j];
		set<int>::const_iterator it;
		for( it=goalitem_pair.first.begin(); it!=goalitem_pair.first.end(); it++)
		{
			goalitem.insert(*it);
		}
		for( it=goalitem_pair.second.begin(); it!=goalitem_pair.second.end();it++)
		{
			goalitem.insert(*it);
		}

		s_set.push_back(goalitem);      //����condition ��consequence�ϲ����� �������s_set
										// s_set�����ǣ�  vector<set<int>> s_set;	
	}									// s_set_pair���ͣ� vector<pair<item,item>>


	num_sen_item = s_set_pair.size();
	
	if(s_set_pair.size()!=s_set.size()){ cout<<"s_set_pair.size()!=s_set.size() in void read_sensitive_item()"<<endl; getchar(); exit(-1); }

	sensitive.close();

	cout<<"\n>>>END of reading sensitive item set from file!\n";
}


/* �⴦���㷨�޸ĵĹؼ����� */
/* ���뵽�ο�����룬 ����ƫ����Ϣ*/
/* �㾦֮�ʣ���������*/

int eval(individual *x)
{
	if(SIZE_threshold <= 2)
	{
		if(AR == FIM_or_AR)
		{
			cal_fitness_two_rules(x);
		}
		else
		{
			cal_fitness_two(x);
		}
	}	
	else
	{
		if( FIM_or_AR == FIM )
			cal_fitness_by_trie_FIM(x);   //����һ��itemset������items�ж��ٸ�����ͨ�á� ����trie������
		else if( FIM_or_AR == AR )
		{
			cal_fitness_by_trie_AR(x);    //����һ��rules������items�ж��ٸ�����ͨ�á� ����trie������
		}
		else
		{
			cout<<"\n������eval(individual *x)----- FIM_or_AR�������ó��� ��ѡ��Ϊ{1,2}��\n";
			getchar();
			exit(-1);
		}	
	}

	//��ע�⡿���Ѿ���֤ cal_fitness_by_trie_AR(x)��cal_fitness_two_rules(x)���н����һ�µģ� ��2-item rules��

	return 0;


    //�޸�void Apriori_Trie::assoc_rule_find( ���� ���������򱣴浽 map<pair<item,item>,vector<int>> ruleset_bodon;

	//const Trie* Trie::is_included() ���Դ����� ������ һ��itemset ��Ƶ�Σ�support���� �÷�Trie::is_included()->counter;  
	//�ο���void Apriori_Trie::assoc_rule_find()�������÷��� main_trie.is_included( condition_part, condition_part.begin() ) ->counter * min_conf
	
	//void Apriori_Trie::assoc_rule_assist(����������α���һ������ �����в�ε�itemsetȡ������

}

bool is_sensitive_rule(set<int> &condition, set<int> &consequence)
{
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		if(s_set_pair[j].first == condition && s_set_pair[j].second == consequence)
		{
			cout<<"\n��ǰ����һ��������";
			return true;
		}
	}
	return false;
}



//��ͨ�á�, ����itemset������item�ж��ٸ�   �����ʣ� �������͵ݹ顿
void cal_fitness_by_trie_FIM( individual *x )   
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;

	set<int> set_transID_delete;       // ����set�������ͣ� ��֤�������ı��Ψһ�����Զ�����, �������ȥ���ظ�����

	//����chromosome���� ��ѡ�е�TID����trans��
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //��Ϊtrans��set<int>���ͣ�set�����Զ��������򣩣�����Ҫ���ǲ������ظ�Ԫ��
		}
		else
		{
			printf("ERROR in function eval() : ");		exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();


	/*----------���ѡ�д�ɾ����ÿ��trans������main_trie_copy���ж�Ӧ����count-------------------*/

	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> goaltrans = thin_DB[*CI];		// ��"����" ���ݿ� ��ȡ����Ӧtransaction; 
													// ���ڡ���α��롱,��Ϊѡ�е�trans�������ˣ��������й��� ���goaltrans.size()>=2; ���ڲ������˵ġ�һ����롱������ѡ�е�����transΪ��trans
		p_main_trie_copy->update_trie_tree(goaltrans.begin(), goaltrans.end());
	}

	//cout<<"\n\n��ӡ���µĿ�������\n\n";
	//p_main_trie_copy->show_content_preorder();//��ӡ��:����һ��basket,�����ṹ�в��ҿ���ƥ�������"candidate_size -��ѡ��"�� ������counterֵ


	/*----------����main_trie_copy(���µ�copy��),����alpha ----------------------------*/

	// ����alpha, ��ÿ������itemsets, �������� �����counterֵ�Ƿ� >= (unsigned)sup_l_actual
	for(vector<item>::const_iterator sen_CI = s_set.begin(); sen_CI!=s_set.end();sen_CI++)  //���ÿ�� ����itemset�Ƿ����ص���<sup_l_actual�� 
	{
		item sen_itemset = *sen_CI;
		item sen_itemset_recoded; 

		basket_recode3(sen_itemset, sen_itemset_recoded);//����
		unsigned long sen_counter = p_main_trie_copy->get_itemset_counter_from_trie(sen_itemset_recoded);		//get_itemset_counter_from_trie()�����Լ�д�ĺ���

		cout<<"sen_counter "<<sen_counter<<endl;
		if(sen_counter >= (unsigned long)sup_l_actual)
			alpha++;
	}


	/*----------�Ա�main_trie��(ԭ��)�� main_trie_copy(���µ�copy��), ����beta, gamma--------*/
	//Trie * p_main_trie = p_apriori->get_pointer_to_Apriori_Trie()->get_pointer_to_main_trie();  //���ָ��main_trie������ָ��

	p_apriori->get_pointer_to_Apriori_Trie()->compare_two_trees(p_main_trie_copy, beta, gamma, sup_u, sup_l_actual);

	beta = beta - (s_set_pair.size()-alpha);  //compare_two_trees()��õ�beta������ԭ��Ƶ�δ��ڵ���sup_u,����С��sup_l_actual��Ƶ����
											  //�����˱����ص������ ������ټ�ȥ �����˵�����������

	/*---------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
	}
	else
	{
		cout<<"\nERROR in dimension in FUNC eval()!\n";
	}

	printf("\n");
	printf("alpha=%d  beta=%d  gamma=%d  |X|= %u  [Objectives]: ", alpha, beta, gamma, actual_num_del );

	for(int i=0;i<dimension;i++)
		printf(" %.2lf  ", x->objective_value[i]);


	cout<<"\n\n";


	// ��main_trie��main_trie_copyͬ���� ʹmain_trie_copy�ָ�����ǰ��ֵ
	p_apriori->get_pointer_to_Apriori_Trie()->reset_tree_to_same(p_main_trie_copy);


}




//��ͨ�á�, ����rule������item�ж��ٸ�   �����ʣ� �������͵ݹ顿

void cal_fitness_by_trie_AR( individual *x )  
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;

	set<int> set_transID_delete;       // ����set�������ͣ� ��֤�������ı��Ψһ�����Զ�����, �������ȥ���ظ�����

	//����chromosome���� ��ѡ�е�TID����trans��
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //set�����Զ��������򣩣�����Ҫ���ǲ������ظ�Ԫ��
		}
		else
		{
			printf("ERROR in function eval() : ");	getchar();	exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();


	/*----------���ѡ�д�ɾ����ÿ��trans������main_trie_copy���ж�Ӧ����count-------------------*/

	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> goaltrans = thin_DB[*CI];		// ��"����" ���ݿ� ��ȡ����Ӧtransaction;
													// ���ڡ���α��롱,��Ϊѡ�е�trans�������ˣ��������й��� ���goaltrans.size()>=2; ���ڲ������˵ġ�һ����롱������ѡ�е�����transΪ��trans
		p_main_trie_copy->update_trie_tree(goaltrans.begin(), goaltrans.end());
	}

	//cout<<"\n\n��ӡ���µĿ�������\n\n";
	//p_main_trie_copy->show_content_preorder();//��ӡ��:����һ��basket,�����ṹ�в��ҿ���ƥ�������"candidate_size -��ѡ��"�� ������counterֵ

	/*----------------------------------------------------------------------------------------------*/
	// ͳ��alpha
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		item sen_condition = s_set_pair[j].first;
		item sen_consequence = s_set_pair[j].second;
		item sen_union = s_set[j];

		item sen_condition_encode, sen_consequence_encode, sen_union_encode;

		basket_recode3(sen_condition, sen_condition_encode);//����
		basket_recode3(sen_consequence, sen_consequence_encode);//����	
		basket_recode3(sen_union, sen_union_encode);//����

		unsigned long sen_new_condition_sup   =  p_main_trie_copy->get_itemset_counter_from_trie(sen_condition_encode) ;   
		unsigned long sen_new_consequence_sup =  p_main_trie_copy->get_itemset_counter_from_trie(sen_consequence_encode)  ;  
		unsigned long sen_new_union_sup =  p_main_trie_copy->get_itemset_counter_from_trie(sen_union_encode) ;

		double sen_new_union_sup_ratio =  (double)( sen_new_union_sup)/(db_size-actual_num_del);
		double sen_new_rule_confidence =  (double)( sen_new_union_sup)/( sen_new_condition_sup );
		double sen_new_rule_lift = (double)(sen_new_union_sup)/(sen_new_condition_sup)*(double)(db_size-actual_num_del)/sen_new_consequence_sup;

		if(sen_new_union_sup_ratio >= MST && sen_new_rule_confidence >=MCT && sen_new_rule_lift >= MLT)
		{	
			alpha++;	cout<<"\nalpha is incereasing for sensitive rules not hidden!   alpha="<<alpha;
		}
	}

	/*----------------------------------------------------------------------------------------------*/
	// ͳ��beta��gamma
			
	map<pair<item,item>,vector<int>>::const_iterator CI = bodon_ruleset.begin();  //��ע�⡿��bodon_rulesetʹ�õ����ڲ�����
	for(; CI != bodon_ruleset.end(); CI++)
	{
			item condition = CI->first.first;			//ע�⣺�������itemʹ������.item����set<int>����. 
			item consequence = CI->first.second;		//ע�⣺�������itemʹ������
			item union_two;

			for(item::const_iterator it_condition = condition.begin() ; it_condition!= condition.end(); it_condition++)
			{
				union_two.insert(*it_condition);
				//cout<< *it_condition << " ";
			}
			//cout<<" -> ";
			for(item::const_iterator it_consequence = consequence.begin() ; it_consequence!= consequence.end(); it_consequence++)
			{
				union_two.insert(*it_consequence);
				//cout<<*it_consequence<<" ";
			}
			//cout<<endl;
			double old_union_sup_ratio = (double)(CI->second[2])/db_size;
			double old_rule_confidence = (double)(CI->second[2])/(CI->second[0]);
			double old_rule_lift = (double)(CI->second[2])/(CI->second[0])*(double)db_size/CI->second[1];

			unsigned long new_condition_sup   =  p_main_trie_copy->get_itemset_counter_from_trie(condition) ;   
			unsigned long new_consequence_sup =  p_main_trie_copy->get_itemset_counter_from_trie(consequence) ;  
			unsigned long new_union_sup = p_main_trie_copy->get_itemset_counter_from_trie(union_two) ;

			double new_union_sup_ratio =  (double)( new_union_sup)/(db_size-actual_num_del);
			double new_rule_confidence =  (double)( new_union_sup)/( new_condition_sup );
			double new_rule_lift = (double)(new_union_sup)/(new_condition_sup)*(double)(db_size-actual_num_del)/new_consequence_sup;

			//cout<<"old_sup="<< old_union_sup_ratio<<" old_conf="<<old_rule_confidence<<endl;
			//cout<<"new_sup="<< new_union_sup_ratio<<" old_conf="<<new_rule_confidence<<endl;

			// ͳ��beta��gamma
			if (  (old_union_sup_ratio>=MST && old_rule_confidence>=MCT && old_rule_lift>=MLT) && (new_union_sup_ratio<MST || new_rule_confidence<MCT || new_rule_lift<MLT)  )
			{
				beta++;  cout<<"\nbeta is incereasing for lost rules!   beta="<<beta;
			}

			if (  (old_union_sup_ratio<MST || old_rule_confidence<MCT || old_rule_lift<MLT) && (new_union_sup_ratio>=MST && new_rule_confidence>=MCT && new_rule_lift>=MLT)  )
			{
				gamma++;  cout<<"\ngamma is incereasing for spurious rules!   gamma="<<gamma;
			}

	}

	beta = beta - (s_set_pair.size()-alpha); 

	
	/*----------------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
	}
	else
	{
		cout<<"\nERROR in dimension in FUNC eval()!\n";
	}

	printf("\n");
	printf("alpha=%d  beta=%d  gamma=%d  |X|= %u  [Objectives]: ", alpha, beta, gamma, actual_num_del );

	for(int i=0;i<dimension;i++)
		printf(" %.2lf  ", x->objective_value[i]);

	cout<<"\n\n";

	// ��main_trie��main_trie_copyͬ���� ʹmain_trie_copy�ָ�����ǰ��ֵ
	p_apriori->get_pointer_to_Apriori_Trie()->reset_tree_to_same(p_main_trie_copy);

}




void cal_fitness_two_rules(individual *x)
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;

	set<int> set_transID_delete;       // ����set�������ͣ� ��֤�������ı��Ψһ�����Զ�����, �������ȥ���ظ�����

	//����chromosome���� ��ѡ�е�TID����trans��
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //��Ϊtrans��set<int>���ͣ�set�����Զ��������򣩣�����Ҫ���ǲ������ظ�Ԫ��
		}
		else
		{
			printf("ERROR in function eval() : ");		exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();


	//=============================================================================================
	//�Ľ������δ����ʡ��ʱ�䣡��������ʱ1�����ң�
	//=============================================================================================

	// ����ѡ�е�trans, ����copy_support_of_items_one2 �� copy_temp_counter_array2

	vector<unsigned long> copy_support_of_items_one2(copy_support_of_items_one);
	vector< vector<unsigned long> > copy_temp_counter_array2(copy_temp_counter_array);


	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> thin_trans = thin_DB[*CI];		// ��"����" ���ݿ� ��ȡ����Ӧtransaction; 
													
		vector<int>::const_iterator it1, it2;		

		for(it1=thin_trans.begin(); it1!=thin_trans.end(); it1++)
			copy_support_of_items_one2[*it1] -= 1;

		//���¡�2-item������
		if(thin_trans.size()>=2)										//����������ǳ���Ҫ��������,basket.size()==0�� �ᵼ������ʱ�쳣���������� thin_trans.size()����>=1
		{																//���ڡ���α��롱�� ��Ϊѡ�е�trans�������ˣ��������й��� ���size>=2; ���ڲ������˵ġ�һ����롱������ѡ�е�����transΪ��trans
			for(it1=thin_trans.begin(); it1!=thin_trans.end()-1; it1++)		// thin_DB���ݿ��ÿһ���Ѿ����򣨰����±��롱��
				for(it2=it1+1; it2!=thin_trans.end(); it2++)
					copy_temp_counter_array2[*it1][*it2-*it1-1]   -= 1;
		}
	}

	//---------------------- ͳ��alpha-----------------------------------
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		item sen_condition = s_set_pair[j].first;
		item sen_consequence = s_set_pair[j].second;

		itemtype sen_inter_code_condition   = copy_new_code[ *( sen_condition.begin() )  ] - 1;// ��Ʒ���-������
		itemtype sen_inter_code_consequence = copy_new_code[ *( sen_consequence.begin()) ] - 1;

		int sen_new_condition_sup   =  copy_support_of_items_one2[ sen_inter_code_condition  ]  ;   
		int sen_new_consequence_sup =  copy_support_of_items_one2[ sen_inter_code_consequence ] ;  
		int sen_new_union_sup ;

		if( sen_inter_code_condition <  sen_inter_code_consequence   )
			sen_new_union_sup = copy_temp_counter_array2[ sen_inter_code_condition ][sen_inter_code_consequence - sen_inter_code_condition - 1] ;
		else if( sen_inter_code_condition > sen_inter_code_consequence   )
			sen_new_union_sup = copy_temp_counter_array2[ sen_inter_code_consequence][sen_inter_code_condition- sen_inter_code_consequence - 1];
		else
		{	cout<<"\n\nERROR in cal_fitness_two_rules(): (sen_inter_code_condition) and (sen_inter_code_consequence)Ӧ�ò���ȣ�\n   ";
			getchar(); exit(-1);
		}

		double sen_new_union_sup_ratio =  (double)( sen_new_union_sup)/(db_size-actual_num_del);
		double sen_new_rule_confidence =  (double)( sen_new_union_sup)/( sen_new_condition_sup );
		double sen_new_rule_lift = (double)(sen_new_union_sup)/(sen_new_condition_sup)*(double)(db_size-actual_num_del)/sen_new_consequence_sup;

		if(sen_new_union_sup_ratio >= MST && sen_new_rule_confidence >=MCT && sen_new_rule_lift >= MLT)
		{	
			alpha++;	cout<<"\nalpha is incereasing for sensitive rules not hidden!   alpha="<<alpha;
		}

	}
	/*------------���´��� �Ա��¾ɽṹ  ͳ��beta��gamma-----------------*/

	//for(unsigned n=0;n<candi_rulesets.size(); n++)
	//{
		//map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		//for(; CI!=candi_rulesets[n].end(); CI++)
		map<pair<item,item>,vector<int>>::const_iterator CI = bodon_ruleset.begin();  //��ע�⡿��bodon_rulesetʹ�õ����ڲ�����
		for(; CI != bodon_ruleset.end(); CI++)
		{
			item condition = CI->first.first;			//item����set<int>����
			item consequence = CI->first.second;

			double old_union_sup_ratio = (double)(CI->second[2])/db_size;
			double old_rule_confidence = (double)(CI->second[2])/(CI->second[0]);
			double old_rule_lift = (double)(CI->second[2])/(CI->second[0])*(double)db_size/CI->second[1];

			//itemtype inter_code_condition   = copy_new_code[ *( condition.begin() )  ] - 1;// ��Ʒ���-������
			//itemtype inter_code_consequence = copy_new_code[ *( consequence.begin()) ] - 1;
			itemtype inter_code_condition   =  *( condition.begin() ); //��ע�⡿��bodon_rulesetʹ�õ������ڲ�����
			itemtype inter_code_consequence =  *( consequence.begin()) ;
			 

			int new_condition_sup   =  copy_support_of_items_one2[ inter_code_condition  ]  ;   
			int new_consequence_sup =  copy_support_of_items_one2[ inter_code_consequence ] ;  
			int new_union_sup ;

			if( inter_code_condition <  inter_code_consequence   )
				new_union_sup = copy_temp_counter_array2[ inter_code_condition ][inter_code_consequence - inter_code_condition - 1] ;
			else if( inter_code_condition > inter_code_consequence   )
				new_union_sup = copy_temp_counter_array2[ inter_code_consequence][inter_code_condition- inter_code_consequence - 1];
			else
			{	cout<<"\n\nERROR in cal_fitness_two_rules(): (inter_code_condition) and (inter_code_consequence)Ӧ�ò���ȣ�\n   ";
				getchar(); exit(-1);
			}

			double new_union_sup_ratio =  (double)( new_union_sup)/(db_size-actual_num_del);
			double new_rule_confidence =  (double)( new_union_sup)/( new_condition_sup );
			double new_rule_lift = (double)(new_union_sup)/(new_condition_sup)*(double)(db_size-actual_num_del)/new_consequence_sup;

			/*
			if(  is_sensitive_rule(condition, consequence)   )
			{
				if(new_union_sup_ratio >= MST && new_rule_confidence >=MCT && new_rule_lift >= MLT)
				{	
					alpha++;	cout<<"\nalpha is incereasing for rules!   alpha="<<alpha;
				}
			}
			*/

			// ͳ��beta��gamma
			if (  (old_union_sup_ratio>=MST && old_rule_confidence>=MCT && old_rule_lift>=MLT) && (new_union_sup_ratio<MST || new_rule_confidence<MCT || new_rule_lift<MLT)  )
			{
				//if( ! is_sensitive_rule(condition, consequence) )
				//{
					beta++;  cout<<"\nbeta is incereasing for lost rules!   beta="<<beta;
				//}
			}

			if (  (old_union_sup_ratio<MST || old_rule_confidence<MCT || old_rule_lift<MLT) && (new_union_sup_ratio>=MST && new_rule_confidence>=MCT && new_rule_lift>=MLT)  )
			{
				gamma++;  cout<<"\ngamma is incereasing for spurious rules!   gamma="<<gamma;
			}

		}//CI
	//}//n

	beta = beta -  (s_set_pair.size()-alpha); 


	//x->fitness = (double)(alpha*50+beta*25+gamma*25)/100.0;			//��Ϊ�ö�Ŀ�����

	/*---------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
	}
	else
	{
		cout<<"\nERROR in dimension in FUNC eval()!\n";
	}

	printf("\n");
	printf("alpha=%d  beta=%d  gamma=%d  |X|= %u  [Objectives]: ", alpha, beta, gamma, actual_num_del );

	for(int i=0;i<dimension;i++)
		printf(" %.2lf  ", x->objective_value[i]);

	cout<<"\n\n";

	return;
}




void cal_fitness_two(individual *x)
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;


	set<int> set_transID_delete;       // ����set�������ͣ� ��֤�������ı��Ψһ�����Զ�����, �������ȥ���ظ�����

	//vector<vector<int>>   affect;

	//����chromosome���� ��ѡ�е�TID����trans��
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //��Ϊtrans��set<int>���ͣ�set�����Զ��������򣩣�����Ҫ���ǲ������ظ�Ԫ��
		}
		else
		{
			printf("ERROR in function eval() : ");
			exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();

	//=============================================================================================
	//�Ľ������δ����ʡ��ʱ�䣡��������ʱ1�����ң�
	//=============================================================================================

	// ���������ݽṹ ���� ��ԭ���Ľṹ �Ƚϣ������õ�alpha,beta,gamma

	vector<unsigned long> copy_support_of_items_one2(copy_support_of_items_one);
	vector< vector<unsigned long> > copy_temp_counter_array2(copy_temp_counter_array);


	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> thin_trans = thin_DB[*CI];		//��"����" ���ݿ� ��ȡ����Ӧtransaction
		vector<int>::const_iterator it1, it2;

		for(it1=thin_trans.begin(); it1!=thin_trans.end(); it1++)
			copy_support_of_items_one2[*it1] -= 1;

		//���¡�2-item������
		if(thin_trans.size()>=2)		//���������ǳ���Ҫ������һ������ҵ���bug������,basket.size()==0�� �ᵼ������ʱ�쳣���������� thin_trans.size()����>=1
		{								//���ڡ���α��롱�� ��Ϊѡ�е�trans�������ˣ��������й��� ���size>=2; ���ڲ������˵ġ�һ����롱������ѡ�е�����transΪ��trans
			for(it1=thin_trans.begin(); it1!=thin_trans.end()-1; it1++)
				for(it2=it1+1; it2!=thin_trans.end(); it2++)
					copy_temp_counter_array2[*it1][*it2-*it1-1]   -= 1;
		}
	}

	/*------------���´��� �Ա� һά �� ��ά �¾ɽṹ-----------------*/

	// "1-item"��"2-item"��alpha����
	for(vector<item>::const_iterator sen_CI = s_set.begin(); sen_CI!=s_set.end();sen_CI++)  //���ÿ�� ����itemset�Ƿ����ص���<sup_l_actual�� 
	{
		item sen_item = *sen_CI;
		if(sen_item.size()== (unsigned int)1)   // "1-item"������
		{
			int internal_code =  copy_new_code[*(sen_item.begin())]-1;		//����Ʒ���롱 ת���� ���ڲ��롱
			if( copy_support_of_items_one2[internal_code] >= (unsigned)sup_l_actual  )
				alpha++;
		}
		
		if(sen_item.size()== (unsigned int)2)   // "2-item"������
		{
			unsigned int code1, code2;
			code1 = copy_new_code[*(sen_item.begin())]-1;					//����Ʒ���롱 ת���� ���ڲ��롱
			code2 = copy_new_code[*(sen_item.end())]-1;

			if(copy_temp_counter_array2[code1][code2-code1-1] >= (unsigned)sup_l_actual  )
			{
				alpha++;
				cout<<"\n Alpha is increasing!   Alpha = "<<alpha;
			}
		}
	}

	// "1-item"��beta��gamma����
	for(unsigned i=0; i<copy_support_of_items_one.size(); i++)
	{
		if( copy_support_of_items_one[i]< (unsigned)sup_u  &&  copy_support_of_items_one2[i]>= (unsigned)sup_l_actual )   // ����Ƶ����--����Ƶ����
		{
			gamma++;  cout<<"\n Gamma is increasing(1-item)!   Gamma = "<<gamma;
		}

		if( copy_support_of_items_one[i] >= (unsigned)sup_u &&  copy_support_of_items_one2[i]< (unsigned)sup_l_actual )  // ��Ƶ����--����Ƶ��(������)
		{
			//�ȼ���ǲ���������
			int flag=0;
			item temp_item;
			temp_item.insert( (int)(copy_new_code_inverse[i]) );

			for(vector<item>::const_iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)				
			{		
				if ( temp_item == *s_setCI) {cout<<"\n ��1-item������"; flag =1; break;}
			}

			if(0 ==flag )  // "������"�� ������
			{
				beta++;  cout<<"\n Beta is increasing(1-item)! Beta = "<<beta;
			}	
		}
	}

	// "2-item"��beta��gamma����
	for(unsigned i=0; i<copy_temp_counter_array.size(); i++)
		for(unsigned j=0; j<copy_temp_counter_array[i].size(); j++)   //copy_temp_counter_array[i].size() == copy_temp_counter_array.size()+1 - (i+1)
		{
			if( copy_temp_counter_array[i][j]<sup_u && copy_temp_counter_array2[i][j]>=sup_l_actual )   // ����Ƶ����--����Ƶ����
			{
				gamma++;  cout<<"\n Gamma is increasing (2-item)!   Gamma = "<<gamma;
			}
			if( copy_temp_counter_array[i][j]>=sup_u && copy_temp_counter_array2[i][j]<sup_l_actual )	// ��Ƶ����--����Ƶ�� (������)
			{
				//�ȼ���ǲ���������	
				item temp_item;
				temp_item.insert( (int)(copy_new_code_inverse[i]) );		// ���ڲ��롱-->����Ʒ��š�
				temp_item.insert( (int)(copy_new_code_inverse[j+i+1]) );    // ע�⣺����j�±꣬ ��Ҫ������ġ��洢�±ꡱת���ɡ��߼��±ꡱ��
																				// �Σ�Apriori_Trie.cpp�ļ���find_candidate_two()���������ע��
				int flag=0;
				for(vector<item>::const_iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)				
				{		
					if ( temp_item == *s_setCI) {cout<<"\n ��2-item������"; flag =1; break;}
				}
				if(0==flag)  // "������"�� ������
				{
					beta++; cout<<"\n Beta is increasing(2-item)! Beta = "<<beta;
				}
			}
		}



	// ע�⣺ ��ȡsensitive item�Ĺ����Ѿ��Ƶ��� state0()�����У�����


	//x->fitness = (double)(alpha*50+beta*25+gamma*25)/100.0;			//��Ϊ�ö�Ŀ�����
	/*---------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
	}
	else
	{
		cout<<"\nERROR in dimension in FUNC eval()!\n";
	}

	printf("\n");
	printf("alpha=%d  beta=%d  gamma=%d  |X|= %u  [Objectives]: ", alpha, beta, gamma, actual_num_del );

	for(int i=0;i<dimension;i++)
		printf(" %.2lf  ", x->objective_value[i]);


	cout<<"\n\n";

	return;
}



void cal_fitness(individual *x)
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;


	set<int> trans;    // ����set�������ͣ� ��֤�������ı��Ψһ�����Զ�����
	DB affect;

	//��ҪӋ��fitness���е�TID����trans��
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			trans.insert(x->bit_string[i]);  //��Ϊtrans��set<int>���ͣ�set�����Զ��������򣩣�����Ҫ���ǲ������ظ�Ԫ��
		}
		else
		{
			printf("ERROR in function eval() : ");
			exit(-1);
		}
	}

	//�����x��Transaction��������item����affect
	/*=========================================
	[3](1,5,7)                [0](1,5,7)
	[5](3,5)        -->       [1](3,5)
	[6](4,6,8)                [2](4,6,8)
	===========================================*/

	for (set<int>::const_iterator transCI=trans.begin(); transCI!=trans.end(); transCI++)
	{
		//affect.push_back(thin_DB[*transCI]);
		affect.push_back(database[*transCI]);
	}

	//*****************************************
	int count = 0;
	collectset affectlarge = large_sets;          // ����large_sets��	large_sets�Ȱ���real_large,Ҳ����pre_large


	/*=========================================
	[0]( [0]( (234),345)),  [1]((245),567),...)
	[1]( [0]( (234,245),212) ),  ...)
	[2]( [0]( (234,245,345), 230 ), ...    )
	===========================================*/
	//collectset affectreallarge = reallarge_sets;
	//collectset affectprelarge = prelarge_sets;

	// affectlarge����������large item sets,����1-item,2-item,3-item,...
	// affectlarge��ÿһ�д洢��n-item��Ƶ����(large item sets)
	// ������Ĳ�ѭ���߼��������ģ�
	//		���affectlarge��ÿһ�е�ÿһ��Ƶ�����Ӧ��������for��, 
	//			ɨ��affect���ݿ�(���ݻ���ѡ��transaction���ɵ�΢����database)��ÿһ��(��Ӧ�������������for)
	//				�����Ƶ����(goalitem)��affect���ݿ��ÿ�н��бȶ�(��Ӧwhileѭ��)�� 
	//				���ĳ��transaction��ĳ�У�������item set,��affectlarge[n][goalitem]--;
	//				�����׼��ɾ����transaction���������ԣ��������������1�����Կ�������������Ƶ���ģ�����large item sets��

	// Ҳ���ǣ�����ÿһ��large item,ɨ��׼��ɾ����transaction���ɵ����ݿ��ÿһ�У����бȶԣ����������Ƶ����͸������ļ���ֵ��

	//====================================================================================================
	//������δ���ǳ���ʱ���൱������ɨ��500��(length=500) �����ݿ⣬�ٽ���ģʽƥ�䣡��������ʱ12�����ң�
	//====================================================================================================
	//ȡ��affect���ݿ���ÿһ��Transaction����ÿһ��Ƶ����ƥ��

	for (unsigned int n=0; n < affectlarge.size(); n++)
	{
		//ȡ��large_sets��ÿһ��item	
		for (itemsetCI tempCI=affectlarge[n].begin(); tempCI!=affectlarge[n].end(); tempCI++)
		{
			item goalitem = tempCI->first;
			//affect�����Ҫ�u��TID������items
			for (DBCI db_CI=affect.begin(); db_CI!=affect.end(); db_CI++)
			{
				transaction goaltrans = *db_CI;				//affect��transaction;  transaction�����Զ����򣬵�database��ÿ��������
				itemCI itemCI = goalitem.begin();			//large_sets��item ; itemΪset<int>, �����Ԫ���Զ���������
				transactionCI transCI = goaltrans.begin();	
				//һ��itemһ��item�Ȍ�
				while (itemCI!=goalitem.end() && transCI!=goaltrans.end())  //  �����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���
																			//  ��set����set<int>���Զ�����Ĺ���,vector<int>û�У�����initial()�����˹�����
				{
					if (*transCI == *itemCI)
					{
						count++;
						itemCI++;
						transCI++;
					}
					else
					{
						transCI++;
					}
				}
				if (count == n+1)
				{
					affectlarge[n][goalitem]--;					//ɾ����length��transactions,ֻ��ʹ����Ƶ������Ƶ�ν���,����������
					//affectreallarge[n][goalitem]--;
					//affectprelarge[n][goalitem]--;
					//cout<< "UPDATE affectedlarge" << endl;
				}
				count = 0;
			}
		}
	}


	// ע�⣺ ��ȡsensitive item�Ĺ����Ѿ��Ƶ��� state0()�����У�����


	// ��Ը��º��large item��ÿһ�е�ÿһ��Ƶ����(���º�ļ����Դ���sup_l)������Ӧ�������forѭ����
	// ��������������������s_set�У�˵��ͨ��ɾ����Щ����ָ����transaction����δ�����������������������������Ƶ���ģ�����sup_l��
	// �� alpha++;

	// ��ע�⡿�� large_sets������reallarge_sets��prelarge_sets


	for (unsigned int n=0; n<affectlarge.size(); n++)											//����n+1Ƶ���
	{
		for (itemsetCI tempCI=affectlarge[n].begin(); tempCI!=affectlarge[n].end(); tempCI++)	//����n+1Ƶ�����ÿһ��Ƶ����
		{
				//�����������
				int sen_flag=0;
				for(vector<set<int>>::iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)				//�������м���(vector)�е�ÿһ��
				{																				//ÿһ����һ��item(��set<int>����)
					//cout<<"\n�ߵ���������if�ж�������\n";									
					if (tempCI->first==*s_setCI)
					{
						cout<<"\nAlpha��Alpha�� "<<" sup_l_actual ="<< sup_l_actual <<"  (tempCI->second)="<<(tempCI->second)<<endl ;

						if((tempCI->second) >= (db_size-length)*MST  )											//����Ҫ�޸ġ����˴�����sup_l��(db_size-length)*MST����
						{
							alpha++;    // alpha��ʾδ���ص��������������
										// �����ɾ����ָ��Transaction�� ��Ӧ��large item��(���º��)��Ȼ����
							sen_flag=1;
						
							cout << "ALPHA is increasing:  "<<"alpha = " << alpha << "   curGen:  "<< gen << endl;;
						}
					}
				}
				if(1== sen_flag)   //���affectlarge���е�n�еĵ�ǰ���ǡ�����������������жϣ�������һ��ѭ��
					continue;

			
				// large_sets��ɾ������ָ����transaction֮ǰ��Ƶ����� affectlarge��ɾ������֮���Ƶ�����

				for (itemsetCI tempCII=large_sets[n].begin(); tempCII!=large_sets[n].end(); tempCII++)
				{   // tempCIIΪԭ����Ƶ�����large_sets��ָ�룻 tempCI Ϊ����ɾ����transactions ����ֵ���Ƶ������ָ��

					if ((tempCI->first==tempCII->first) && ((tempCII->second)>=sup_u)   && (  (tempCI->second)< sup_l_actual  )  )
					{		
						// ĳ��ԭ����Ƶ��(real_large)�ģ� ɾ��ָ��transactions֮�� ��ɲ�Ƶ���ģ�ע�⣺Ӧ���ǲ���������������������
						// missing large items
					
						//�˴������Ƿ�Ϊ��������ж�
						//-------------------------���´�����Cheng Peng���-----------------------------------

						//cout<<"\n ĳ��ԭ����Ƶ���ģ� ɾ��ָ��transactions֮�� ��ɲ�Ƶ���ģ�\n";
						//getchar();

						int flag=0;
						for(vector<set<int>>::iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)
						{
							if (tempCI->first==*s_setCI)  //˵�������Ƶ����Ϊ��Ƶ��������������
							{
								flag=1;break;
							}
						}

						if( 0==flag ) //flag==0˵���������������ɾ��ָ��transaction֮����Ƶ����ɲ�Ƶ��
						{
							beta++; 
							cout<<"Beta! Beta! "<<" sup_l_actual ="<< sup_l_actual <<" (tempCI->second)="<<(tempCI->second)<< "  (tempCII->second)= "<<(tempCII->second)<< endl ;
							cout<<"\nBETA is increasing:  " << " beta == "<< beta << "   curGen:  "<< gen << "\n";
						}

						//---------------------------------------------------------------------------------------	  
					}

					if(   (tempCI->first==tempCII->first) && ((tempCII->second)<sup_u) && ((tempCI->second)>=sup_l_actual  )  )
					{				
								//ʵ���ϼ�ʹpre_large �������item ��ɾ��trans�� ��Ȼ>=sup_l
								// tempCI������Ƶ�������У� ���Ѿ�˵����Ƶ�����ˣ� ��˲�����(tempCII->second)<sup_u)
								// ������Ҫ���¿���! �ѵ����������д��� ע�⣺ԭ����large_sets����real_large��pre_large.
								// ��ˣ�ԭ��С��sup_uͬʱ����sup_l��Ƶ���������ɾ��length��transactions֮���Դ���sup_l
								// ��ˣ�ɾ���ı���length����̫�࣬�����Ȼ����Ƶ��λ��[sup_l,sup_u]�����Ƶ����̫�࣬ gammaֵ����

								
								// ĳ��ԭ���ǲ�Ƶ����(������real_large,����pre_large)�� ɾ��ָ��transactions֮�� ���Ƶ����
								// ���������Ҳ����������� ��ǰ���continue����ʾ��
								// ������������ɾ��transaction������Ƶ���ģ�(tempCI->second)>=sup_l)��,
								// ��alpha++������������������������ж�
								// �����ִ�е���� ˵����ǰ�� (*tempCI) �ǡ������еġ�
					
								// generate new frequent items which doesn't exist before; 
								// But why the condition is ((tempCII->second)<sup_u)  && ((tempCI->second)>=sup_l)
						gamma++;

						
						cout<<"Gamma! Gamma! "<<" sup_l_actual ="<< sup_l_actual <<" (tempCI->second)="<<(tempCI->second)<< "  (tempCII->second)= "<<(tempCII->second)<< endl ;
						cout << "\n\n GAMMA is increasing:  gamma = " << gamma << "   curGen:  "<< gen << endl;
					}
				}
		}
	}


	//x->fitness = (double)(alpha*50+beta*25+gamma*25)/100.0;			//��Ϊ�ö�Ŀ�����

	x->objective_value[0] = (double)alpha;   //ע�����ǵ�Ŀ��������double�ͣ� ��alpha, beta, gamma��int��
	x->objective_value[1] = (double)beta;
	x->objective_value[2] = (double)gamma;


	printf("\n");
	printf("alpha=%d  beta=%d  gamma=%d     [Objectives]: ", alpha, beta, gamma);

	for(int i=0;i<3;i++)
		printf(" %.2lf  ", x->objective_value[i]);

	cout<<"\n\n";

	return;
}



void init_random_shuffle()
{
	cout<<"\n>>>����my_random_shuffle()\n";

	//for (int i = 0; i < prelarge_trans_array_size; i++) 
	//	v1.push_back(i); 

	for(unsigned int i=0;i<s_set.size(); i++)
	{
		 vector<int> v2_row;
		 for( int j=0;j<sen_frequency_array[i]; j++)
			 v2_row.push_back(j);

		 v2.push_back(v2_row);
		 v2_row.clear();
	}
	cout<<"\n>>>�˳�my_random_shuffle()\n";
} 


/* create a random new individual and allocate memory for it,
   returns a pointer to the new individual */
individual *new_individual()
{
     individual *return_ind;
     int result;

     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->bit_string = (int *) malloc(sizeof(int) * length);  // ע��lengthΪ�����α��롱�ܳ���

	 // ��ȻΪĿ����������double�͵Ŀռ�
	 return_ind->objective_value = (double *) malloc(sizeof(double) * dimension);


	 /*
	 //ǰ��α���
	random_shuffle(v1.begin(), v1.end());

	cout<<"\nǰ��α���:\n";
	 for(int i=0;i<length1;i++)
	 {
		 return_ind->bit_string[i] = prelarge_trans_array[ v1[i] ];
		 //cout<<v1[i]<<" ";
	 }
	 */

	 // ���α���
	for(unsigned int i=0;i<s_set.size(); i++)
	{
		random_shuffle( v2[i].begin(),v2[i].end() );
	}

	cout<<"\n���α��룺\n";
	 int offset=0;
	 for ( unsigned int i=0; i < s_set_pair.size(); i++)
	 {
		 cout<<"\n��"<<i<<"�������\n";
		 for( int j=0; j < encoding_len[i]; j++ )
		 {
				
				//int pos;
				
				//pos = irand( sen_frequency_array[i]  );  // ����û��Ҫ����������ظ��������ݿ��¼�� (ͨ��<set>���Ϳ���ʵ�����ظ�)
				//pos = irand( sen_trans_vec[i].size() );

				return_ind->bit_string[ length1 + offset + j] = sen_trans_array[i][ v2[i][j] ] ;   
				
				//cout<<v2[i][j]<<" ";
												
		 }
		 offset = offset + encoding_len[i];	
	 }
	 cout<<endl;

	 if(offset!=length2){cout<<"\n ERROR! in *new_individual():  offset!=length2 "; getchar(); exit(-1);}

     //return_ind->length = length1+length2;
	 return_ind->length = length;

     /* evaluating the objective functions */

     result = eval(return_ind);

	 cout<<"\n ������һ����һ���¸���";
     return (return_ind);
}


/* copy an individual and return the pointer to it */
individual *copy_individual(individual *ind)
{
     individual *return_ind;
     int i;

     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->bit_string = (int *) malloc(sizeof(int) * length);
     return_ind->objective_value = (double *) malloc(sizeof(double) * dimension);
     
     for (i = 0; i < length; i++)
		return_ind->bit_string[i] = ind->bit_string[i];

     for (i = 0; i < dimension; i++)
		return_ind->objective_value[i] = ind->objective_value[i];

     return_ind->length = ind->length;

     return(return_ind);
}


/* ע�⣺����������������ʵ�������

		1) һ����write_output_file(), д�����š�Ŀ�������������������ļ�

		2����һ����write_output_obj()��д�����š�Ŀ���������ļ�������������������
		
*/
/* Writes the index, objective values and bit string of
   all individuals in global_population to 'out_filename'. */
void write_output_file()
{
     int j, current_id;
     FILE *fp_out;
     individual *temp;

	 char outfile_name[128];     //ע�Ȿ�ز����ļ��е� outputfile �����Ѿ�������
	 	 
	 sprintf(outfile_name,"%s_MST_%.3f_MCT_%.2f_MLT_%.2f_output.dat",datafrom, MST,MCT,MLT);
     
     fp_out = fopen(outfile_name, "w");  
     assert(fp_out != NULL);

	 /*******************************************************************/
	 //д����ز���������ļ�



	 fprintf(fp_out, "Length	    = %d \n", length); 
	 fprintf(fp_out, "MaxGen	    = %d \n", maxgen); 
	 fprintf(fp_out, "datafrom	    = %s \n", datafrom); 
	 fprintf(fp_out, "sensfile	    = %s \n", sensfile); 

	 fprintf(fp_out, "FIM(1)AR(2)   = %.3f \n",FIM_or_AR );  
	 fprintf(fp_out, "MST	        = %.3f \n",MST );   // �����������
	 fprintf(fp_out, "MCT	        = %.3f \n",MCT );	// �����������
	 fprintf(fp_out, "MLT	        = %.3f \n",MLT );	// �����������

	 fprintf(fp_out, "sup_u		    = %d \n", sup_u); 
	 fprintf(fp_out, "sup_l		    = %d \n", sup_l); 
	 fprintf(fp_out, "db_size	    = %d \n", db_size); 

	 char tempStr[128];
	 if(0==mutation_type) strcpy(tempStr, "no mutation");
	 else if(1== mutation_type) strcpy(tempStr,"one bit mutation");
	 else if(2== mutation_type) strcpy(tempStr,"independent bit mutation");
	 else if(3== mutation_type) strcpy(tempStr,"shuffle_random mutation");
	 else strcpy(tempStr,"INVALID mutation type");
	 	 
	 fprintf(fp_out, "mutation_type	= %s \n", tempStr); 

	 if(0==recombination_type) strcpy(tempStr, "no recombination");
	 else if(1== recombination_type) strcpy(tempStr,"one point crossover");
	 else if(2== recombination_type) strcpy(tempStr,"uniform crossover");
	 else if(3== recombination_type) strcpy(tempStr,"shuffle_random crossover");
	 else strcpy(tempStr,"INVALID recombination type");
	 	 
	 fprintf(fp_out, "recombination_type	= %s \n", tempStr); 

	 fprintf(fp_out, "mutat_prob	= %.3f \n",mutat_prob ); 
	 fprintf(fp_out, "recom_prob	= %.3f \n",recom_prob ); 
	 fprintf(fp_out, "bit_turn_prob	= %.3f \n\n",bit_turn_prob ); 

	 /*******************************************************************/


     current_id = get_first();

     while (current_id != -1)
     {       
			temp = get_individual(current_id);
			fprintf(fp_out, "%d \n", current_id); /* write index */

			cout<<endl;
			for (j = 0; j < dimension; j++)
			{
				printf( "%.2lf ",   get_objective_value(current_id, j)    ); //��Ϊ��Ӧ������ʱ�����˶�Ŀ����С���������ڻ�ԭ
				fprintf(fp_out, "%.2lf ",  get_objective_value(current_id, j)  );
			}
            fprintf(fp_out, "\n");

			for (j = 0; j < temp->length; j++)
			{
				fprintf(fp_out, "%d  ", temp->bit_string[j]);
			}
			fprintf(fp_out, "\n\n\n");

			current_id = get_next(current_id);
     }

     fclose(fp_out);

	 printf("\nWriting outcome to the file \"%s\" is complete!",outfile);

	 getchar();
}


/**********| addition for PPDM  |*******/

void write_output_obj()
{
     int j, current_id;
     FILE *fp_out;
     individual *temp;

	 char output_obj_file[128];

	 //sprintf(output_obj_file, "%s_obj",outfile);
	 //sprintf(output_obj_file, "%s_obj.dat","PPDM_dt_output");
	 	 
	 sprintf(output_obj_file,"%s_MST_%.3f_MCT_%.2f_MLT_%.2f_obj.dat",datafrom, MST,MCT,MLT);
     

     fp_out = fopen(output_obj_file, "w");
     assert(fp_out != NULL);


     current_id = get_first();

     while (current_id != -1)
     {       
			temp = get_individual(current_id);
			fprintf(fp_out, "%d ", current_id); /* write index */

			cout<<endl;
			for (j = 0; j < dimension; j++)
			{
				printf( "%.2lf ",   get_objective_value(current_id, j)    ); //��Ϊ��Ӧ������ʱ�����˶�Ŀ����С���������ڻ�ԭ
				fprintf(fp_out, "%.2lf ",  get_objective_value(current_id, j)  );
			}
            fprintf(fp_out, "\n");

			current_id = get_next(current_id);
     }



	 /*******************************************************************/
	 //д����ز���������ļ�

	 fprintf(fp_out, "\n\n\n");

	 fprintf(fp_out, "Length	    = %d \n", length); 
	 fprintf(fp_out, "MaxGen	    = %d \n", maxgen); 
	 fprintf(fp_out, "datafrom	    = %s \n", datafrom); 
	 fprintf(fp_out, "sensfile	    = %s \n", sensfile); 

	 fprintf(fp_out, "FIM(1)AR(2)   = %.3f \n",FIM_or_AR );  
	 fprintf(fp_out, "MST	        = %.3f \n",MST );   // �����������
	 fprintf(fp_out, "MCT	        = %.3f \n",MCT );	// �����������
	 fprintf(fp_out, "MLT	        = %.3f \n",MLT );	// �����������

	 fprintf(fp_out, "sup_u		    = %d \n", sup_u); 
	 fprintf(fp_out, "sup_l		    = %d \n", sup_l); 
	 fprintf(fp_out, "db_size	    = %d \n", db_size); 

	 char tempStr[128];
	 if(0==mutation_type) strcpy(tempStr, "no mutation");
	 else if(1== mutation_type) strcpy(tempStr,"one bit mutation");
	 else if(2== mutation_type) strcpy(tempStr,"independent bit mutation");
	 else if(3== mutation_type) strcpy(tempStr,"shuffle_random mutation");
	 else strcpy(tempStr,"INVALID mutation type");
	 	 
	 fprintf(fp_out, "mutation_type	= %s \n", tempStr); 

	 if(0==recombination_type) strcpy(tempStr, "no recombination");
	 else if(1== recombination_type) strcpy(tempStr,"one point crossover");
	 else if(2== recombination_type) strcpy(tempStr,"uniform crossover");
	 else if(3== recombination_type) strcpy(tempStr,"shuffle_random crossover");
	 else strcpy(tempStr,"INVALID recombination type");
	 	 
	 fprintf(fp_out, "recombination_type	= %s \n", tempStr); 

	 fprintf(fp_out, "mutat_prob	= %.3f \n",mutat_prob ); 
	 fprintf(fp_out, "recom_prob	= %.3f \n",recom_prob ); 
	 fprintf(fp_out, "bit_turn_prob	= %.3f \n\n",bit_turn_prob ); 

	 /*******************************************************************/

     fclose(fp_out);

	 printf("\nWriting objective vector to the file \"%s\" is complete!",output_obj_file);

	 getchar();
}



//�� ��һ����Ⱥ��Ŀ��ռ�����д���ļ�
void write_ini_obj()
{
	   int j, current_id;
     FILE *fp_out;
     individual *temp;

	 char output_obj_file[128];

	 sprintf(output_obj_file,"%s_MST_%.3f_MCT_%.2f_MLT_%.2f_ini.dat",datafrom, MST,MCT,MLT);
     
     fp_out = fopen(output_obj_file, "w");
     assert(fp_out != NULL);


     current_id = get_first();

     while (current_id != -1)
     {       
			temp = get_individual(current_id);
			fprintf(fp_out, "%d ", current_id); /* write index */

			cout<<endl;
			for (j = 0; j < dimension; j++)
			{
				printf( "%.2lf ",   get_objective_value(current_id, j)    ); //��Ϊ��Ӧ������ʱ�����˶�Ŀ����С���������ڻ�ԭ
				fprintf(fp_out, "%.2lf ",  get_objective_value(current_id, j)  );
			}
            fprintf(fp_out, "\n");

			current_id = get_next(current_id);
     }


	 fclose(fp_out);

	 printf("\nWriting objective vectors in the first generation to the file \"%s\" is complete!",output_obj_file);


}
