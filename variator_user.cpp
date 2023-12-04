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

//int db_size; /* 改在 apriori.cpp文件中定义*/


/*=================================Add for PPDM ===================================*/

	int length; /* length of the binary string */  // MAX_del_trans是底线， 决定了sup_l,但真正的删除笔数是length(即编码长度)
	int maxgen; /* maximum number of generations (stop criterion) */


	//------------------------------------------------------------------------------[void cal_length()函数]
	int sum1_prelarge =0,  sum2_sensitive = 0;
	int length1;	// length1的值代表了“二段编码”的分段位置； 前半段长度为length1，前半段在prelarge_trans_array[]中搜索
	int length2;	// 后半段长度为length2, 后半段编码在 sen_trans_array[] 中搜索
	//------------------------------------------------------------------------------


	//----------------------------------------------
	//			保存涉及sensitive item sets的trans
	//----------------------------------------------

		
	int sen_trans_array_size=0;				// sen_trans_array[]数组长度
 											// 从sen_trans_set转存过来，保存了包含敏感项的数据库记录在 数据库中的序列号，序列号范围[0,database.size()-1];
	int ** sen_trans_array;					// 在filter_sen_trans()中申请内存空间,保存每个sen item涉及的trans编号


	int * sen_frequency_array=NULL;			// 保存了每个敏感项的出现频次；也是在filter_sen_trans()中申请内存空间

	//----------------------------------------------
	//			保存涉及prelarge item sets的trans
	//----------------------------------------------
	set<int> prelarge_trans_set;				// pre_large项涉及的trans集合， 使用set<int>类型可以防止 插入重复的元素

	int prelarge_trans_array_size=0;			// prelarge_trans_array[]数组长度
	int * prelarge_trans_array=NULL;			// 从prelarge_trans_set转存过来，保存了包含敏感项的数据库记录在 数据库中的序列号，序列号范围[0,database.size()-1];
												// 在filter_prelarge_trans()中申请内存空间

	int * prelarge_frequency_array=NULL;		// 保存了每个prelarge项的出现频次；也是在filter_prelarge_trans()中申请内存空间

	
	//----------------------------------------------
	//			计算overlap for each large item set
	//----------------------------------------------

	double sel_overlap_percent=0.1;		// 确定在高出sup_l的百分之多少比例的trans用来计算overlap,因为高频频繁集对低频项overlap计算影响较大，
										// 容易模糊低频频繁集之间的overlap差别。且超过sup_l过多的高频频繁集不适合做敏感项，因为采用它们删除的trans太多

	int *encoding_len;					// 动态申请数组， 用于创建每个sensitive item对应的染色体中的编码长度。请参考PPT中“多段编码”部分。


	vector<int> v1;					// 用于存放打乱随机排列的0~prelarge_trans_array_size-1 序列， 请参照void init_random_shuffle()

	vector<vector<int>> v2;			// 每一行用于存放打乱随机排列的0~sen_trans_vec[i].size()-1序列， 请参照void init_random_shuffle()

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

	 gen = 1;			// gen为全局变量， state0()为第一代

     initial_population = (int *) malloc(alpha * sizeof(int)); 
     if (initial_population == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }

     /*=============| added for PPDM |======================*/

     result = read_local_parameters();

	 apriori_bodon();			// 先由apriori_bodon()计算得到所有“关联规则”写入文件(并保存于bodon_ruleset)， 再从中选择“敏感规则”
								// 获取“1-item”和“2-item(矩阵)”频繁项数据结构；获取 复制树main_trie_copy。

	 cout<<"\nPlease define the sensitive rules in the file.  When finished press [Enter] to continue! \n ";
	 getchar();

	 read_sensitive_item();		// 打开敏感项文件，读取到vector<item>类型的 s_set_pair中

	 apriori_large_sets();		// 调用旧apriori算法 [如果已调了apriori_bodon,旧apriori算法只计算sup_u和sup_l]
								// 进一步调cal_support()计算sup_u和sup_l，计算sup_l-->需知道MAX_del_trans-->需知道"每个敏感项的support"-->调用filter_sen_trans();
								// 注意它间接 调用了filter_sen_trans()

	 //if(AR == FIM_or_AR && 2 == SIZE_threshold)
	 //	apriori_gen_rules();	// 自己写的产生关联规则!!!（只能产生2-item规则） apriori_bodon()已经产生了关联规则并保存

								
	 //filter_sen_trans();		// 找出包含敏感项的数据库记录编号；并得到各敏感项的频次
								// 改为在old_apriori.cpp文件 cal_support()-->cal_MAX_delTrans()函数中调用

	 //filter_prelarge_trans();	// 找出包含prelarge项的数据库记录【涉及prelarge】
								// 删除item的策略 就不需要 此函数。

	 //cal_length();				// 计算length(染色体编码长度),根据敏感项频次和sup_l 计算；read_local_parameters()读来的length值被替换 
								// 【涉及prelarge】

	 write_stat();				//  向文件写入统计信息【涉及prelarge】

	 //cal_overlap();			// 计算每个频繁项的overlap degree，写入文件
								// (对于一个数据集 计算一次就够了！)

	 //init_random_shuffle();		// 借用STL中random_shuffle生成不重复的随机数序列，初始化。在函数new_individual(),crossover(), mutation()中用到不重复的随机序列
								// 【涉及prelarge】

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
          initial_population[i] = add_individual(new_individual());  // new_individual()为新个体申请内存，随机初始化
												//【涉及prelarge】	 // 将新个体地址 插入种群 global_population.individual_array[]
																	 // add_individual（）返回值为 个体在种群中的位置（即身份证号）
          if(initial_population[i] == -1)
               return(1);
     } 

     result = write_ini(initial_population);	//将第一代群体的身份证号（即个体在种群中的插入位置）及 各维目标取值 写入ini_file文件

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
		 	write_ini_obj();   // add for PPDM， 只写入第一代群体个体编号、目标向量

     /*=================| added for PPDM |=========================*/

     gen++;		// 进化代数加1
	 cout<<"\n\n"<<"gen = "<<gen<<endl;


     result = variate(parent_identities, offspring_identities);  //【涉及prelarge】
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
        write_output_file();  // 写入个体编号、目标向量、决策向量

		write_output_obj();   // add for PPDM， 只写入个体编号、目标向量
     }


	 /*-----------释放“多段编码”涉及的动态申请的额外内存--------------*/
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
     
   result = read_arc();   //删除掉了种群中的非存档个体

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
// 根据profitWeightRatios数组内容对下标进行排序

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



//删掉很多knapsack相关代码
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
     fscanf(fp, "%lf", &SM);				//读入safety margin
	 cout<<str<<" == "<<SM<<endl;

	 fscanf(fp, "%s", str);
     assert(strcmp(str, "A01") == 0);
     fscanf(fp, "%lf", &A01);				//读入The proportion of 0's and 1's which are blocked
	 cout<<str<<" == "<<A01<<endl;

	 /************* add for PPDM ***************************/

	 // 从param_file中读入数据源，敏感项文件，和suppor阈值

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
//       功能： 计算每个频繁项的overlap degree， 写入文件 
//		 先计算出一个频繁项涉及到的所有transactions，也包含其他频繁项的总频次，
//		 coverage(large_item_i) = 总频次/large_item_i涉及到的transaction总数
//****************************************************************************

void cal_overlap()
{
	//int count = 0;
	//int sum_item = 0;

	////double sel_percent=0.1;	//【改为全局变量】确定在高出sup_l的百分之多少比例的trans用来计算overlap,因为高频频繁集对overlap影响大，
	//							//容易模糊低频频繁集之间的overlap差别。不必选，且超过sup_l过多的高频频繁集不适合做敏感项，否则删除的trans太多

	////map<item,double> item_overlap;
	////typedef map<item, double>::const_iterator itemsetCI_double;
	////typedef map<item, double>::value_type newdata_double;

	//cout<<"\n>>>Begin to calculate overlap degree for each large item set\n\n";

	//ofstream outfile;
	//char filename[80];
	//sprintf(filename, "%s_%.3f_overlap.txt",datafrom, MST);

	//outfile.open(filename, ios::out);



	//for (unsigned int i=0; i<large_sets.size(); i++)								// 统计频繁项数量，包括reallarge和prelarge
	//{
	//	for (itemsetCI tempCI=large_sets[i].begin(); tempCI!=large_sets[i].end(); tempCI++)	 
	//	{
	//			if(tempCI->second< sup_u + (int)(db_size* sel_overlap_percent ) )
	//			{				  //原来是sup_l
	//				item goalitem = tempCI->first;
	//				sum_item = 0;     // 用来统计 overlap degree;

	//				vector<int> Titem;		//用来记录 支持 当前频繁项 的数据库记录 序号

	//				for (unsigned int j=0; j<database.size();j++)					// 找出支持当前频繁项的trans,将其序号存入Titem
	//				{
	//					transaction goaltrans = database[j];				
	//					itemCI itemCI = goalitem.begin();							//  item为set<int>, 里面的元素自动升序排列
	//					transactionCI transCI = goaltrans.begin();	
	//				
	//					while (itemCI!=goalitem.end() && transCI!=goaltrans.end())  //  这样比对的前提是：goalitem 和 goaltrans已经自动排好序
	//																				//  即set类型set<int>有自动排序的功能,vector<int>没有，我们initial()做了人工排序。
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
	//					if (count == goalitem.size())								// if成立表示匹配
	//					{
	//						Titem.push_back(j);				 
	//					}
	//					count = 0;
	//				}

	//				cout<<"找出了Titem   Titem.size()="<<Titem.size()<<endl;

	//				for(unsigned int j=0;j<Titem.size();j++)						// 为Titem中每个trans, 寻找匹配的large item, 每找到1个，计数器sum_item加1
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

	//										while (itemCII!=goalitem2.end() && transCII!=goaltrans2.end())  //  这样比对的前提是：goalitem 和 goaltrans已经自动排好序
	//																			//  即set类型set<int>有自动排序的功能,vector<int>没有，我们initial()做了人工排序。
	//										{
	//											if (*transCII == *itemCII)
	//											{
	//												count2++;
	//												itemCII++;
	//												transCII++;
	//											}
	//											else{transCII++;}
	//										}
	//										if (count2 == goalitem2.size())								// if成立表示匹配
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
	//				//将当前频繁项和它的overlap值，插入item_overlap映射

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

	//				Titem.clear();  //清空Titem

	//				cout<<"...";
	//			}

	//	}
	//}


	//outfile.close();

	//cout<<"\n>>>End of calculating overlap degree. Writing into file! \n\n";
}



//*******************************************************************
//       向文件写出统计信息。 包括： 
//		 database_size, MST, sup_u, sup_l,datafile, sen_file, 
//*******************************************************************

void write_stat()
{
	ofstream outfile;
	char filename[80];

	cout<<"\n>>> Begin to write the statistics information into file!\n";

	sprintf(filename,"%s_MST_%.3f_MCT_%.2f_MLT_%.2f_statistic.dat",datafrom, MST,MCT,MLT);
	outfile.open(filename, ios::out);

	outfile<<"多段编码--删除items\n\n";

	outfile<<"datafile:                 "<<datafrom<<endl;
	outfile<<"database size:            "<<db_size<<endl;
	outfile<<"Frequent Itemset(1) or Association Rule(2):"<<FIM_or_AR<<endl;;   //【FIM与AR选择开关】关联规则相关
	outfile<<"Minimum support threshold (MST):      "<<MST<<endl;
	outfile<<"Minimum confidence threshold (MCT):   "<<MCT<<endl;				// 关联规则相关
	outfile<<"Minimum lift threshold (MLT):         "<<MLT<<endl;				// 关联规则相关

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
	//新设计的交叉算子（ 洗牌方式交叉 shuffle random）, 保证子代基因无重复

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

    /*  //删除item不会减少trans数量， 故没有prelarge itemsets
	outfile<<"\nprelarge item sets LIST: \n";
	for (unsigned int i=0; i<prelarge_sets.size(); i++)												//对于n+1频繁项集
	{

		outfile<<"\nprelarge"<<i+1<<"-item sets LIST (num="<< prelarge_sets[i].size()<<") :\n";
		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	//对于n+1频繁项集的每一个频繁项
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
	for (unsigned int i=0; i<large_sets.size(); i++)									//统计频繁项数量，包括reallarge和prelarge
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
          result_ids[i] = add_individual(copy_individual(get_individual(selected[i])));  //返回值为新个体插入种群的位置（身份证号）

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
			   else if (recombination_type == 3)   //新设计的交叉算子（ 洗牌方式交叉 shuffle random）, 保证子代基因无重复
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
                    result = one_bit_mutation(get_individual(result_ids[i]));  //one_bit_mutation()搜索里有限， 推荐 indep_bit_mutation()
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



// 针对PPDM问题, mutation操作需要进行强有力的修改！

// 怎样生成的染色体中的编号不重复？ (通过<set>类型可以实现)
/* flip one bit at random position */
int one_bit_mutation(individual *ind)
{
     int position1, position2;

     if(ind == NULL)
          return (1);
     /* flip bit at position */

	position1 = irand(length1);
	ind->bit_string[position1] = prelarge_trans_array[ irand(prelarge_trans_array_size)] ;  

	 //ind->bit_string[length1+position2] = sen_trans_array[ irand(sen_trans_array_size)] ;  // 怎样保证与染色体中其他的编号不重复？

	 int offset=0;
	 for(int j=0;j<num_sen_item;j++)
	 {
			position2 = irand(encoding_len[j]);
			ind->bit_string[ length1+offset+position2 ] = sen_trans_array[j][ irand( sen_frequency_array[j]) ] ;	// 怎样保证与染色体中其他的编号不重复？
																										// 怎样保证与染色体中其他的编号不重复	
			offset=offset+encoding_len[j];
	 }
											

     return (0);
}


// 怎样生成的染色体中的编号不重复？ (通过<set>类型可以实现)
/* flip all bits with certain probability */
int indep_bit_mutation(individual *ind)
{

     double probability;
     /* absolute probability */
     probability = bit_turn_prob;  //param_file中为0.01,是否太小了， 实验确定（先改成0.1）
     if(ind == NULL)
          return (1);

	 //前半段：for prelarge_item 编码
	 random_shuffle(v1.begin(), v1.end());

     for(int i = 0; i < length1; i++)
     {
          if(drand(1) < probability)
          {
			  ind->bit_string[i] = prelarge_trans_array[ v1[i] ] ;			// 怎样保证与染色体中其他的编号不重复？random_shuffle()函数													
          }
     }

	 int offset=0;
	 for(int j=0;j<num_sen_item;j++)
	 {
		for(int i = 0; i < encoding_len[j]; i++)
		{
			if(drand(1) < probability)
			{
				ind->bit_string[ length1+offset+i ] = sen_trans_array[j][ irand( sen_frequency_array[j]) ] ;	// 怎样保证与染色体中其他的编号不重复？																									// 怎样保证与染色体中其他的编号不重复
			}
		}
		offset=offset+encoding_len[j];
	 }
     
     return (0);
}


// 怎样生成的染色体中的编号不重复？ (通过<set>类型可以实现)
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
          ind2->bit_string[i] = ind1->bit_string[i];		// 怎样生成的染色体中的编号不重复？
          ind1->bit_string[i] = bit_string_ind2[i];   		// 怎样生成的染色体中的编号不重复？		  
	 }

	 position2 = new int(num_sen_item);

	 for(int j=0;j<num_sen_item;j++)
		position2[j] = irand(encoding_len[j]);

	 int offset=0;
     for(int j=0;j<num_sen_item;j++) 
	 {
		for(i = length1+ offset ; i < length1+ offset+ position2[j]; i++) 
		{
			ind2->bit_string[i] = ind1->bit_string[i];		// 怎样生成的染色体中的编号不重复？
			ind1->bit_string[i] = bit_string_ind2[i];   		// 怎样生成的染色体中的编号不重复？
		}  

		offset = offset+encoding_len[j];
	 }

     free(bit_string_ind2);
	 free(position2);

     return(0);
}


// 怎样生成的染色体中的编号不重复？
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

     for(i = 0; i < length; i++) {							// 对于"二段编码"， uniform_crossover似乎不需要修改
          choose = irand(2);
          if(choose == 1) /* switch around bits */
          { 
               ind2->bit_string[i] = ind1->bit_string[i];			// 怎样生成的染色体中的编号不重复？ (通过<set>类型可以实现)
               ind1->bit_string[i] = bit_string_ind2[i];			// 怎样生成的染色体中的编号不重复？ (通过<set>类型可以实现)
          } /* else leave bit as is */   
     }  

     free(bit_string_ind2);

     return (0);
}



// 生成的子代染色体中的gene取值 不重复
// 新的交叉算子！！！！！！
// 新的交叉算子！！！！！！
int shuffle_crossover(individual *ind1, individual *ind2)
{
     int  i;

	 set<int> pre_combine_gene_set;
	 vector<int> pre_combine_gene_vec;

	 cout<<"\nshuffle_crossover\n";

	 //“前半段”Pre_large编码合并,从set转入vector
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



	 //“后半段”中的每个section中的gene合并
	 int offset=0;
     for(int j=0;j<num_sen_item;j++) 
	 {
		set<int> sen_combine_gene_set;
		vector<int> sen_combine_gene_vec;

		for(i = length1+ offset ; i < length1+ offset + encoding_len[j]; i++) 
		{
			sen_combine_gene_set.insert(ind1->bit_string[i]);			//将两个父体“后半段”第j个基因片段中的gene合并
			sen_combine_gene_set.insert(ind2->bit_string[i]);
		}  

		set<int>::const_iterator setCI2= sen_combine_gene_set.begin();	//将set导入vector
		while(setCI2 != sen_combine_gene_set.end())
		{
			sen_combine_gene_vec.push_back(*setCI2);
			setCI2++;
		}

		random_shuffle(sen_combine_gene_vec.begin(), sen_combine_gene_vec.end());	// 重新对合并的gene进行洗牌 【保证不重复的 关键所在！】

		for(i = length1+ offset ; i < length1+ offset + encoding_len[j]; i++) 
		{
			ind1->bit_string[i] = sen_combine_gene_vec[ i-(length1+ offset) ];
		}

		random_shuffle(sen_combine_gene_vec.begin(), sen_combine_gene_vec.end());	// 再次对合并的gene进行洗牌 【保证不重复的 关键所在！】

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



// 生成的子代染色体中的gene取值 不重复
// 新的变异算子！！！！！！
// 新的变异算子！！！！！！
int shuffle_mutation(individual *ind)
{
	cout<<"shuffle_mutation"<<endl;
	//前半段编码
	random_shuffle(v1.begin(), v1.end());

	for(int i=0;i<length1;i++)
	{
		ind->bit_string[i] = prelarge_trans_array[ v1[i] ];
	}

	//后半段编码
		
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




// 这个随机数发生器还有待改进！
// 这个随机数发生器还有待改进！
// 这个随机数发生器还有待改进！
/* Generate a random integer. */
//随机函数满足需要
int irand(int range)
{
     int j;

	 //Paramfile文件中已经有了种子，

	 //srand((unsigned)time(NULL));		//以当前时间作为随机种子

     j=(int) ((double) range * (double) rand() / (RAND_MAX + 1.0)); // RAND_MAX是32767

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
//       计算length(染色体编码长度),根据敏感项频次和sup_l 
//		 length = ( SIGMA ( frequency_i - sup_l) )*2
//*******************************************************************

void cal_length()
{
	//int sum1_prelarge =0,  sum2_sensitive = 0;
	//int length1=0,length2=0;							// 改用全局变量，因为length1的值代表了“二段编码”的分段位置

	cout<<">>>Begin to calculate length!\n";

	encoding_len = new int[ s_set.size() ];				// 用于存放“多段编码”中每段编码的长度；一段编码对应一个sensitive item.

	//先计算： 每个sensitive item对应的编码长度
	for(unsigned int i=0;i<s_set.size();i++)
	{
		if(sen_frequency_array[i]>sup_u)
		{
			encoding_len[i] = sen_frequency_array[i]-sup_l + 1;
			sum2_sensitive = sum2_sensitive + encoding_len[i] ;
		}
		else
		{
			cout<<"\nERROR in void cal_length(): 敏感项"<<i<<"不是频繁项"; getchar(); exit(-1);
		}
	}

	length2 = sum2_sensitive;

	if( MAX_del_trans < (unsigned)length2)  {
		cout<<"MAX_del_trans取值 过小 ， 敏感项编码长度超出！\n";
		cout<<"MAX_del_trans = "<<MAX_del_trans <<"    length2 = "<<length2<<endl;

		write_stat();

		getchar();

		length2 = MAX_del_trans ;  

		exit(-1);
	}


	//前半段编码长度 //删除item策略 就不需要此部分了！
	/*
	for (unsigned int i=0; i<prelarge_sets.size(); i++)												//对于n+1频繁项集
	{
		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	//对于n+1频繁项集的每一个频繁项
		{
			if(tempCI->second>=sup_l && tempCI->second<sup_u)
				sum1_prelarge = sum1_prelarge + (tempCI->second - sup_l) + 1;
			else{
				cout<< "prelarge项频次 应位于 [sup_l,sup_u)区间。 出错！！"; getchar(); exit(-1);
			}
		}
	}
	*/


	/*
	if (sum1*1.5 < (double)prelarge_trans_array_size)			length1 = (int)(sum1*1.5);
	else if (sum1*1.3 < (double)prelarge_trans_array_size)  	length1 = (int)(sum1*1.3);
	else if (sum1*1.1 < (double)prelarge_trans_array_size)		length1 = (int)(sum1*1.1);
	else {
		cout<<"\nMST 设置过小\n"; 
		getchar(); exit(-1);
	}
	*/

	/*
	if (sum1_prelarge  <= (double)prelarge_trans_array_size && sum1_prelarge<= MAX_del_trans - length2 )			
		length1 = (int)(sum1_prelarge);
	else {

		cout<<"\n pre_large段的编码长度 越界， 允许的最大值：" ;

		int temp_min =  ( (MAX_del_trans - length2)< prelarge_trans_array_size) ? (MAX_del_trans - length2): prelarge_trans_array_size;

		cout<<temp_min<<endl;

		length1 = temp_min;

		//getchar(); 
		
		//exit(-1);
	}
	*/
	sum1_prelarge = 0;
	length1 = 0;

	length = length1 + length2;   // 染色体编码总长度， 3个都是全局变量

	cout<<"\nprelarge项超阈值sup_l 频次总计 sum1="<<sum1_prelarge<<"  length1="<<length1<<endl;
	cout<<"敏感项超阈值sup_l 频次总计     sum2="<<sum2_sensitive<< "  length2="<<length2<<endl;
	cout<<"编码总长度  length = length1+length2 ="<<length<<endl;
	cout<<"\nactual_sup_l = "<<(db_size*MST)<<"\n";
	cout<<"original sup_l = "<<sup_l<<endl;

	sup_l_actual = (int)(db_size*MST);


	getchar();

	cout<<"\n>>>END of calculating length!\n";
}




//*************************************************************************************
//       功能：过滤出包含pre_large_set中item的记录; 计算每个pre_large项频次
//
//		 结果: 保存 “包含pre_large项记录的记录号”的prelarge_trans_array[]和prelarge_trans_set,
//*************************************************************************************

/*
       只需针对 pre_large 1-itemset, 而不必针对pre_large 2-item, 
	   因为包含pre_large 2-item 的transactions 一定也包含构成2-item的1-item.

*/
void filter_prelarge_trans()	//使用了“瘦身”数据库
{
	//apriori 算法 得到的 结果是: large_sets, reallarge_sets, prelarge_sets

	int count =0 ;

	cout<<"\n\nBEGIN:  void filter_prelarge_trans()\n";

	cout<<"\nprelarge_sets.size()="<<prelarge_sets.size()<<"\n";

	for (unsigned int i=0; i<prelarge_sets.size(); i++)												//对于n+1频繁项集	
	{
		if( prelarge_trans_set.size() >= db_size/5 ) break;


		cout<<"\n--------prelarge_sets["<<i<<"].size = "<<prelarge_sets[i].size()<<"\n";

		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	//对于n+1频繁项集的每一个频繁项
		{
				if( prelarge_trans_set.size() >= db_size/4 ) break;

				item goalitem = tempCI->first;

				for(unsigned int j=0;j<thin_DB.size();j++)// thin_DB 瘦身数据库
				{
					transaction goaltrans = thin_DB[j];   // thin_DB 瘦身数据库

					itemCI itemCI = goalitem.begin();										//item为set<int>, 里面的元素自动升序排列
					transactionCI transCI = goaltrans.begin();	

					while (itemCI!=goalitem.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
					{
						if (*transCI == (copy_new_code[*itemCI]-1) ){  // 【注意】：copy_new_code[]“译码”
							count++;
							itemCI++;
							transCI++;
						}
						else { transCI++;}
					}
					if (count == goalitem.size() )
					{
						//cout<<"匹配了一个 pre_large 项 "<< j  ;
						prelarge_trans_set.insert(j);										//使用set类型可以防止 插入重复记录号
						//cout<<"   当前prelarge trans 总数为: "<< prelarge_trans_set.size()<<endl;
					}
					count = 0;
				}

		}
	}

	//将prelarge记录 由集合prelarge_trans_set转存到数组prelarge_trans_array[]
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
	if(k!=prelarge_trans_array_size)		//检查计数是否出错！
	{
		cout<<endl<<"涉及prelarge项 的 transaction计数错误！！！"; getchar(); exit(-1);
	}

	//getchar();

	// 显示所有prelarge项，频次， 涉及敏感项的 数据库记录号
	cout<<"\nprelarge项数量:  "<<  prelarge_trans_set.size()<<endl;
	for (unsigned int i=0; i<prelarge_sets.size(); i++)										
	{
		for (itemsetCI tempCI=prelarge_sets[i].begin(); tempCI!=prelarge_sets[i].end(); tempCI++)	
		{
			cout<< "prelarge项：<";
			item goalitem = (tempCI->first) ;
			itemCI itemCI = goalitem.begin();
			while(itemCI!=goalitem.end())
			{
				cout<<" "<<*itemCI<<" ";
				itemCI++;
			}
			cout<<">    frquency频次："<< (tempCI->second)<<endl;
			goalitem.clear();
		}
	}
		
	/*
	cout<<"\n显示包含prelarge项的数据库记录号："<<endl<<"【总数为】："<< prelarge_trans_set.size() <<endl;
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
//       功能：过滤出包含敏感项的 transactions；计算每个敏感项频次  
//
//		 结果:  保存 “包含敏感项记录的记录号” sen_trans_array[]和sen_trans_set, 
//				以及保存了各敏感项出现频次的数组sen_frequency_array[]
//**************************************************************************************


void filter_sen_trans()
{													
	int count=0;

	vector<set<int>> sen_trans_vec;			// 敏感项涉及的trans集合， 使用set<int>类型可以防止 插入重复的元素
											// 每个敏感项对应一个sensitive transaction集合

	sen_frequency_array=(int *)new int[s_set_pair.size()];								//sen_frequency_array[]保存每个敏感项的出现频次

	cout<<"\n=======BEGIN to filter out transactions containing sensitive rules or itemsets=======\n";

	//模式匹配， 找出敏感记录，记录号存入集合sen_trans_set 
	//计算每个敏感项的出现频次，存入sen_frequency_array[]
	for(unsigned int j=0;j<s_set_pair.size();j++)				
	{																				//对于敏感集合(vector)中的每一行；每一行是一个item(即set<int>类型)	
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

		//s_set.push_back(goalitem);      //“将condition 和consequence合并”， 结果存入s_set
										// s_set类型是：  vector<set<int>> s_set;	
		item goalitem_encoding;
		basket_recode3(goalitem, goalitem_encoding);		//将goalitem改为使用内部码表示（trie树及瘦身数据库都使用内幕码）

		//for(unsigned int i=0;i<database.size(); i++)				
		for(unsigned int i=0;i<thin_DB.size(); i++)				//使用“瘦身”数据库(thin_DB中每行经过了排序处理)
		{
			//transaction goaltrans = database[i];
			transaction goaltrans = thin_DB[i];					//使用“瘦身”数据库			

			itemCI itemCI = goalitem_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
			transactionCI transCI = goaltrans.begin();	

			while (itemCI!=goalitem_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )	//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
				//if ( *transCI  == *itemCI )//原数据库
				if( *transCI  == *itemCI )   //尽管形式一样， 但这里两边使用的是内部码
				{
					count++;
					itemCI++;
					transCI++;
				}
				else { transCI++;}
			}
			if (count == goalitem_encoding.size() && goalitem.size()!=0)
			{
				//cout<<"敏感项匹配了一个Trans:  "<< i  ;
				sen_trans_set.insert(i);											//使用set类型可以防止 插入重复记录号
				
				//cout<<"   当前sensitive trans 总数为: "<< sen_trans_set.size()<<endl;
				sen_frequency_count++;
			}
			count = 0;
		}

		sen_trans_vec.push_back(sen_trans_set);			// 将支持 第j个敏感项 的 sensitive transactions ID集合 存入sen_trans_vec

		sen_frequency_array[j]=sen_trans_set.size();	// 将支持 第j个敏感项 的 transaction数量 存入全局数组 sen_frequency_array[ ]

		cout<<"\n\n>>>第"<<j<<"个敏感项的频次为："<<sen_frequency_count<< "  当前sen_trans_set.size()="<<sen_trans_set.size()<<"\n\n";

		sen_trans_set.clear();
	}

	if(s_set_pair.size()!=s_set.size()) {cout<<"\n 保存敏感项的两个容器 大小不一致！\n"; getchar();exit(-1);}

	for(unsigned int i=0;i<s_set_pair.size();i++)		// 检验机制
	{
		if(sen_trans_vec[i].size()!=sen_frequency_array[i])		{cout<<"\nERROR 2 in void filter_sen_trans()"; getchar();exit(-1);}
		if(sen_frequency_array[i]<sup_u)	
		{
			cout<<"sen_frequency_array[i]="<<sen_frequency_array[i]<<"  sup_u="<<sup_u;
			cout<<"\nERROR 3 in void filter_sen_trans():敏感项"<<i<<"不是频繁项. 请检查sensitive data file是否包含空行"; 
			getchar();exit(-1);
		}
	}

	//将敏感记录 由集合sen_trans_set转存到数组sen_trans_array[]
	int k=0;

	//int **sen_trans_array;			//多加的这个局部变量把握害苦了！ 仔细，不留下任何隐患！
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
			sen_trans_array[i][j]=*sen_trans_setCI;			//  得到一个重要数据结构： 二维动态数组（每行不等长）sen_trans_array[i][j]
			sen_trans_setCI++;								//  每行对应一个敏感规则， 每行内容是对应敏感规则的sensitive transactions ID
			j++;
		}
		if( j!= sen_frequency_array[i] )  {cout<<"\nsen_trans_array[][]每行的计数出错！！！\n "; getchar(); exit(-1);}
	}



	// 显示所有敏感项，频次， 涉及敏感项的 数据库记录号

	cout<<"\n敏感项数量:  "<<  s_set_pair.size()<<endl;
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
		cout<<  "  频次为： "<< sen_frequency_array[j]<<endl;
	}


	int total_sen_trans_num=0;
	for(unsigned int i=0;i<s_set_pair.size();i++)
	{
		cout<<"\n\n第"<<i<<"个敏感项涉及的trans数量："<<sen_trans_vec[i].size();
		cout<<"\n第"<<i<<"个敏感项涉及的trans如下：\n";

		for(unsigned int j=0;j<sen_trans_vec[i].size();j++)
		{		
			if(j%20==0) cout<<endl;
			cout<<" "<<sen_trans_array[i][j];
			total_sen_trans_num++;
		}
	}
	cout<<"\n\nSensitive trans总数为："<<total_sen_trans_num<<endl;

	cout<<"\n=======END of sensitive_trans filtering!=================\n";

	//getchar();
}


//****************************************************************
//			打开敏感项文件，读取到vector<item>类型的 s_set中
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
		cout << "无法打开dataset\n";	 

	while (!sensitive.eof())
	{
		getline(sensitive, line);					//读取一行， 注意文件里使用的是“商品编号”表示物品
		stringstream ss(line);	

		int consequence_flag = 0;

		while (ss >> token)
		{
			if(token.compare("->") == 0 )
			{
				consequence_flag = 1;
				cout<<"\n读到了 敏感规则 A->B  \n";
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
	[1]( （245,247）-> (277) )
	[2]( （281）-> (190)  )
	===========================================*/

	//进一步找到包含敏感itemsets的其他更多itemsets.如：用户指定AB敏感，则ABC,ABD都敏感
	if( FIM_or_AR == FIM )		
	{
		set<set<int>> all_include_condition;

		vector<pair<item,item>>::const_iterator sen_CI;
		for( sen_CI = s_set_pair.begin(); sen_CI != s_set_pair.end(); sen_CI++ )
		{
			condition.clear();
			consequence.clear();
			basket_recode3(sen_CI->first, condition); //编码
			if(sen_CI->first.size()!=condition.size())  {cout<<"read_sensitive_item()：有敏感项不是频繁项"; getchar(); exit(-1);}
				
			p_apriori->get_pointer_to_Apriori_Trie()->find_all_include_itemset(condition, all_include_condition, sup_u);

			//现在得到的all_include_condition是使用内部编码的， 需要解码！

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

	//if(1==gen)			// 在第一代验证敏感项是否读入成功
	{
		cout<<"\n验证sensative items sets read from file! \n";
		
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
				cout<<"敏感项的设置不是 关联规则， 但 FIM_or_AR == AR\n";
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

	//将s_set_pair中的condition与consequence合并， 存入s_set
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

		s_set.push_back(goalitem);      //“将condition 和consequence合并”， 结果存入s_set
										// s_set类型是：  vector<set<int>> s_set;	
	}									// s_set_pair类型： vector<pair<item,item>>


	num_sen_item = s_set_pair.size();
	
	if(s_set_pair.size()!=s_set.size()){ cout<<"s_set_pair.size()!=s_set.size() in void read_sensitive_item()"<<endl; getchar(); exit(-1); }

	sensitive.close();

	cout<<"\n>>>END of reading sensitive item set from file!\n";
}


/* 这处是算法修改的关键所在 */
/* 加入到参考点距离， 加入偏好信息*/
/* 点睛之笔！！！！！*/

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
			cal_fitness_by_trie_FIM(x);   //不管一个itemset包含的items有多少个，【通用】 利用trie树计算
		else if( FIM_or_AR == AR )
		{
			cal_fitness_by_trie_AR(x);    //不管一个rules包含的items有多少个，【通用】 利用trie树计算
		}
		else
		{
			cout<<"\n函数：eval(individual *x)----- FIM_or_AR参数设置出错， 候选项为{1,2}！\n";
			getchar();
			exit(-1);
		}	
	}

	//【注意】：已经验证 cal_fitness_by_trie_AR(x)和cal_fitness_two_rules(x)运行结果是一致的， 在2-item rules上

	return 0;


    //修改void Apriori_Trie::assoc_rule_find( ）， 将关联规则保存到 map<pair<item,item>,vector<int>> ruleset_bodon;

	//const Trie* Trie::is_included() 可以从树中 检索出 一个itemset 的频次（support）， 用法Trie::is_included()->counter;  
	//参考：void Apriori_Trie::assoc_rule_find()函数中用法： main_trie.is_included( condition_part, condition_part.begin() ) ->counter * min_conf
	
	//void Apriori_Trie::assoc_rule_assist(）告诉你如何遍历一棵树， 把所有层次的itemset取出来！

}

bool is_sensitive_rule(set<int> &condition, set<int> &consequence)
{
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		if(s_set_pair[j].first == condition && s_set_pair[j].second == consequence)
		{
			cout<<"\n当前项是一个敏感项";
			return true;
		}
	}
	return false;
}



//【通用】, 不管itemset包含的item有多少个   【精彩， 利用树和递归】
void cal_fitness_by_trie_FIM( individual *x )   
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;

	set<int> set_transID_delete;       // 采用set集合类型， 保证读进来的编号唯一，且自动排序, 另外可以去除重复基因！

	//根据chromosome编码 将选中的TID存入trans中
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //因为trans是set<int>类型，set可以自动排序（升序），更重要的是不会有重复元素
		}
		else
		{
			printf("ERROR in function eval() : ");		exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();


	/*----------针对选中待删除的每个trans，更新main_trie_copy树中对应结点的count-------------------*/

	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> goaltrans = thin_DB[*CI];		// 从"瘦身" 数据库 中取出对应transaction; 
													// 对于“多段编码”,因为选中的trans经过过滤，包含敏感规则； 因此goaltrans.size()>=2; 对于不带过滤的“一般编码”，可能选中的瘦身trans为空trans
		p_main_trie_copy->update_trie_tree(goaltrans.begin(), goaltrans.end());
	}

	//cout<<"\n\n打印更新的拷贝树！\n\n";
	//p_main_trie_copy->show_content_preorder();//打印树:对于一条basket,在树结构中查找可以匹配的所有"candidate_size -候选项"， 更新其counter值


	/*----------检索main_trie_copy(更新的copy树),决定alpha ----------------------------*/

	// 计算alpha, 对每个敏感itemsets, 检索树， 检查其counter值是否 >= (unsigned)sup_l_actual
	for(vector<item>::const_iterator sen_CI = s_set.begin(); sen_CI!=s_set.end();sen_CI++)  //检查每个 敏感itemset是否被隐藏掉（<sup_l_actual） 
	{
		item sen_itemset = *sen_CI;
		item sen_itemset_recoded; 

		basket_recode3(sen_itemset, sen_itemset_recoded);//编码
		unsigned long sen_counter = p_main_trie_copy->get_itemset_counter_from_trie(sen_itemset_recoded);		//get_itemset_counter_from_trie()是我自己写的函数

		cout<<"sen_counter "<<sen_counter<<endl;
		if(sen_counter >= (unsigned long)sup_l_actual)
			alpha++;
	}


	/*----------对比main_trie树(原树)和 main_trie_copy(更新的copy树), 决定beta, gamma--------*/
	//Trie * p_main_trie = p_apriori->get_pointer_to_Apriori_Trie()->get_pointer_to_main_trie();  //获得指向main_trie根结点的指针

	p_apriori->get_pointer_to_Apriori_Trie()->compare_two_trees(p_main_trie_copy, beta, gamma, sup_u, sup_l_actual);

	beta = beta - (s_set_pair.size()-alpha);  //compare_two_trees()求得的beta是所有原来频次大于等于sup_u,现在小于sup_l_actual的频繁项
											  //包含了被隐藏的敏感项， 因此需再减去 隐藏了的敏感项数量

	/*---------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
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


	// 将main_trie与main_trie_copy同步， 使main_trie_copy恢复更新前的值
	p_apriori->get_pointer_to_Apriori_Trie()->reset_tree_to_same(p_main_trie_copy);


}




//【通用】, 不管rule包含的item有多少个   【精彩， 利用树和递归】

void cal_fitness_by_trie_AR( individual *x )  
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;

	set<int> set_transID_delete;       // 采用set集合类型， 保证读进来的编号唯一，且自动排序, 另外可以去除重复基因！

	//根据chromosome编码 将选中的TID存入trans中
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //set可以自动排序（升序），更重要的是不会有重复元素
		}
		else
		{
			printf("ERROR in function eval() : ");	getchar();	exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();


	/*----------针对选中待删除的每个trans，更新main_trie_copy树中对应结点的count-------------------*/

	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> goaltrans = thin_DB[*CI];		// 从"瘦身" 数据库 中取出对应transaction;
													// 对于“多段编码”,因为选中的trans经过过滤，包含敏感规则； 因此goaltrans.size()>=2; 对于不带过滤的“一般编码”，可能选中的瘦身trans为空trans
		p_main_trie_copy->update_trie_tree(goaltrans.begin(), goaltrans.end());
	}

	//cout<<"\n\n打印更新的拷贝树！\n\n";
	//p_main_trie_copy->show_content_preorder();//打印树:对于一条basket,在树结构中查找可以匹配的所有"candidate_size -候选项"， 更新其counter值

	/*----------------------------------------------------------------------------------------------*/
	// 统计alpha
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		item sen_condition = s_set_pair[j].first;
		item sen_consequence = s_set_pair[j].second;
		item sen_union = s_set[j];

		item sen_condition_encode, sen_consequence_encode, sen_union_encode;

		basket_recode3(sen_condition, sen_condition_encode);//编码
		basket_recode3(sen_consequence, sen_consequence_encode);//编码	
		basket_recode3(sen_union, sen_union_encode);//编码

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
	// 统计beta和gamma
			
	map<pair<item,item>,vector<int>>::const_iterator CI = bodon_ruleset.begin();  //【注意】：bodon_ruleset使用的是内部编码
	for(; CI != bodon_ruleset.end(); CI++)
	{
			item condition = CI->first.first;			//注意：这里读出item使用内码.item就是set<int>类型. 
			item consequence = CI->first.second;		//注意：这里读出item使用内码
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

			// 统计beta和gamma
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
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
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

	// 将main_trie与main_trie_copy同步， 使main_trie_copy恢复更新前的值
	p_apriori->get_pointer_to_Apriori_Trie()->reset_tree_to_same(p_main_trie_copy);

}




void cal_fitness_two_rules(individual *x)
{
	int alpha = 0;
	int beta = 0;
	int gamma = 0;

	set<int> set_transID_delete;       // 采用set集合类型， 保证读进来的编号唯一，且自动排序, 另外可以去除重复基因！

	//根据chromosome编码 将选中的TID存入trans中
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //因为trans是set<int>类型，set可以自动排序（升序），更重要的是不会有重复元素
		}
		else
		{
			printf("ERROR in function eval() : ");		exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();


	//=============================================================================================
	//改进后的这段代码节省了时间！！！（耗时1秒左右）
	//=============================================================================================

	// 根据选中的trans, 更新copy_support_of_items_one2 和 copy_temp_counter_array2

	vector<unsigned long> copy_support_of_items_one2(copy_support_of_items_one);
	vector< vector<unsigned long> > copy_temp_counter_array2(copy_temp_counter_array);


	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> thin_trans = thin_DB[*CI];		// 从"瘦身" 数据库 中取出对应transaction; 
													
		vector<int>::const_iterator it1, it2;		

		for(it1=thin_trans.begin(); it1!=thin_trans.end(); it1++)
			copy_support_of_items_one2[*it1] -= 1;

		//更新“2-item”矩阵
		if(thin_trans.size()>=2)										//这个条件【非常重要】！否则,basket.size()==0， 会导致运行时异常！！！！！ thin_trans.size()必须>=1
		{																//对于“多段编码”， 因为选中的trans经过过滤，包含敏感规则； 因此size>=2; 对于不带过滤的“一般编码”，可能选中的瘦身trans为空trans
			for(it1=thin_trans.begin(); it1!=thin_trans.end()-1; it1++)		// thin_DB数据库的每一行已经排序（按“新编码”）
				for(it2=it1+1; it2!=thin_trans.end(); it2++)
					copy_temp_counter_array2[*it1][*it2-*it1-1]   -= 1;
		}
	}

	//---------------------- 统计alpha-----------------------------------
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		item sen_condition = s_set_pair[j].first;
		item sen_consequence = s_set_pair[j].second;

		itemtype sen_inter_code_condition   = copy_new_code[ *( sen_condition.begin() )  ] - 1;// 商品编号-》内码
		itemtype sen_inter_code_consequence = copy_new_code[ *( sen_consequence.begin()) ] - 1;

		int sen_new_condition_sup   =  copy_support_of_items_one2[ sen_inter_code_condition  ]  ;   
		int sen_new_consequence_sup =  copy_support_of_items_one2[ sen_inter_code_consequence ] ;  
		int sen_new_union_sup ;

		if( sen_inter_code_condition <  sen_inter_code_consequence   )
			sen_new_union_sup = copy_temp_counter_array2[ sen_inter_code_condition ][sen_inter_code_consequence - sen_inter_code_condition - 1] ;
		else if( sen_inter_code_condition > sen_inter_code_consequence   )
			sen_new_union_sup = copy_temp_counter_array2[ sen_inter_code_consequence][sen_inter_code_condition- sen_inter_code_consequence - 1];
		else
		{	cout<<"\n\nERROR in cal_fitness_two_rules(): (sen_inter_code_condition) and (sen_inter_code_consequence)应该不相等！\n   ";
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
	/*------------以下代码 对比新旧结构  统计beta和gamma-----------------*/

	//for(unsigned n=0;n<candi_rulesets.size(); n++)
	//{
		//map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		//for(; CI!=candi_rulesets[n].end(); CI++)
		map<pair<item,item>,vector<int>>::const_iterator CI = bodon_ruleset.begin();  //【注意】：bodon_ruleset使用的是内部编码
		for(; CI != bodon_ruleset.end(); CI++)
		{
			item condition = CI->first.first;			//item就是set<int>类型
			item consequence = CI->first.second;

			double old_union_sup_ratio = (double)(CI->second[2])/db_size;
			double old_rule_confidence = (double)(CI->second[2])/(CI->second[0]);
			double old_rule_lift = (double)(CI->second[2])/(CI->second[0])*(double)db_size/CI->second[1];

			//itemtype inter_code_condition   = copy_new_code[ *( condition.begin() )  ] - 1;// 商品编号-》内码
			//itemtype inter_code_consequence = copy_new_code[ *( consequence.begin()) ] - 1;
			itemtype inter_code_condition   =  *( condition.begin() ); //【注意】：bodon_ruleset使用的已是内部编码
			itemtype inter_code_consequence =  *( consequence.begin()) ;
			 

			int new_condition_sup   =  copy_support_of_items_one2[ inter_code_condition  ]  ;   
			int new_consequence_sup =  copy_support_of_items_one2[ inter_code_consequence ] ;  
			int new_union_sup ;

			if( inter_code_condition <  inter_code_consequence   )
				new_union_sup = copy_temp_counter_array2[ inter_code_condition ][inter_code_consequence - inter_code_condition - 1] ;
			else if( inter_code_condition > inter_code_consequence   )
				new_union_sup = copy_temp_counter_array2[ inter_code_consequence][inter_code_condition- inter_code_consequence - 1];
			else
			{	cout<<"\n\nERROR in cal_fitness_two_rules(): (inter_code_condition) and (inter_code_consequence)应该不相等！\n   ";
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

			// 统计beta和gamma
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


	//x->fitness = (double)(alpha*50+beta*25+gamma*25)/100.0;			//改为用多目标求解

	/*---------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
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


	set<int> set_transID_delete;       // 采用set集合类型， 保证读进来的编号唯一，且自动排序, 另外可以去除重复基因！

	//vector<vector<int>>   affect;

	//根据chromosome编码 将选中的TID存入trans中
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			set_transID_delete.insert(x->bit_string[i]);  //因为trans是set<int>类型，set可以自动排序（升序），更重要的是不会有重复元素
		}
		else
		{
			printf("ERROR in function eval() : ");
			exit(-1);
		}
	}

	unsigned actual_num_del = set_transID_delete.size();

	//=============================================================================================
	//改进后的这段代码节省了时间！！！（耗时1秒左右）
	//=============================================================================================

	// 这两个数据结构 用于 与原来的结构 比较，进而得到alpha,beta,gamma

	vector<unsigned long> copy_support_of_items_one2(copy_support_of_items_one);
	vector< vector<unsigned long> > copy_temp_counter_array2(copy_temp_counter_array);


	for (set<int>::const_iterator CI=set_transID_delete.begin();  CI!=set_transID_delete.end();  CI++)
	{
		vector<int> thin_trans = thin_DB[*CI];		//从"瘦身" 数据库 中取出对应transaction
		vector<int>::const_iterator it1, it2;

		for(it1=thin_trans.begin(); it1!=thin_trans.end(); it1++)
			copy_support_of_items_one2[*it1] -= 1;

		//更新“2-item”矩阵
		if(thin_trans.size()>=2)		//此条件【非常重要】！找一整天才找到此bug！否则,basket.size()==0， 会导致运行时异常！！！！！ thin_trans.size()必须>=1
		{								//对于“多段编码”， 因为选中的trans经过过滤，包含敏感规则； 因此size>=2; 对于不带过滤的“一般编码”，可能选中的瘦身trans为空trans
			for(it1=thin_trans.begin(); it1!=thin_trans.end()-1; it1++)
				for(it2=it1+1; it2!=thin_trans.end(); it2++)
					copy_temp_counter_array2[*it1][*it2-*it1-1]   -= 1;
		}
	}

	/*------------以下代码 对比 一维 和 二维 新旧结构-----------------*/

	// "1-item"及"2-item"的alpha计算
	for(vector<item>::const_iterator sen_CI = s_set.begin(); sen_CI!=s_set.end();sen_CI++)  //检查每个 敏感itemset是否被隐藏掉（<sup_l_actual） 
	{
		item sen_item = *sen_CI;
		if(sen_item.size()== (unsigned int)1)   // "1-item"敏感项
		{
			int internal_code =  copy_new_code[*(sen_item.begin())]-1;		//“商品编码” 转化成 “内部码”
			if( copy_support_of_items_one2[internal_code] >= (unsigned)sup_l_actual  )
				alpha++;
		}
		
		if(sen_item.size()== (unsigned int)2)   // "2-item"敏感项
		{
			unsigned int code1, code2;
			code1 = copy_new_code[*(sen_item.begin())]-1;					//“商品编码” 转化成 “内部码”
			code2 = copy_new_code[*(sen_item.end())]-1;

			if(copy_temp_counter_array2[code1][code2-code1-1] >= (unsigned)sup_l_actual  )
			{
				alpha++;
				cout<<"\n Alpha is increasing!   Alpha = "<<alpha;
			}
		}
	}

	// "1-item"的beta和gamma计算
	for(unsigned i=0; i<copy_support_of_items_one.size(); i++)
	{
		if( copy_support_of_items_one[i]< (unsigned)sup_u  &&  copy_support_of_items_one2[i]>= (unsigned)sup_l_actual )   // “不频繁”--》“频繁”
		{
			gamma++;  cout<<"\n Gamma is increasing(1-item)!   Gamma = "<<gamma;
		}

		if( copy_support_of_items_one[i] >= (unsigned)sup_u &&  copy_support_of_items_one2[i]< (unsigned)sup_l_actual )  // “频繁”--》不频繁(即隐藏)
		{
			//先检查是不是敏感项
			int flag=0;
			item temp_item;
			temp_item.insert( (int)(copy_new_code_inverse[i]) );

			for(vector<item>::const_iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)				
			{		
				if ( temp_item == *s_setCI) {cout<<"\n 是1-item敏感项"; flag =1; break;}
			}

			if(0 ==flag )  // "非敏感"项 被隐藏
			{
				beta++;  cout<<"\n Beta is increasing(1-item)! Beta = "<<beta;
			}	
		}
	}

	// "2-item"的beta和gamma计算
	for(unsigned i=0; i<copy_temp_counter_array.size(); i++)
		for(unsigned j=0; j<copy_temp_counter_array[i].size(); j++)   //copy_temp_counter_array[i].size() == copy_temp_counter_array.size()+1 - (i+1)
		{
			if( copy_temp_counter_array[i][j]<sup_u && copy_temp_counter_array2[i][j]>=sup_l_actual )   // “不频繁”--》“频繁”
			{
				gamma++;  cout<<"\n Gamma is increasing (2-item)!   Gamma = "<<gamma;
			}
			if( copy_temp_counter_array[i][j]>=sup_u && copy_temp_counter_array2[i][j]<sup_l_actual )	// “频繁”--》不频繁 (即隐藏)
			{
				//先检查是不是敏感项	
				item temp_item;
				temp_item.insert( (int)(copy_new_code_inverse[i]) );		// “内部码”-->“商品编号”
				temp_item.insert( (int)(copy_new_code_inverse[j+i+1]) );    // 注意：对于j下标， 需要将矩阵的“存储下标”转化成“逻辑下标”；
																				// 参：Apriori_Trie.cpp文件中find_candidate_two()函数下面的注释
				int flag=0;
				for(vector<item>::const_iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)				
				{		
					if ( temp_item == *s_setCI) {cout<<"\n 是2-item敏感项"; flag =1; break;}
				}
				if(0==flag)  // "非敏感"项 被隐藏
				{
					beta++; cout<<"\n Beta is increasing(2-item)! Beta = "<<beta;
				}
			}
		}



	// 注意： 读取sensitive item的工作已经移到了 state0()函数中！！！


	//x->fitness = (double)(alpha*50+beta*25+gamma*25)/100.0;			//改为用多目标求解
	/*---------------------------------------------------------------------------------------*/
	if(4 == dimension )
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
		x->objective_value[1] = (double)beta;
		x->objective_value[2] = (double)gamma;
		x->objective_value[3] = (double)actual_num_del ;
	}
	else if(3 == dimension)
	{
		x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
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


	set<int> trans;    // 采用set集合类型， 保证读进来的编号唯一，且自动排序
	DB affect;

	//將要計算fitness所有的TID存入trans中
	for (int i=0; i<length; i++)
	{
		if ( (x->bit_string[i] >= 0)  &&  (x->bit_string[i]<  db_size)  )
		{
			trans.insert(x->bit_string[i]);  //因为trans是set<int>类型，set可以自动排序（升序），更重要的是不会有重复元素
		}
		else
		{
			printf("ERROR in function eval() : ");
			exit(-1);
		}
	}

	//將所選的Transaction所包含的item存入affect
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
	collectset affectlarge = large_sets;          // 复制large_sets；	large_sets既包含real_large,也包含pre_large


	/*=========================================
	[0]( [0]( (234),345)),  [1]((245),567),...)
	[1]( [0]( (234,245),212) ),  ...)
	[2]( [0]( (234,245,345), 230 ), ...    )
	===========================================*/
	//collectset affectreallarge = reallarge_sets;
	//collectset affectprelarge = prelarge_sets;

	// affectlarge保存了所有large item sets,包括1-item,2-item,3-item,...
	// affectlarge的每一行存储了n-item的频繁项(large item sets)
	// 下面的四层循环逻辑是这样的：
	//		针对affectlarge的每一行的每一个频繁项（对应最外两层for）, 
	//			扫描affect数据库(根据基因选出transaction构成的微缩版database)的每一行(对应从外往里，第三层for)
	//				将这个频繁项(goalitem)与affect数据库的每行进行比对(对应while循环)， 
	//				如果某个transaction（某行）包含了item set,则affectlarge[n][goalitem]--;
	//				即如果准备删除的transaction包含敏感性，则敏感项计数减1（可以看出，敏感项是频繁的，即是large item sets）

	// 也就是，对于每一个large item,扫描准备删除的transaction构成的数据库的每一行，进行比对，若包含这个频繁项，就更新它的计数值。

	//====================================================================================================
	//下面这段代码非常耗时，相当于重新扫描500行(length=500) 的数据库，再进行模式匹配！！！（耗时12秒左右）
	//====================================================================================================
	//取出affect数据库中每一個Transaction，与每一个频繁项匹配

	for (unsigned int n=0; n < affectlarge.size(); n++)
	{
		//取出large_sets中每一個item	
		for (itemsetCI tempCI=affectlarge[n].begin(); tempCI!=affectlarge[n].end(); tempCI++)
		{
			item goalitem = tempCI->first;
			//affect存的是要評估TID的所有items
			for (DBCI db_CI=affect.begin(); db_CI!=affect.end(); db_CI++)
			{
				transaction goaltrans = *db_CI;				//affect的transaction;  transaction不会自动排序，但database中每行排了序
				itemCI itemCI = goalitem.begin();			//large_sets的item ; item为set<int>, 里面的元素自动升序排列
				transactionCI transCI = goaltrans.begin();	
				//一個item一個item比對
				while (itemCI!=goalitem.end() && transCI!=goaltrans.end())  //  这样比对的前提是：goalitem 和 goaltrans已经自动排好序
																			//  即set类型set<int>有自动排序的功能,vector<int>没有，我们initial()做了人工排序。
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
					affectlarge[n][goalitem]--;					//删除了length个transactions,只会使部分频繁集的频次降低,而不是增加
					//affectreallarge[n][goalitem]--;
					//affectprelarge[n][goalitem]--;
					//cout<< "UPDATE affectedlarge" << endl;
				}
				count = 0;
			}
		}
	}


	// 注意： 读取sensitive item的工作已经移到了 state0()函数中！！！


	// 针对更新后的large item表每一行的每一个频繁项(更新后的计数仍大于sup_l)，（对应外层两个for循环）
	// 如果其存在于敏感项容器s_set中，说明通过删除这些基因指定的transaction后仍未能隐藏这个敏感项，即这个敏感项仍是频繁的（大于sup_l）
	// 则 alpha++;

	// 【注意】： large_sets包含了reallarge_sets和prelarge_sets


	for (unsigned int n=0; n<affectlarge.size(); n++)											//对于n+1频繁项集
	{
		for (itemsetCI tempCI=affectlarge[n].begin(); tempCI!=affectlarge[n].end(); tempCI++)	//对于n+1频繁项集的每一个频繁项
		{
				//如果是敏感项
				int sen_flag=0;
				for(vector<set<int>>::iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)				//对于敏感集合(vector)中的每一行
				{																				//每一行是一个item(即set<int>类型)
					//cout<<"\n走到这里，下面的if判断有问题\n";									
					if (tempCI->first==*s_setCI)
					{
						cout<<"\nAlpha！Alpha！ "<<" sup_l_actual ="<< sup_l_actual <<"  (tempCI->second)="<<(tempCI->second)<<endl ;

						if((tempCI->second) >= (db_size-length)*MST  )											//【重要修改】：此处所有sup_l用(db_size-length)*MST代替
						{
							alpha++;    // alpha表示未隐藏的敏感项的数量，
										// 即如果删除了指定Transaction后， 对应的large item表(更新后的)仍然包含
							sen_flag=1;
						
							cout << "ALPHA is increasing:  "<<"alpha = " << alpha << "   curGen:  "<< gen << endl;;
						}
					}
				}
				if(1== sen_flag)   //如果affectlarge表中第n行的当前项是“敏感项”，跳过下面判断，进入下一次循环
					continue;

			
				// large_sets是删除基因指定的transaction之前的频繁项表； affectlarge是删除它们之后的频繁项表

				for (itemsetCI tempCII=large_sets[n].begin(); tempCII!=large_sets[n].end(); tempCII++)
				{   // tempCII为原来的频繁项表large_sets的指针； tempCI 为根据删除的transactions 更新值后的频繁项表的指针

					if ((tempCI->first==tempCII->first) && ((tempCII->second)>=sup_u)   && (  (tempCI->second)< sup_l_actual  )  )
					{		
						// 某项原来是频繁(real_large)的， 删除指定transactions之后， 变成不频繁的（注意：应考虑部分敏感项就是这种情况）
						// missing large items
					
						//此处少了是否为敏感项的判断
						//-------------------------以下代码由Cheng Peng添加-----------------------------------

						//cout<<"\n 某项原来是频繁的， 删除指定transactions之后， 变成不频繁的！\n";
						//getchar();

						int flag=0;
						for(vector<set<int>>::iterator s_setCI=s_set.begin();s_setCI!=s_set.end(); s_setCI++)
						{
							if (tempCI->first==*s_setCI)  //说明这个由频繁变为不频繁的项是敏感项
							{
								flag=1;break;
							}
						}

						if( 0==flag ) //flag==0说明这个项不是敏感项，但删除指定transaction之后，由频繁变成不频繁
						{
							beta++; 
							cout<<"Beta! Beta! "<<" sup_l_actual ="<< sup_l_actual <<" (tempCI->second)="<<(tempCI->second)<< "  (tempCII->second)= "<<(tempCII->second)<< endl ;
							cout<<"\nBETA is increasing:  " << " beta == "<< beta << "   curGen:  "<< gen << "\n";
						}

						//---------------------------------------------------------------------------------------	  
					}

					if(   (tempCI->first==tempCII->first) && ((tempCII->second)<sup_u) && ((tempCI->second)>=sup_l_actual  )  )
					{				
								//实际上即使pre_large 集合里的item 在删除trans后， 仍然>=sup_l
								// tempCI所属的频繁集表中， 就已经说明是频繁项了， 因此不可能(tempCII->second)<sup_u)
								// 这里需要重新考虑! 难道上面条件有错？？ 注意：原来的large_sets包含real_large和pre_large.
								// 因此，原来小于sup_u同时大于sup_l的频繁项，可能在删除length个transactions之后，仍大于sup_l
								// 因此，删除的笔数length不能太多，否则必然导致频次位于[sup_l,sup_u]区间的频繁项太多， gamma值过大。

								
								// 某项原来是不频繁的(即不是real_large,而是pre_large)， 删除指定transactions之后， 变成频繁的
								// 敏感项可能也有这种情况， 但前面的continue语句表示：
								// 如果是敏感项，且删除transaction后仍是频繁的（(tempCI->second)>=sup_l)）,
								// 则alpha++后，跳过后面另外两种情况的判断
								// 因此能执行到这里， 说明当前项 (*tempCI) 是“非敏感的”
					
								// generate new frequent items which doesn't exist before; 
								// But why the condition is ((tempCII->second)<sup_u)  && ((tempCI->second)>=sup_l)
						gamma++;

						
						cout<<"Gamma! Gamma! "<<" sup_l_actual ="<< sup_l_actual <<" (tempCI->second)="<<(tempCI->second)<< "  (tempCII->second)= "<<(tempCII->second)<< endl ;
						cout << "\n\n GAMMA is increasing:  gamma = " << gamma << "   curGen:  "<< gen << endl;
					}
				}
		}
	}


	//x->fitness = (double)(alpha*50+beta*25+gamma*25)/100.0;			//改为用多目标求解

	x->objective_value[0] = (double)alpha;   //注意我们的目标向量是double型， 而alpha, beta, gamma是int型
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
	cout<<"\n>>>进入my_random_shuffle()\n";

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
	cout<<"\n>>>退出my_random_shuffle()\n";
} 


/* create a random new individual and allocate memory for it,
   returns a pointer to the new individual */
individual *new_individual()
{
     individual *return_ind;
     int result;

     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->bit_string = (int *) malloc(sizeof(int) * length);  // 注意length为“二段编码”总长度

	 // 仍然为目标向量申请double型的空间
	 return_ind->objective_value = (double *) malloc(sizeof(double) * dimension);


	 /*
	 //前半段编码
	random_shuffle(v1.begin(), v1.end());

	cout<<"\n前半段编码:\n";
	 for(int i=0;i<length1;i++)
	 {
		 return_ind->bit_string[i] = prelarge_trans_array[ v1[i] ];
		 //cout<<v1[i]<<" ";
	 }
	 */

	 // 后半段编码
	for(unsigned int i=0;i<s_set.size(); i++)
	{
		random_shuffle( v2[i].begin(),v2[i].end() );
	}

	cout<<"\n后半段编码：\n";
	 int offset=0;
	 for ( unsigned int i=0; i < s_set_pair.size(); i++)
	 {
		 cout<<"\n第"<<i<<"个敏感项：\n";
		 for( int j=0; j < encoding_len[i]; j++ )
		 {
				
				//int pos;
				
				//pos = irand( sen_frequency_array[i]  );  // 这里没有要求产生【无重复】的数据库记录号 (通过<set>类型可以实现无重复)
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

	 cout<<"\n 评价了一个第一代新个体";
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


/* 注意：这里有两个函数都实现输出：

		1) 一个是write_output_file(), 写个体编号、目标向量、决策向量到文件

		2）另一个是write_output_obj()，写个体编号、目标向量到文件，不包括决策向量。
		
*/
/* Writes the index, objective values and bit string of
   all individuals in global_population to 'out_filename'. */
void write_output_file()
{
     int j, current_id;
     FILE *fp_out;
     individual *temp;

	 char outfile_name[128];     //注意本地参数文件中的 outputfile 参数已经废弃掉
	 	 
	 sprintf(outfile_name,"%s_MST_%.3f_MCT_%.2f_MLT_%.2f_output.dat",datafrom, MST,MCT,MLT);
     
     fp_out = fopen(outfile_name, "w");  
     assert(fp_out != NULL);

	 /*******************************************************************/
	 //写入相关参数到结果文件



	 fprintf(fp_out, "Length	    = %d \n", length); 
	 fprintf(fp_out, "MaxGen	    = %d \n", maxgen); 
	 fprintf(fp_out, "datafrom	    = %s \n", datafrom); 
	 fprintf(fp_out, "sensfile	    = %s \n", sensfile); 

	 fprintf(fp_out, "FIM(1)AR(2)   = %.3f \n",FIM_or_AR );  
	 fprintf(fp_out, "MST	        = %.3f \n",MST );   // 关联规则相关
	 fprintf(fp_out, "MCT	        = %.3f \n",MCT );	// 关联规则相关
	 fprintf(fp_out, "MLT	        = %.3f \n",MLT );	// 关联规则相关

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
				printf( "%.2lf ",   get_objective_value(current_id, j)    ); //因为适应度评价时，做了多目标最小化处理，现在还原
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
				printf( "%.2lf ",   get_objective_value(current_id, j)    ); //因为适应度评价时，做了多目标最小化处理，现在还原
				fprintf(fp_out, "%.2lf ",  get_objective_value(current_id, j)  );
			}
            fprintf(fp_out, "\n");

			current_id = get_next(current_id);
     }



	 /*******************************************************************/
	 //写入相关参数到结果文件

	 fprintf(fp_out, "\n\n\n");

	 fprintf(fp_out, "Length	    = %d \n", length); 
	 fprintf(fp_out, "MaxGen	    = %d \n", maxgen); 
	 fprintf(fp_out, "datafrom	    = %s \n", datafrom); 
	 fprintf(fp_out, "sensfile	    = %s \n", sensfile); 

	 fprintf(fp_out, "FIM(1)AR(2)   = %.3f \n",FIM_or_AR );  
	 fprintf(fp_out, "MST	        = %.3f \n",MST );   // 关联规则相关
	 fprintf(fp_out, "MCT	        = %.3f \n",MCT );	// 关联规则相关
	 fprintf(fp_out, "MLT	        = %.3f \n",MLT );	// 关联规则相关

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



//将 第一代种群的目标空间向量写入文件
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
				printf( "%.2lf ",   get_objective_value(current_id, j)    ); //因为适应度评价时，做了多目标最小化处理，现在还原
				fprintf(fp_out, "%.2lf ",  get_objective_value(current_id, j)  );
			}
            fprintf(fp_out, "\n");

			current_id = get_next(current_id);
     }


	 fclose(fp_out);

	 printf("\nWriting objective vectors in the first generation to the file \"%s\" is complete!",output_obj_file);


}
