
/*	apriori算法的输出就是 collectset类型的频繁集表 large_sets  

	GA算法的关键在于个体适应度评价cal_fitness()：
	
		根据chromosom指定的transactions组建一个小database, 
	
		对于apriori算法传来的频繁集表large_sets中每一行上的每一项，与小database的每一行（每一个trans）比对，
		
		若trans包含频繁项item,则对应item的计数减1， 这样， 依据小database, 频繁集表collectset得到了更新。

		接下来，两种可能

		1）	更新后的频繁集表中仍然包含敏感项：
		
			定义敏感项集 sensitemset（类型为vector<item>， item为set<int>）, 其值从文件“sensfile”中读取 , 
		
			对于敏感项集中的每一个item, 与 更新后的频繁集表large_sets中每一行的每一个频繁项 进行比对， 看其是否包含敏感项

			若频繁集表中某行某项 等于 敏感item 且 其计数>= sup_l（为删除这些trans后仍是频繁的）， 则说明仍有频繁的敏感item，alpha加1.

		
		2） 更新后的频繁集表中不含敏感项，又有两种情况：
		
			2.1） 更新前 某非敏感项 是频繁的， 更新后变的不频繁（这种情况，我们只希望仅发生在 敏感且频繁 的项上 ） ， beta加1

				  因此注意：判定时应考虑 排除 敏感项 由频繁 变成 不频繁 的情况， 因为这是应该的， 我们目的本来即是此。 

			2.2） 更新前 某非敏感项 是不频繁的， 更新后变的频繁。gamma加1.
		
			以上2.1）和 2.2）是通过 对比 更新前和更新后的频繁集表collectset来 达到判定的

	

	备注：频繁集表large_sets 第0行存有所有的频繁1-item set， 第1行存有所有的频繁2-item set，  第2行存有所有的频繁2-item set， 以此类推... ...

				==========the storage structure of collectset============
				[0](  [0]( (234),345 ),      [1]( (245),567 ),		...)
				[1](  [0]( (234,245),212 ),  [1]( (232,235), 290 ),	...	)
				=========================================================
*/

/////////////////////////////////////////////////////////////

#ifndef _apriori_
#define _apriori_

#include <iostream>
#include <fstream>

#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <ctime>

#include <utility>      // std::pair, std::make_pair     关联规则产生用到


using namespace std;
/*----------Added for PPDM_dt --------*/

#include "variator_user.h"

#include "apriori_trie\common.hpp"

#include "apriori_trie\Apriori_main.hpp"			// 调用bodon 的apriori快速算法！

//#define low_sup_thres 0.07


extern double MCT_discount;

extern double MST_discount;


extern char datafrom[128];			// datafrom是数据源文件， 它的值从PPDM_dt_param文件中读取
extern char sensfile[128];			// 包含敏感项的文件名,它的值从PPDM_dt_param文件中读取

#define outputfile "out.txt"


//-----------------------------------------------------------------------

typedef set<int> item;						//用来记录一个item的资料形态  							
typedef item::const_iterator itemCI;		//item专用的iter


typedef map<item, int> itemset;				//用来记录item的集合
typedef itemset::const_iterator itemsetCI;	//itemset专用的iter
typedef map<item, int>::value_type newdata;	//新增itemset用:new pair
 
typedef vector<itemset> collectset;			//用来记录多个itemset的集合

		/*----------the storage structure of collectset----------
		[0](  [0]( (234),345 ),      [1]( (245),567 ),		...)
		[1](  [0]( (234,245),212 ),  [1]( (232,235), 290 ),	...	)
		---------------------------------------------------------*/
typedef vector<int> transaction;			//用来记录数据库中一行
typedef transaction::const_iterator transactionCI;

typedef vector<transaction> DB;				//用来记录整个database
typedef DB::const_iterator DBCI;	

//-----------------------------------------------------------------

extern int sup_u ;
extern int sup_l ;
extern unsigned long db_size; 

extern int sup_l_actual;;


//注意：.cpp文件定义了还要在.h文件声明：
extern collectset candi_sets;
extern collectset large_sets;
extern collectset reallarge_sets;
extern collectset prelarge_sets;

extern DB database;


/* ----------------association rules----------------------*/

#include<iomanip>

#define FIM 1
#define AR 2


extern unsigned int SIZE_threshold;

//#define MST 0.022
extern double MST;			// 【MST 最小支持度】，它的值从PPDM_dt_param文件中读取
extern double MCT;			// 【MCT 最小置信度】，它的值也从PPDM_dt_param文件中读取
extern double MLT;			// 【MLT 最小LIFT】,它的值也从PPDM_dt_param文件中读取	


extern unsigned int FIM_or_AR;			// 【程序切换开关】,取值为1表示频繁模式挖掘FIM, 取值为2表示惯量规则挖掘


//#define  ALPHA 3

extern int	num_strong_rules;

extern int  num_strong_rules_by_trie;

extern int  num_itemset_by_trie;

/*===========================Trie 树========================*/


extern vector<vector<int>>  thin_DB ;

extern vector<itemtype> copy_new_code_inverse;

extern vector<itemtype> copy_new_code;

extern vector<unsigned long> copy_support_of_items_one;

extern vector< vector<unsigned long> > copy_temp_counter_array;

//--------------------------------------------------------------------------


extern map<pair<item,item>,vector<int>> bodon_ruleset;  //在Apriori_Trie::assoc_rule_find()中存入Bodon apriori发现的关联规则

//extern map<pair<item,item>,vector<int>> bodon_ruleset_strong;

typedef struct A_RULE_NODE
{
	//pair<item,item> a_rule;
	unsigned int antecedent_supp;
	unsigned int consequent_supp;
	unsigned int union_supp;

	double conf;

}T_RULE_UNKNOWN;

extern map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_strong;				//保存rule的最小可能Supp,最小可能Conf

extern map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_NBRS;					//保存rule的最大可能Supp,最大可能Conf



extern map<pair<item,item>, T_RULE_UNKNOWN>  bodon_ruleset_all_MIN;					//保存所有rule的最小support和confidence

extern map<pair<item,item>, T_RULE_UNKNOWN>  bodon_ruleset_all_MAX;					//保存所有rule的最大support和confidence


typedef map<pair<item,item>, T_RULE_UNKNOWN>::value_type  UNKNOWN_insert_type;



extern map<item, int>  bodon_itemset;

extern map<item, int>  bodon_itemset_frequent;


/*----------------------------------------------------------*/


typedef map<pair<item,item>,vector<int>> ruleset;		// 第一个item是condition(X表示),第二个item是consequence(Y表示), 
														// vector<int>包含sup(X), sup(Y), sup(XUY)

typedef map<pair<item,item>,vector<int>>::value_type   ruleset_insert_type;  // 用于map的insert成员函数，需用此类型参数

typedef map<item, int>::value_type   bodon_itemset_insert_type; 


//------------------------------------------------------------

typedef vector<ruleset> collect_ruleset;

extern collect_ruleset  candi_rulesets;
//------------------------------------------------------------



// 敏感项集合
typedef vector< set<int> > sensitemset;	// 【注意】： 这里采用新的数据结构

typedef sensitemset::const_iterator snesCI;

extern vector<set<int>> s_set;					// 敏感项集合(由s_set_pair转换来)，它的内容在 filter_sen_trans() 函数中 由s_set_pair导入

extern vector<pair<item,item>> s_set_pair;		// 敏感项集合， 包含了从文件读来的敏感项item,

/*
   例如： 存放敏感的文件内容格式如下：
   1->2
   5->3
   则s_set_pair存放格式为：
   vector[0]: ((1),(2))
   vector[1]: ((5),(3))
   则s_set存放格式为：
   vector[0]: (1,2)
   vector[1]: (3,5)
*/


extern int num_sen_item;

extern	int * sen_frequency_array;		// 保存了每个敏感项的出现频次；也是在filter_sen_trans()中申请内存空间


/*
extern 	int ** sen_trans_array;			// 在filter_sen_trans()中申请内存空间,保存每个sen item涉及的trans编号
extern 	int * prelarge_trans_array;		// 从prelarge_trans_set转存过来，保存了包含敏感项的数据库记录在 数据库中的序列号
extern 	int *encoding_len;	
*/

extern unsigned long MAX_del_trans;					//定义位于apriori.cpp   ; MAX_del_trans是底线， 决定了sup_l,但真正的删除笔数是length(即编码长度)

extern double Avg_Dist_MST;					// It defines the average distance of rules set to MST

extern double MAX_del_trans_perc;		//定义位于apriori.cpp

extern double  ALPHA;


//----------------------------function declaration-----------------------------


void initial();	// load data from file database into memory vector<vector<int>>


void cal_support(); // compute sup_l and sup_u


void gen_candione();  // count the frequency of one item


void gen_large(itemset inputset);  //count and select large one item sets


void join(int candi_n);  // form large candidate sets n-1 generate candidate sets n


void scanDB(int candi_n);  // scan database to calculate the frequency of each item sets in candidate sets n 


void output(itemset outputset, int n, int type); // output the content of item sets (one row in large sets "collectset") into outputfile


void load_one_candidate();  // in order to improve efficiency, directly load large 1-item candidate sets, but no computing from database

	
void load_two_candidate();	// in order to improve efficiency, directly load large 2-item candidate sets, but no computing from database



void apriori_large_sets();  // generate all large item sets


void apriori_bodon();		// 调用bodon的apriori算法， 返回其计算得到的 1-频繁， 2-频繁集合； 构建“瘦身”数据库 

void verify_update_trie_tree(); 

/* ----------------association rules----------------------*/

void apriori_gen_rules();   // 由频繁集 产生 rules或准rules

void output_candi_rules();

void output_real_rules();

void print_candi_rules();

void print_one_rule(item condition, item consequence, int condition_sup, int consequence_sup, int union_sup, double lift);


/*-------------------------------------------------------*/


/**------------------------Trie--------------------------*/


void basket_recode3( const set<int>& original_basket, set<int>& new_basket );  //set -> set  编码

void basket_decode3( const set<int>& encoding_basket, set<int>& decoding_basket ); //set -> set  解码

extern vector< vector<unsigned long> > copy_temp_counter_array;

extern Trie *p_main_trie_copy;


#endif