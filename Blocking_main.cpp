

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "Blocking_main.h"

using namespace std;


//-------------------BBRH-Block  IJCAI 2019---------------------

typedef struct NOTDE3{
	unsigned long A_MST;
	unsigned long A_MCT;
}A_MST_MCT;


typedef struct NODE_BORDER_RULE{
	unsigned A_MST;				// support - db_size*MST	[maintain the original value during updating the weight]
	unsigned A_MCT;				// support - antecedent_support*MCT   [maintain the original value during updating the weight]
	unsigned N;					// N = min{A_MST, A_MCT}
	unsigned supp_update;
	unsigned ante_supp_update;
	unsigned supp_original;
	double  weight;
	int		side;			// 1: removed item is on the left side;  2: on the right hand
	int		supp_in_candidate_trans;	// The support of this border rule in candidate transactions (but NOT in the DB)
										// Candidate transactions are those which contain the current sensitive rule
	set<unsigned long> related_DB_ID;
	set<unsigned long> related_DB_ID_only_ante;     //Only the antecedent of the rule is contained in the transaction

}BORDER_RULE_INFO;

typedef struct NODE_BORDER_ITEMSET{
	unsigned supp_original;
	unsigned supp_update;
	unsigned N;				// N = supp_original - db_size*MST
	double	 weight;
	int		supp_in_candidate_trans;    // The support of this border itemset in candidate transactions (but NOT in the DB)
										// Candidate transactions are those which contain the current sensitive rule
	set<unsigned long> related_DB_ID;

}BORDER_ITEMSET_INFO;

	//set< set<int> >  collection_common_itemset;	
	//set< set<int> >  collection_common_itemset_NBRS;
	

	set< pair<item,item> > collection_common_rule_IR;

	//map< pair<item,item>, A_MST_MCT > collection_common_rule_strong;
	map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong; //【NEW】

	set< pair<item,item> > collection_common_rule_IL;		

	//map< pair<item,item>, A_MST_MCT > collection_common_rule_NBRS;
	map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_NBRS;	//【NEW】



typedef struct NODE_REL
{
	//unsigned long trans_ID;
	double		fk;
	double		weight_missing_rules;
	double		weight_ghost_rules;
	double		weight_missing_itemset;

	set< pair<item,item> >  related_border_rules_PBRS;
	set< pair<item,item> >  related_border_rules_PBRS_only_ante;
	set< pair<item,item> >  related_border_rules_NBRS;
	set< item >				related_border_itemsets;

	int PBRS_num;
	int NBRS_num;
	int Border_Itemset_num;

	int related_PBRS_min_dist;
	int related_NBRS_min_dist;
	int related_Border_Itemset_min_dist;
	//set< pair<item,item> >	related_rules;
	//set< item >				related_itemsets;

}T_RELATED;

	map<unsigned long, T_RELATED>  supporting_trans_related_map_1;	// store fully supporting transactions for blocking 1s
	map<unsigned long, T_RELATED>  supporting_trans_related_map_0;	// store partially supporting transactions for blocking 0s


	vector<unsigned long>			copy_support_of_items_one2;			//后面会更新copy_support_of_items_one2[]
	vector< vector<unsigned long> > copy_temp_counter_array2;			//后面会更新copy_temp_counter_array2[][]

	unsigned long N0, N1;						// store the current values for N0 and N1
	vector<int> N0_vec, N1_vec;		// store the the number of blocking 0's and blocking 1's for each sensitive rule

	double SM;			//	Safety margin
	double A01;			//	The proportion of trans which change 0's into ?.  0 < A01 < 1

						//	A01/(1-A01) == N0/N1

//-------------------------------------------


//===========================================
double *db_fk;   //db_fk长度同db_size(main函数申请空间),实际只保存每个supp trans的fk值，方便从trans ID(对应数组下标)直接存取到fk
//===========================================


set<set<int>> s_set_update;					// 敏感项集合， 为DSS2007算法准备，此算法将其中itemset一个一个隐藏,直到为空。

set<set<int>> deleted_non_sensitive_all;	// 全局变量： 存放了所有被隐藏的non-sensitivie itemset


//本方法的借鉴意义就在于此数据结构： 
//	1) 为每个sensitive trans 保存了它所支持的所有non-sensitive frequent itemset list
//  2) 也为每个sensitive trans 保存了它所支持的所有  sensitive frequent itemset list
//当修改或删除此trans时， 需要更新当前trans所支持的sensitive itemset或non-sensitive itemset的support, 
//因为有此数据结构，就不需要重新计算寻找当前trans支持的itemset; 以空间换时间，提高运算速度

typedef struct NODE
{
	unsigned long trans_ID;

	set< pair<item,item> >  supp_all_rule_set;
	set< pair<item,item> >  supp_only_antecedent_rule_set;

	//set< pair<item,item> >  supp_strong_rule_set;
	//set< pair<item,item> >  supp_NBRS_rule_set;

	double        weight;								
	//set<set<int>> sup_itemsets;						//存放本trans所支持的所有non-sensitive itemsets
	//set<set<int>> sup_itemsets_sen;					//存放本trans所支持的所有sensitive itemsets
}T_LEN;



vector<T_LEN> supporting_trans_vec;							// 保存了所有supporting transaction的序号ID和对应长度

vector<T_LEN> supporting_trans_vec_NBRS;					// 保存了所有“部分支持antecedent且不完全支持consequent”的transaction的序号ID和对应长度

vector<unsigned long> all_supp_trans_ID;


//-------------| 计算程序运行时间 |----------------
//-------------------------------------------------
clock_t start_clock, start_clock2, end_clock;
double total_run_time, total_run_time_after_apriori;
//-------------------------------------------------



bool sort_by_fk(const T_LEN& n1, const T_LEN& n2)
{
	return (n1.weight < n2.weight);						//目前是【升序】排序：“<”导致【升序】， “>”导致【降序】
	//这里应为<
}

bool sort_by_fk_NBRS(const T_LEN& n1, const T_LEN& n2)
{
	return (n1.weight > n2.weight);						//目前是【升序】排序：“<”导致【升序】， “>”导致【降序】
	//return (n1.weight < n2.weight);
	//这里应为>
}


vector<vector<T_LEN>>  supporting_trans_vec2;		// supporting_trans_vec2[][]是二维向量，用来分别保存每个sensitive rule的supporting trans ID和对应trans length


unsigned long * support_sen_rule_array;				// 动态创建的数组， 每个元素表示对应敏感规则的出现频次（支持度）


//针对TKDE04-2.b算法, 将sensitive rules 按照size和support排序
typedef struct NODE2
{
	unsigned long sen_rules_ID;
	unsigned long itemset_size;
	unsigned long itemset_sup;						
}ITEMSET_NODE;


vector<ITEMSET_NODE> sen_itemset_sort_vec;

//qsort()中的compare()比较函数返回“负”或“正”			用法：qsort(s,100,sizeof(s[0]),cmp);
//sort()中的compare()比较函数返回1或0，即“真”或“假”。	用法：sort(a,a+5,cmp);    


bool sort_by_itemset_size_sup(const ITEMSET_NODE& n1, const ITEMSET_NODE& n2)     //“<”导致【升序】， “>”导致【降序】
{
	if(n1.itemset_size == n2.itemset_size)
	{
		return n1.itemset_sup > n2.itemset_sup;
	}
	else
		return n1.itemset_size > n2.itemset_size;

}



int * N_iteration_array;

int main(int argc, char *argv[])
{
	start_clock = clock();			// 【开始计时】
	cout<<"\n第1次计时开始！！！"<<endl;

	strcpy_s(paramfile,FILE_NAME_LENGTH, "BA_param.txt");  //【注意】：这里将 问题参数文件名固话了，其实可以是argv[1]

    read_local_parameters();

	apriori_bodon();			// 先由apriori_bodon()计算得到所有“关联规则”写入文件(并保存于bodon_ruleset)， 再从中选择“敏感规则”
								// 获取“1-item”和“2-item(矩阵)”频繁项数据结构；获取复制树main_trie_copy；获取“瘦身”DB

	sup_u = (int)(db_size*MST);

	//initial();		//【注意】： 尽管已经得到“瘦身"数据库thin_DB。参： get_thin_DB(thin_DB, datafrom);
					//			 initial()函数得到“非瘦身”“内存”数据库database, 以下情况用到它：
					//          【TKDE04】1.a,1.b,2.a,2.b需要根据"原始数据库"每行trans大小排序；也为将来修改后的数据库的输出到文件做准备


	cout<<"\n\n\nPlease specify sensitive rules in the file.  When finished press [Enter] to continue! \n ";
	getchar();

	read_sensitive_item();		// 从敏感项文件，读取rules到 vector<pair<item,item>>类型的s_set_pair 和 vector<item>类型的s_set中

	start_clock2 = clock();		// 【计时】： 从文件中 读完敏感项后， 又一次开始计时！！！
	cout<<"\n第2次计时开始[读取完敏感项后]！！！"<<endl;
	//-----------------------------------------------

	//sort_sen_itemset();			//根据sensitive itemset的size和support对s_set_pair做降序排列.注意：一定是先排序，后找supporting trans

	//功能： s_set==>s_set_update
	//注意： s_set_update是set类型, s_set是vector类型； s_set_update使用的是内码
	for(vector<set<int>>::iterator it=s_set.begin(); it!=s_set.end(); it++)
	{
		item union_item, union_item_encoding;
		union_item = *it;								 
		basket_recode3(union_item, union_item_encoding);
		s_set_update.insert(union_item_encoding);      
	}

	copy_support_of_items_one2 = copy_support_of_items_one;		//后面会更新copy_support_of_items_one2[]
	copy_temp_counter_array2 = copy_temp_counter_array;	//后面会更新copy_temp_counter_array2[][]


	//db_fk = new double[db_size];						//db_fk长度同db_size,实际只保存每个supp trans的fk值，方便从trans ID(对应数组下标)直接存取到fk
	//for(unsigned int i=0;i<db_size;i++)  db_fk[i] = 0.0;	

	//filter_sup_trans3();			// 只扫描一遍数据库； 找到所有sensitive rules的supporting trans集合,保存到all_supp_trans_ID
	//cal_fk_ratio();				// 得到到每个trans对应的fk值

	hide_rules_process();


	getchar();
	return 0;
}



//根据sensitive itemset的size和support对s_set_pair做降序排列；【算法2.b的特色之一】
void sort_sen_itemset()
{
	cout<<"\nBEGIN to sort sensitive itemset!\n";

	for(unsigned int i=0;i<s_set_pair.size();i++)				
	{
		item union_item;
		item antecedent, consequent, antecedent_encoding, consequent_encoding;

		antecedent = s_set_pair[i].first;
		consequent = s_set_pair[i].second;

		basket_recode3(antecedent, antecedent_encoding);							//转为内码
		basket_recode3(consequent, consequent_encoding);			

		set<int>::const_iterator it;
		for( it= s_set_pair[i].first.begin(); it!= s_set_pair[i].first.end(); it++)
		{
			union_item.insert(*it);
		}
		for( it= s_set_pair[i].second.begin(); it!= s_set_pair[i].second.end(); it++)
		{
			union_item.insert(*it);			//合并antecedent和consequent
		}	

		//vector<int> rule_3_attribute = bodon_ruleset[s_set_pair[i]];				// 【错误写法】，导致找不到rule；因为bodon_ruleset使用内码
		vector<int> rule_3_attribute = bodon_ruleset[pair<item,item>(antecedent_encoding, consequent_encoding)];
		
		if(rule_3_attribute.empty()) {cout<<"Error--sort_sen_itemset()! bodon_ruleset中找不到此敏感规则\n"; getchar(); exit(-1);}

		unsigned long union_sup = rule_3_attribute[2];								// 整个rule的支持度
		unsigned long antecedent_sup = rule_3_attribute[0];							// rule左边antecedent的支持度

		ITEMSET_NODE tmp;
		tmp.sen_rules_ID = i;
		tmp.itemset_size = union_item.size();
		tmp.itemset_sup = union_sup;

		sen_itemset_sort_vec.push_back(tmp);
	}

	// 利用sort()库函数（C++支持）或qsort()库函数(C支持) 排序
	// "升序"排序：函数sort_by_len()“<”代表升序， “>”代表降序
	sort(sen_itemset_sort_vec.begin(), sen_itemset_sort_vec.end(), sort_by_itemset_size_sup);	


	//验证排序
	for(unsigned int i=0;i<s_set_pair.size();i++)				
	{
		cout<<"  size: "<<sen_itemset_sort_vec[i].itemset_size;
		cout<<"  sup:  "<<sen_itemset_sort_vec[i].itemset_sup;
		cout<<"  ID:   "<<sen_itemset_sort_vec[i].sen_rules_ID<<endl;
	}

	//重新安排s_set_pair中rule的存放顺序
	vector<pair<item,item>> s_set_pair_copy;

	for(unsigned i=0;i<sen_itemset_sort_vec.size();i++)
	{
		unsigned ID = sen_itemset_sort_vec[i].sen_rules_ID;
		s_set_pair_copy.push_back( s_set_pair[ID] );
	}

	for(unsigned i=0; i<s_set_pair_copy.size();i++)
	{
		s_set_pair[i]=s_set_pair_copy[i];
	}

	//验证按照size和union_sup重新安排“存放顺序”的s_set_pair[]
	for(unsigned i=0; i<s_set_pair.size(); i++)
	{
		display_rule(s_set_pair[i]);
	}

	transfer_s_set_pair_2_s_set();

	cout<<"\nEND to sort sensitive itemset!\n";

	getchar();
}


void transfer_s_set_pair_2_s_set()
{
	s_set.clear();

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
	}	

}


//计算t包含的sensitive itemset数量; 注意形参goaltrans接收的是thin_DB中的transaction
unsigned long cal_bk(vector<int> & goaltrans,  set<set<int>> &itemset_collection)
{
	unsigned num_bk = 0;

	for(set<set<int>>::iterator it=s_set_update.begin(); it!=s_set_update.end(); it++)				
	{
			item union_item, union_item_encoding;

			/*
			union_item = *it;								//【注意】：这里用的是s_set,其保存的是每个rules的antecedent和consequent合并itemset
			basket_recode3(union_item, union_item_encoding);
			*/

			union_item_encoding = *it;

			transactionCI transCI = goaltrans.begin();	
			itemCI itemCI = union_item_encoding.begin();
			int count = 0;

			while (itemCI!=union_item_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )								//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
				if( *transCI  == *itemCI )														//尽管形式一样， 但这里两边使用的是内部码
				{
					count++;
					itemCI++;
					transCI++;
				}
				else { transCI++;}
			}
			if (count == union_item_encoding.size() && union_item_encoding.size()!=0)
			{
				num_bk++;
				itemset_collection.insert(union_item_encoding);
			}
	}

	return num_bk;

}



//获得t包含的所有itemset的数量
unsigned long cal_ak(vector<int> & goaltrans)
{
	unsigned long num_ak = 0;

	//提示： map<item, int>  bodon_itemset;
	//       typedef map<item, int>::value_type   bodon_itemset_insert_type; 

	for(map<item,int>::iterator it=bodon_itemset.begin(); it!=bodon_itemset.end(); it++)
	{
		set<int> itemset_encode;
		itemset_encode	= it->first;  //bodon_ruleset中保存的是内码itemset

		transactionCI transCI = goaltrans.begin();	
		itemCI itemCI = itemset_encode.begin();
		int count = 0;

		while (itemCI!=itemset_encode.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
		{
			//if ( copy_new_code_inverse[*transCI] == *itemCI )								//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
			if( *transCI  == *itemCI )														//尽管形式一样， 但这里两边使用的是内部码
			{
					count++;
					itemCI++;
					transCI++;
			}
			else { transCI++;}
		}
			
		if (count == itemset_encode.size() && itemset_encode.size()!=0)
		{
			num_ak++;
		}	

	}

	//漏洞： 找到t包含的所有频繁itemset后， 没有排除敏感itemset

	return num_ak;
}


//通过trie树获得t包含的所有itemset的数量
unsigned long cal_ak_by_trie(vector<int> & goaltrans,  set<set<int>> &itemset_collection )
{
	set<int> itemset;

	p_main_trie_copy->find_sup_itemset_trie_tree(goaltrans.begin(), goaltrans.end(), itemset, itemset_collection);


	set<set<int>>::iterator it=itemset_collection.begin(); 
	while( it!=itemset_collection.end()  )
	{
		if( is_belong_to_sensitive(*it) )
		{
			set<set<int>>::iterator tmp_it;
			tmp_it=itemset_collection.erase(it);

			//if(tmp_it == itemset_collection.end() ) {cout<<"cal_ak_by_trie()中删除的是最后一个item"<<endl; break;}

			it = tmp_it;

			//cout<<"排除sensitive itemset 发生"<<endl;
		}
		else
		{
			it++;
		}
	}


	//注意itemset_collection里的itemset使用的是内码
	return itemset_collection.size();
}


//计算每个snestive trans的fk ratio; fk=bk/(ak+1)
void cal_fk_ratio()
{
	cout<<"\nBEGIN to calculate fk ratio for each supporting trans!\n";

	for(unsigned i = 0; i<all_supp_trans_ID.size(); i++)
	{
		unsigned long ID = all_supp_trans_ID[i] ;

		vector<int> t =  thin_DB[ ID ];

		//unsigned long ak = cal_ak(t);			//获得t包含的所有itemset的数量 [时间复杂度高]
		set<set<int>> itemset_collection;
		unsigned long ak = cal_ak_by_trie(t, itemset_collection);	//通过trie树获得t包含的所有non-sensitive itemset的数量 [时间复杂度低]
		//supporting_trans_vec[i].sup_itemsets = itemset_collection;	//注意itemset_collection里的itemset使用的是内码
		
		itemset_collection.clear();
		unsigned long bk = cal_bk(t, itemset_collection);			//获得t包含的sensitive itemset的数量
		//supporting_trans_vec[i].sup_itemsets_sen = itemset_collection;

		double fk;

		//if(ak > 0) fk= (double)bk/ak;
		//else fk = (double)bk; 
		fk = (double)bk/(1+ak);

		db_fk[ID] = fk;

		cout<<".";

		
	}

	printf_file_fk();

	cout<<"\nEND to calculate fk ratio for each supporting trans!\n";

	display_sup_trans();
}


void printf_file_fk()
{

	ofstream outfile;
	char filename[160];

	cout<<"\n>>> BEGIN to printf_file_fk() \n";

	sprintf_s(filename,160,"%s_MST_%.3f_MCT_%.2f_App_Intell14_TID_fk_outcome.dat",datafrom, MST,MCT);
	outfile.open(filename, ios::out);

	outfile<<"Trans ID\t\tfk\n\n";	
	for(unsigned i = 0; i<all_supp_trans_ID.size(); i++)
	{
		unsigned long ID = all_supp_trans_ID[i] ;
	
		outfile<<ID<<"\t\t"<<db_fk[ID]<<endl;

	}

	outfile.close();

	cout<<"\n>>> END to printf_file_fk()\n";

}


void filter_sup_trans3()
{

	cout<<"\nBEGIN to filter out supporting transactions for all sensitive rule [ver 3]!\n";

	for(unsigned long i=0;i<thin_DB.size(); i++)				//使用“瘦身”数据库(thin_DB中每行经过了排序处理)
	{
		transaction goaltrans = thin_DB[i];					 		

		for(unsigned int j=0;j<s_set_pair.size();j++)				
		{
			item union_item, union_item_encoding;

			union_item = s_set[j];								//【注意】：这里用的是s_set,其保存的是每个rules的antecedent和consequent合并itemset
			basket_recode3(union_item, union_item_encoding);

			transactionCI transCI = goaltrans.begin();	
			itemCI itemCI = union_item_encoding.begin();
			int count = 0;

			while (itemCI!=union_item_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )								//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
				if( *transCI  == *itemCI )														//尽管形式一样， 但这里两边使用的是内部码
				{
					count++;
					itemCI++;
					transCI++;
				}
				else { transCI++;}
			}
			if (count == union_item_encoding.size() && union_item_encoding.size()!=0)
			{
				all_supp_trans_ID.push_back(i); //i为当前trans ID

				break;		//只要当前trans包含1个sensitive itemset,就不必判断它是否还包含其他sensitive itemset
			}
		}
	}

	cout<<"\nEND to filter out supporting transactions for all sensitive rule [ver 3]!\n";

}


//只扫描一遍数据库； 找到每个sensitive rules的supporting trans集合,分别保存
void filter_sup_trans2()
{
	cout<<"\nBEGIN to filter out supporting transactions for each sensitive rule!\n";

	for(unsigned int j=0;j<s_set_pair.size();j++)
	{
		vector<T_LEN> sup_vec;
		supporting_trans_vec2.push_back(sup_vec);	
	}


	for(unsigned long i=0;i<thin_DB.size(); i++)				//使用“瘦身”数据库(thin_DB中每行经过了排序处理)
	{
		transaction goaltrans = thin_DB[i];					 		

		for(unsigned int j=0;j<s_set_pair.size();j++)				
		{
			item union_item, union_item_encoding;

			union_item = s_set[j];								//【注意】：这里用的是s_set,其保存的是每个rules的antecedent和consequent合并itemset
			basket_recode3(union_item, union_item_encoding);

			transactionCI transCI = goaltrans.begin();	
			itemCI itemCI = union_item_encoding.begin();
			int count = 0;

			while (itemCI!=union_item_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )								//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
				if( *transCI  == *itemCI )														//尽管形式一样， 但这里两边使用的是内部码
				{
					count++;
					itemCI++;
					transCI++;
				}
				else { transCI++;}
			}
			if (count == union_item_encoding.size() && union_item_encoding.size()!=0)
			{
				T_LEN tmp;
				tmp.trans_ID = i;
				//tmp.trans_len = goaltrans.size();   //【注意】:根据1.a,1.b,2.a,2.b算法的本意，这里不是瘦身DB的trans长度；而是原始DB的trans长度
				//tmp.trans_len = database[i].size();   //         这里使用了原数据库中trans的长度， 而非“瘦身数据库”

				supporting_trans_vec2[j].push_back( tmp );					//supporting_trans_vec2[][]是二维向量，用来分别保存每个sensitive rule的supporting tras ID和对应trans length
			}

		}

	}

	/*
	// 改到hide_rules_process()中排序处理

	//按照trans_len大小排序
	for(unsigned int j=0;j<s_set_pair.size();j++)				
	{
		// 利用sort()库函数（C++支持）或qsort()库函数(C支持) 排序
		// "升序"排序：函数sort_by_len()“<”代表升序， “>”代表降序
		sort(supporting_trans_vec2[j].begin(), supporting_trans_vec2[j].end(), sort_by_len);		

	}

	//验证
	for(unsigned int j=0;j<s_set_pair.size();j++)				
	{

		cout<<"\n--------------------------------------------------------\n\n";

		display_rule( s_set_pair[j]);

		for(unsigned i = 0; i<supporting_trans_vec2[j].size(); i++)
		{
			cout<<"  "<<supporting_trans_vec2[j][i].trans_ID <<" ------ "<<supporting_trans_vec2[j][i].trans_len<<endl;
		}
	}
	*/

		
	cout<<"\nEND to filter out supporting transactions for each sensitive rule!\n";

}



bool is_rule_contain_x(pair<item,item> a_rule, int x )
{
		item union_item,union_item_encoding;

		set<int>::const_iterator it;
		for( it= a_rule.first.begin(); it!= a_rule.first.end(); it++)
		{
			union_item.insert(*it);
		}
		for( it= a_rule.second.begin(); it!= a_rule.second.end(); it++)
		{
			union_item.insert(*it);			//合并antecedent和consequent
		}	
		basket_recode3(union_item, union_item_encoding);

		it = union_item_encoding.begin();
		for(; it!=union_item_encoding.end(); it++ )
		{
			if( *it == x ) return true;
		}

		return false;

}



bool TNODE_point_less(const T_LEN node, const unsigned ID)
{
	return node.trans_ID < ID;
}


void remove_t(unsigned long k, unsigned long selected_ID)		//k是当前 敏感规则序号， selected_ID是要删除的trans_ID
{

	if( k<0 || k>=s_set_pair.size())  {cout<<"Error in remove_t! k越界"; getchar(); exit(-1);}
	  
	vector<T_LEN>::const_iterator it = lower_bound(	supporting_trans_vec2[k].begin(), supporting_trans_vec2[k].end(), selected_ID , TNODE_point_less);
	  
	  //注意：使用lower_bound()前提是向量有序！
	  //lower_bound()返回一个 iterator, 它指向在[first,last)标记的有序序列中可以插入value，而不会破坏容器顺序的第一个位置，
	  //而这个位置标记了一个大于等于value 的值。调用lower_bound之前必须确定序列为有序序列，否则调用出错。
	  //参考： “编程”文件夹--》STL--》lower_bound.docx

	if( it != supporting_trans_vec2[k].end() &&  (*it).trans_ID == selected_ID )
	{
		//找到了！从supporting_trans_vec2[k]集合中删除
		;
		vector<T_LEN>::const_iterator tmp_it;
		tmp_it = supporting_trans_vec2[k].erase(it);			//erase返回被删除元素的下一个元素的指针(iterator)	
	}


}



bool is_belong_to_sensitive(set<int> itemset)
{

	for(set<set<int>>::iterator it=s_set_update.begin(); it!=s_set_update.end(); it++)
	{
		if(*it == itemset) return true;
	}

	return false;
}



bool is_itemset_contain_item(set<int> itemset, int item)
{
	for(set<int>::iterator it=itemset.begin(); it!=itemset.end(); it++ )
	{
		if(item == *it)
			return true;
	}

	return false;
}



void combine_item(item &union_encoding, item &antecedent_encoding, item &consequent_encoding)
{
		set<int>::const_iterator it;
		for( it= antecedent_encoding.begin(); it!= antecedent_encoding.end(); it++)
		{
			union_encoding.insert(*it);
		}
		for( it= consequent_encoding.begin(); it!= consequent_encoding.end(); it++)
		{
			union_encoding.insert(*it);			
		}

}

bool is_belong_to(const int an_item, const set<int> & an_itemset)
{
	set<int>::iterator it2 = an_itemset.begin();
	for(; it2 != an_itemset.end(); it2++)
	{
		if( an_item == *it2 ) return true;
	}

	return false;
}


//overload
bool is_belong_to(const set<int> & items, const set<int> & an_itemset)
{
	set<int>::iterator it1 = items.begin();

	for(; it1!=items.end(); it1++)
	{
		set<int>::iterator it2 = an_itemset.begin();

		for(; it2 != an_itemset.end(); it2++)
			if( *it1 == *it2 ) return true;
	}

	return false;
}



bool vec_is_belong_to_set( vector<int> added_items, const set<int> & an_itemset)
{
	vector<int>::iterator it = added_items.begin();
	for(; it != added_items.end(); it++)
	{
		if ( an_itemset.end() != an_itemset.find(*it) ) return true;
	}
	return false;
}



double func_weight(double weight)
{
	return weight*weight;
}



void	find_R_common( set<int> &IR )
{
		
	cout<<"\n\n[BEGIN] void	find_R_common( set<int> &IR )	\n";
		
	//set< set<int> >  collection_common_itemset;		【改用全局】
	//set< pair<item,item> > collection_common_rule;	【改用全局】
	set<int>::iterator it = IR.begin();
	for(; it!= IR.end(); it++)							//Actually, there is only one item in IR. 
	{													//It is enough to blocking one item for blocing 1s to ? for a candidate transaction.
		int cur_item = *it;

		//TYPE： map<pair<item,item>,vector<int>> bodon_ruleset_strong, 使用内码
		//TYPE： map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MIN;使用内码

		map<pair<item,item>,T_RULE_UNKNOWN>::iterator  it_bodon = bodon_ruleset_all_MIN.begin();
		for(; it_bodon!=bodon_ruleset_all_MIN.end(); it_bodon++)
		{
			item union_encoding, antecedent_encoding, consequent_encoding;

			antecedent_encoding = it_bodon->first.first;
			consequent_encoding = it_bodon->first.second;
			combine_item(union_encoding, antecedent_encoding, consequent_encoding);

			unsigned	MIN_supp		= it_bodon->second.union_supp;
			double		MIN_conf		= it_bodon->second.conf;
			unsigned	MIN_antecedent	= it_bodon->second.antecedent_supp;

			unsigned	supp2 = MIN_supp;

			BORDER_RULE_INFO tmp;
			tmp.supp_original			= MIN_supp;
			tmp.supp_update				= MIN_supp;  
			tmp.ante_supp_update		= MIN_antecedent;
			tmp.supp_in_candidate_trans = 0;


			if( is_belong_to(cur_item, consequent_encoding) )								//1st screening: having common items
			{
					//collection_common_itemset.insert(union_encoding);
					collection_common_rule_IR.insert( pair<item,item>(antecedent_encoding,consequent_encoding ));


					if( MIN_supp/(double)db_size >= MST  &&  MIN_conf >= MCT )				//2nd screening: strong rules
					{
						unsigned num_above_MST = (unsigned)(MIN_supp - db_size*MST);
						unsigned num_above_MCT = (unsigned)(MIN_supp - MIN_antecedent*MCT);

						//if(num_above_MST < N1  ||  num_above_MCT < N1)
						if(num_above_MST <= N1  ||  num_above_MCT <= N1)					//3rd screening: possible after N1 modifications
						{
							//A_MST_MCT tmp;
							//BORDER_RULE_INFO tmp;
							tmp.A_MST = num_above_MST;    tmp.A_MCT = num_above_MCT;

							//map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong; //【NEW】
							if( num_above_MST <= num_above_MCT)
								tmp.N = num_above_MST;									// 'N' denotes how many removing operations at minimum  
							else														//  to cause the border rule lost
								tmp.N = num_above_MCT;

							tmp.weight = (double)(tmp.supp_original - tmp.supp_update + 1)/(tmp.N + 1);
							tmp.weight = func_weight(tmp.weight);						//adjust the weight change rate through math function

							tmp.side  = 2;
							collection_common_rule_strong[it_bodon->first] = tmp ;
						}
					}
			}
			else if( is_belong_to(cur_item, antecedent_encoding) )
			{

					collection_common_rule_IR.insert( pair<item,item>(antecedent_encoding,consequent_encoding ));

					if( MIN_supp/(double)db_size >= MST  &&  MIN_conf >= MCT )
					{
						unsigned num_above_MST = (unsigned)(MIN_supp - db_size*MST);
						unsigned num_above_MCT = (unsigned)((double)(MIN_supp - MIN_antecedent*MCT)/(1-MCT));

						//if(num_above_MST < N1  ||  num_above_MCT < N1)
						if(num_above_MST <= N1  ||  num_above_MCT <= N1)
						{
							//A_MST_MCT tmp;
							//BORDER_RULE_INFO tmp;
							tmp.A_MST = num_above_MST;    tmp.A_MCT = num_above_MCT;


							//map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong; //【NEW】
							if( num_above_MST <= num_above_MCT)
								tmp.N = num_above_MST;									// 'N' denotes how many blocking 1sto '?' operations at minimum  
							else														//  to cause the border rule lost
								tmp.N = num_above_MCT;

							tmp.weight = (double)(tmp.supp_original - tmp.supp_update  + 1)/(tmp.N + 1);
							tmp.weight = func_weight(tmp.weight);						//adjust the weight change rate through math function

							tmp.side  = 1;
							collection_common_rule_strong[it_bodon->first] = tmp ;
						}
					}			
			}//else if


		}

	}

	cout<<"\n   -- bodon_ruleset_all_MIN.size() == "<<bodon_ruleset_all_MIN.size();
	cout<<"\n   -- collection_common_rule_IR.size() == "<<collection_common_rule_IR.size();
	cout<<"\n   -- collection_common_rule_strong.size() == "<<collection_common_rule_strong.size();
		
	cout<<"\n\n[END] void	find_R_common( set<int> &IR )	\n";

}



void	find_R_common_NBRS( set<int> &IL )  //【做了重要修改！】
{
		
	cout<<"\n\n[BEGIN] void	find_R_common_NBRS( set<int> &IL )	\n";
		
	//set< set<int> >  collection_common_itemset;			【改用全局】
	//set< pair<item,item> > collection_common_rule_NBRS;	【改用全局】
	set<int>::iterator it = IL.begin();
	//for(; it!= IL.end(); it++)
	{
		int cur_item = *it;

		//类型： map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_strong_MIN; 使用内码
		//类型： map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MAX;使用内码

		map<pair<item,item>, T_RULE_UNKNOWN>::iterator  it_bodon = bodon_ruleset_all_MAX.begin();
		for(; it_bodon != bodon_ruleset_all_MAX.end(); it_bodon++)
		{						

			item union_encoding, antecedent_encoding, consequent_encoding;						//[NOTE]: for blocking 1s to ?, only one item needs to be removed
																								//for bolcking 0s to ?, all empty item 0s in the sensitive rule need to be added
			antecedent_encoding = it_bodon->first.first;					
			consequent_encoding = it_bodon->first.second;
			combine_item(union_encoding, antecedent_encoding, consequent_encoding);

			//if( is_belong_to(cur_item, union_encoding) )										//1st screening：share the common item
			if( is_belong_to( IL, union_encoding )  )	
			{
					pair<item,item> a_rule	= it_bodon->first;
					unsigned MAX_supp		= it_bodon->second.union_supp;
					double   MAX_conf		= it_bodon->second.conf;
					unsigned MAX_antecedent = it_bodon->second.antecedent_supp;

					unsigned supp2			= MAX_supp;

					collection_common_rule_IL.insert(a_rule);

					//T_RULE_UNKNOWN & a_record_min = bodon_ruleset_all_MIN[a_rule];
					//if(  ( a_record_min.conf >= MCT                       &&  MAX_supp/(double)db_size < MST)   ||
					//	 ( a_record_min.union_supp/(double)db_size >= MST &&  MAX_conf < MCT )	  )																			
				
					//if(  ( MAX_conf >= MCT                       &&  MAX_supp/(double)db_size < MST)   ||
					//	   ( MAX_supp/(double)db_size >= MST       &&  MAX_conf < MCT )	  )		

					if(  MAX_supp/(double)db_size < MST   ||   MAX_conf < MCT   )				//2nd screening																				
					{					
						double num_supp_below_MST = 0;
						double num_conf_below_MCT = 0;

						if(  MAX_supp/(double)db_size < MST   )
							num_supp_below_MST = (db_size*MST - MAX_supp);

						if(  MAX_conf < MCT    )
						{
							//if ( is_belong_to(cur_item, consequent_encoding)   )
							if ( !is_belong_to(IL, antecedent_encoding)   )
								num_conf_below_MCT =   MAX_antecedent*MCT - MAX_supp;
							else
								num_conf_below_MCT =  (MAX_antecedent*MCT - MAX_supp)/(1-MCT) ;
						}

						//if( is_belong_to(cur_item, consequent_encoding)  &&  MAX_supp/(double)db_size < MST  ) 
						//	num_conf_below_MCT = ( MAX_antecedent*MCT - MAX_supp );
						//else 
						//	num_conf_below_MCT = ( (MAX_antecedent*MCT - MAX_supp)/(1-MCT) );

						if(  num_supp_below_MST <= (double)N0  && num_conf_below_MCT <= (double)N0 )
						  																		//3rd screening：添加N0个?，有变成ghost rules的可能
						{
							//A_MST_MCT tmp;
							BORDER_RULE_INFO tmp;

							tmp.supp_original			= MAX_supp;
							tmp.supp_update				= MAX_supp;  
							tmp.ante_supp_update		= MAX_antecedent;
							tmp.supp_in_candidate_trans = 0;

							tmp.A_MST = (int)num_supp_below_MST;   tmp.A_MCT = (int)num_conf_below_MCT;	//【注意】：取整后可能为0，
							 																			//因此在cal_P_for_each_trans_NBRS（）中对应值应+1

							if( num_supp_below_MST <= num_conf_below_MCT)
								tmp.N = (unsigned)num_conf_below_MCT;								// 'N' denotes how many blocking 0s to '?' operations at minimum  
							else																	//  to cause the border rule lost
								tmp.N = (unsigned)num_supp_below_MST;

							tmp.weight = (double)(tmp.supp_update - tmp.supp_original + 1)/(tmp.N + 1);
							tmp.weight = func_weight(tmp.weight);									//adjust the weight change rate through math function


							collection_common_rule_NBRS[pair<item,item>(antecedent_encoding,consequent_encoding)] = tmp; //满足此条件方是一个Negative Border rule
						
						}
					}
			}

		}
	}
	cout<<"\n   -- bodon_ruleset_all_MAX.size() =="<<bodon_ruleset_all_MAX.size();
	cout<<"\n   -- collection_common_rule_IL.size() =="<<collection_common_rule_IL.size();
	cout<<"\n   -- collection_common_rule_NBRS.size() =="<<collection_common_rule_NBRS.size();

	cout<<"\n\n[END] void	find_R_common_NBRS( set<int> &IL )	\n";

}




//Filter out the transactions which support the current sensitive rule
//In the following steps, some supporting transactions will be selected to block 1s to ?


void filter_sup_trans( const item  & goalitem_encoding)   //过滤出敏感项的每个supporting trans的ID和len, 结果存入结构体向量
{
	cout<<"\n\n[BEGIN]: void filter_sup_trans( const item  goalitem) \n";
	int count = 0;
	//item goalitem_encoding;

	map< pair<item,item>, BORDER_RULE_INFO >::iterator	ruleCI;		// Used to traverse "collection_common_rule_strong"
	map< set<int>,BORDER_ITEMSET_INFO >::iterator		itemsetCI;	// 

	//basket_recode3(goalitem, goalitem_encoding);					//将goalitem改为使用内部码表示（trie树及瘦身数据库都使用内部码）

	//for(unsigned int i=0;i<database.size(); i++)				
	for(unsigned long i=0;i<thin_DB.size(); i++)					//使用“瘦身”数据库(thin_DB中每行经过了排序处理)
	{
		//transaction goaltrans = database[i];
		transaction goaltrans = thin_DB[i];							//使用“瘦身”数据库			

		itemCI itemCI = goalitem_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
		transactionCI transCI = goaltrans.begin();	

		while (itemCI!=goalitem_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
		{
			//if ( copy_new_code_inverse[*transCI] == *itemCI )	//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
			if( *transCI  == *itemCI )   //尽管形式一样， 但这里两边使用的是内部码
			{
				count++;
				itemCI++;
				transCI++;
			}
			else { transCI++;}
		}

		int flag = 0;

		if (count == goalitem_encoding.size() && goalitem_encoding.size()!=0)
		{
			//T_LEN tmp;
			//tmp.trans_ID = i;		//i为trans在数据库中的记录序号
			//tmp.fk = db_fk[i];
			//supporting_trans_vec.push_back( tmp );

			T_RELATED tmp2;
			tmp2.fk = 0;
			tmp2.Border_Itemset_num = 0;
			tmp2.PBRS_num = 0;
			tmp2.NBRS_num = 0;

			//map<unsigned long, T_RELATED>  supporting_trans_related_map;   MAP: i-->T_RELATED
			supporting_trans_related_map_1[i]=tmp2;

			flag = 1;
		}

		count = 0;

		//check the count of border rules in the transactions which support the current sensitive rule

		if(flag)
		{
			//map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong;		
			//map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_NBRS;

			ruleCI = collection_common_rule_strong.begin();						// For each supporting transaction, traverse all positive border rules
			for(; ruleCI != collection_common_rule_strong.end(); ruleCI++)
			{
				if( is_trans_contain_rule(goaltrans, ruleCI->first) )
					ruleCI->second.supp_in_candidate_trans += 1;
			}

			ruleCI = collection_common_rule_NBRS.begin();
			for(; ruleCI != collection_common_rule_NBRS.end(); ruleCI++)			//For each supporting transaction, traverse all negative border rules
			{
				if( is_trans_only_contain_antecedent(goaltrans, ruleCI->first) )	//only support the antecedent, not the consequent 
					ruleCI->second.supp_in_candidate_trans += 1; 
			}

			//【NOTE】:
			// ghost rules come from two ways:
			// (1)  Assume a negative border rule'support is greater than MST but its confidence is less than MCT.
			//		Block 1 to ? in the transactions which fully support the current sensitive rule, 
			//		and this will decrease the sensitive rule's support and confidence.
			//		If the current transaction fully support a negative border rule's antecedent but not the consequent,
			//		and the blocked item belongs to the antecedent part of the border rule,
			//		then the minimum support of the border rule's antecedent will decrease and its maximum confidence will increase, 
			//
			// (2)	Assume a negative border rule's support is less than MST, or its confidence is less than MCT.
			//		Block 0 to ? on transactions which partially support the sensitive rule'antecedent but not support the consequent,
			//		and this will make it fully support sensitive rule's antecedent but not the consequent. 
			//		Thus, the maximum support of the sensitive rule's antecedent and its minimum confidence will decrease. 
			//		If this blocked transaction becomes to fully support a negative border rule but it previously does not support it, 
			//		then the border rule's maximum support and maximum confidence will rise simultaneously.


			//itemsetCI = collection_border_itemset.begin();
			//for(; itemsetCI != collection_border_itemset.end(); itemsetCI++)	// For each supporting transaction, traverse all border itemsets
			//{
			//	if( is_trans_contain_itemset(goaltrans, itemsetCI->first) )
			//		itemsetCI->second.supp_in_candidate_trans += 1;
			//}
		}
	}
	cout<<"   -- supporting_trans_vec.size() == "<<supporting_trans_vec.size();
		
	cout<<"\n[END]:  void filter_sup_trans( const item  goalitem) \n";
}



//从数据库中过滤出 所有部分支持condition且不完全支持consequent的transactions

void filter_sup_trans_NBRS( set<int> & antecedent_encoding, set<int> & consequent_encoding )   //过滤出敏感项的每个supporting trans的ID和len, 结果存入结构体向量
{
		cout<<"\n\n[BEGIN]: void filter_sup_trans_NBRS( ) \n";

		map< pair<item,item>, BORDER_RULE_INFO >::iterator	ruleCI;	
	
		for(unsigned int i=0;i<thin_DB.size(); i++)										//使用“瘦身”数据库(thin_DB中每行经过了排序处理)
		{
				transaction & goaltrans = thin_DB[i];														 
		
				item::iterator CI1 = antecedent_encoding.begin();											//item为set<int>, 里面的元素自动升序排列
				vector<int>::const_iterator transCI = goaltrans.begin();	

				int count1= 0;
				int flag1 = 0;
				while ( CI1 != antecedent_encoding.end() && transCI!=goaltrans.end())	//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
				{
					//if ( copy_new_code_inverse[*transCI] == *itemCI )					//可能会出错！【注意】：copy_new_code_inverse[*transCI]译码后得到的商品编号不一定是升序排列的！
					if( abs(*transCI)  == *CI1 )										//尽管形式一样， 但这里两边使用的是内部码
					{
						count1++;
						CI1++;
						transCI++;
					}
					else { transCI++;}
				}
				if (count1 == antecedent_encoding.size() && antecedent_encoding.size()!=0)
				{
					flag1 = 1;
				}

				//检查当前trans是否匹配consequent_encoding, 标志flag2==1表示匹配
				int count2 = 0;
				int flag2 = 0;

				item::iterator  CI2 = consequent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
				transCI = goaltrans.begin();	

				while ( CI2 != consequent_encoding.end() && transCI!=goaltrans.end())					//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
				{										
					if( abs(*transCI)  == *CI2 )					 
					{
						count2++;
						CI2++;
						transCI++;
					}
					else { transCI++;}
				}
				if (count2 == consequent_encoding.size() && consequent_encoding.size()!=0)
				{
					flag2 = 1;
				}

				int flag = 0;

				if(0==flag1 && 0==flag2)															// 不完全支持antecedent, 同时也不完全支持consequent
				{
					//T_LEN tmp;
					//tmp.trans_ID = i;		//i为trans在数据库中的记录序号
					//tmp.fk = db_fk[i];
					//supporting_trans_vec_NBRS.push_back( tmp );


					T_RELATED tmp2;
					tmp2.fk = 0;
					tmp2.Border_Itemset_num = 0;
					tmp2.PBRS_num = 0;
					tmp2.NBRS_num = 0;

					//map<unsigned long, T_RELATED>  supporting_trans_related_map;   MAP: i-->T_RELATED
					supporting_trans_related_map_0[i]=tmp2;

					flag = 1;
				}

				if(flag)
				{
					ruleCI = collection_common_rule_NBRS.begin();
					for(; ruleCI != collection_common_rule_NBRS.end(); ruleCI++)			//For each supporting transaction, traverse all negative border rules
					{
						if( is_trans_adding_antecedent_contains_rule(goaltrans, antecedent_encoding, ruleCI->first) )	//only support the antecedent, not the consequent 
							ruleCI->second.supp_in_candidate_trans += 1; 
					}
				
				}

		}

		cout<<"  --supporting_trans_vec_NBRS.size()=="<<supporting_trans_vec_NBRS.size();
		
	cout<<"\n[END]:  void filter_sup_trans_NBRS(  ) \n";
}



bool is_trans_contain_rule(transaction & goaltrans, pair<item,item>   a_rule)		//要求参数goaltrans已经采用内部编码
{
	item union_encoding;
	int count = 0;

	combine_item(union_encoding, a_rule.first, a_rule.second);

	itemCI itemCI = union_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transactionCI transCI = goaltrans.begin();	

	while (itemCI!=union_encoding.end() && transCI!=goaltrans.end())			//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *transCI  == *itemCI )												//尽管形式一样， 但这里两边使用的是内部码
		{
			count++;
			itemCI++;
			transCI++;
		}
		else { transCI++;}
	}
	if (count == union_encoding.size() && union_encoding.size()!=0)
		return true;
	else
		return false;
}



int is_trans_contain_rule_L_or_R_del(transaction &goaltrans, pair<item,item> a_rule, int removed_item)		//要求参数goaltrans已经采用内部编码
{
	item antecedent_encoding = a_rule.first;

	itemCI itemCI1			= antecedent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transactionCI transCI	= goaltrans.begin();	

	int count1 = 0;
	int flag1 = 0;
	while (itemCI1!=antecedent_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *transCI  == *itemCI1 )														//尽管形式一样， 但这里两边使用的是内部码
		{
			itemCI1++;
			transCI++;
			
			count1++;
		}
		else { transCI++;}
	}
	if (count1 == antecedent_encoding.size() && antecedent_encoding.size()!=0)
		flag1 = 1;

	//----------------------------------------------------------------------------------

	item consequent_encoding = a_rule.second;

	itemCI itemCI2	= consequent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transCI			= goaltrans.begin();	

	int count2 = 0;
	int flag2 = 0;
	while (itemCI2!=consequent_encoding.end() && transCI!=goaltrans.end())																		
	{
		if( *transCI  == *itemCI2 )														
		{
			itemCI2++;
			transCI++;

			count2++;
		}
		else { transCI++;}
	}
	if (count2 == consequent_encoding.size() && consequent_encoding.size()!=0)
		flag2 = 1;

	//----------------------------------------------------------
	if( 1 == flag1 && 1 == flag2 )
	{
		if( is_belong_to(removed_item, antecedent_encoding)  )
			return 1;
		else if ( is_belong_to(removed_item, consequent_encoding)  )
			return 2;
		else
		{	cout<<"\nFunction--is_trans_contain_rule_L_or_R(): removed item is not in antecedent or consequent.";  //【说明】：被删item并非common item in IR,IR至少包含两个item
			return 0;
		}

	}
	else if(  1 == flag1 && 0 == flag2 )
	{
		if( is_belong_to(removed_item, antecedent_encoding)  )
			return 3;
		else
			return 0;
	}
	else if(  0 == flag1 && 1 == flag2 )
	{
		if( is_belong_to(removed_item, consequent_encoding)  )
			return 4;
		else
			return 0;
	}
	else
		return 0;

}



int is_trans_contain_rule_L_or_R_add(transaction &goaltrans, pair<item,item> a_rule, vector<int> &added_items)		//要求参数goaltrans已经采用内部编码
{
	item antecedent_encoding = a_rule.first;

	itemCI itemCI1 = antecedent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transactionCI transCI = goaltrans.begin();	

	int count1 = 0;
	int flag1 = 0;
	while (itemCI1!=antecedent_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *transCI  == *itemCI1 )														//尽管形式一样， 但这里两边使用的是内部码
		{
			itemCI1++;
			transCI++;
			
			count1++;
		}
		else { transCI++;}
	}
	if (count1 == antecedent_encoding.size() && antecedent_encoding.size()!=0)
		flag1 = 1;

	//----------------------------------------------------------------------------------

	item consequent_encoding = a_rule.second;

	itemCI itemCI2 = consequent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transCI = goaltrans.begin();	

	int count2 = 0;
	int flag2 = 0;
	while (itemCI2!=consequent_encoding.end() && transCI!=goaltrans.end())																		
	{
		if( *transCI  == *itemCI2 )														
		{
			itemCI2++;
			transCI++;

			count2++;
		}
		else { transCI++;}
	}
	if (count2 == consequent_encoding.size() && consequent_encoding.size()!=0)
		flag2 = 1;

	//----------------------------------------------------------
	if( 1 == flag1 && 1 == flag2 )
	{
		if( vec_is_belong_to_set(added_items, antecedent_encoding)  )
			return 1;
		else if ( vec_is_belong_to_set(added_items, consequent_encoding)  )
			return 2;
		else
		{	cout<<"\nFunction--is_trans_contain_rule_L_or_R(): added item is not in antecedent or consequent.";  //【说明】：被添加的item并非common item in IL,IL至少包含两个item
			return 0;
		}

	}
	else if(  1 == flag1 && 0 == flag2 )
	{
		if( vec_is_belong_to_set(added_items, antecedent_encoding)  )
			return 3;
		else
			return 0;
	}
	else if(  0 == flag1 && 1 == flag2 )
	{
		if( vec_is_belong_to_set(added_items, consequent_encoding)  )
			return 4;
		else
			return 0;
	}
	else
		return 0;

}



bool is_trans_only_contain_antecedent(transaction & goaltrans, pair<item,item> a_rule)		//要求参数goaltrans已经采用内部编码
{
	item antecedent_encoding = a_rule.first;

	itemCI itemCI1 = antecedent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transactionCI transCI = goaltrans.begin();	

	int count1 = 0;
	int flag1 = 0;
	while (itemCI1!=antecedent_encoding.end() && transCI!=goaltrans.end())				//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *transCI  == *itemCI1 )														//尽管形式一样， 但这里两边使用的是内部码
		{
			itemCI1++;
			transCI++;
			
			count1++;
		}
		else { transCI++;}
	}
	if (count1 == antecedent_encoding.size() && antecedent_encoding.size()!=0)
		flag1 = 1;
	else
		return false;

	//----------------------------------------------------------------------------------

	item consequent_encoding = a_rule.second;

	itemCI itemCI2 = consequent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	transCI = goaltrans.begin();	

	int count2 = 0;
	int flag2 = 0;
	while (itemCI2!=consequent_encoding.end() && transCI!=goaltrans.end())																		
	{
		if( *transCI  == *itemCI2 )														
		{
			itemCI2++;
			transCI++;

			count2++;
		}
		else { transCI++;}
	}
	if (count2 == consequent_encoding.size() && consequent_encoding.size()!=0)
		flag2 = 1;

	//----------------------------------------------------------
	if( 1 == flag1 && 0 == flag2 )
		return true;
	else
		return false;

}




bool is_trans_adding_antecedent_contains_rule(transaction  origin_trans, item  senstive_rule_antecedent, pair<item,item> a_rule )
{

	set<int> combine;

	set<int>::iterator it_itemset = senstive_rule_antecedent.begin();
	while( it_itemset != senstive_rule_antecedent.end() )
	{
		combine.insert(*it_itemset);
		it_itemset++;
	}

	transactionCI it_trans	= origin_trans.begin();
	while( it_trans != origin_trans.end() )
	{
		combine.insert(*it_trans);
		it_trans++;
	}

	//----------------------------------------------
	item  antecedent_encoding		= a_rule.first;

	itemCI				itemCI1		= antecedent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	set<int>::iterator	combineCI	= combine.begin();	

	int	count1	= 0;
	int	flag1	= 0;

	while (itemCI1!=antecedent_encoding.end() && combineCI!=combine.end())					//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *combineCI  == *itemCI1 )														//尽管形式一样， 但这里两边使用的是内部码
		{
			itemCI1++;
			combineCI++;
			count1++;
		}
		else { combineCI++;}
	}
	if (count1 == antecedent_encoding.size() && antecedent_encoding.size()!=0)
		flag1 = 1;
	else
		return false;

	//--------------------------------------------------

	item consequent_encoding = a_rule.second;

	itemCI itemCI2	= consequent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	combineCI		= combine.begin();	

	int count2 = 0;
	int flag2 = 0;
	while (itemCI2!=consequent_encoding.end() && combineCI!=combine.end())																		
	{
		if( *combineCI  == *itemCI2 )														
		{
			itemCI2++;
			combineCI++;
			count2++;
		}
		else { combineCI++;}
	}
	if (count2 == consequent_encoding.size() && consequent_encoding.size()!=0)
		flag2 = 1;

	//----------------------------------------------------------
	if( 1 == flag1 && 1 == flag2 )
		return true;
	else
		return false;


}




bool is_trans_contain_rule_NBRS( transaction &goaltrans, pair<item,item> a_rule,  set<int> antecedent_encoding)
{

	item union_encoding;
	combine_item(union_encoding, a_rule.first, a_rule.second);

	set<int> goaltrans_set;
	for(transactionCI it = goaltrans.begin(); it!=goaltrans.end(); it++)
		goaltrans_set.insert(*it);

	for(set<int>::iterator it2 = antecedent_encoding.begin(); it2!=antecedent_encoding.end(); it2++)
		goaltrans_set.insert(*it2);

	itemCI itemCI = union_encoding.begin();										 
	set<int>::iterator transCI_set = goaltrans_set.begin();	

	int count = 0;
	while (itemCI!=union_encoding.end() && transCI_set!=goaltrans_set.end())			//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *transCI_set  == *itemCI )													//尽管形式一样， 但这里两边使用的是内部码
		{
			count++;
			itemCI++;
			transCI_set++;
		}
		else { transCI_set++;}
	}
	if (count == union_encoding.size() && union_encoding.size()!=0)
		return true;
	else
		return false;

}


bool is_trans_only_contain_antecedent_NBRS(transaction goaltrans, pair<item,item> a_rule,  set<int> add_antecedent_encoding)		//要求参数goaltrans已经采用内部编码
{

	set<int> goaltrans_set;
	for(transactionCI it = goaltrans.begin(); it!=goaltrans.end(); it++)
		goaltrans_set.insert(*it);

	for(set<int>::iterator it2 = add_antecedent_encoding.begin(); it2!=add_antecedent_encoding.end(); it2++)
		goaltrans_set.insert(*it2);

	//------------------------------------------------
	item antecedent_encoding = a_rule.first;

	itemCI it1 = antecedent_encoding.begin();										//item为set<int>, 里面的元素自动升序排列
	set<int>::iterator trans_set_CI = goaltrans_set.begin();	

	int count1 = 0;
	int flag1 = 0;
	while ( it1!=antecedent_encoding.end() && trans_set_CI!=goaltrans_set.end())			//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *trans_set_CI  == *it1 )												//尽管形式一样， 但这里两边使用的是内部码
		{
			count1++;
			it1++;
			trans_set_CI++;
		}
		else { trans_set_CI++;}
	}
	if (count1 == antecedent_encoding.size() && antecedent_encoding.size()!=0)
		flag1 =1;
	else 
		return false;

	//-------------------------------------------------
	item consequent_encoding = a_rule.second;

	itemCI it2 = consequent_encoding.begin();
	trans_set_CI = goaltrans_set.begin();

	int count2 = 0;
	int flag2 = 0;
	while ( it2!=consequent_encoding.end() && trans_set_CI!=goaltrans_set.end())			//这样比对的前提是：goalitem 和 goaltrans已经自动排好序															
	{
		if( *trans_set_CI  == *it2 )														//尽管形式一样， 但这里两边使用的是内部码
		{
			it2++;
			trans_set_CI++;

			count2++;
		}
		else { trans_set_CI++;}
	}
	if (count2 == consequent_encoding.size() && consequent_encoding.size()!=0)
		flag2 =1;
	//--------------------------------------------------

	if( 1==flag1 && 0==flag2  ) 
		return true;
	else 
		return false;

}




void cal_P_for_each_trans()							//每个trans的权重P对应fk
{
	//提示： vector<T_LEN> supporting_trans_vec;	
	cout<<"\n\n[BEGIN] void cal_P_for_each_trans()	\n";

	cout<<"  supporting_trans_vec.size() = "<<supporting_trans_vec.size()<<endl;
	cout<<"  collection_common_rule_strong.size() = "<<collection_common_rule_strong.size()<<endl;
	cout<<"  please wait...\n";


	//map< pair<item,item>, A_MST_MCT >  supp_rule_set;
	map< pair<item,item>, BORDER_RULE_INFO >  supp_rule_set;			//【NEW】


	vector<T_LEN>::iterator it = supporting_trans_vec.begin();
	for(; it != supporting_trans_vec.end(); it++)
	{
		double A_MCT_max = 0, A_MCT_min = (double)db_size;
		double trans_weight = 0;

		unsigned long DB_ID = it->trans_ID;
		transaction &goaltrans = thin_DB[DB_ID];

		//提示：set< pair<item,item> > collection_common_rule;						//collection_common_rule保存与当前sen_rule存在common item的其他rules
		//set< pair<item,item> >::iterator it2 = collection_common_rule_strong.begin();
		//map< pair<item,item>, A_MST_MCT >::iterator it2 = collection_common_rule_strong.begin();

		map< pair<item,item>, BORDER_RULE_INFO >::iterator it2 = collection_common_rule_strong.begin();		//【NEW】
		

		for(; it2 != collection_common_rule_strong.end(); it2++)
		{
			if( is_trans_contain_rule(goaltrans, it2->first) )
			{
				//it->supp_all_rule_set.insert(*it2);									//store all rules（Strong or Not Strong） which might be influenced by modification.
				//if( bodon_ruleset_all_MIN[*it2].union_supp/(double)db_size >= MST && bodon_ruleset_all_MIN[*it2].conf >= MCT)
						supp_rule_set.insert(*it2);										//"supp_rule_set" only store related strong rules
			}
			//else if( is_trans_only_contain_antecedent(goaltrans, *it2)  )
			//	it->supp_only_antecedent_rule_set.insert(*it2);

		}

		//cout<<" supp_rule_set.size() == "<<supp_rule_set.size()<<endl;
		//it->supp_strong_rule_set = supp_rule_set;			//保存当前trans完全支持的strong rules

		for(it2 = supp_rule_set.begin(); it2!=supp_rule_set.end(); it2++)
		{
			unsigned num_above_MST = it2->second.A_MST;
			unsigned num_above_MCT = it2->second.A_MCT;
																					//加1是为了防止A_MCT_R == 0;
			double A_MCT_R;
			if( num_above_MST<num_above_MCT )
			     A_MCT_R = (double)num_above_MST + 1;	
			else
				 A_MCT_R = (double)num_above_MCT + 1;	

			if (A_MCT_R > A_MCT_max)	A_MCT_max = A_MCT_R;
			if (A_MCT_R < A_MCT_min)	A_MCT_min = A_MCT_R;

			trans_weight += 1.0/A_MCT_R; 
		}

		//必须考虑supp_rule_set,size()为0的情况
		if(supp_rule_set.size()==0)									//说明当前trans不包含任何满足3层过滤的rule, 这种情况副作用最小！
		{
			//cout<<"\n\nERROR! 至少当前要隐藏的sen rule应在supp_rule_set中"; getchar(); exit(-1);
			trans_weight = 0;
		}

		//trans_weight = trans_weight* A_MCT_max/A_MCT_min;

		it->weight = trans_weight;				// trans_weight为当前trans的权重， 保存于fk, fk这个名称从其他程序(DSS07)移植过来，实际表示的是P
	
		//it->weight = (double)(supp_rule_set.size());

		supp_rule_set.clear();

		//cout<<".";
	}

	cout<<"\n[END] void cal_P_for_each_trans()	\n";
}




void cal_P_for_each_trans_NBRS( set<int> & antecedent_encoding)										//为每个NBRS的supporting trans计算权重P，  P对应fk
{
	//提示： vector<T_LEN> supporting_trans_vec_NBRS;	

	cout<<"\n\n[BEGIN]:  void cal_P_for_each_trans_NBRS()	\n";


	//map< pair<item,item>, A_MST_MCT >  supp_rule_set_NBRS;	
	map< pair<item,item>, BORDER_RULE_INFO >  supp_rule_set_NBRS;


	cout<<"  supporting_trans_vec_NBRS.size() = "<<supporting_trans_vec_NBRS.size()<<endl;
	cout<<"  collection_common_rule_NBRS.size() = "<<collection_common_rule_NBRS.size()<<endl;
	//int i=0;
	cout<<"  please wait...\n";

	vector<T_LEN>::iterator it = supporting_trans_vec_NBRS.begin();

	for(; it != supporting_trans_vec_NBRS.end(); it++)
	{
		double B_MCT_max = 0, B_MCT_min = (double)db_size;
		double trans_weight = 0;

		unsigned long DB_ID = it->trans_ID;
		transaction & goaltrans = thin_DB[DB_ID];

		//提示：set< pair<item,item> > collection_common_rule_NBRS;							//【以前】collection_common_rule_NBRS保存与当前sen_rule的antecedent存在common item的其他Negative Boder Rules Set
		//set< pair<item,item> >::iterator it2 = collection_common_rule_NBRS.begin();			//【现在】collection_common_rule_NBRS保存与当前sen_rule的antecedent存在common item的【所有rules】
		//map< pair<item,item>, A_MST_MCT >::iterator it2 = collection_common_rule_NBRS.begin();
		map< pair<item,item>, BORDER_RULE_INFO >::iterator it2 = collection_common_rule_NBRS.begin();

		
		for(; it2 != collection_common_rule_NBRS.end(); it2++)
		{
			//if(is_trans_contain_rule_NBRS(goaltrans, *it2, antecedent_encoding))			//【注意】: is_trans_contain_rule_NBRS()函数需要先补齐IL，再判断是否包含rule
			if(is_trans_contain_rule_NBRS(goaltrans, it2->first, antecedent_encoding))	
			{
				//it->supp_all_rule_set.insert(*it2);
				//if(  ( bodon_ruleset_all_MIN[*it2].conf > MCT  &&  bodon_ruleset_all_MAX[*it2].union_supp/(double)db_size < MST )  ||
				//	 ( bodon_ruleset_all_MAX[*it2].conf < MCT  &&  bodon_ruleset_all_MIN[*it2].union_supp/(double)db_size > MST )  
				//  )
						supp_rule_set_NBRS.insert(*it2);											//supp_rule_set保存了 从collection_common_rule选出的本trans支持的rules
			}
			//else if( is_trans_only_contain_antecedent_NBRS(goaltrans, *it2, antecedent_encoding)  )
			//	it->supp_only_antecedent_rule_set.insert(*it2);
		}

		//cout<<" it->supp_all_rule_set.size()== "<<it->supp_all_rule_set.size();
		//cout<<" supp_rule_set_NBRS.size() == "<<supp_rule_set_NBRS.size()<<endl;
		//it->supp_NBRS_rule_set = supp_rule_set_NBRS;

		it2 = supp_rule_set_NBRS.begin();
		for(; it2!=supp_rule_set_NBRS.end(); it2++)
		{
			//unsigned int MAX_union_sup, MIN_union_sup, MAX_antecedent_sup;
			//double MAX_conf, MIN_conf;

			//T_RULE_UNKNOWN & a_record = bodon_ruleset_all_MAX[*it2];

			//MAX_union_sup = a_record.union_supp;
			//MAX_antecedent_sup = a_record.antecedent_supp;
			//MAX_conf = a_record.conf;

			//T_RULE_UNKNOWN & a_record2 = bodon_ruleset_all_MIN[*it2];
			//MIN_union_sup = a_record2.union_supp;
			//MIN_conf = a_record2.conf;

			//double B_MCT_R;
			//if( (double)MIN_union_sup/db_size >= MST  && MAX_conf<MCT )		
			//{
			//	B_MCT_R =  (MCT - MAX_conf)*MAX_antecedent_sup  + 1;				//conf小于MCT    加1是为了防止A_MCT_R == 0;
			//}
			//else if( MIN_conf >= MCT && (double)MAX_union_sup/db_size < MST  )
			//{
			//	B_MCT_R =  MST * db_size - (double)MAX_union_sup  + 1;				//supp小于MST
			//}
			//else
			//{
			//	cout<<"\n Error in cal_P_for_each_trans_NBRS():  NBRS rule的support和conf都大于阈值！";	getchar(); exit(-1);
			//}

			double B_MCT_R;
			unsigned num_below_MST = it2->second.A_MST;
			unsigned num_below_MCT = it2->second.A_MCT;

			if( num_below_MST >  num_below_MCT )
				B_MCT_R = (double)(  num_below_MST ) + 1 ;				//加1是为了防止A_MCT_R == 0;
			else
				B_MCT_R = (double)(  num_below_MCT ) + 1;

			if (B_MCT_R > B_MCT_max)	B_MCT_max = B_MCT_R;
			if (B_MCT_R < B_MCT_min)	B_MCT_min = B_MCT_R;

			trans_weight += 1.0/B_MCT_R; 
		}

		if(supp_rule_set_NBRS.size() == 0)				//应该考虑supp_rule_set_NBRS.size()==0的情况，即当前trans插入item后不包含任何Negative-border rule的情况
			trans_weight = 0;
		else
		{
			//trans_weight = trans_weight* B_MCT_max/B_MCT_min;
			;
		}

		it->weight = trans_weight;				// trans_weight为当前trans的权重， 保存于fk, fk这个名称从其他程序(DSS07)移植过来，实际表示的是P

		//it->weight = (double)(supp_rule_set_NBRS.size());

		supp_rule_set_NBRS.clear();

		//cout<<".";

		//cout<<"."<<" "<<i;
		//i++;
	}

	cout<<"\n[END]:  void cal_P_for_each_trans_NBRS()	\n";
}




void modify_trans_thin_DB_by_adding( vector<int> & trans_vec, item antecedent, vector<int> & added_items_vec, vector<int> & t_copy )
{
	set<int> trans_set;										//【警惕】:插入新?时，务必保证新插入的负值应放在 对应正值按排序应存放的位置，这样才保证与item匹配时不出错！！

	added_items_vec.clear();

	vector<int>::iterator it_vec = trans_vec.begin();
	for(; it_vec!=trans_vec.end() ; it_vec++)				//将记录trans从vector<int>转入set<int>
	{
		trans_set.insert(*it_vec);						
	}

	item::iterator l_it = antecedent.begin();				
	for(; l_it != antecedent.end(); l_it++ )
	{
		if( trans_set.end() == trans_set.find(*l_it) )			//find()函数返回值== end()表示 在当前trans中找不到*l_it
		{
			added_items_vec.push_back(*l_it);					//找出待添加的items, 放入added_items_set
			trans_set.insert(*l_it);	//实际更新trans			//【版本1】【不实际更新trans】【严格要求sensitive rules之间无comon items】
																//【采用】：【版本2】【实际更新trans】
		}
	}

	sort(added_items_vec.begin(), added_items_vec.end());

	trans_vec.clear();
	set<int>::iterator it_set = trans_set.begin();				//将新的trans由set<int>导回vector<int>
	for(; it_set!= trans_set.end(); it_set++)
	{
		trans_vec.push_back(*it_set);
	}

	t_copy = trans_vec;

	vector<int>::iterator it_add = added_items_vec.begin();
	for( ; it_add != added_items_vec.end(); it_add++ )
	{
		for(vector<int>::iterator it_trans = trans_vec.begin(); it_trans != trans_vec.end(); it_trans++   )
		{
			if( *it_trans == *it_add)
			{
				*it_trans = -(*it_add);
				break;
			}
		}
	}

}



void compute_N0_N1(item & antecedent_encoding, item & consequent_encoding)
{
	cout<<"\n\nBEGIN to compute N0 and N1\n";
	
	pair<item,item> cur_sen_rule = pair<item,item>(antecedent_encoding, consequent_encoding);
	
	unsigned  union_sup = bodon_ruleset_all_MIN[ cur_sen_rule ].union_supp;	// 整个rule的支持度
	unsigned  antecedent_sup = bodon_ruleset_all_MIN[ cur_sen_rule ].antecedent_supp; // rule左边antecedent的支持度	

	assert(A01 != 0);

	double N0_double =  ( (double)union_sup - (MCT - MCT*SM)*antecedent_sup ) / ( (MCT - MCT*SM) + (1.0/A01 -1) )  ;
	N0 =  (unsigned)( N0_double ) + 1;

	N0_vec.push_back(N0);

	double N1_double =  (1-A01)/A01 * (N0_double);
	N1 =  (unsigned)( N1_double ) + 1;

	N1_vec.push_back(N1);


	printf("\n N0 = %u     N1 = %u", N0, N1);

	cout<<"\n\nEND to compute N0 and N1\n";
}



//////////////////////////////////////////////////////////
void hide_rules_process(  )
{
	int alpha, beta, gamma;

	unsigned long modified_item_num=0;

	set<int> set_modified_trans_ID;						//【新增】： 统计被修改trans数量

	alpha = beta = gamma = 0;

	N_iteration_array = new int[s_set_pair.size()];

	for(unsigned int i=0;i<s_set_pair.size();i++)				
	{
		item union_item, union_encoding;
		item antecedent, consequent, antecedent_encoding, consequent_encoding;

		set<int>::const_iterator it;
		for( it= s_set_pair[i].first.begin(); it!= s_set_pair[i].first.end(); it++)
		{
			union_item.insert(*it);
			antecedent.insert(*it);
		}
		for( it= s_set_pair[i].second.begin(); it!= s_set_pair[i].second.end(); it++)
		{
			union_item.insert(*it);			
			consequent.insert(*it);
		}

		cout<<"\n\n\n------------------------------------------------------------\n";
		cout<<"\nTo hide the "<<i+1<<" sensitive rule"<<endl;
		display_rule2( pair<item,item>(antecedent, consequent ) );

		basket_recode3(antecedent, antecedent_encoding);
		basket_recode3(consequent, consequent_encoding);				// 将“条件”和“结论”转化为内码表示，因为bodon_ruleset使用内码
		basket_recode3(union_item, union_encoding);

		compute_N0_N1( antecedent_encoding, consequent_encoding  );


		//----------------------------------------------------------------------------------------------


		int selected_item ;

		item::const_iterator it2 = consequent_encoding.begin();			//注意：使用“内码”版本的consequent
		selected_item = *it2;

		for(it2++; it2!= consequent_encoding.end(); it2++)				//选支持度最高的item； 此方法与当前trans无关
		{
			if( copy_support_of_items_one[*it2] > copy_support_of_items_one[selected_item] )	
				selected_item = *it2;   
		}
		set<int> temp_IR;
		temp_IR.insert(selected_item);
		//-------------- 从consequent中选择一个item删除【注意：此法选择的待删item不与当前trans相关】----

		//find_R_common(antecedent_encoding);	【这里犯了严重错误！！！！！】【应根据consequent找包含common item的rules】
		find_R_common(temp_IR);
		find_R_common_NBRS( antecedent_encoding );


		filter_sup_trans( union_encoding );													// 扫描数据库(scan DB)找出当前sen_rule的supporting trans
		filter_sup_trans_NBRS(antecedent_encoding, consequent_encoding );


		cal_P_for_each_trans();															//为每个supp_trans计算P(权重)值，对应fk
		cal_P_for_each_trans_NBRS(antecedent_encoding);									//为每个NBRS_supp_trans计算P(权重)值，对应fk


		collection_common_rule_strong.clear();											//【清空】collection_common_rule保存与当前sen_rule存在common item的其他rules
		collection_common_rule_NBRS.clear();


		//利用sort()库函数（C++支持）或qsort()库函数(C支持) 排序; 按权重对supp. trans排序
		//sort函数头文件#include <algorithm>    ---  qsort函数头文件#include<stdlib.h>

		cout<<"\n\n[BEGIN]: sort!";
		sort(supporting_trans_vec.begin(), supporting_trans_vec.end(), sort_by_fk);		//【升序】排序：函数sort_by_len()“<”代表升序， “>”代表降序												
		sort(supporting_trans_vec_NBRS.begin(), supporting_trans_vec_NBRS.end(), sort_by_fk_NBRS);
		cout<<"\n\n[END]: sort!\n";

		//display_sup_trans();															//屏幕显示，验证排序结果
		//display_sup_trans_NBRS();


		//[NEW!!!]
		//[NEW!!!] make changes from sorting to finding the minimum in each iteration of database modification
		//[NEW!!!] weights of border rules and priority of related transactions will be DYNAMICALLY UPDATED
		//-----------------------------------------------------------------------------------------------------
				
		//【Blocking N1 1'】


		cout<<"\n\n 【Blocking N1 1'】  N1 == "<<N1<<endl;
		cout<<"\n   -- collection_common_rule_IR.size() == "<<collection_common_rule_IR.size()<<endl;
		for(unsigned j=0; j<N1; j++)
		{
				//-----------------------------------------------------------------------------
				//在内存数据库databas中真实删除selected_item对应项
				//int flag = 0;
				//vector<int>::const_iterator it3 = database[ supporting_trans_vec[j].trans_ID ].begin();
				//while( it3 != database[ supporting_trans_vec[j].trans_ID ].end() )
				//{
				//	if( copy_new_code[*it3] -1  == selected_item  )    // 将*it3转化为内码
				//	{
				//		vector<int>::const_iterator tmp_it;
				//		tmp_it = database[ supporting_trans_vec[j].trans_ID ].erase(it3);  //erase返回被删除元素的下一个元素的指针(iterator)
				//		it3 = tmp_it;
				//		flag = 1;
				//	}
				//	else
				//		it3++;
				//}
				//------------------------------------------------------------------------------------

				unsigned long DB_ID = supporting_trans_vec[j].trans_ID;
				vector<int> &t = thin_DB[ DB_ID ];	 
				vector<int> t_copy = thin_DB[ DB_ID ];	

				int flag = 0;
				vector<int>::iterator it3 = t.begin();
				while( it3 != t.end() )
				{
					if( *it3  == selected_item  )    
					{
						*it3 = -selected_item;					//【注意】：这里将item符号取反，变成对称的负数，表示此item在当前trans中为Unknown

						flag = 1;
						break;
					}
					else
						it3++;
				}
				//--------------------------------------------------------------------------------
			
				if(flag)	//加了这个标志判断可防止重复删除同一trans的同一item, 防止下面的重复更新
				{
					modified_item_num++;
					set_modified_trans_ID.insert(DB_ID);				//【新增】： 统计被修改trans数量
					//--------------------------------------------

					cout<<".";
					//cout<<"\nj=="<<j<<"  collection_common_rule_IR.size() == "<<collection_common_rule_IR.size()<<endl;

					set< pair<item,item> >::iterator it4 = collection_common_rule_IR.begin();
					for(; it4 != collection_common_rule_IR.end(); it4++)
					{
						//类型： map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MIN;
						//T_RULE_UNKNOWN a_record;

						pair<item,item> related_rule = *it4;
						T_RULE_UNKNOWN & a_record = bodon_ruleset_all_MIN[*it4];

						int result = is_trans_contain_rule_L_or_R_del(t_copy, *it4, selected_item);

						if(  1 == result  )
						{
							a_record.antecedent_supp -= 1;
							a_record.union_supp -= 1;
							a_record.conf = (double)(a_record.union_supp)/a_record.antecedent_supp; 

							//cout<<"\n之后:\n";
							//cout<<" bodon_ruleset_all_MIN[*it4].antecedent_supp = "<<bodon_ruleset_all_MIN[*it4].antecedent_supp;
							//cout<<" bodon_ruleset_all_MIN[*it4].union_supp = "<<bodon_ruleset_all_MIN[*it4].union_supp;
							//cout<<" bodon_ruleset_all_MIN[*it4].conf = "<<bodon_ruleset_all_MIN[*it4].conf ;
						}
						else if( 2 == result )
						{	
							a_record.consequent_supp -= 1;
							a_record.union_supp -= 1;
							a_record.conf = (double)(a_record.union_supp)/a_record.antecedent_supp; 
						}
						else if( 3 == result )
						{
							T_RULE_UNKNOWN & a_record_MAX = bodon_ruleset_all_MAX[*it4];
							a_record_MAX.antecedent_supp -= 1;
							a_record_MAX.conf = (double)(a_record_MAX.union_supp)/a_record_MAX.antecedent_supp;
						}
					}
				}
				else
				{
					cout<<"发生重复删除item!  Sensitive rules之间有common items!"<<endl;		//getchar();
				}
		}

		//------------------------------------------------------------------------------------------
		//【Block N0 0'】

		cout<<"\n\n 【Block N0 0'】  N0 == "<<N0<<endl;
		cout<<"\n   -- collection_common_rule_IL.size() == "<<collection_common_rule_IL.size()<<endl;

		unsigned N0_iteration;

		if( supporting_trans_vec_NBRS.size() < N0 )  //如果可用的partial support trans不够， 
		{
			cout<<"   \n【Warning】:The partially supporting trans is not sufficient!!\n";
			cout<<"    N0(required num) = "<<N0<<"    supporting_trans_vec_NBRS.size() = "<<supporting_trans_vec_NBRS.size()<<endl;

			N0_iteration = supporting_trans_vec_NBRS.size();
		}
		else
			N0_iteration = N0;


		for(unsigned j=0; j<N0_iteration; j++)
		{
				unsigned long DB_ID = supporting_trans_vec_NBRS[j].trans_ID;
				vector<int> &t = thin_DB[ DB_ID ];	 

				//vector<int> t_copy = thin_DB[ DB_ID ];	
				vector<int> t_copy;

				vector<int> added_items_vec_encode;

				modify_trans_thin_DB_by_adding(thin_DB[ DB_ID ], antecedent_encoding, added_items_vec_encode, t_copy);	//update database

				if(added_items_vec_encode.empty()) cout<<"!!!!!!!!!!!!!!!!!!!!!!!\n";

				if( !added_items_vec_encode.empty())										// 可防止重复更新
				{

					set_modified_trans_ID.insert(DB_ID);

					modified_item_num = modified_item_num + added_items_vec_encode.size();						//update data side effects
					//-------------------------------------------------------------------------
					cout<<".";
					//cout<<"\nj=="<<j<<"  collection_common_rule_IL.size() == "<<collection_common_rule_IL.size()<<endl;

					set<pair<item,item>>::iterator it5 = collection_common_rule_IL.begin();						//update knowledge side effects
					for( ; it5 != collection_common_rule_IL.end(); it5++  )
					{
						//类型： map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MAX;
						//T_RULE_UNKNOWN a_record;

						pair<item,item> related_rule = *it5;
						T_RULE_UNKNOWN & a_record = bodon_ruleset_all_MAX[*it5];

						int result = is_trans_contain_rule_L_or_R_add(t_copy, *it5, added_items_vec_encode) ;
						if(  1 == result )
						{
							a_record.antecedent_supp += 1;
							a_record.union_supp += 1;
							a_record.conf = (double)(a_record.union_supp)/a_record.antecedent_supp; 
						}
						else if( 2 == result )
						{	
							a_record.consequent_supp += 1;
							a_record.union_supp += 1;
							a_record.conf = (double)(a_record.union_supp)/a_record.antecedent_supp; 
						}
						else if( 3 == result )
						{
							T_RULE_UNKNOWN & a_record_MIN = bodon_ruleset_all_MIN[*it5];
							a_record_MIN.antecedent_supp += 1;
							a_record_MIN.conf = (double)(a_record_MIN.union_supp)/a_record_MIN.antecedent_supp;
						}	
					}

				}

				t_copy.clear();	//It can be omitted because t_copt is a local variable
				//-------------------------------------------------------------
		}

		collection_common_rule_IR.clear();
		collection_common_rule_IL.clear();

		supporting_trans_vec.clear();							// 【清空】supporting_trans_vec是全局结构体向量
		vector<T_LEN>().swap(supporting_trans_vec);				// 这样才会真正释放物理空间

		supporting_trans_vec_NBRS.clear();
		vector<T_LEN>().swap(supporting_trans_vec_NBRS);		//Free the memory

	}

	//----------------------------------------------------------------------------------------
	
	 write_DB_to_file();										//output transactions from the memory to files

	//---------------------- 统计alpha--------------------------------------------------------
	//类型：vector<pair<item,item>> s_set_pair;	
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		item antecedent_encoding, consequent_encoding;

		basket_recode3(s_set_pair[j].first, antecedent_encoding);
		basket_recode3(s_set_pair[j].second, consequent_encoding);				// 将“条件”和“结论”转化为内码表示，因为bodon_ruleset使用内码

		pair<set<int>, set<int>> a_rule = pair<set<int>, set<int>>(antecedent_encoding, consequent_encoding );

		unsigned sen_new_union_sup = bodon_ruleset_all_MIN[a_rule].union_supp;
		double sen_new_rule_confidence = bodon_ruleset_all_MIN[a_rule].conf;

		if(sen_new_union_sup/(double)db_size >= MST && sen_new_rule_confidence >=MCT   )
		{	
			alpha++;	cout<<"\nalpha is incereasing for sensitive rules not hidden!   alpha="<<alpha;
			
			display_rule(s_set_pair[j]);

			cout<<"union_sup_ratio = "<<sen_new_union_sup/(double)db_size<<" confidence = "<<sen_new_rule_confidence<<endl;
			//getchar();
		}

	}


	/*-------------------- 统计beta ----------------------*/

	map<pair<item,item>, T_RULE_UNKNOWN>::iterator it = bodon_ruleset_strong.begin();
	for(; it != bodon_ruleset_strong.end(); it++)
	{
		unsigned min_union_sup = bodon_ruleset_all_MIN[it->first].union_supp;
		double min_conf = bodon_ruleset_all_MIN[it->first].conf;

		if( min_union_sup/(double)db_size < MST || min_conf < MCT )
			beta ++;
	}

	beta = beta -  (s_set_pair.size()-alpha); 


	/*-------------------- 统计gamma ----------------------*/

	it = bodon_ruleset_NBRS.begin();
	for(; it != bodon_ruleset_NBRS.end(); it++)								//Rules in "bodon_ruleset_NBRS" are originally below MST or MCT
	{
		unsigned max_union_sup = bodon_ruleset_all_MAX[it->first].union_supp;
		double max_conf = bodon_ruleset_all_MAX[it->first].conf;

		if( max_union_sup/(double)db_size >= MST && max_conf >= MCT )
			gamma ++;
	}


	cout<<endl;
	cout<<" alpha = "<<alpha<<"    beta = "<<beta<<"    gamma = "<<gamma<<"   #of modified items = "<<modified_item_num<<endl;


	write_performance(alpha, beta, gamma, set_modified_trans_ID.size(), modified_item_num);

}




void display_sup_trans() 
{
	cout<<"\n\n\n";
	for(unsigned i = 0; i<supporting_trans_vec.size(); i++)
	{
		cout<<"  "<<supporting_trans_vec[i].trans_ID <<" ------   weight = "<<supporting_trans_vec[i].weight<<endl;
		//cout<<"  bk = "<<supporting_trans_vec[i].bk<<"  ak = "<<supporting_trans_vec[i].ak<<endl;
	}

}

void display_sup_trans_NBRS() 
{
	cout<<"\n\n\n";
	for(unsigned i = 0; i<supporting_trans_vec_NBRS.size(); i++)
	{
		cout<<"  "<<supporting_trans_vec_NBRS[i].trans_ID <<" ------   weight(NBRS) = "<<supporting_trans_vec_NBRS[i].weight<<endl;
	}

}

void display_rule(pair<item,item> a_rule)
{
	cout<<"\n ( ";

	item::const_iterator it = a_rule.first.begin();
	for(; it != a_rule.first.end()  ; it++)
	{
		cout<<*it<<" ";
	}
	cout<<") -> ( ";
	for(it = a_rule.second.begin(); it!=a_rule.second.end(); it++)
	{
		cout<<*it<<" ";
	}
	cout<<")"<<endl;


	//item union_item;
	//item antecedent, consequent, antecedent_encoding, consequent_encoding;

	//antecedent = a_rule.first;
	//consequent = a_rule.second;

	//display_rule( pair<item,item>(antecedent,consequent) );

	//basket_recode3(antecedent, antecedent_encoding);							//转为内码
	//basket_recode3(consequent, consequent_encoding);			

	//for( it= a_rule.first.begin(); it!= a_rule.first.end(); it++)
	//{
	//		union_item.insert(*it);
	//}
	//for( it= a_rule.second.begin(); it!= a_rule.second.end(); it++)
	//{
	//		union_item.insert(*it);			//合并antecedent和consequent
	//}	

	////vector<int> rule_3_attribute = bodon_ruleset[s_set_pair[i]];				// 【错误写法】，导致找不到rule；因为bodon_ruleset使用内码
	//vector<int> rule_3_attribute = bodon_ruleset[pair<item,item>(antecedent_encoding, consequent_encoding)];
	//	
	//if(rule_3_attribute.empty()) {cout<<"Error--display_rule()! bodon_ruleset中找不到此敏感规则\n"; getchar(); exit(-1);}

	//unsigned long union_sup = rule_3_attribute[2];								// 整个rule的支持度
	//unsigned long antecedent_sup = rule_3_attribute[0];							// rule左边antecedent的支持度

	//cout<<"   union_sup = "<<union_sup<<"   antecedent_sup = "<<antecedent_sup<<endl;
}



void display_rule2(pair<item,item> a_rule)
{
	cout<<"\n------------------------------------------------------------\n";
	cout<<"\n ( ";

	item::const_iterator it = a_rule.first.begin();
	for(; it != a_rule.first.end()  ; it++)
	{
		cout<<*it<<" ";
	}
	cout<<") -> ( ";
	for(it = a_rule.second.begin(); it!=a_rule.second.end(); it++)
	{
		cout<<*it<<" ";
	}
	cout<<")"<<endl;

}


void write_DB_to_file()
{


}


void write_performance(int alpha, int beta, int gamma, int modified_trans_num, int modified_items_num)
{
	ofstream outfile;
	char filename[256];

	cout<<"\n>>> BEGIN to write the outcome of App_Intell_2014 into file!\n";

	sprintf_s(filename,256,"%s_MST_%.4f_MCT_%.2f_rules_%d_SM_%.2f_A01_%.2f_[New_BA]_outcome.dat",datafrom, MST,MCT, s_set_pair.size(), SM, A01 );
	outfile.open(filename, ios::out);

	double alpha_perc = (double)alpha/(s_set_pair.size());
	double beta_perc = (double)beta/( num_strong_rules_by_trie-s_set_pair.size() );
	double gamma_perc = (double)gamma/num_strong_rules_by_trie;
	double trans_modify_perc = (double)modified_trans_num/db_size;

	outfile<<"alpha \tbeta \tgamma \tmodified_trans_num \tremoved_items_num\n";

	outfile<<"(alpha,beta,gamma):\t\t"<<alpha<<", "<<beta<<", "<<gamma<<endl;
	outfile<<"(modified_trans_num, modified_items_num):  "<<modified_trans_num<<", "<<modified_items_num<<endl;
	outfile<<"(alpha,beta,gamma, trans_modiy)(percent %):\t"<<100*alpha_perc<<", ";
	outfile<<100*beta_perc <<", "<<100*gamma_perc << ", " << 100*trans_modify_perc << endl;
	outfile<<"(#item_modify) "<<modified_items_num<<"\t"<<" (#trans_modify) "<<modified_trans_num<<"("<< trans_modify_perc <<"%)"<<endl;

	outfile<<"\n\n/*-----------------------------------------------------------------------------*/\n\n";

	outfile<<"算法BA 2007  删除items \n\n";

	outfile<<"datafile:                 "<<datafrom<<endl;
	outfile<<"database size:            "<<db_size<<endl;
	outfile<<"Frequent Itemset(1) or Association Rule(2):"<<FIM_or_AR<<endl;;   //【FIM与AR选择开关】关联规则相关
	outfile<<"Minimum support threshold (MST):      "<<MST<<endl;
	outfile<<"Minimum confidence threshold (MCT):   "<<MCT<<endl;				// 关联规则相关
	//outfile<<"Minimum lift threshold (MLT):         "<<MLT<<endl;				// 关联规则相关
	outfile<<"MST discount:             "<<MST_discount<<endl;
	outfile<<"MCT discount:             "<<MCT_discount<<endl;				

	outfile<<"Safety Margin (SM) :      "<<SM<<endl;
	outfile<<"The proportion of blocking 0's to 1's (A01):   "<<A01<<endl;


	outfile<<endl;
 
	outfile<<"Number of strong rules calculated by trie tree     =   "<<num_strong_rules_by_trie<<endl; 
	outfile<<"Number of frequent itemset calculated by trie tree =   "<<bodon_itemset_frequent.size()<<endl;
	outfile<<"Number of all itemset calculated by trie tree =   "<<bodon_itemset.size()<<endl;

	outfile<<"\n bodon_itemset.size() == "<<bodon_itemset.size();								//【2014-Aug-4 19:41 更正】
	outfile<<"\n bodon_itemset_frequent.size() == "<<bodon_itemset_frequent.size();				//【2014-Sep-2 18:35 更正】
	outfile<<"\n num_itemset_by_trie == "<< num_itemset_by_trie << endl;

	outfile<<"\n bodon_ruleset_all_MIN.size() == "<<bodon_ruleset_all_MIN.size()<<endl;		//bodon_ruleset包含pre-strong rules;因此统计必须用num_strong_rules_by_trie
	outfile<<"\n bodon_ruleset_all_MAX.size() == "<<bodon_ruleset_all_MAX.size()<<endl;		//bodon_ruleset包含pre-strong rules;因此统计必须用num_strong_rules_by_trie
	outfile<<"\n bodon_ruleset_strong.size()  == "<<bodon_ruleset_strong.size()<<endl;		//bodon_ruleset包含pre-strong rules;因此统计必须用num_strong_rules_by_trie
	outfile<<"\n bodon_ruleset_NBRS.size()    == "<<bodon_ruleset_NBRS.size()<<endl;		//bodon_ruleset包含pre-strong rules;因此统计必须用num_strong_rules_by_trie
				
	outfile<<" num_strong_rules_by_trie == "<<num_strong_rules_by_trie<<endl;

	outfile<<" \nbodon_ruleset_all包含pre-strong rules;因此统计必须用num_strong_rules_by_trie\n";
	outfile<<" \n由于复制了一棵树;，num_itemset_by_trie取值是bodon_itemset_frequent.size()的2倍 \n";



	outfile<<"\n/*-----------------------------------------------------------------------------*/\n\n";

	outfile<<"sensitive file: \t"<<sensfile<<endl;
	outfile<<"Count of sensitive item sets:"<<s_set.size()<<endl;
	outfile<<"sensitive item sets LIST:"<<endl;

	int i=0;
	for(vector<pair<item,item>>::const_iterator sen_CI=s_set_pair.begin();sen_CI!=s_set_pair.end();sen_CI++,i++)				
	{
		item antecedent, consequent;

		antecedent = sen_CI->first;
		consequent = sen_CI->second;

		set<int>::const_iterator CI1,CI2;

		outfile<<"\n";
		CI1=antecedent.begin();
		for(;CI1!=antecedent.end();CI1++)
		{
			outfile<< *CI1 <<" ";
		}

		if(AR==FIM_or_AR)
		{
			outfile<<"-> ";

			CI2=consequent.begin();
			for(;CI2!=consequent.end();CI2++)
			{
				outfile<< *CI2 <<" ";
			}
		}

		item antecedent_encoding, consequent_encoding;

		basket_recode3(antecedent, antecedent_encoding);	//转为内码
		basket_recode3(consequent, consequent_encoding);	

		pair<set<int>, set<int>>  a_rule = pair<set<int>, set<int>>(antecedent_encoding, consequent_encoding );
		
		//if(rule_3_attribute.empty()) {cout<<"Error in write_performance()! 没有在bodon_ruleset中找到此敏感规则"<<endl; getchar(); exit(-1);}

		unsigned long union_sup = bodon_ruleset_strong[a_rule].union_supp;					// 整个rule的支持度
		double conf =  bodon_ruleset_strong[a_rule].conf;

		outfile<<"\n   union_sup = "<<union_sup/(double)db_size ;
		outfile<<"\n   conf      = "<<conf<<endl;


		outfile<<"\n\n   N0 = " << N0_vec[i] << "       N1 = "<< N1_vec[i] <<endl;

	}

	//-----------------------------------------------------------------------------
	end_clock = clock();		//【结束计时】
	double total_run_time = (double)(end_clock - start_clock)/CLOCKS_PER_SEC;
	double total_run_time_after_apriori = (double)(end_clock - start_clock2)/CLOCKS_PER_SEC;

	outfile<<"total_run_time: \t"<<total_run_time<<endl;
	outfile<<"total_run_time_after_apriori"<<total_run_time_after_apriori<<endl;
	//------------------------------------------------------------------------------
								
	outfile.close();

	cout<<"\n>>> END to write the outcome of App_Intell14 into file!\n";

}


