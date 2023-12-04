#include"old_apriori.h"





collectset candi_sets;
collectset large_sets;
collectset reallarge_sets;
collectset prelarge_sets;

/*==========================================================
			 Association Rule	(AR)
===========================================================*/

collect_ruleset  candi_rulesets;

//-------------------BR方-----------------------------------

double MCT_discount = 0;

double MST_discount = 0;

//----------------------------------------------------------


unsigned int SIZE_threshold;

double MST;					// 【MST 最小支持度】,它的值从PPDM_dt_param文件中读取
double MCT;					// 【MCT 最小置信度】，它的值也从PPDM_dt_param文件中读取
double MLT;					// 【MLT 最小LIFT】,它的值也从PPDM_dt_param文件中读取	

unsigned int FIM_or_AR;			// 【程序切换开关】,取值为1表示频繁模式挖掘FIM, 取值为2表示关联规则挖掘


							//【注意】：	read_local_parameters()函数也做了相应修改！！！ 
							//				state0()也要修改

							//				load_one_candidate()函数, load_two_candidate()函数也修改,

							//				修改up_sup_threshold 为MST

							//				void write_output_file()函数， void write_output_obj()函数，void write_stat()函数

							//				本地参数文件 PPDM_dt_param.txt 需修改

							//				修改cal_support()函数中的1.5为ALPHA, ALPHA为符号常量


int	num_strong_rules = 0 ;

int num_strong_rules_by_trie = 0;

int num_itemset_by_trie = 0;

/*--------------------------------------------------------*/

/*===========================Trie 树========================*/

vector<itemtype> copy_new_code_inverse;

vector<itemtype> copy_new_code;

vector<unsigned long> copy_support_of_items_one;

vector< vector<unsigned long> > copy_temp_counter_array;

vector<vector<int>>  thin_DB ;



map<pair<item,item>,vector<int>> bodon_ruleset;   //在Apriori_Trie::assoc_rule_find()中存入Bodon apriori发现的关联规则

//map<pair<item,item>,vector<int>> bodon_ruleset_strong; //【暂时不用】




map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_strong;			//保存所有原先数据库中strong rules

map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_NBRS;			//保存所有原先数据库中的NBRS rules



map<pair<item,item>, T_RULE_UNKNOWN>  bodon_ruleset_all_MIN;					//保存所有rule的最小support和confidence

map<pair<item,item>, T_RULE_UNKNOWN>  bodon_ruleset_all_MAX;					//保存所有rule的最大support和confidence




map<item, int>  bodon_itemset;

map<item, int>  bodon_itemset_frequent;


Trie *p_main_trie_copy;

/*----------------------------------------------------------*/

DB database;

unsigned long MAX_del_trans;			// MAX_del_trans是底线， 决定了sup_l,但真正的删除笔数是length(即编码长度)

double Avg_Dist_MST;					// It defines the average distance of rules set to MST


double MAX_del_trans_perc = 0.3;

double  ALPHA=1.5;

/*----------------------------------------------------------*/


char datafrom[128];			// datafrom是数据源文件， 它的值从PPDM_dt_param文件中读取
char sensfile[128];			// 包含敏感项的文件名,它的值从PPDM_dt_param文件中读取


int sup_u = 0;
int sup_l = 0;
int sup_l_actual;

unsigned long db_size;



	
vector<set<int>> s_set;					// 敏感项集合,  区别在于它将每个敏感项的condition和consequene合并到一个set
										//				它的内容在 filter_sen_trans() 函数中 由s_set_pair导入
vector<pair<item,item>> s_set_pair;		// 敏感项集合， 包含了从文件读来的敏感项item, 它的每一个元素是pair类型，对应一个敏感项
										//				(第一个item表示condition,第二个item表示consequence)
int num_sen_item;

/****************************************************************************/


// 输入：
// 输出： 内存数据库database，  大小db_size 

void initial()										//将外存文件中的数据库---》读取进内存， 生成"非瘦身"“内存”数据库database
{													//相应的， 还有一个瘦身版数据库"thin_DB"
	cout<<"\n>>>Begin to Initialize!\n";
	cout << "Loading database..." << endl;

	ifstream dataset;
	string line;
	string token;
	int data = 0;
	int total_item_count = 0;
	int max_trans_size = 0;
	int trans_size = 0;
	int trans_count = 0;// 统计数据库中trans的数量 
	double ats;
	transaction trans;
	//********************
	dataset.open(datafrom, ios::in);

	if (!dataset.is_open())
		cout << "無法開啟dataset\n";

	// create memory database

	while (!dataset.eof())			// 把数据从文件读到一个个vector中，这些vector<int>又存储到database中（类型是vector< vector<int> >）
	{

		getline(dataset, line);
		//切成token,塞进ss
		stringstream ss(line);
		//将item一笔一笔读进token作判断
		trans_size = 0;
		while (ss >> token)
		{
			//将读进来的字串token，转换成int 型 data
			data = atoi(token.c_str());
			trans.push_back(data);			// trans是vector类型的，并无自动排序功能
											// 但经检查，dataset的每一行(每个Transaction)的数据都是按升序排放的
			total_item_count++;
			trans_size++;
		}
		if(trans_size > max_trans_size)
		{
			max_trans_size = trans_size;
		}

		// 为了保证后续的比对正确（频繁集与transaction比对；  敏感项与transaction的比对）
		// 建议这里对trans进行排序。
		sort(trans.begin(),trans.end(),less<int>());  //sort()是C++提供的库函数；是qsort()的改进版

		database.push_back(trans);
		trans.clear();

		trans_count++;
	}
	ats = total_item_count/database.size ();
	cout << "Max Transaction size :"<< max_trans_size  << endl;
	cout << "# of items :"<< total_item_count  << endl;
	cout << "Average Transaction size :"<< ats  << endl;
	cout << "The count of transactions in database : " << trans_count << endl;

	db_size = database.size();  // added by Cheng Peng

	dataset.close();

	cout<<">>>END of initializing!\n\n";
}

void cal_db_size()
{
	ifstream dataset;
	string line;
	unsigned long trans_count=0;
 
	dataset.open(datafrom, ios::in);
	if (!dataset.is_open())
		cout << "無法開啟dataset\n";

	while (!dataset.eof())			 
	{
		getline(dataset, line);
		trans_count++;
	}

	db_size = trans_count;

}



// 输出 MAX_del_trans

void cal_MAX_delTrans()
{
	unsigned long dif_sen_frequency=0;	// 统计敏感项频次超出sup_u的累计和

	filter_sen_trans();			// 【需扫描数据库】找出包含敏感项的数据库记录编号；并得到各敏感项的频次
								// 考虑到计算sup_l需要用到sen_frequency_array[]数组，因此在这里调用
								// 调用前提是已经读取sensitive item到s_set中, 已经读数据库记录到内存中(database)

	int num_sen = s_set.size();

	for( int j=0;j<num_sen;j++)				
	{															
		dif_sen_frequency = (sen_frequency_array[j] - sup_u ) + dif_sen_frequency;
	}

	Avg_Dist_MST = ( (double)dif_sen_frequency /(double)num_sen )/db_size;	// 计算敏感项到MST的平均距离

	if ( dif_sen_frequency * ALPHA > db_size*MAX_del_trans_perc )
	{
		cout<<"\n需要删除的trans比例 超过 MAX_del_trans_perc, 请检查每个sensitive item 超出sup_u的频次\n";
		//cout<<"重新确定合适的MAX_del_trans_perc\n";
		//getchar();	exit(-1);
		MAX_del_trans = (unsigned long)(db_size*MAX_del_trans_perc);
	}
	else
		MAX_del_trans = (int)(dif_sen_frequency * ALPHA) ;



//	MAX_del_trans = db_size*MAX_del_trans_perc;
}


//计算support
//输出： sup_u, sup_l
void cal_support()
{
	cout<<"\n>>>Begin to Calculate Support!\n\n";
	cout << "database.size: " << db_size << endl;
	sup_u = (int)(db_size * MST);
	cout << "sup_u = " << sup_u << endl;

	cal_MAX_delTrans();		// 由它再进一步调用 filter_sen_trans();	为每个sensitive rule检索出sensitive transactions 

	//sup_l = (int)((db_size-MAX_del_trans) * MST);		//sup_l的作用是实现保存pre_large集合，避免重新扫描数据库

	cout<<"\n因为是删除item, sup_u==sup_l\n"; 

	sup_l = sup_u ;									

	if(sup_l<0) {cout<<"ERROR in cal_support(): sup_l<0 "<<endl; getchar(); exit(-1);}

	cout << "sup_l = " << sup_l << endl;
	cout << "sup_u = " << sup_u << endl;

	cout<<"\n>>>END of Calculating Support!\n\n";
}


//产生candidate_1
void gen_candione()
{
	cout << "\n>>>Generate candidate_1..." << endl;
	itemset candione;
	int data = 0;	//将item转成int，以利于set内自动做排序
	item temp;	// item即set<int>, 有自动排序功能， 但vector类型不具有自动排序功能

	for (DBCI databaseCI=database.begin(); databaseCI!=database.end(); databaseCI++)
	{
		transaction trans_temp = *databaseCI;	
		for (transactionCI trans_tempCI=trans_temp.begin();trans_tempCI!=trans_temp.end(); trans_tempCI++)
		{
			temp.insert(*trans_tempCI);   //  "temp" only contains one int number

			if (!candione.count(temp))
				candione.insert(newdata(temp, 1) );
			else
				candione[temp]++;
			temp.erase(temp.begin(),temp.end());
		}
		//cout<<".";
	}
	candi_sets.push_back(candione);
	cout << ">>>\nEND:  Generate candidate_1 complete  候选集-1产生!\n"<< endl;
}


//产生large set_n (存在large_sets[n-1])
//传入candi_sets_n (存在candi_sets[n-1])
void gen_large(itemset inputset)
{
	cout << endl<<"Prepare to Gen_large..." << endl;
	itemset large_temp;
	itemset reallarge_temp;
	itemset prelarge_temp;
	itemsetCI iter = inputset.begin();
	while(iter != inputset.end())
	{
		//大于sup_u叫做reallarge
		if ((iter->second) >= sup_u)
		{
			large_temp.insert(newdata(iter->first,iter->second));
			reallarge_temp.insert(newdata(iter->first,iter->second));
		}
		//否则， 大于sup_l叫做prelarge
		else if (iter->second >= sup_l)
		{
			large_temp.insert(newdata(iter->first,iter->second));
			prelarge_temp.insert(newdata(iter->first,iter->second));
		}
		iter++;
	}
	//塞进large_n集合 (存入large[n-1])
	large_sets.push_back(large_temp);
	//塞进reallarge_n集合 (存入reallarge[n-1])
	reallarge_sets.push_back(reallarge_temp);
	//塞进prelarge_n集合 (存入prelarge[n-1])
	prelarge_sets.push_back(prelarge_temp);

	cout << "\nGenerate large/prelarge/reallarge complete  !" << endl;
}


//large_n-1 join 出 candiset_n item sets; 注意： 是"large" n-1 set
void join(int candi_n)
{
	cout << "\n Join candidate " << candi_n << endl;

	itemset temp = large_sets[candi_n-2];
	itemsetCI tempCI;
	itemset temp2 = large_sets[candi_n-2];
	itemsetCI temp2CI;

	item old;
	itemCI oldCI;
	item old2;
	itemCI old2CI;

	item new_item;
	itemset aft_join;
	int a = 0;
	int count_join = 0;

	// temp and temp2 are two  map<item, int>
	// temp->first is item; temp->second is count;
	// create any combination of two sets of items, which's size is function input parameter candi_n

	for (tempCI=temp.begin(); tempCI!=temp.end(); tempCI++)	//频繁集里两两itemset的组合。如果结果集合元素个数是candi_n,则加入aft_join
	{
		old = tempCI->first;
		for (temp2CI=temp2.begin(); temp2CI!=temp2.end(); temp2CI++)	
		{
			old2 = temp2CI->first;
			
			for (oldCI=old.begin(); oldCI!=old.end(); oldCI++)
			{
				new_item.insert(*oldCI);
			}
			for (old2CI=old2.begin(); old2CI!=old2.end(); old2CI++)
			{
				new_item.insert(*old2CI);				//merge two items 
			}
			if (new_item.size() == candi_n)
			{
				count_join++;
				aft_join.insert(newdata(new_item, 0));
			}
			new_item.erase(new_item.begin(),new_item.end());
		}
		
		cout<<".";
	}
	candi_sets.push_back(aft_join);

	cout<<endl<<"count_join== "<< count_join<< endl;

	cout<<"\n END joining!   candidate_n = "<< candi_n<<endl;
}


//calculate support for each item in candi_sets[candi_n-1]
//扫描的是读入内存的database
void scanDB(int candi_n)
{
	cout << "Scan DB: caculate the count..." << endl;
	int count = 0;

	for (itemsetCI tempCI=candi_sets[candi_n-1].begin(); tempCI!=candi_sets[candi_n-1].end(); tempCI++)
	{
		item goalitem = tempCI->first;

		// for each item scan the memory database once
		for (DBCI db_CI=database.begin(); db_CI!=database.end(); db_CI++)
		{
			transaction goaltrans = *db_CI;
			itemCI itemCI = goalitem.begin();          //////// the same name for "itemCI"
			transactionCI transCI = goaltrans.begin();
			
			while (itemCI!=goalitem.end() && transCI!=goaltrans.end())
			{
				// Assume the vector<int> (for goaltrans) and the set<int> (for goalitem) are both sorted automatically
				// vector类型不会自动排序； set类型会自动排序
				// 对于vector类型的goaltrans, 因为dtatset中每一行的item是升序的，或者将database读入内存时， 每一行进行了排序。
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

			if (count == candi_n)
			{
				candi_sets[candi_n-1][goalitem]++;
			}
			count = 0;
		}

		cout<<".";
	}
	cout<<endl<<"Scan DB complete!"<<endl;
}

//
void output(itemset outputset, int n, int type)
{
	//cout << "output" << endl;
	ofstream outfile;
	itemsetCI iter;
	item temp;
	itemCI tempCI;
	string state;

	if (type == 0)
		state = "candidate ";
	else if (type == 1)
		state = "large ";
	else if (type == 2)
		state = "prelarge ";

	outfile.open(outputfile, ios::app);  //将candidate, large和prelarge set放在了同一个文件中
	outfile <<  "Output " << state << n << " ..." << endl;
	if (outputset.size()!=0)
	{
		iter = outputset.begin();	
		while (iter != outputset.end())
		{
			temp = iter->first;
			outfile << "< ";
			for (tempCI=temp.begin(); tempCI!=temp.end(); tempCI++)
			{
				outfile << *tempCI << " ";
			}
			outfile << "> : " << iter->second << " times" << endl;
			iter++;
		}
	}
	outfile.close();
}


// 输出频繁项到外部文件
void output_LARGEset(itemset outputset, int n, int type)
{
	//cout << "output" << endl;
	ofstream outfile;
	itemsetCI iter;
	item temp;
	itemCI tempCI;
	char state[32];
	char filename[80];

	if (type == 1)
		strcpy(state, "large");
	else if (type == 2)
		strcpy(state, "prelarge");
	else if (type == 3)
		strcpy(state, "reallarge");
	else
	{	
		cout<< "Error in output_largeset() function!!!!!\n"; exit(-1);
	}

	sprintf(filename, "%s_MST_%.3f_MCT_%.2f_MLT_%.2f_large_sets.dat ",datafrom, MST,MCT,MLT);
	outfile.open(filename, ios::app);

	outfile <<"\n------------------------------------------------------------\n";
	outfile <<  "Output "<< state <<" "<< n << "-frequent item " << " ...\n";
	outfile <<  "MST: "<<MST<<"  sup_u： "<<sup_u<<"  sup_l： "<<sup_l;
	outfile <<"\n------------------------------------------------------------\n";

	cout<<endl<<"开始：写入文件 "<<filename<< "  "<< state << endl;

	if (outputset.size()!=0)
	{
		iter = outputset.begin();	
		while (iter != outputset.end())
		{
			temp = iter->first;
			//outfile << "< ";
			for (tempCI=temp.begin(); tempCI!=temp.end(); tempCI++)
			{
				outfile << *tempCI <<" ";
			}
			//outfile << "> : " << iter->second << " times" << endl;

			outfile << ": " << iter->second << endl;

			iter++;
		}
	}

	outfile <<  "\nEND --- Output " << state << n << "- frequent item " << " ...\n\n" ;

	outfile.close();

	cout<<endl<<"结束： 写入文件  "<<filename<< "  "<< state << endl;

}


//载入one candidate itemset
void load_one_candidate()
{
	cout << "Loading one candidate..." << endl;

	ifstream dataset;
	string line;
	string token;
	int data = 0;
	item iter;
	int count = 0;
	itemset large_temp;
	itemset prelarge_temp;

	//itemset candidate_temp;  // added by Cheng Peng
	char candi_filename[128];

	char filename[128];
	sprintf(filename, "%s_MST_%.2f_MCT_%.2f_MLT_%.2f_large_sets.dat ",datafrom, MST,MCT,MLT);
	ofstream outfile;
	outfile.open(filename,ios::out);	// 因为是第一次写large_sets文件，先清除保存上次运行结果的文件内容
	outfile.close();				

	sprintf(candi_filename, "%s_candi_1",datafrom);
	dataset.open(candi_filename, ios::in);

	if (!dataset.is_open())
		cout << "无法打开 one candidate file\n";
	while (!dataset.eof())
	{
		//读进一行
		getline(dataset, line);		//切成token,塞进ss

		stringstream ss(line);

		ss >> token;				//读入第一个值
		data = atoi(token.c_str());
		iter.insert(data);  

		ss >> token;				//读入第二个值
		count = atoi(token.c_str());

		if(count >= sup_l)			//包括real_large和pre_large
		{
			large_temp.insert(newdata(iter, count) );
			if(count < sup_u)
				prelarge_temp.insert( newdata(iter, count) );
		}
		//candidate_temp.insert(newdata(iter, count) );  // added by Cheng Peng

		iter.clear ();
	}

	large_sets.push_back(large_temp);
	output_LARGEset(large_temp,1,1);		//将“1-large” 项写入文件
	large_temp.clear();

	prelarge_sets.push_back(prelarge_temp);
	output_LARGEset(prelarge_temp,1,2);		//将“1-prelarge” 项写入文件
	prelarge_temp.clear();

	cout<<"\n Complete loading one candidate!!\n";

	//candi_sets.push_back(candidate_temp);		// added by Cheng Peng
	//candidate_temp.clear();					

	dataset.close();



	cout<<"\n当前 1-item prelarge项数量:  "<< prelarge_sets.size()<<endl;
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




}

//载入two candidate itemset
void load_two_candidate()
{
	cout << "\nLoading two candidate..." << endl;

	ifstream dataset;
	string line;
	string token;
	int data = 0;
	int data2 = 0;
	item iter;
	int count = 0;

	itemset large_temp;
	itemset prelarge_temp;
	//itemset candidate_temp;  // added by Cheng Peng

	char candi_filename[80];

	//sprintf(candi_filename, "%s_%.3f_candi_2",datafrom, MST); 【将来采用加了MST的文件名！！！】

	sprintf(candi_filename, "%s_candi_2",datafrom);
	dataset.open(candi_filename, ios::in);

	if (!dataset.is_open())
		cout << "無法開啟 two candidate file\n";
	while (!dataset.eof())
	{

		getline(dataset, line);

		stringstream ss(line);

		ss >> token;
		//读入第一个值
		data = atoi(token.c_str());
		iter.insert(data);  
		ss >> token;
		//读入第二个值
		data2 = atoi(token.c_str());
		iter.insert(data2);

		ss >> token;
		//读入第三个值
		count = atoi(token.c_str());

		//新增到itemset中
		if(count >= sup_l)
		{
			large_temp.insert( newdata(iter, count) );
			if( count < sup_u)
				prelarge_temp.insert( newdata(iter, count) );
		}

		//candidate_temp.insert(newdata(iter, count) );  // added by Cheng Peng

		iter.clear ();
	}
	large_sets.push_back(large_temp);
	output_LARGEset(large_temp,2,1);		//将“2-large” 项写入文件
	large_temp.clear();

	prelarge_sets.push_back(prelarge_temp);
	output_LARGEset(prelarge_temp,2,2);		//将“2-prelarge” 项写入文件
	prelarge_temp.clear();

	//candi_sets.push_back(candidate_temp);		// added by Cheng Peng
	//candidate_temp.clear();					// added by Cheng Peng

	dataset.close();

	cout<<"\nComplete loading two-item candidate!\n";

}


void get_trie_data()
{
	cout<<"\n\nBEGIN:  void get_trie_data()\n\n";

	// "1-item"数据结构
	copy_new_code_inverse = p_apriori->get_pointer_to_input_output_manager()->get_new_code_inverse()  ;
	copy_support_of_items_one = p_apriori->get_pointer_to_input_output_manager()->get_support_of_items_one();

	cout<<"\n1-item数据结构\n";
	//for(unsigned n=0; n<copy_new_code_inverse.size(); n++)
	//{
	//	cout<<copy_new_code_inverse[n];
	//	cout<<"  ("<<copy_support_of_items_one[n]<<")"<<endl;
	//}

	copy_new_code = p_apriori->get_pointer_to_input_output_manager()->get_new_code()  ;


	cout<<"\n\n";
	//for(unsigned n=0; n<copy_new_code.size(); n++)
	//	cout<<copy_new_code[n]<<" ";

	copy_temp_counter_array;
	//【注】：“2-item” 数据结构 已经通过全局变量 copy_temp_counter_array , 位于Apriori_Trie::delete_infrequent_two()函数
	//【注】：“2-item” 数据结构 已经通过全局变量 copy_temp_counter_array , 位于Apriori_Trie::delete_infrequent_two()函数

	cout<<"\n2-item数据结构\n";

	//for(unsigned i=0; i<copy_temp_counter_array.size(); i++)
	//{
	//	//for(unsigned j=0;j<copy_temp_counter_array.size() - (i+1); j++)  //上界小1；缺了一个
	//	for(unsigned j=0;j<copy_temp_counter_array[i].size(); j++)
	//	{
	//		cout<<copy_temp_counter_array[i][j]<<" ";
	//	}
	//	cout<<endl;
	//}

	cout<<"\nEND:  void get_trie_data()\n\n";

}


void transfer_trie_to_largesets()
{

	cout<<"\n\nBEGIN: void transfer_trie_to_largesets()! \n";


	// 传导1-item large sets，  注意 “内部编码”还原成了“商品编码”
	itemset large_temp;
	itemset prelarge_temp;
	for(unsigned n=0; n<copy_new_code_inverse.size(); n++)
	{	
		item temp_item;

		if( copy_support_of_items_one[n] >=  (unsigned)sup_l )		//这个条件是否为废话？？copy_support_of_items_one[n]所有元素都大于sup_l
		{
			temp_item.insert( copy_new_code_inverse[(int)n]  );				//  译码,还原成“原商品编号”
			large_temp.insert( newdata(temp_item, (int)copy_support_of_items_one[n])  );

			if( copy_support_of_items_one[n] >= (unsigned)sup_l  &&  copy_support_of_items_one[n] < (unsigned)sup_u )
			{
				prelarge_temp.insert( newdata(temp_item, copy_support_of_items_one[n])  );
			}
		}

	}
	large_sets.push_back(large_temp);
	large_temp.clear();
	prelarge_sets.push_back(prelarge_temp);
	prelarge_temp.clear();

	// 传导2-item large sets， 注意 “内部编码”还原成了“商品编码”
	for(unsigned i=0; i<copy_temp_counter_array.size(); i++)
	{
		cout<<"\n"<<i<<" copy_temp_counter_array[i].size() = "<<copy_temp_counter_array[i].size();
		cout<<"   "<<"copy_temp_counter_array.size()-(i+1) = "<<copy_temp_counter_array.size()-(i+1)<<"\n";

		//上面的两个cout证明copy_temp_counter_array.size()-(i+1) 比 copy_temp_counter_array[i].size() 小1
		//因为copy_temp_counter_array.size()【即二维矩阵的行数】，已经比1-频繁项的数量小1了【相当于已经减1了】

		for(unsigned j=0;j<copy_temp_counter_array[i].size(); j++)  
		{
			item temp_item;
			temp_item.insert( copy_new_code_inverse[ (int)i ]  );				//   译码 存入 larger_sets, 还原成“原商品编号”
			temp_item.insert( copy_new_code_inverse[ (int)(j+(i+1)) ]  );
			int temp_count= (int)copy_temp_counter_array[i][j]  ;

			if(temp_count>=sup_l)
				large_temp.insert( newdata(temp_item, temp_count ) );
			if(temp_count>=sup_l && temp_count<sup_u) 
				prelarge_temp.insert( newdata(temp_item, temp_count ) );
		}
	}
	large_sets.push_back(large_temp);
	large_temp.clear();
	prelarge_sets.push_back(prelarge_temp);
	prelarge_temp.clear();

	cout<<"\nEND: void transfer_trie_to_largesets()! \n\n";
}



void basket_recode2( const set<int>& original_basket, vector<int>& new_basket )  // set -> vector  编码 (且滤去非频繁item)
{
   new_basket.clear();
   for( set<int>::iterator it_basket = original_basket.begin(); it_basket != original_basket.end(); it_basket++ )
     if( copy_new_code[*it_basket] ) 
		 new_basket.push_back( copy_new_code[*it_basket]-1 );	//如果当前 商品（ID编号为*it_basket)是频繁的（new_code[*it_basket]取值非0） 
																//则将其“新编号”(对应“1-频繁”商品在 所有“1-频繁”商品中 按频次升序排列时 的序号<从0开始计数>)存入new_basket向量
   sort( new_basket.begin(), new_basket.end() );				
}																//【注意】：copy_new_code[]和copy_new_code_inverse[]只对频繁项进行编解码


void basket_recode3( const set<int>& original_basket, set<int>& new_basket ) //set -> set  编码
{
   new_basket.clear();
   for( set<int>::iterator it_basket = original_basket.begin(); it_basket != original_basket.end(); it_basket++ )
   {
     if( copy_new_code[*it_basket] ) 
		 new_basket.insert( copy_new_code[*it_basket]-1 );	    //set集合可自动排序	 	
   }
}	


void basket_decode3( const set<int>& encoding_basket, set<int>& decoding_basket ) //set -> set  解码
{
   decoding_basket.clear();
   for( set<int>::iterator it_basket = encoding_basket.begin(); it_basket != encoding_basket.end(); it_basket++ )
   {
		 decoding_basket.insert( copy_new_code_inverse[*it_basket] );	    //set集合可自动排序	 	
   }
}	





int read_in_a_line_thin_DB( set<int>& basket, ifstream & basket_file2 )		//我写的 函数重载: 需通过实参 指定文件对象
{																						
   if( basket_file2.eof() ) return 0;
   char          c;
   unsigned int      pos;

   basket.clear();
   do														   //  读取一行（对应一个basket,也称transaction）
   {
      int item = 0;
      pos = 0;
      basket_file2.get(c);
      while(basket_file2.good() && (c >= '0') && (c <= '9'))	   //  读取一个item（对应一个“商品编号”）
      {
         item *= 10;
         item += int(c)-int('0');
         basket_file2.get(c);
         pos++;
      }
      if( pos ) basket.insert( item );
   }														   //  结果（一行数据）存于形式参数 basket中（引用类型）
   while( !basket_file2.eof() && c != '\n' );
   return 1;
}


//注意： get_thin_DB（）要用到针对size_threshold==2时的 copy_new_code[]数据结构
void get_thin_DB( vector<vector<int>>  & thin_DB, char basket_filename[128] ) 		// 我新增加的成员	
{

	cout<<"\n\nBEGIN to get_thin_DB......\n";

    set<int> basket;
    vector<int> basket_v;	

	ifstream temp_file_obj;
	temp_file_obj.open(basket_filename);

	while(  read_in_a_line_thin_DB( basket, temp_file_obj ) )	//扫描数据库每一行， 瘦身！
	{
		  basket_recode2(basket, basket_v);
		  thin_DB.push_back(basket_v);
		  // read_in_a_line()函数和basket_recode()函数 已经对basket和basket_v向量 进行了清空 
	}

	temp_file_obj.close();


	//输出瘦身数据库--验证
	/*
	for(unsigned int i=0; i< thin_DB.size(); i++)
	{
		cout<<"thin_DB 第"<<i<<"行：";
		for(unsigned int j=0; j<thin_DB[i].size(); j++)
		{
			cout <<thin_DB[i][j]<<" ";
		}
		cout<<"\n";
	}
	*/

	db_size=thin_DB.size();

	cout<<"\n--The thin DB size is "<<db_size<<"\n";

	cout<<"\nEND to get_thin_DB \n\n";
}




void apriori_bodon()
{

	// 调用bodon 的apriori快速算法！
	char outcomefile[128];
	sprintf(outcomefile, "%s_MST_%.4f_MCT_%.2f_large_sets.dat ",datafrom, MST,MCT);

	 															
	//double prelarge_MST = (double)(1-MAX_del_trans_perc)*MST;		 // 注意：  删除item或增加item不会导致db_size增加，absolute support不会变 

	//Apriori_main(datafrom, outcomefile, MST, MCT, SIZE_threshold);   // 调用Bodon_apriori算法 找到频繁集 后， 已经将数据集文件关闭

	Apriori_main(datafrom, outcomefile, MST*MST_discount, MCT, SIZE_threshold); 	 // 为了找到NBRS, 需要将MST降低，已找到supp<MST 的itemset
																									   
	//cout<<"\nThe size of database is "<<db_size<<"\n\n";


	//if(2==SIZE_threshold)
	//{		
		
	     get_trie_data();				// 获取"1-item"数据结构	copy_new_code_inverse[], copy_support_of_items_one[], copy_new_code[]
										// 获取"2-item"数据结构 通过全局变量copy_temp_counter_array[][]得到
										// 【注意】： trie树 还未获取！

		 //transfer_trie_to_largesets();	// 因为要检索出prelarge项， 需事先计算sup_l和sup_u，而计算sup_l-->要知道MAX_del_trans-->要知道"敏感项是哪几个？"

		 get_thin_DB(thin_DB, datafrom);  // 获得 “瘦身” 的 “译码” 数据库
										  // 注意： get_thin_DB（）要用到针对size_threshold==2时的 copy_new_code[]数据结构

	     cout<<"\n已经获取了瘦身数据库！"<<endl;


	//}
	
	    // 修改void Apriori_Trie::assoc_rule_find( ）， 将关联规则保存到 map<pair<item,item>,vector<int>> bodon_ruleset;

		// 获取树结构拷贝
		p_main_trie_copy = new Trie(db_size);

		p_apriori->get_pointer_to_Apriori_Trie()->copy_Trie_tree(p_main_trie_copy);  //拷贝树, 真正在内存中构造出另一颗一模一样的树，占用不同内存空间

		p_apriori->get_pointer_to_Apriori_Trie()->write_content_to_file(*(p_apriori->get_pointer_to_input_output_manager()), p_main_trie_copy);  
		// 调用了我自己写的重载函数,将拷贝树main_trie_copy的结点信息写入 与原树main_trie相同的文件

		p_apriori->get_pointer_to_input_output_manager()->close(); // 原来的关闭操作在Apriori.cpp文件void Apriori::APRIORI_alg()函数中倒数第2行

		cout<<"\n\n打印拷贝树！\n\n";
		//p_main_trie_copy->show_content_preorder();//打印树

		//getchar();
		//verify_update_trie_tree();    //验证Trie::update_trie_tree()函数是否正确


		// 在cal_fitness（）函数中 【复制一棵树】， 根据选中的trans, 【扫描树】， 减小"trans对应结点"的count
		// 【对比两颗树】的结点count值：对比原“树” 及 更新后的“树”
		// 更进一步， 如果调用MFC把这棵【树绘制】出来更好。


}


// 假设删除thin_DB中所有trans, 更新树main_trie_copy, 验证 所有树结点的counter被所有trans更新后是否为0
void verify_update_trie_tree()
{

		for(unsigned int j=0;j<thin_DB.size();j++)// thin_DB 瘦身数据库
		{
			vector<int> goaltrans = thin_DB[j];   
			p_main_trie_copy->update_trie_tree(goaltrans.begin(), goaltrans.end());
		}
		cout<<"\n\n再次打印拷贝树！\n\n";
		p_main_trie_copy->show_content_preorder();//打印树

}



/* 此函数被state0()调用 */
void apriori_large_sets()
{

	if(AR == FIM_or_AR)
		cout<<"\n隐藏 关联规则 Rule hiding!\n";
	else
		cout<<"\n隐藏 频繁项 Frequent itemset hiding!\n";
	cout<<"\nSIZE_threshold == "<< SIZE_threshold<<endl;

	//cal_db_size();


	initial();		//【注意】： 尽管已经得到“瘦身"数据库thin_DB。参： get_thin_DB(thin_DB, datafrom);
					//			 initial()函数得到“非瘦身”“内存”数据库database, 以下情况用到它：
					//          【TKDE04】1.a,1.b,2.a,2.b需要根据原始数据库每行trans大小排序；算法为将来修改后的数据库的输出到文件做准备



	cal_support();		//  cal_support()主要计算出sup_u和sup_l.  计算sup_l-->需知道MAX_del_trans-->需知道"每个敏感项的support"-->调用filter_sen_trans();
						// 【注意】： 采用删除item 的策略 就不用 计算sup_u和sup_l， 不用 获取 “占用庞大内存”的database,
						//			  不用事先计算 过滤 sensitive trans, 而是在thin DB中 过滤

	transfer_trie_to_largesets();	// 因为要统计prelarge项， 需事先计算出sup_l和sup_u，因此本函数调用置于cal_support()之后
									// 因本函数传递"1-item"频繁项和"2-item"频繁项， 因此也需置于apriori_bodon之后


	// 原来apriori_bodon()的调用位置！！


	//int n = 1;	// frequent set大小	  //C1()

	//【节省时间】从文件载入large itemsets
	//load_one_candidate();
	//load_two_candidate();

/*
	gen_candione();
	// 输出--测试用
	output(candi_sets[0],1,0);
	// L1();
	gen_large(candi_sets[0]);
	// 输出--测试用

	cout << "Output large_1..." << endl;
	cout << "Output large_1..." << endl;
	//output(reallarge_sets[0],1,1);
	output_LARGEset(large_sets[0],1,1);			//记录1-item large set
	output_LARGEset(prelarge_sets[0], 1, 2);	//记录1-item prelarge set
	output_LARGEset(reallarge_sets[0], 1, 3);	//记录1-item reallarge set


	//cout << "Output prelarge_1..." << endl;
	//output(prelarge_sets[0],1,2);
	// ---------------------------------
	// 当large n的size不等于0
	// 继续寻找large n+1
	// ---------------------------------
	// 因large n实际存储于large_sets[n-1]

	
	while(large_sets[n-1].size()!=0)
	{
		++n;
							// 产生新的candidate set
		join(n);			// generate candidate sets of size n by bruce force search

		scanDB(n);			// compute support for each item in candi_sets[n-1] by bruce force search
							// 扫描数据库， 计算每个itemset在数据库 一行中 出现的次数(频次)
		
		//output(candi_sets[n-1],n,0);					// 输出--测试用
									
		gen_large(candi_sets[n-1]);						// 利用threshold判断，找出large set 
		
		cout << "Output large "<< n << " ..." << endl;
		//output(reallarge_sets[n-1],n,1);				// 输出--测试用
			
		output_LARGEset(large_sets[n-1], n, 1);			//记录n-item large set
		output_LARGEset(prelarge_sets[n-1], n, 2);		//记录n-item prelarge set
		output_LARGEset(reallarge_sets[n-1], n, 3);		//记录n-item reallarge set

	}
*/

	//至此， 已产生了所有的large item sets
	// At present, we have found all large item sets of all sizes

}





/*=======================================================================================================
				以下代码是有关从frequent item sets 产生 association rules
=========================================================================================================*/


/*An association rule is valid if its 【confidence, support and lift】 are greater than or equal to corresponding threshold values.*/
/* 以下代码 由频繁集 产生 候选rules  */
void apriori_gen_rules()
{
	cout<<"Begin to generate association rules!"<<endl;
	cout<<"开始产生关联规则"<<endl;

	item condition_part, consequence_part;
	for(unsigned int n= 0+1 ; n<large_sets.size(); n++)			//从“2-频繁”项开始推导 关联规则
	{
		map<pair<item,item>,vector<int>>  rule_n;

		for (itemsetCI tempCI=large_sets[n].begin(); tempCI!=large_sets[n].end(); tempCI++)
		{
			item goalitem = tempCI->first; 

			itemCI goalitemCI = goalitem.begin();
			for(; goalitemCI != goalitem.end(); goalitemCI++ )
			{
				int con_value = *goalitemCI;

				/*方法1*/
				consequence_part.insert( con_value );	//rule的结果部分
				goalitem.erase(goalitemCI);				//rule的条件部分

				int condition_sup = large_sets[n-1][goalitem];  //condition  支持度
				int consequence_sup = large_sets[0][consequence_part];  // consequence 支持度
				int union_sup = (tempCI->second) ;

				vector<int> sup_vec;
				sup_vec.push_back(condition_sup);
				sup_vec.push_back(consequence_sup);
				sup_vec.push_back(union_sup);

				rule_n.insert(  ruleset_insert_type(pair<item,item>(goalitem,consequence_part), sup_vec)  ); //将候选规则加入rule_n
				//rule_n.insert(  pair<pair<item,item>,vector<int>>(pair<item,item>(goalitem,consequence_part), sup_vec)  );  //写成这样也对！！

				goalitemCI = (goalitem.insert(con_value)).first;  //【小心】 此处insert()的返回值是pair<iterator,bool>, first为指向新插入元素的迭代器;
																  //		 item为set<int>类型， 可以自动排序， 因此插入的con_value仍在原来的位置
				consequence_part.clear();

				/*方法2*/
			
			}
		}

		candi_rulesets.push_back(rule_n);			// candi_rulesets是全局变量，类型vector<  map<pair<item,item>,vector<int>>  >

		rule_n.clear();
	}

	cout<<"产生关联规则 结束";

	output_candi_rules();

	output_real_rules();

	print_candi_rules();

	//getchar();
}


//文件输出 所有候选规则
void output_candi_rules()
{
	ofstream outfile;

	char filename[80];

	sprintf(filename, "%s_MST_%.3f_MCT_%.2f_MLT_%.2f_rules.dat ",datafrom, MST,MCT,MLT);
	outfile.open(filename, ios::out);

	outfile<<"\n---------------------------------------------------------------------------\n";
	outfile <<  "\tOutput " << "- Association Rules Candidates" << " ...";
	outfile<<"\n---------------------------------------------------------------------------\n";
	
	cout<<endl<<"开始：写入文件 "<<filename<< "  "<< " -- association rules"<< endl;

	for(unsigned n=0;n<candi_rulesets.size(); n++)
	{
		map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		for(; CI!=candi_rulesets[n].end(); CI++)
		{
			 item condition,consequence;     //item就是set<int>类型
			 item::const_iterator it_condition, it_consequence;

			 condition = CI->first.first;
			 consequence = CI->first.second;

			 it_condition = condition.begin();
			 it_consequence = consequence.begin();

			 double union_sup =  (double)(CI->second[2])/db_size;
			 double rule_confidence = (double)(CI->second[2])/(CI->second[0]);
			 double rule_lift = (double)(CI->second[2])/(CI->second[0])*(double)db_size/CI->second[1];

			 outfile<<" (";
			 for(; it_condition!=condition.end();it_condition++)
			 {
				outfile<<" "<<*it_condition;			 
			 }

			 outfile<<" ) -> (";
			 for(; it_consequence != consequence.end(); it_consequence++)
			 {
				outfile<<" "<<*it_consequence<<" ) \t" ;
			 }

			 outfile.setf(ios::fixed);   
			 //setiosflags(ios::fixed)

			 outfile<<" condition:"<<CI->second[0]<<" consequence:"<<CI->second[1]<<setprecision(3)<<" union_sup:"<<union_sup;
			 outfile<<" confidence:"<<rule_confidence;
			 outfile<<" lift:"<< rule_lift <<endl;
		
		}//CI
	
	}//n

	outfile.close();
		
	cout<<endl<<"结束：写入文件 "<<filename<< "  "<< " -- association rules"<< endl;

}


//文件输出 大于MST,MCT,MLT的规则
void output_real_rules()
{
	ofstream outfile;

	char filename[80];

	sprintf(filename, "%s_MST_%.3f_MCT_%.2f_MLT_%.2f_rules.dat ",datafrom, MST,MCT,MLT);
	outfile.open(filename, ios::app);

	//outfile.setf(ios::fixed);
	outfile<<"\n---------------------------------------------------------------------------\n";
	outfile <<  "\tOutput " << "- Association Rules ( >MST and >MCT and >MLT )" << " ...\n";
	outfile <<  "\t       " << "- MST:"<<MST<<"  MCT:"<<MCT<<"  MLT:"<<MLT;
	outfile<<"\n---------------------------------------------------------------------------\n";
	cout<<endl<<"开始：写入文件 "<<filename<< "  "<< " -- association rules"<< endl;

	num_strong_rules = 0;		// 全局变量， 表示strong rules的数量。 
								// 注意：candi_rulesets里除了strong rules, 还有pre_large itemsets生成的或置信度小于MCT的rules

	for(unsigned n=0;n<candi_rulesets.size(); n++)
	{
		map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		for(; CI!=candi_rulesets[n].end(); CI++)
		{
			 item condition,consequence;     //item就是set<int>类型
			 item::const_iterator it_condition, it_consequence;

			 condition = CI->first.first;
			 consequence = CI->first.second;

			 it_condition = condition.begin();
			 it_consequence = consequence.begin();

			 double union_sup =  (double)(CI->second[2])/db_size;
			 double rule_confidence = (double)(CI->second[2])/(CI->second[0]);
			 double rule_lift = (double)(CI->second[2])/(CI->second[0])*(double)db_size/CI->second[1];

			 if( union_sup>=MST && rule_confidence>=MCT && rule_lift>=MLT)
			 {
				num_strong_rules++;

				outfile<<" (";
				for(; it_condition!=condition.end();it_condition++)
				{
					outfile<<" "<<*it_condition;			 
				}

				outfile<<" ) -> (";
				for(; it_consequence != consequence.end(); it_consequence++)
				{
					outfile<<" "<<*it_consequence;
				}

				outfile<<" ) \t" ;

				outfile.setf(ios::fixed);   
				//setiosflags(ios::fixed)

				outfile<<" condition:"<<CI->second[0]<<" consequence:"<<CI->second[1]<<setprecision(3)<<" union_sup:"<<union_sup;
				outfile<<" confidence:"<<rule_confidence;
				outfile<<" lift:"<< rule_lift <<endl;
			 }
		
		}//CI
	
	}//n

	outfile.close();
		
	cout<<endl<<"结束：写入文件 "<<filename<< "  "<< " -- association rules"<< endl;
	
}


// 屏幕打印 候选规则
void print_candi_rules()
{

	for(unsigned int n= 0; n<candi_rulesets.size(); n++)
	{
		map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		for(; CI != candi_rulesets[n].end(); CI++)
		{
			item condition, consequence;
			vector<int> sup_vec;

			condition = CI->first.first;
			consequence = CI->first.second;
			sup_vec = CI->second;

			double rule_confidence =  (double)(sup_vec[2])/(sup_vec[0]);
			double lift =  (double)sup_vec[2] /(double)(sup_vec[0]) ;     // lift(X->Y) = (sup(XUY)/db_size)/[(sup(X)/db_size)*(sup(Y)/db_size)]
			lift = lift*(double)db_size/(double)sup_vec[1];
			
			if(rule_confidence > MCT)
				print_one_rule(condition, consequence, sup_vec[0], sup_vec[1], sup_vec[2], lift );
		
		}
	}



}

// 屏幕打印 一条规则
void print_one_rule(item condition, item consequence, int condition_sup, int consequence_sup, int union_sup, double lift)
{
	
	item::const_iterator it_temp1, it_temp2;

	cout.setf(ios::fixed);   

	cout<<" (";
	for(it_temp1 = condition.begin(); it_temp1 != condition.end(); it_temp1++ )
	{
		cout<<" "<< *it_temp1;	
	}

	cout<<" ) -> (";
	for(it_temp2 = consequence.begin(); it_temp2 != consequence.end(); it_temp2++ )
	{
		cout<<" "<< *it_temp2;	
	}
	cout<<" )   condi:"<<condition_sup<<"  conse:"<<consequence_sup<<"  union:"<<(double)union_sup/db_size;
	std::setprecision(4);
	cout<<"  confidence:"<<(double)union_sup/condition_sup<<"  lift:"<<lift<<"\n";

}