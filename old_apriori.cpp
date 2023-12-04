#include"old_apriori.h"





collectset candi_sets;
collectset large_sets;
collectset reallarge_sets;
collectset prelarge_sets;

/*==========================================================
			 Association Rule	(AR)
===========================================================*/

collect_ruleset  candi_rulesets;

//-------------------BR��-----------------------------------

double MCT_discount = 0;

double MST_discount = 0;

//----------------------------------------------------------


unsigned int SIZE_threshold;

double MST;					// ��MST ��С֧�ֶȡ�,����ֵ��PPDM_dt_param�ļ��ж�ȡ
double MCT;					// ��MCT ��С���Ŷȡ�������ֵҲ��PPDM_dt_param�ļ��ж�ȡ
double MLT;					// ��MLT ��СLIFT��,����ֵҲ��PPDM_dt_param�ļ��ж�ȡ	

unsigned int FIM_or_AR;			// �������л����ء�,ȡֵΪ1��ʾƵ��ģʽ�ھ�FIM, ȡֵΪ2��ʾ���������ھ�


							//��ע�⡿��	read_local_parameters()����Ҳ������Ӧ�޸ģ����� 
							//				state0()ҲҪ�޸�

							//				load_one_candidate()����, load_two_candidate()����Ҳ�޸�,

							//				�޸�up_sup_threshold ΪMST

							//				void write_output_file()������ void write_output_obj()������void write_stat()����

							//				���ز����ļ� PPDM_dt_param.txt ���޸�

							//				�޸�cal_support()�����е�1.5ΪALPHA, ALPHAΪ���ų���


int	num_strong_rules = 0 ;

int num_strong_rules_by_trie = 0;

int num_itemset_by_trie = 0;

/*--------------------------------------------------------*/

/*===========================Trie ��========================*/

vector<itemtype> copy_new_code_inverse;

vector<itemtype> copy_new_code;

vector<unsigned long> copy_support_of_items_one;

vector< vector<unsigned long> > copy_temp_counter_array;

vector<vector<int>>  thin_DB ;



map<pair<item,item>,vector<int>> bodon_ruleset;   //��Apriori_Trie::assoc_rule_find()�д���Bodon apriori���ֵĹ�������

//map<pair<item,item>,vector<int>> bodon_ruleset_strong; //����ʱ���á�




map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_strong;			//��������ԭ�����ݿ���strong rules

map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_NBRS;			//��������ԭ�����ݿ��е�NBRS rules



map<pair<item,item>, T_RULE_UNKNOWN>  bodon_ruleset_all_MIN;					//��������rule����Сsupport��confidence

map<pair<item,item>, T_RULE_UNKNOWN>  bodon_ruleset_all_MAX;					//��������rule�����support��confidence




map<item, int>  bodon_itemset;

map<item, int>  bodon_itemset_frequent;


Trie *p_main_trie_copy;

/*----------------------------------------------------------*/

DB database;

unsigned long MAX_del_trans;			// MAX_del_trans�ǵ��ߣ� ������sup_l,��������ɾ��������length(�����볤��)

double Avg_Dist_MST;					// It defines the average distance of rules set to MST


double MAX_del_trans_perc = 0.3;

double  ALPHA=1.5;

/*----------------------------------------------------------*/


char datafrom[128];			// datafrom������Դ�ļ��� ����ֵ��PPDM_dt_param�ļ��ж�ȡ
char sensfile[128];			// ������������ļ���,����ֵ��PPDM_dt_param�ļ��ж�ȡ


int sup_u = 0;
int sup_l = 0;
int sup_l_actual;

unsigned long db_size;



	
vector<set<int>> s_set;					// �������,  ������������ÿ���������condition��consequene�ϲ���һ��set
										//				���������� filter_sen_trans() ������ ��s_set_pair����
vector<pair<item,item>> s_set_pair;		// ������ϣ� �����˴��ļ�������������item, ����ÿһ��Ԫ����pair���ͣ���Ӧһ��������
										//				(��һ��item��ʾcondition,�ڶ���item��ʾconsequence)
int num_sen_item;

/****************************************************************************/


// ���룺
// ����� �ڴ����ݿ�database��  ��Сdb_size 

void initial()										//������ļ��е����ݿ�---����ȡ���ڴ棬 ����"������"���ڴ桱���ݿ�database
{													//��Ӧ�ģ� ����һ����������ݿ�"thin_DB"
	cout<<"\n>>>Begin to Initialize!\n";
	cout << "Loading database..." << endl;

	ifstream dataset;
	string line;
	string token;
	int data = 0;
	int total_item_count = 0;
	int max_trans_size = 0;
	int trans_size = 0;
	int trans_count = 0;// ͳ�����ݿ���trans������ 
	double ats;
	transaction trans;
	//********************
	dataset.open(datafrom, ios::in);

	if (!dataset.is_open())
		cout << "�o���_��dataset\n";

	// create memory database

	while (!dataset.eof())			// �����ݴ��ļ�����һ����vector�У���Щvector<int>�ִ洢��database�У�������vector< vector<int> >��
	{

		getline(dataset, line);
		//�г�token,����ss
		stringstream ss(line);
		//��itemһ��һ�ʶ���token���ж�
		trans_size = 0;
		while (ss >> token)
		{
			//�����������ִ�token��ת����int �� data
			data = atoi(token.c_str());
			trans.push_back(data);			// trans��vector���͵ģ������Զ�������
											// ������飬dataset��ÿһ��(ÿ��Transaction)�����ݶ��ǰ������ŷŵ�
			total_item_count++;
			trans_size++;
		}
		if(trans_size > max_trans_size)
		{
			max_trans_size = trans_size;
		}

		// Ϊ�˱�֤�����ıȶ���ȷ��Ƶ������transaction�ȶԣ�  ��������transaction�ıȶԣ�
		// ���������trans��������
		sort(trans.begin(),trans.end(),less<int>());  //sort()��C++�ṩ�Ŀ⺯������qsort()�ĸĽ���

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
		cout << "�o���_��dataset\n";

	while (!dataset.eof())			 
	{
		getline(dataset, line);
		trans_count++;
	}

	db_size = trans_count;

}



// ��� MAX_del_trans

void cal_MAX_delTrans()
{
	unsigned long dif_sen_frequency=0;	// ͳ��������Ƶ�γ���sup_u���ۼƺ�

	filter_sen_trans();			// ����ɨ�����ݿ⡿�ҳ���������������ݿ��¼��ţ����õ����������Ƶ��
								// ���ǵ�����sup_l��Ҫ�õ�sen_frequency_array[]���飬������������
								// ����ǰ�����Ѿ���ȡsensitive item��s_set��, �Ѿ������ݿ��¼���ڴ���(database)

	int num_sen = s_set.size();

	for( int j=0;j<num_sen;j++)				
	{															
		dif_sen_frequency = (sen_frequency_array[j] - sup_u ) + dif_sen_frequency;
	}

	Avg_Dist_MST = ( (double)dif_sen_frequency /(double)num_sen )/db_size;	// ���������MST��ƽ������

	if ( dif_sen_frequency * ALPHA > db_size*MAX_del_trans_perc )
	{
		cout<<"\n��Ҫɾ����trans���� ���� MAX_del_trans_perc, ����ÿ��sensitive item ����sup_u��Ƶ��\n";
		//cout<<"����ȷ�����ʵ�MAX_del_trans_perc\n";
		//getchar();	exit(-1);
		MAX_del_trans = (unsigned long)(db_size*MAX_del_trans_perc);
	}
	else
		MAX_del_trans = (int)(dif_sen_frequency * ALPHA) ;



//	MAX_del_trans = db_size*MAX_del_trans_perc;
}


//����support
//����� sup_u, sup_l
void cal_support()
{
	cout<<"\n>>>Begin to Calculate Support!\n\n";
	cout << "database.size: " << db_size << endl;
	sup_u = (int)(db_size * MST);
	cout << "sup_u = " << sup_u << endl;

	cal_MAX_delTrans();		// �����ٽ�һ������ filter_sen_trans();	Ϊÿ��sensitive rule������sensitive transactions 

	//sup_l = (int)((db_size-MAX_del_trans) * MST);		//sup_l��������ʵ�ֱ���pre_large���ϣ���������ɨ�����ݿ�

	cout<<"\n��Ϊ��ɾ��item, sup_u==sup_l\n"; 

	sup_l = sup_u ;									

	if(sup_l<0) {cout<<"ERROR in cal_support(): sup_l<0 "<<endl; getchar(); exit(-1);}

	cout << "sup_l = " << sup_l << endl;
	cout << "sup_u = " << sup_u << endl;

	cout<<"\n>>>END of Calculating Support!\n\n";
}


//����candidate_1
void gen_candione()
{
	cout << "\n>>>Generate candidate_1..." << endl;
	itemset candione;
	int data = 0;	//��itemת��int��������set���Զ�������
	item temp;	// item��set<int>, ���Զ������ܣ� ��vector���Ͳ������Զ�������

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
	cout << ">>>\nEND:  Generate candidate_1 complete  ��ѡ��-1����!\n"<< endl;
}


//����large set_n (����large_sets[n-1])
//����candi_sets_n (����candi_sets[n-1])
void gen_large(itemset inputset)
{
	cout << endl<<"Prepare to Gen_large..." << endl;
	itemset large_temp;
	itemset reallarge_temp;
	itemset prelarge_temp;
	itemsetCI iter = inputset.begin();
	while(iter != inputset.end())
	{
		//����sup_u����reallarge
		if ((iter->second) >= sup_u)
		{
			large_temp.insert(newdata(iter->first,iter->second));
			reallarge_temp.insert(newdata(iter->first,iter->second));
		}
		//���� ����sup_l����prelarge
		else if (iter->second >= sup_l)
		{
			large_temp.insert(newdata(iter->first,iter->second));
			prelarge_temp.insert(newdata(iter->first,iter->second));
		}
		iter++;
	}
	//����large_n���� (����large[n-1])
	large_sets.push_back(large_temp);
	//����reallarge_n���� (����reallarge[n-1])
	reallarge_sets.push_back(reallarge_temp);
	//����prelarge_n���� (����prelarge[n-1])
	prelarge_sets.push_back(prelarge_temp);

	cout << "\nGenerate large/prelarge/reallarge complete  !" << endl;
}


//large_n-1 join �� candiset_n item sets; ע�⣺ ��"large" n-1 set
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

	for (tempCI=temp.begin(); tempCI!=temp.end(); tempCI++)	//Ƶ����������itemset����ϡ�����������Ԫ�ظ�����candi_n,�����aft_join
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
//ɨ����Ƕ����ڴ��database
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
				// vector���Ͳ����Զ����� set���ͻ��Զ�����
				// ����vector���͵�goaltrans, ��Ϊdtatset��ÿһ�е�item������ģ����߽�database�����ڴ�ʱ�� ÿһ�н���������
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

	outfile.open(outputfile, ios::app);  //��candidate, large��prelarge set������ͬһ���ļ���
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


// ���Ƶ����ⲿ�ļ�
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
	outfile <<  "MST: "<<MST<<"  sup_u�� "<<sup_u<<"  sup_l�� "<<sup_l;
	outfile <<"\n------------------------------------------------------------\n";

	cout<<endl<<"��ʼ��д���ļ� "<<filename<< "  "<< state << endl;

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

	cout<<endl<<"������ д���ļ�  "<<filename<< "  "<< state << endl;

}


//����one candidate itemset
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
	outfile.open(filename,ios::out);	// ��Ϊ�ǵ�һ��дlarge_sets�ļ�������������ϴ����н�����ļ�����
	outfile.close();				

	sprintf(candi_filename, "%s_candi_1",datafrom);
	dataset.open(candi_filename, ios::in);

	if (!dataset.is_open())
		cout << "�޷��� one candidate file\n";
	while (!dataset.eof())
	{
		//����һ��
		getline(dataset, line);		//�г�token,����ss

		stringstream ss(line);

		ss >> token;				//�����һ��ֵ
		data = atoi(token.c_str());
		iter.insert(data);  

		ss >> token;				//����ڶ���ֵ
		count = atoi(token.c_str());

		if(count >= sup_l)			//����real_large��pre_large
		{
			large_temp.insert(newdata(iter, count) );
			if(count < sup_u)
				prelarge_temp.insert( newdata(iter, count) );
		}
		//candidate_temp.insert(newdata(iter, count) );  // added by Cheng Peng

		iter.clear ();
	}

	large_sets.push_back(large_temp);
	output_LARGEset(large_temp,1,1);		//����1-large�� ��д���ļ�
	large_temp.clear();

	prelarge_sets.push_back(prelarge_temp);
	output_LARGEset(prelarge_temp,1,2);		//����1-prelarge�� ��д���ļ�
	prelarge_temp.clear();

	cout<<"\n Complete loading one candidate!!\n";

	//candi_sets.push_back(candidate_temp);		// added by Cheng Peng
	//candidate_temp.clear();					

	dataset.close();



	cout<<"\n��ǰ 1-item prelarge������:  "<< prelarge_sets.size()<<endl;
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




}

//����two candidate itemset
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

	//sprintf(candi_filename, "%s_%.3f_candi_2",datafrom, MST); ���������ü���MST���ļ�����������

	sprintf(candi_filename, "%s_candi_2",datafrom);
	dataset.open(candi_filename, ios::in);

	if (!dataset.is_open())
		cout << "�o���_�� two candidate file\n";
	while (!dataset.eof())
	{

		getline(dataset, line);

		stringstream ss(line);

		ss >> token;
		//�����һ��ֵ
		data = atoi(token.c_str());
		iter.insert(data);  
		ss >> token;
		//����ڶ���ֵ
		data2 = atoi(token.c_str());
		iter.insert(data2);

		ss >> token;
		//���������ֵ
		count = atoi(token.c_str());

		//������itemset��
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
	output_LARGEset(large_temp,2,1);		//����2-large�� ��д���ļ�
	large_temp.clear();

	prelarge_sets.push_back(prelarge_temp);
	output_LARGEset(prelarge_temp,2,2);		//����2-prelarge�� ��д���ļ�
	prelarge_temp.clear();

	//candi_sets.push_back(candidate_temp);		// added by Cheng Peng
	//candidate_temp.clear();					// added by Cheng Peng

	dataset.close();

	cout<<"\nComplete loading two-item candidate!\n";

}


void get_trie_data()
{
	cout<<"\n\nBEGIN:  void get_trie_data()\n\n";

	// "1-item"���ݽṹ
	copy_new_code_inverse = p_apriori->get_pointer_to_input_output_manager()->get_new_code_inverse()  ;
	copy_support_of_items_one = p_apriori->get_pointer_to_input_output_manager()->get_support_of_items_one();

	cout<<"\n1-item���ݽṹ\n";
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
	//��ע������2-item�� ���ݽṹ �Ѿ�ͨ��ȫ�ֱ��� copy_temp_counter_array , λ��Apriori_Trie::delete_infrequent_two()����
	//��ע������2-item�� ���ݽṹ �Ѿ�ͨ��ȫ�ֱ��� copy_temp_counter_array , λ��Apriori_Trie::delete_infrequent_two()����

	cout<<"\n2-item���ݽṹ\n";

	//for(unsigned i=0; i<copy_temp_counter_array.size(); i++)
	//{
	//	//for(unsigned j=0;j<copy_temp_counter_array.size() - (i+1); j++)  //�Ͻ�С1��ȱ��һ��
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


	// ����1-item large sets��  ע�� ���ڲ����롱��ԭ���ˡ���Ʒ���롱
	itemset large_temp;
	itemset prelarge_temp;
	for(unsigned n=0; n<copy_new_code_inverse.size(); n++)
	{	
		item temp_item;

		if( copy_support_of_items_one[n] >=  (unsigned)sup_l )		//��������Ƿ�Ϊ�ϻ�����copy_support_of_items_one[n]����Ԫ�ض�����sup_l
		{
			temp_item.insert( copy_new_code_inverse[(int)n]  );				//  ����,��ԭ�ɡ�ԭ��Ʒ��š�
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

	// ����2-item large sets�� ע�� ���ڲ����롱��ԭ���ˡ���Ʒ���롱
	for(unsigned i=0; i<copy_temp_counter_array.size(); i++)
	{
		cout<<"\n"<<i<<" copy_temp_counter_array[i].size() = "<<copy_temp_counter_array[i].size();
		cout<<"   "<<"copy_temp_counter_array.size()-(i+1) = "<<copy_temp_counter_array.size()-(i+1)<<"\n";

		//���������cout֤��copy_temp_counter_array.size()-(i+1) �� copy_temp_counter_array[i].size() С1
		//��Ϊcopy_temp_counter_array.size()������ά��������������Ѿ���1-Ƶ���������С1�ˡ��൱���Ѿ���1�ˡ�

		for(unsigned j=0;j<copy_temp_counter_array[i].size(); j++)  
		{
			item temp_item;
			temp_item.insert( copy_new_code_inverse[ (int)i ]  );				//   ���� ���� larger_sets, ��ԭ�ɡ�ԭ��Ʒ��š�
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



void basket_recode2( const set<int>& original_basket, vector<int>& new_basket )  // set -> vector  ���� (����ȥ��Ƶ��item)
{
   new_basket.clear();
   for( set<int>::iterator it_basket = original_basket.begin(); it_basket != original_basket.end(); it_basket++ )
     if( copy_new_code[*it_basket] ) 
		 new_basket.push_back( copy_new_code[*it_basket]-1 );	//�����ǰ ��Ʒ��ID���Ϊ*it_basket)��Ƶ���ģ�new_code[*it_basket]ȡֵ��0�� 
																//���䡰�±�š�(��Ӧ��1-Ƶ������Ʒ�� ���С�1-Ƶ������Ʒ�� ��Ƶ����������ʱ �����<��0��ʼ����>)����new_basket����
   sort( new_basket.begin(), new_basket.end() );				
}																//��ע�⡿��copy_new_code[]��copy_new_code_inverse[]ֻ��Ƶ������б����


void basket_recode3( const set<int>& original_basket, set<int>& new_basket ) //set -> set  ����
{
   new_basket.clear();
   for( set<int>::iterator it_basket = original_basket.begin(); it_basket != original_basket.end(); it_basket++ )
   {
     if( copy_new_code[*it_basket] ) 
		 new_basket.insert( copy_new_code[*it_basket]-1 );	    //set���Ͽ��Զ�����	 	
   }
}	


void basket_decode3( const set<int>& encoding_basket, set<int>& decoding_basket ) //set -> set  ����
{
   decoding_basket.clear();
   for( set<int>::iterator it_basket = encoding_basket.begin(); it_basket != encoding_basket.end(); it_basket++ )
   {
		 decoding_basket.insert( copy_new_code_inverse[*it_basket] );	    //set���Ͽ��Զ�����	 	
   }
}	





int read_in_a_line_thin_DB( set<int>& basket, ifstream & basket_file2 )		//��д�� ��������: ��ͨ��ʵ�� ָ���ļ�����
{																						
   if( basket_file2.eof() ) return 0;
   char          c;
   unsigned int      pos;

   basket.clear();
   do														   //  ��ȡһ�У���Ӧһ��basket,Ҳ��transaction��
   {
      int item = 0;
      pos = 0;
      basket_file2.get(c);
      while(basket_file2.good() && (c >= '0') && (c <= '9'))	   //  ��ȡһ��item����Ӧһ������Ʒ��š���
      {
         item *= 10;
         item += int(c)-int('0');
         basket_file2.get(c);
         pos++;
      }
      if( pos ) basket.insert( item );
   }														   //  �����һ�����ݣ�������ʽ���� basket�У��������ͣ�
   while( !basket_file2.eof() && c != '\n' );
   return 1;
}


//ע�⣺ get_thin_DB����Ҫ�õ����size_threshold==2ʱ�� copy_new_code[]���ݽṹ
void get_thin_DB( vector<vector<int>>  & thin_DB, char basket_filename[128] ) 		// �������ӵĳ�Ա	
{

	cout<<"\n\nBEGIN to get_thin_DB......\n";

    set<int> basket;
    vector<int> basket_v;	

	ifstream temp_file_obj;
	temp_file_obj.open(basket_filename);

	while(  read_in_a_line_thin_DB( basket, temp_file_obj ) )	//ɨ�����ݿ�ÿһ�У� ����
	{
		  basket_recode2(basket, basket_v);
		  thin_DB.push_back(basket_v);
		  // read_in_a_line()������basket_recode()���� �Ѿ���basket��basket_v���� ��������� 
	}

	temp_file_obj.close();


	//����������ݿ�--��֤
	/*
	for(unsigned int i=0; i< thin_DB.size(); i++)
	{
		cout<<"thin_DB ��"<<i<<"�У�";
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

	// ����bodon ��apriori�����㷨��
	char outcomefile[128];
	sprintf(outcomefile, "%s_MST_%.4f_MCT_%.2f_large_sets.dat ",datafrom, MST,MCT);

	 															
	//double prelarge_MST = (double)(1-MAX_del_trans_perc)*MST;		 // ע�⣺  ɾ��item������item���ᵼ��db_size���ӣ�absolute support����� 

	//Apriori_main(datafrom, outcomefile, MST, MCT, SIZE_threshold);   // ����Bodon_apriori�㷨 �ҵ�Ƶ���� �� �Ѿ������ݼ��ļ��ر�

	Apriori_main(datafrom, outcomefile, MST*MST_discount, MCT, SIZE_threshold); 	 // Ϊ���ҵ�NBRS, ��Ҫ��MST���ͣ����ҵ�supp<MST ��itemset
																									   
	//cout<<"\nThe size of database is "<<db_size<<"\n\n";


	//if(2==SIZE_threshold)
	//{		
		
	     get_trie_data();				// ��ȡ"1-item"���ݽṹ	copy_new_code_inverse[], copy_support_of_items_one[], copy_new_code[]
										// ��ȡ"2-item"���ݽṹ ͨ��ȫ�ֱ���copy_temp_counter_array[][]�õ�
										// ��ע�⡿�� trie�� ��δ��ȡ��

		 //transfer_trie_to_largesets();	// ��ΪҪ������prelarge� �����ȼ���sup_l��sup_u��������sup_l-->Ҫ֪��MAX_del_trans-->Ҫ֪��"���������ļ�����"

		 get_thin_DB(thin_DB, datafrom);  // ��� ������ �� �����롱 ���ݿ�
										  // ע�⣺ get_thin_DB����Ҫ�õ����size_threshold==2ʱ�� copy_new_code[]���ݽṹ

	     cout<<"\n�Ѿ���ȡ���������ݿ⣡"<<endl;


	//}
	
	    // �޸�void Apriori_Trie::assoc_rule_find( ���� ���������򱣴浽 map<pair<item,item>,vector<int>> bodon_ruleset;

		// ��ȡ���ṹ����
		p_main_trie_copy = new Trie(db_size);

		p_apriori->get_pointer_to_Apriori_Trie()->copy_Trie_tree(p_main_trie_copy);  //������, �������ڴ��й������һ��һģһ��������ռ�ò�ͬ�ڴ�ռ�

		p_apriori->get_pointer_to_Apriori_Trie()->write_content_to_file(*(p_apriori->get_pointer_to_input_output_manager()), p_main_trie_copy);  
		// ���������Լ�д�����غ���,��������main_trie_copy�Ľ����Ϣд�� ��ԭ��main_trie��ͬ���ļ�

		p_apriori->get_pointer_to_input_output_manager()->close(); // ԭ���Ĺرղ�����Apriori.cpp�ļ�void Apriori::APRIORI_alg()�����е�����2��

		cout<<"\n\n��ӡ��������\n\n";
		//p_main_trie_copy->show_content_preorder();//��ӡ��

		//getchar();
		//verify_update_trie_tree();    //��֤Trie::update_trie_tree()�����Ƿ���ȷ


		// ��cal_fitness���������� ������һ�������� ����ѡ�е�trans, ��ɨ�������� ��С"trans��Ӧ���"��count
		// ���Ա����������Ľ��countֵ���Ա�ԭ������ �� ���º�ġ�����
		// ����һ���� �������MFC����á������ơ��������á�


}


// ����ɾ��thin_DB������trans, ������main_trie_copy, ��֤ ����������counter������trans���º��Ƿ�Ϊ0
void verify_update_trie_tree()
{

		for(unsigned int j=0;j<thin_DB.size();j++)// thin_DB �������ݿ�
		{
			vector<int> goaltrans = thin_DB[j];   
			p_main_trie_copy->update_trie_tree(goaltrans.begin(), goaltrans.end());
		}
		cout<<"\n\n�ٴδ�ӡ��������\n\n";
		p_main_trie_copy->show_content_preorder();//��ӡ��

}



/* �˺�����state0()���� */
void apriori_large_sets()
{

	if(AR == FIM_or_AR)
		cout<<"\n���� �������� Rule hiding!\n";
	else
		cout<<"\n���� Ƶ���� Frequent itemset hiding!\n";
	cout<<"\nSIZE_threshold == "<< SIZE_threshold<<endl;

	//cal_db_size();


	initial();		//��ע�⡿�� �����Ѿ��õ�������"���ݿ�thin_DB���Σ� get_thin_DB(thin_DB, datafrom);
					//			 initial()�����õ������������ڴ桱���ݿ�database, ��������õ�����
					//          ��TKDE04��1.a,1.b,2.a,2.b��Ҫ����ԭʼ���ݿ�ÿ��trans��С�����㷨Ϊ�����޸ĺ�����ݿ��������ļ���׼��



	cal_support();		//  cal_support()��Ҫ�����sup_u��sup_l.  ����sup_l-->��֪��MAX_del_trans-->��֪��"ÿ���������support"-->����filter_sen_trans();
						// ��ע�⡿�� ����ɾ��item �Ĳ��� �Ͳ��� ����sup_u��sup_l�� ���� ��ȡ ��ռ���Ӵ��ڴ桱��database,
						//			  �������ȼ��� ���� sensitive trans, ������thin DB�� ����

	transfer_trie_to_largesets();	// ��ΪҪͳ��prelarge� �����ȼ����sup_l��sup_u����˱�������������cal_support()֮��
									// �򱾺�������"1-item"Ƶ�����"2-item"Ƶ��� ���Ҳ������apriori_bodon֮��


	// ԭ��apriori_bodon()�ĵ���λ�ã���


	//int n = 1;	// frequent set��С	  //C1()

	//����ʡʱ�䡿���ļ�����large itemsets
	//load_one_candidate();
	//load_two_candidate();

/*
	gen_candione();
	// ���--������
	output(candi_sets[0],1,0);
	// L1();
	gen_large(candi_sets[0]);
	// ���--������

	cout << "Output large_1..." << endl;
	cout << "Output large_1..." << endl;
	//output(reallarge_sets[0],1,1);
	output_LARGEset(large_sets[0],1,1);			//��¼1-item large set
	output_LARGEset(prelarge_sets[0], 1, 2);	//��¼1-item prelarge set
	output_LARGEset(reallarge_sets[0], 1, 3);	//��¼1-item reallarge set


	//cout << "Output prelarge_1..." << endl;
	//output(prelarge_sets[0],1,2);
	// ---------------------------------
	// ��large n��size������0
	// ����Ѱ��large n+1
	// ---------------------------------
	// ��large nʵ�ʴ洢��large_sets[n-1]

	
	while(large_sets[n-1].size()!=0)
	{
		++n;
							// �����µ�candidate set
		join(n);			// generate candidate sets of size n by bruce force search

		scanDB(n);			// compute support for each item in candi_sets[n-1] by bruce force search
							// ɨ�����ݿ⣬ ����ÿ��itemset�����ݿ� һ���� ���ֵĴ���(Ƶ��)
		
		//output(candi_sets[n-1],n,0);					// ���--������
									
		gen_large(candi_sets[n-1]);						// ����threshold�жϣ��ҳ�large set 
		
		cout << "Output large "<< n << " ..." << endl;
		//output(reallarge_sets[n-1],n,1);				// ���--������
			
		output_LARGEset(large_sets[n-1], n, 1);			//��¼n-item large set
		output_LARGEset(prelarge_sets[n-1], n, 2);		//��¼n-item prelarge set
		output_LARGEset(reallarge_sets[n-1], n, 3);		//��¼n-item reallarge set

	}
*/

	//���ˣ� �Ѳ��������е�large item sets
	// At present, we have found all large item sets of all sizes

}





/*=======================================================================================================
				���´������йش�frequent item sets ���� association rules
=========================================================================================================*/


/*An association rule is valid if its ��confidence, support and lift�� are greater than or equal to corresponding threshold values.*/
/* ���´��� ��Ƶ���� ���� ��ѡrules  */
void apriori_gen_rules()
{
	cout<<"Begin to generate association rules!"<<endl;
	cout<<"��ʼ������������"<<endl;

	item condition_part, consequence_part;
	for(unsigned int n= 0+1 ; n<large_sets.size(); n++)			//�ӡ�2-Ƶ�����ʼ�Ƶ� ��������
	{
		map<pair<item,item>,vector<int>>  rule_n;

		for (itemsetCI tempCI=large_sets[n].begin(); tempCI!=large_sets[n].end(); tempCI++)
		{
			item goalitem = tempCI->first; 

			itemCI goalitemCI = goalitem.begin();
			for(; goalitemCI != goalitem.end(); goalitemCI++ )
			{
				int con_value = *goalitemCI;

				/*����1*/
				consequence_part.insert( con_value );	//rule�Ľ������
				goalitem.erase(goalitemCI);				//rule����������

				int condition_sup = large_sets[n-1][goalitem];  //condition  ֧�ֶ�
				int consequence_sup = large_sets[0][consequence_part];  // consequence ֧�ֶ�
				int union_sup = (tempCI->second) ;

				vector<int> sup_vec;
				sup_vec.push_back(condition_sup);
				sup_vec.push_back(consequence_sup);
				sup_vec.push_back(union_sup);

				rule_n.insert(  ruleset_insert_type(pair<item,item>(goalitem,consequence_part), sup_vec)  ); //����ѡ�������rule_n
				//rule_n.insert(  pair<pair<item,item>,vector<int>>(pair<item,item>(goalitem,consequence_part), sup_vec)  );  //д������Ҳ�ԣ���

				goalitemCI = (goalitem.insert(con_value)).first;  //��С�ġ� �˴�insert()�ķ���ֵ��pair<iterator,bool>, firstΪָ���²���Ԫ�صĵ�����;
																  //		 itemΪset<int>���ͣ� �����Զ����� ��˲����con_value����ԭ����λ��
				consequence_part.clear();

				/*����2*/
			
			}
		}

		candi_rulesets.push_back(rule_n);			// candi_rulesets��ȫ�ֱ���������vector<  map<pair<item,item>,vector<int>>  >

		rule_n.clear();
	}

	cout<<"������������ ����";

	output_candi_rules();

	output_real_rules();

	print_candi_rules();

	//getchar();
}


//�ļ���� ���к�ѡ����
void output_candi_rules()
{
	ofstream outfile;

	char filename[80];

	sprintf(filename, "%s_MST_%.3f_MCT_%.2f_MLT_%.2f_rules.dat ",datafrom, MST,MCT,MLT);
	outfile.open(filename, ios::out);

	outfile<<"\n---------------------------------------------------------------------------\n";
	outfile <<  "\tOutput " << "- Association Rules Candidates" << " ...";
	outfile<<"\n---------------------------------------------------------------------------\n";
	
	cout<<endl<<"��ʼ��д���ļ� "<<filename<< "  "<< " -- association rules"<< endl;

	for(unsigned n=0;n<candi_rulesets.size(); n++)
	{
		map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		for(; CI!=candi_rulesets[n].end(); CI++)
		{
			 item condition,consequence;     //item����set<int>����
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
		
	cout<<endl<<"������д���ļ� "<<filename<< "  "<< " -- association rules"<< endl;

}


//�ļ���� ����MST,MCT,MLT�Ĺ���
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
	cout<<endl<<"��ʼ��д���ļ� "<<filename<< "  "<< " -- association rules"<< endl;

	num_strong_rules = 0;		// ȫ�ֱ����� ��ʾstrong rules�������� 
								// ע�⣺candi_rulesets�����strong rules, ����pre_large itemsets���ɵĻ����Ŷ�С��MCT��rules

	for(unsigned n=0;n<candi_rulesets.size(); n++)
	{
		map<pair<item,item>,vector<int>>::const_iterator CI = candi_rulesets[n].begin();
		for(; CI!=candi_rulesets[n].end(); CI++)
		{
			 item condition,consequence;     //item����set<int>����
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
		
	cout<<endl<<"������д���ļ� "<<filename<< "  "<< " -- association rules"<< endl;
	
}


// ��Ļ��ӡ ��ѡ����
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

// ��Ļ��ӡ һ������
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