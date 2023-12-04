

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
	map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong; //��NEW��

	set< pair<item,item> > collection_common_rule_IL;		

	//map< pair<item,item>, A_MST_MCT > collection_common_rule_NBRS;
	map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_NBRS;	//��NEW��



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


	vector<unsigned long>			copy_support_of_items_one2;			//��������copy_support_of_items_one2[]
	vector< vector<unsigned long> > copy_temp_counter_array2;			//��������copy_temp_counter_array2[][]

	unsigned long N0, N1;						// store the current values for N0 and N1
	vector<int> N0_vec, N1_vec;		// store the the number of blocking 0's and blocking 1's for each sensitive rule

	double SM;			//	Safety margin
	double A01;			//	The proportion of trans which change 0's into ?.  0 < A01 < 1

						//	A01/(1-A01) == N0/N1

//-------------------------------------------


//===========================================
double *db_fk;   //db_fk����ͬdb_size(main��������ռ�),ʵ��ֻ����ÿ��supp trans��fkֵ�������trans ID(��Ӧ�����±�)ֱ�Ӵ�ȡ��fk
//===========================================


set<set<int>> s_set_update;					// ������ϣ� ΪDSS2007�㷨׼�������㷨������itemsetһ��һ������,ֱ��Ϊ�ա�

set<set<int>> deleted_non_sensitive_all;	// ȫ�ֱ����� ��������б����ص�non-sensitivie itemset


//�������Ľ����������ڴ����ݽṹ�� 
//	1) Ϊÿ��sensitive trans ����������֧�ֵ�����non-sensitive frequent itemset list
//  2) ҲΪÿ��sensitive trans ����������֧�ֵ�����  sensitive frequent itemset list
//���޸Ļ�ɾ����transʱ�� ��Ҫ���µ�ǰtrans��֧�ֵ�sensitive itemset��non-sensitive itemset��support, 
//��Ϊ�д����ݽṹ���Ͳ���Ҫ���¼���Ѱ�ҵ�ǰtrans֧�ֵ�itemset; �Կռ任ʱ�䣬��������ٶ�

typedef struct NODE
{
	unsigned long trans_ID;

	set< pair<item,item> >  supp_all_rule_set;
	set< pair<item,item> >  supp_only_antecedent_rule_set;

	//set< pair<item,item> >  supp_strong_rule_set;
	//set< pair<item,item> >  supp_NBRS_rule_set;

	double        weight;								
	//set<set<int>> sup_itemsets;						//��ű�trans��֧�ֵ�����non-sensitive itemsets
	//set<set<int>> sup_itemsets_sen;					//��ű�trans��֧�ֵ�����sensitive itemsets
}T_LEN;



vector<T_LEN> supporting_trans_vec;							// ����������supporting transaction�����ID�Ͷ�Ӧ����

vector<T_LEN> supporting_trans_vec_NBRS;					// ���������С�����֧��antecedent�Ҳ���ȫ֧��consequent����transaction�����ID�Ͷ�Ӧ����

vector<unsigned long> all_supp_trans_ID;


//-------------| �����������ʱ�� |----------------
//-------------------------------------------------
clock_t start_clock, start_clock2, end_clock;
double total_run_time, total_run_time_after_apriori;
//-------------------------------------------------



bool sort_by_fk(const T_LEN& n1, const T_LEN& n2)
{
	return (n1.weight < n2.weight);						//Ŀǰ�ǡ��������򣺡�<�����¡����򡿣� ��>�����¡�����
	//����ӦΪ<
}

bool sort_by_fk_NBRS(const T_LEN& n1, const T_LEN& n2)
{
	return (n1.weight > n2.weight);						//Ŀǰ�ǡ��������򣺡�<�����¡����򡿣� ��>�����¡�����
	//return (n1.weight < n2.weight);
	//����ӦΪ>
}


vector<vector<T_LEN>>  supporting_trans_vec2;		// supporting_trans_vec2[][]�Ƕ�ά�����������ֱ𱣴�ÿ��sensitive rule��supporting trans ID�Ͷ�Ӧtrans length


unsigned long * support_sen_rule_array;				// ��̬���������飬 ÿ��Ԫ�ر�ʾ��Ӧ���й���ĳ���Ƶ�Σ�֧�ֶȣ�


//���TKDE04-2.b�㷨, ��sensitive rules ����size��support����
typedef struct NODE2
{
	unsigned long sen_rules_ID;
	unsigned long itemset_size;
	unsigned long itemset_sup;						
}ITEMSET_NODE;


vector<ITEMSET_NODE> sen_itemset_sort_vec;

//qsort()�е�compare()�ȽϺ������ء�����������			�÷���qsort(s,100,sizeof(s[0]),cmp);
//sort()�е�compare()�ȽϺ�������1��0�������桱�򡰼١���	�÷���sort(a,a+5,cmp);    


bool sort_by_itemset_size_sup(const ITEMSET_NODE& n1, const ITEMSET_NODE& n2)     //��<�����¡����򡿣� ��>�����¡�����
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
	start_clock = clock();			// ����ʼ��ʱ��
	cout<<"\n��1�μ�ʱ��ʼ������"<<endl;

	strcpy_s(paramfile,FILE_NAME_LENGTH, "BA_param.txt");  //��ע�⡿�����ｫ ��������ļ����̻��ˣ���ʵ������argv[1]

    read_local_parameters();

	apriori_bodon();			// ����apriori_bodon()����õ����С���������д���ļ�(��������bodon_ruleset)�� �ٴ���ѡ�����й���
								// ��ȡ��1-item���͡�2-item(����)��Ƶ�������ݽṹ����ȡ������main_trie_copy����ȡ������DB

	sup_u = (int)(db_size*MST);

	//initial();		//��ע�⡿�� �����Ѿ��õ�������"���ݿ�thin_DB���Σ� get_thin_DB(thin_DB, datafrom);
					//			 initial()�����õ������������ڴ桱���ݿ�database, ��������õ�����
					//          ��TKDE04��1.a,1.b,2.a,2.b��Ҫ����"ԭʼ���ݿ�"ÿ��trans��С����ҲΪ�����޸ĺ�����ݿ��������ļ���׼��


	cout<<"\n\n\nPlease specify sensitive rules in the file.  When finished press [Enter] to continue! \n ";
	getchar();

	read_sensitive_item();		// ���������ļ�����ȡrules�� vector<pair<item,item>>���͵�s_set_pair �� vector<item>���͵�s_set��

	start_clock2 = clock();		// ����ʱ���� ���ļ��� ����������� ��һ�ο�ʼ��ʱ������
	cout<<"\n��2�μ�ʱ��ʼ[��ȡ���������]������"<<endl;
	//-----------------------------------------------

	//sort_sen_itemset();			//����sensitive itemset��size��support��s_set_pair����������.ע�⣺һ���������򣬺���supporting trans

	//���ܣ� s_set==>s_set_update
	//ע�⣺ s_set_update��set����, s_set��vector���ͣ� s_set_updateʹ�õ�������
	for(vector<set<int>>::iterator it=s_set.begin(); it!=s_set.end(); it++)
	{
		item union_item, union_item_encoding;
		union_item = *it;								 
		basket_recode3(union_item, union_item_encoding);
		s_set_update.insert(union_item_encoding);      
	}

	copy_support_of_items_one2 = copy_support_of_items_one;		//��������copy_support_of_items_one2[]
	copy_temp_counter_array2 = copy_temp_counter_array;	//��������copy_temp_counter_array2[][]


	//db_fk = new double[db_size];						//db_fk����ͬdb_size,ʵ��ֻ����ÿ��supp trans��fkֵ�������trans ID(��Ӧ�����±�)ֱ�Ӵ�ȡ��fk
	//for(unsigned int i=0;i<db_size;i++)  db_fk[i] = 0.0;	

	//filter_sup_trans3();			// ֻɨ��һ�����ݿ⣻ �ҵ�����sensitive rules��supporting trans����,���浽all_supp_trans_ID
	//cal_fk_ratio();				// �õ���ÿ��trans��Ӧ��fkֵ

	hide_rules_process();


	getchar();
	return 0;
}



//����sensitive itemset��size��support��s_set_pair���������У����㷨2.b����ɫ֮һ��
void sort_sen_itemset()
{
	cout<<"\nBEGIN to sort sensitive itemset!\n";

	for(unsigned int i=0;i<s_set_pair.size();i++)				
	{
		item union_item;
		item antecedent, consequent, antecedent_encoding, consequent_encoding;

		antecedent = s_set_pair[i].first;
		consequent = s_set_pair[i].second;

		basket_recode3(antecedent, antecedent_encoding);							//תΪ����
		basket_recode3(consequent, consequent_encoding);			

		set<int>::const_iterator it;
		for( it= s_set_pair[i].first.begin(); it!= s_set_pair[i].first.end(); it++)
		{
			union_item.insert(*it);
		}
		for( it= s_set_pair[i].second.begin(); it!= s_set_pair[i].second.end(); it++)
		{
			union_item.insert(*it);			//�ϲ�antecedent��consequent
		}	

		//vector<int> rule_3_attribute = bodon_ruleset[s_set_pair[i]];				// ������д�����������Ҳ���rule����Ϊbodon_rulesetʹ������
		vector<int> rule_3_attribute = bodon_ruleset[pair<item,item>(antecedent_encoding, consequent_encoding)];
		
		if(rule_3_attribute.empty()) {cout<<"Error--sort_sen_itemset()! bodon_ruleset���Ҳ��������й���\n"; getchar(); exit(-1);}

		unsigned long union_sup = rule_3_attribute[2];								// ����rule��֧�ֶ�
		unsigned long antecedent_sup = rule_3_attribute[0];							// rule���antecedent��֧�ֶ�

		ITEMSET_NODE tmp;
		tmp.sen_rules_ID = i;
		tmp.itemset_size = union_item.size();
		tmp.itemset_sup = union_sup;

		sen_itemset_sort_vec.push_back(tmp);
	}

	// ����sort()�⺯����C++֧�֣���qsort()�⺯��(C֧��) ����
	// "����"���򣺺���sort_by_len()��<���������� ��>��������
	sort(sen_itemset_sort_vec.begin(), sen_itemset_sort_vec.end(), sort_by_itemset_size_sup);	


	//��֤����
	for(unsigned int i=0;i<s_set_pair.size();i++)				
	{
		cout<<"  size: "<<sen_itemset_sort_vec[i].itemset_size;
		cout<<"  sup:  "<<sen_itemset_sort_vec[i].itemset_sup;
		cout<<"  ID:   "<<sen_itemset_sort_vec[i].sen_rules_ID<<endl;
	}

	//���°���s_set_pair��rule�Ĵ��˳��
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

	//��֤����size��union_sup���°��š����˳�򡱵�s_set_pair[]
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
	}	

}


//����t������sensitive itemset����; ע���β�goaltrans���յ���thin_DB�е�transaction
unsigned long cal_bk(vector<int> & goaltrans,  set<set<int>> &itemset_collection)
{
	unsigned num_bk = 0;

	for(set<set<int>>::iterator it=s_set_update.begin(); it!=s_set_update.end(); it++)				
	{
			item union_item, union_item_encoding;

			/*
			union_item = *it;								//��ע�⡿�������õ���s_set,�䱣�����ÿ��rules��antecedent��consequent�ϲ�itemset
			basket_recode3(union_item, union_item_encoding);
			*/

			union_item_encoding = *it;

			transactionCI transCI = goaltrans.begin();	
			itemCI itemCI = union_item_encoding.begin();
			int count = 0;

			while (itemCI!=union_item_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )								//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
				if( *transCI  == *itemCI )														//������ʽһ���� ����������ʹ�õ����ڲ���
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



//���t����������itemset������
unsigned long cal_ak(vector<int> & goaltrans)
{
	unsigned long num_ak = 0;

	//��ʾ�� map<item, int>  bodon_itemset;
	//       typedef map<item, int>::value_type   bodon_itemset_insert_type; 

	for(map<item,int>::iterator it=bodon_itemset.begin(); it!=bodon_itemset.end(); it++)
	{
		set<int> itemset_encode;
		itemset_encode	= it->first;  //bodon_ruleset�б����������itemset

		transactionCI transCI = goaltrans.begin();	
		itemCI itemCI = itemset_encode.begin();
		int count = 0;

		while (itemCI!=itemset_encode.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
		{
			//if ( copy_new_code_inverse[*transCI] == *itemCI )								//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
			if( *transCI  == *itemCI )														//������ʽһ���� ����������ʹ�õ����ڲ���
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

	//©���� �ҵ�t����������Ƶ��itemset�� û���ų�����itemset

	return num_ak;
}


//ͨ��trie�����t����������itemset������
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

			//if(tmp_it == itemset_collection.end() ) {cout<<"cal_ak_by_trie()��ɾ���������һ��item"<<endl; break;}

			it = tmp_it;

			//cout<<"�ų�sensitive itemset ����"<<endl;
		}
		else
		{
			it++;
		}
	}


	//ע��itemset_collection���itemsetʹ�õ�������
	return itemset_collection.size();
}


//����ÿ��snestive trans��fk ratio; fk=bk/(ak+1)
void cal_fk_ratio()
{
	cout<<"\nBEGIN to calculate fk ratio for each supporting trans!\n";

	for(unsigned i = 0; i<all_supp_trans_ID.size(); i++)
	{
		unsigned long ID = all_supp_trans_ID[i] ;

		vector<int> t =  thin_DB[ ID ];

		//unsigned long ak = cal_ak(t);			//���t����������itemset������ [ʱ�临�Ӷȸ�]
		set<set<int>> itemset_collection;
		unsigned long ak = cal_ak_by_trie(t, itemset_collection);	//ͨ��trie�����t����������non-sensitive itemset������ [ʱ�临�Ӷȵ�]
		//supporting_trans_vec[i].sup_itemsets = itemset_collection;	//ע��itemset_collection���itemsetʹ�õ�������
		
		itemset_collection.clear();
		unsigned long bk = cal_bk(t, itemset_collection);			//���t������sensitive itemset������
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

	for(unsigned long i=0;i<thin_DB.size(); i++)				//ʹ�á��������ݿ�(thin_DB��ÿ�о�����������)
	{
		transaction goaltrans = thin_DB[i];					 		

		for(unsigned int j=0;j<s_set_pair.size();j++)				
		{
			item union_item, union_item_encoding;

			union_item = s_set[j];								//��ע�⡿�������õ���s_set,�䱣�����ÿ��rules��antecedent��consequent�ϲ�itemset
			basket_recode3(union_item, union_item_encoding);

			transactionCI transCI = goaltrans.begin();	
			itemCI itemCI = union_item_encoding.begin();
			int count = 0;

			while (itemCI!=union_item_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )								//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
				if( *transCI  == *itemCI )														//������ʽһ���� ����������ʹ�õ����ڲ���
				{
					count++;
					itemCI++;
					transCI++;
				}
				else { transCI++;}
			}
			if (count == union_item_encoding.size() && union_item_encoding.size()!=0)
			{
				all_supp_trans_ID.push_back(i); //iΪ��ǰtrans ID

				break;		//ֻҪ��ǰtrans����1��sensitive itemset,�Ͳ����ж����Ƿ񻹰�������sensitive itemset
			}
		}
	}

	cout<<"\nEND to filter out supporting transactions for all sensitive rule [ver 3]!\n";

}


//ֻɨ��һ�����ݿ⣻ �ҵ�ÿ��sensitive rules��supporting trans����,�ֱ𱣴�
void filter_sup_trans2()
{
	cout<<"\nBEGIN to filter out supporting transactions for each sensitive rule!\n";

	for(unsigned int j=0;j<s_set_pair.size();j++)
	{
		vector<T_LEN> sup_vec;
		supporting_trans_vec2.push_back(sup_vec);	
	}


	for(unsigned long i=0;i<thin_DB.size(); i++)				//ʹ�á��������ݿ�(thin_DB��ÿ�о�����������)
	{
		transaction goaltrans = thin_DB[i];					 		

		for(unsigned int j=0;j<s_set_pair.size();j++)				
		{
			item union_item, union_item_encoding;

			union_item = s_set[j];								//��ע�⡿�������õ���s_set,�䱣�����ÿ��rules��antecedent��consequent�ϲ�itemset
			basket_recode3(union_item, union_item_encoding);

			transactionCI transCI = goaltrans.begin();	
			itemCI itemCI = union_item_encoding.begin();
			int count = 0;

			while (itemCI!=union_item_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
			{
				//if ( copy_new_code_inverse[*transCI] == *itemCI )								//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
				if( *transCI  == *itemCI )														//������ʽһ���� ����������ʹ�õ����ڲ���
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
				//tmp.trans_len = goaltrans.size();   //��ע�⡿:����1.a,1.b,2.a,2.b�㷨�ı��⣬���ﲻ������DB��trans���ȣ�����ԭʼDB��trans����
				//tmp.trans_len = database[i].size();   //         ����ʹ����ԭ���ݿ���trans�ĳ��ȣ� ���ǡ��������ݿ⡱

				supporting_trans_vec2[j].push_back( tmp );					//supporting_trans_vec2[][]�Ƕ�ά�����������ֱ𱣴�ÿ��sensitive rule��supporting tras ID�Ͷ�Ӧtrans length
			}

		}

	}

	/*
	// �ĵ�hide_rules_process()��������

	//����trans_len��С����
	for(unsigned int j=0;j<s_set_pair.size();j++)				
	{
		// ����sort()�⺯����C++֧�֣���qsort()�⺯��(C֧��) ����
		// "����"���򣺺���sort_by_len()��<���������� ��>��������
		sort(supporting_trans_vec2[j].begin(), supporting_trans_vec2[j].end(), sort_by_len);		

	}

	//��֤
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
			union_item.insert(*it);			//�ϲ�antecedent��consequent
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


void remove_t(unsigned long k, unsigned long selected_ID)		//k�ǵ�ǰ ���й�����ţ� selected_ID��Ҫɾ����trans_ID
{

	if( k<0 || k>=s_set_pair.size())  {cout<<"Error in remove_t! kԽ��"; getchar(); exit(-1);}
	  
	vector<T_LEN>::const_iterator it = lower_bound(	supporting_trans_vec2[k].begin(), supporting_trans_vec2[k].end(), selected_ID , TNODE_point_less);
	  
	  //ע�⣺ʹ��lower_bound()ǰ������������
	  //lower_bound()����һ�� iterator, ��ָ����[first,last)��ǵ����������п��Բ���value���������ƻ�����˳��ĵ�һ��λ�ã�
	  //�����λ�ñ����һ�����ڵ���value ��ֵ������lower_bound֮ǰ����ȷ������Ϊ�������У�������ó���
	  //�ο��� ����̡��ļ���--��STL--��lower_bound.docx

	if( it != supporting_trans_vec2[k].end() &&  (*it).trans_ID == selected_ID )
	{
		//�ҵ��ˣ���supporting_trans_vec2[k]������ɾ��
		;
		vector<T_LEN>::const_iterator tmp_it;
		tmp_it = supporting_trans_vec2[k].erase(it);			//erase���ر�ɾ��Ԫ�ص���һ��Ԫ�ص�ָ��(iterator)	
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
		
	//set< set<int> >  collection_common_itemset;		������ȫ�֡�
	//set< pair<item,item> > collection_common_rule;	������ȫ�֡�
	set<int>::iterator it = IR.begin();
	for(; it!= IR.end(); it++)							//Actually, there is only one item in IR. 
	{													//It is enough to blocking one item for blocing 1s to ? for a candidate transaction.
		int cur_item = *it;

		//TYPE�� map<pair<item,item>,vector<int>> bodon_ruleset_strong, ʹ������
		//TYPE�� map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MIN;ʹ������

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

							//map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong; //��NEW��
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


							//map< pair<item,item>, BORDER_RULE_INFO > collection_common_rule_strong; //��NEW��
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



void	find_R_common_NBRS( set<int> &IL )  //��������Ҫ�޸ģ���
{
		
	cout<<"\n\n[BEGIN] void	find_R_common_NBRS( set<int> &IL )	\n";
		
	//set< set<int> >  collection_common_itemset;			������ȫ�֡�
	//set< pair<item,item> > collection_common_rule_NBRS;	������ȫ�֡�
	set<int>::iterator it = IL.begin();
	//for(; it!= IL.end(); it++)
	{
		int cur_item = *it;

		//���ͣ� map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_strong_MIN; ʹ������
		//���ͣ� map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MAX;ʹ������

		map<pair<item,item>, T_RULE_UNKNOWN>::iterator  it_bodon = bodon_ruleset_all_MAX.begin();
		for(; it_bodon != bodon_ruleset_all_MAX.end(); it_bodon++)
		{						

			item union_encoding, antecedent_encoding, consequent_encoding;						//[NOTE]: for blocking 1s to ?, only one item needs to be removed
																								//for bolcking 0s to ?, all empty item 0s in the sensitive rule need to be added
			antecedent_encoding = it_bodon->first.first;					
			consequent_encoding = it_bodon->first.second;
			combine_item(union_encoding, antecedent_encoding, consequent_encoding);

			//if( is_belong_to(cur_item, union_encoding) )										//1st screening��share the common item
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
						  																		//3rd screening�����N0��?���б��ghost rules�Ŀ���
						{
							//A_MST_MCT tmp;
							BORDER_RULE_INFO tmp;

							tmp.supp_original			= MAX_supp;
							tmp.supp_update				= MAX_supp;  
							tmp.ante_supp_update		= MAX_antecedent;
							tmp.supp_in_candidate_trans = 0;

							tmp.A_MST = (int)num_supp_below_MST;   tmp.A_MCT = (int)num_conf_below_MCT;	//��ע�⡿��ȡ�������Ϊ0��
							 																			//�����cal_P_for_each_trans_NBRS�����ж�ӦֵӦ+1

							if( num_supp_below_MST <= num_conf_below_MCT)
								tmp.N = (unsigned)num_conf_below_MCT;								// 'N' denotes how many blocking 0s to '?' operations at minimum  
							else																	//  to cause the border rule lost
								tmp.N = (unsigned)num_supp_below_MST;

							tmp.weight = (double)(tmp.supp_update - tmp.supp_original + 1)/(tmp.N + 1);
							tmp.weight = func_weight(tmp.weight);									//adjust the weight change rate through math function


							collection_common_rule_NBRS[pair<item,item>(antecedent_encoding,consequent_encoding)] = tmp; //�������������һ��Negative Border rule
						
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


void filter_sup_trans( const item  & goalitem_encoding)   //���˳��������ÿ��supporting trans��ID��len, �������ṹ������
{
	cout<<"\n\n[BEGIN]: void filter_sup_trans( const item  goalitem) \n";
	int count = 0;
	//item goalitem_encoding;

	map< pair<item,item>, BORDER_RULE_INFO >::iterator	ruleCI;		// Used to traverse "collection_common_rule_strong"
	map< set<int>,BORDER_ITEMSET_INFO >::iterator		itemsetCI;	// 

	//basket_recode3(goalitem, goalitem_encoding);					//��goalitem��Ϊʹ���ڲ����ʾ��trie�����������ݿⶼʹ���ڲ��룩

	//for(unsigned int i=0;i<database.size(); i++)				
	for(unsigned long i=0;i<thin_DB.size(); i++)					//ʹ�á��������ݿ�(thin_DB��ÿ�о�����������)
	{
		//transaction goaltrans = database[i];
		transaction goaltrans = thin_DB[i];							//ʹ�á��������ݿ�			

		itemCI itemCI = goalitem_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
		transactionCI transCI = goaltrans.begin();	

		while (itemCI!=goalitem_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
		{
			//if ( copy_new_code_inverse[*transCI] == *itemCI )	//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
			if( *transCI  == *itemCI )   //������ʽһ���� ����������ʹ�õ����ڲ���
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
			//tmp.trans_ID = i;		//iΪtrans�����ݿ��еļ�¼���
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

			//��NOTE��:
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



//�����ݿ��й��˳� ���в���֧��condition�Ҳ���ȫ֧��consequent��transactions

void filter_sup_trans_NBRS( set<int> & antecedent_encoding, set<int> & consequent_encoding )   //���˳��������ÿ��supporting trans��ID��len, �������ṹ������
{
		cout<<"\n\n[BEGIN]: void filter_sup_trans_NBRS( ) \n";

		map< pair<item,item>, BORDER_RULE_INFO >::iterator	ruleCI;	
	
		for(unsigned int i=0;i<thin_DB.size(); i++)										//ʹ�á��������ݿ�(thin_DB��ÿ�о�����������)
		{
				transaction & goaltrans = thin_DB[i];														 
		
				item::iterator CI1 = antecedent_encoding.begin();											//itemΪset<int>, �����Ԫ���Զ���������
				vector<int>::const_iterator transCI = goaltrans.begin();	

				int count1= 0;
				int flag1 = 0;
				while ( CI1 != antecedent_encoding.end() && transCI!=goaltrans.end())	//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
				{
					//if ( copy_new_code_inverse[*transCI] == *itemCI )					//���ܻ������ע�⡿��copy_new_code_inverse[*transCI]�����õ�����Ʒ��Ų�һ�����������еģ�
					if( abs(*transCI)  == *CI1 )										//������ʽһ���� ����������ʹ�õ����ڲ���
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

				//��鵱ǰtrans�Ƿ�ƥ��consequent_encoding, ��־flag2==1��ʾƥ��
				int count2 = 0;
				int flag2 = 0;

				item::iterator  CI2 = consequent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
				transCI = goaltrans.begin();	

				while ( CI2 != consequent_encoding.end() && transCI!=goaltrans.end())					//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
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

				if(0==flag1 && 0==flag2)															// ����ȫ֧��antecedent, ͬʱҲ����ȫ֧��consequent
				{
					//T_LEN tmp;
					//tmp.trans_ID = i;		//iΪtrans�����ݿ��еļ�¼���
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



bool is_trans_contain_rule(transaction & goaltrans, pair<item,item>   a_rule)		//Ҫ�����goaltrans�Ѿ������ڲ�����
{
	item union_encoding;
	int count = 0;

	combine_item(union_encoding, a_rule.first, a_rule.second);

	itemCI itemCI = union_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
	transactionCI transCI = goaltrans.begin();	

	while (itemCI!=union_encoding.end() && transCI!=goaltrans.end())			//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *transCI  == *itemCI )												//������ʽһ���� ����������ʹ�õ����ڲ���
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



int is_trans_contain_rule_L_or_R_del(transaction &goaltrans, pair<item,item> a_rule, int removed_item)		//Ҫ�����goaltrans�Ѿ������ڲ�����
{
	item antecedent_encoding = a_rule.first;

	itemCI itemCI1			= antecedent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
	transactionCI transCI	= goaltrans.begin();	

	int count1 = 0;
	int flag1 = 0;
	while (itemCI1!=antecedent_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *transCI  == *itemCI1 )														//������ʽһ���� ����������ʹ�õ����ڲ���
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

	itemCI itemCI2	= consequent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
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
		{	cout<<"\nFunction--is_trans_contain_rule_L_or_R(): removed item is not in antecedent or consequent.";  //��˵��������ɾitem����common item in IR,IR���ٰ�������item
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



int is_trans_contain_rule_L_or_R_add(transaction &goaltrans, pair<item,item> a_rule, vector<int> &added_items)		//Ҫ�����goaltrans�Ѿ������ڲ�����
{
	item antecedent_encoding = a_rule.first;

	itemCI itemCI1 = antecedent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
	transactionCI transCI = goaltrans.begin();	

	int count1 = 0;
	int flag1 = 0;
	while (itemCI1!=antecedent_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *transCI  == *itemCI1 )														//������ʽһ���� ����������ʹ�õ����ڲ���
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

	itemCI itemCI2 = consequent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
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
		{	cout<<"\nFunction--is_trans_contain_rule_L_or_R(): added item is not in antecedent or consequent.";  //��˵����������ӵ�item����common item in IL,IL���ٰ�������item
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



bool is_trans_only_contain_antecedent(transaction & goaltrans, pair<item,item> a_rule)		//Ҫ�����goaltrans�Ѿ������ڲ�����
{
	item antecedent_encoding = a_rule.first;

	itemCI itemCI1 = antecedent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
	transactionCI transCI = goaltrans.begin();	

	int count1 = 0;
	int flag1 = 0;
	while (itemCI1!=antecedent_encoding.end() && transCI!=goaltrans.end())				//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *transCI  == *itemCI1 )														//������ʽһ���� ����������ʹ�õ����ڲ���
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

	itemCI itemCI2 = consequent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
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

	itemCI				itemCI1		= antecedent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
	set<int>::iterator	combineCI	= combine.begin();	

	int	count1	= 0;
	int	flag1	= 0;

	while (itemCI1!=antecedent_encoding.end() && combineCI!=combine.end())					//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *combineCI  == *itemCI1 )														//������ʽһ���� ����������ʹ�õ����ڲ���
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

	itemCI itemCI2	= consequent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
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
	while (itemCI!=union_encoding.end() && transCI_set!=goaltrans_set.end())			//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *transCI_set  == *itemCI )													//������ʽһ���� ����������ʹ�õ����ڲ���
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


bool is_trans_only_contain_antecedent_NBRS(transaction goaltrans, pair<item,item> a_rule,  set<int> add_antecedent_encoding)		//Ҫ�����goaltrans�Ѿ������ڲ�����
{

	set<int> goaltrans_set;
	for(transactionCI it = goaltrans.begin(); it!=goaltrans.end(); it++)
		goaltrans_set.insert(*it);

	for(set<int>::iterator it2 = add_antecedent_encoding.begin(); it2!=add_antecedent_encoding.end(); it2++)
		goaltrans_set.insert(*it2);

	//------------------------------------------------
	item antecedent_encoding = a_rule.first;

	itemCI it1 = antecedent_encoding.begin();										//itemΪset<int>, �����Ԫ���Զ���������
	set<int>::iterator trans_set_CI = goaltrans_set.begin();	

	int count1 = 0;
	int flag1 = 0;
	while ( it1!=antecedent_encoding.end() && trans_set_CI!=goaltrans_set.end())			//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *trans_set_CI  == *it1 )												//������ʽһ���� ����������ʹ�õ����ڲ���
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
	while ( it2!=consequent_encoding.end() && trans_set_CI!=goaltrans_set.end())			//�����ȶԵ�ǰ���ǣ�goalitem �� goaltrans�Ѿ��Զ��ź���															
	{
		if( *trans_set_CI  == *it2 )														//������ʽһ���� ����������ʹ�õ����ڲ���
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




void cal_P_for_each_trans()							//ÿ��trans��Ȩ��P��Ӧfk
{
	//��ʾ�� vector<T_LEN> supporting_trans_vec;	
	cout<<"\n\n[BEGIN] void cal_P_for_each_trans()	\n";

	cout<<"  supporting_trans_vec.size() = "<<supporting_trans_vec.size()<<endl;
	cout<<"  collection_common_rule_strong.size() = "<<collection_common_rule_strong.size()<<endl;
	cout<<"  please wait...\n";


	//map< pair<item,item>, A_MST_MCT >  supp_rule_set;
	map< pair<item,item>, BORDER_RULE_INFO >  supp_rule_set;			//��NEW��


	vector<T_LEN>::iterator it = supporting_trans_vec.begin();
	for(; it != supporting_trans_vec.end(); it++)
	{
		double A_MCT_max = 0, A_MCT_min = (double)db_size;
		double trans_weight = 0;

		unsigned long DB_ID = it->trans_ID;
		transaction &goaltrans = thin_DB[DB_ID];

		//��ʾ��set< pair<item,item> > collection_common_rule;						//collection_common_rule�����뵱ǰsen_rule����common item������rules
		//set< pair<item,item> >::iterator it2 = collection_common_rule_strong.begin();
		//map< pair<item,item>, A_MST_MCT >::iterator it2 = collection_common_rule_strong.begin();

		map< pair<item,item>, BORDER_RULE_INFO >::iterator it2 = collection_common_rule_strong.begin();		//��NEW��
		

		for(; it2 != collection_common_rule_strong.end(); it2++)
		{
			if( is_trans_contain_rule(goaltrans, it2->first) )
			{
				//it->supp_all_rule_set.insert(*it2);									//store all rules��Strong or Not Strong�� which might be influenced by modification.
				//if( bodon_ruleset_all_MIN[*it2].union_supp/(double)db_size >= MST && bodon_ruleset_all_MIN[*it2].conf >= MCT)
						supp_rule_set.insert(*it2);										//"supp_rule_set" only store related strong rules
			}
			//else if( is_trans_only_contain_antecedent(goaltrans, *it2)  )
			//	it->supp_only_antecedent_rule_set.insert(*it2);

		}

		//cout<<" supp_rule_set.size() == "<<supp_rule_set.size()<<endl;
		//it->supp_strong_rule_set = supp_rule_set;			//���浱ǰtrans��ȫ֧�ֵ�strong rules

		for(it2 = supp_rule_set.begin(); it2!=supp_rule_set.end(); it2++)
		{
			unsigned num_above_MST = it2->second.A_MST;
			unsigned num_above_MCT = it2->second.A_MCT;
																					//��1��Ϊ�˷�ֹA_MCT_R == 0;
			double A_MCT_R;
			if( num_above_MST<num_above_MCT )
			     A_MCT_R = (double)num_above_MST + 1;	
			else
				 A_MCT_R = (double)num_above_MCT + 1;	

			if (A_MCT_R > A_MCT_max)	A_MCT_max = A_MCT_R;
			if (A_MCT_R < A_MCT_min)	A_MCT_min = A_MCT_R;

			trans_weight += 1.0/A_MCT_R; 
		}

		//���뿼��supp_rule_set,size()Ϊ0�����
		if(supp_rule_set.size()==0)									//˵����ǰtrans�������κ�����3����˵�rule, ���������������С��
		{
			//cout<<"\n\nERROR! ���ٵ�ǰҪ���ص�sen ruleӦ��supp_rule_set��"; getchar(); exit(-1);
			trans_weight = 0;
		}

		//trans_weight = trans_weight* A_MCT_max/A_MCT_min;

		it->weight = trans_weight;				// trans_weightΪ��ǰtrans��Ȩ�أ� ������fk, fk������ƴ���������(DSS07)��ֲ������ʵ�ʱ�ʾ����P
	
		//it->weight = (double)(supp_rule_set.size());

		supp_rule_set.clear();

		//cout<<".";
	}

	cout<<"\n[END] void cal_P_for_each_trans()	\n";
}




void cal_P_for_each_trans_NBRS( set<int> & antecedent_encoding)										//Ϊÿ��NBRS��supporting trans����Ȩ��P��  P��Ӧfk
{
	//��ʾ�� vector<T_LEN> supporting_trans_vec_NBRS;	

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

		//��ʾ��set< pair<item,item> > collection_common_rule_NBRS;							//����ǰ��collection_common_rule_NBRS�����뵱ǰsen_rule��antecedent����common item������Negative Boder Rules Set
		//set< pair<item,item> >::iterator it2 = collection_common_rule_NBRS.begin();			//�����ڡ�collection_common_rule_NBRS�����뵱ǰsen_rule��antecedent����common item�ġ�����rules��
		//map< pair<item,item>, A_MST_MCT >::iterator it2 = collection_common_rule_NBRS.begin();
		map< pair<item,item>, BORDER_RULE_INFO >::iterator it2 = collection_common_rule_NBRS.begin();

		
		for(; it2 != collection_common_rule_NBRS.end(); it2++)
		{
			//if(is_trans_contain_rule_NBRS(goaltrans, *it2, antecedent_encoding))			//��ע�⡿: is_trans_contain_rule_NBRS()������Ҫ�Ȳ���IL�����ж��Ƿ����rule
			if(is_trans_contain_rule_NBRS(goaltrans, it2->first, antecedent_encoding))	
			{
				//it->supp_all_rule_set.insert(*it2);
				//if(  ( bodon_ruleset_all_MIN[*it2].conf > MCT  &&  bodon_ruleset_all_MAX[*it2].union_supp/(double)db_size < MST )  ||
				//	 ( bodon_ruleset_all_MAX[*it2].conf < MCT  &&  bodon_ruleset_all_MIN[*it2].union_supp/(double)db_size > MST )  
				//  )
						supp_rule_set_NBRS.insert(*it2);											//supp_rule_set������ ��collection_common_ruleѡ���ı�trans֧�ֵ�rules
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
			//	B_MCT_R =  (MCT - MAX_conf)*MAX_antecedent_sup  + 1;				//confС��MCT    ��1��Ϊ�˷�ֹA_MCT_R == 0;
			//}
			//else if( MIN_conf >= MCT && (double)MAX_union_sup/db_size < MST  )
			//{
			//	B_MCT_R =  MST * db_size - (double)MAX_union_sup  + 1;				//suppС��MST
			//}
			//else
			//{
			//	cout<<"\n Error in cal_P_for_each_trans_NBRS():  NBRS rule��support��conf��������ֵ��";	getchar(); exit(-1);
			//}

			double B_MCT_R;
			unsigned num_below_MST = it2->second.A_MST;
			unsigned num_below_MCT = it2->second.A_MCT;

			if( num_below_MST >  num_below_MCT )
				B_MCT_R = (double)(  num_below_MST ) + 1 ;				//��1��Ϊ�˷�ֹA_MCT_R == 0;
			else
				B_MCT_R = (double)(  num_below_MCT ) + 1;

			if (B_MCT_R > B_MCT_max)	B_MCT_max = B_MCT_R;
			if (B_MCT_R < B_MCT_min)	B_MCT_min = B_MCT_R;

			trans_weight += 1.0/B_MCT_R; 
		}

		if(supp_rule_set_NBRS.size() == 0)				//Ӧ�ÿ���supp_rule_set_NBRS.size()==0�����������ǰtrans����item�󲻰����κ�Negative-border rule�����
			trans_weight = 0;
		else
		{
			//trans_weight = trans_weight* B_MCT_max/B_MCT_min;
			;
		}

		it->weight = trans_weight;				// trans_weightΪ��ǰtrans��Ȩ�أ� ������fk, fk������ƴ���������(DSS07)��ֲ������ʵ�ʱ�ʾ����P

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
	set<int> trans_set;										//�����衿:������?ʱ����ر�֤�²���ĸ�ֵӦ���� ��Ӧ��ֵ������Ӧ��ŵ�λ�ã������ű�֤��itemƥ��ʱ��������

	added_items_vec.clear();

	vector<int>::iterator it_vec = trans_vec.begin();
	for(; it_vec!=trans_vec.end() ; it_vec++)				//����¼trans��vector<int>ת��set<int>
	{
		trans_set.insert(*it_vec);						
	}

	item::iterator l_it = antecedent.begin();				
	for(; l_it != antecedent.end(); l_it++ )
	{
		if( trans_set.end() == trans_set.find(*l_it) )			//find()��������ֵ== end()��ʾ �ڵ�ǰtrans���Ҳ���*l_it
		{
			added_items_vec.push_back(*l_it);					//�ҳ�����ӵ�items, ����added_items_set
			trans_set.insert(*l_it);	//ʵ�ʸ���trans			//���汾1������ʵ�ʸ���trans�����ϸ�Ҫ��sensitive rules֮����comon items��
																//�����á������汾2����ʵ�ʸ���trans��
		}
	}

	sort(added_items_vec.begin(), added_items_vec.end());

	trans_vec.clear();
	set<int>::iterator it_set = trans_set.begin();				//���µ�trans��set<int>����vector<int>
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
	
	unsigned  union_sup = bodon_ruleset_all_MIN[ cur_sen_rule ].union_supp;	// ����rule��֧�ֶ�
	unsigned  antecedent_sup = bodon_ruleset_all_MIN[ cur_sen_rule ].antecedent_supp; // rule���antecedent��֧�ֶ�	

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

	set<int> set_modified_trans_ID;						//���������� ͳ�Ʊ��޸�trans����

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
		basket_recode3(consequent, consequent_encoding);				// �����������͡����ۡ�ת��Ϊ�����ʾ����Ϊbodon_rulesetʹ������
		basket_recode3(union_item, union_encoding);

		compute_N0_N1( antecedent_encoding, consequent_encoding  );


		//----------------------------------------------------------------------------------------------


		int selected_item ;

		item::const_iterator it2 = consequent_encoding.begin();			//ע�⣺ʹ�á����롱�汾��consequent
		selected_item = *it2;

		for(it2++; it2!= consequent_encoding.end(); it2++)				//ѡ֧�ֶ���ߵ�item�� �˷����뵱ǰtrans�޹�
		{
			if( copy_support_of_items_one[*it2] > copy_support_of_items_one[selected_item] )	
				selected_item = *it2;   
		}
		set<int> temp_IR;
		temp_IR.insert(selected_item);
		//-------------- ��consequent��ѡ��һ��itemɾ����ע�⣺�˷�ѡ��Ĵ�ɾitem���뵱ǰtrans��ء�----

		//find_R_common(antecedent_encoding);	�����ﷸ�����ش��󣡣�����������Ӧ����consequent�Ұ���common item��rules��
		find_R_common(temp_IR);
		find_R_common_NBRS( antecedent_encoding );


		filter_sup_trans( union_encoding );													// ɨ�����ݿ�(scan DB)�ҳ���ǰsen_rule��supporting trans
		filter_sup_trans_NBRS(antecedent_encoding, consequent_encoding );


		cal_P_for_each_trans();															//Ϊÿ��supp_trans����P(Ȩ��)ֵ����Ӧfk
		cal_P_for_each_trans_NBRS(antecedent_encoding);									//Ϊÿ��NBRS_supp_trans����P(Ȩ��)ֵ����Ӧfk


		collection_common_rule_strong.clear();											//����ա�collection_common_rule�����뵱ǰsen_rule����common item������rules
		collection_common_rule_NBRS.clear();


		//����sort()�⺯����C++֧�֣���qsort()�⺯��(C֧��) ����; ��Ȩ�ض�supp. trans����
		//sort����ͷ�ļ�#include <algorithm>    ---  qsort����ͷ�ļ�#include<stdlib.h>

		cout<<"\n\n[BEGIN]: sort!";
		sort(supporting_trans_vec.begin(), supporting_trans_vec.end(), sort_by_fk);		//���������򣺺���sort_by_len()��<���������� ��>��������												
		sort(supporting_trans_vec_NBRS.begin(), supporting_trans_vec_NBRS.end(), sort_by_fk_NBRS);
		cout<<"\n\n[END]: sort!\n";

		//display_sup_trans();															//��Ļ��ʾ����֤������
		//display_sup_trans_NBRS();


		//[NEW!!!]
		//[NEW!!!] make changes from sorting to finding the minimum in each iteration of database modification
		//[NEW!!!] weights of border rules and priority of related transactions will be DYNAMICALLY UPDATED
		//-----------------------------------------------------------------------------------------------------
				
		//��Blocking N1 1'��


		cout<<"\n\n ��Blocking N1 1'��  N1 == "<<N1<<endl;
		cout<<"\n   -- collection_common_rule_IR.size() == "<<collection_common_rule_IR.size()<<endl;
		for(unsigned j=0; j<N1; j++)
		{
				//-----------------------------------------------------------------------------
				//���ڴ����ݿ�databas����ʵɾ��selected_item��Ӧ��
				//int flag = 0;
				//vector<int>::const_iterator it3 = database[ supporting_trans_vec[j].trans_ID ].begin();
				//while( it3 != database[ supporting_trans_vec[j].trans_ID ].end() )
				//{
				//	if( copy_new_code[*it3] -1  == selected_item  )    // ��*it3ת��Ϊ����
				//	{
				//		vector<int>::const_iterator tmp_it;
				//		tmp_it = database[ supporting_trans_vec[j].trans_ID ].erase(it3);  //erase���ر�ɾ��Ԫ�ص���һ��Ԫ�ص�ָ��(iterator)
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
						*it3 = -selected_item;					//��ע�⡿�����ｫitem����ȡ������ɶԳƵĸ�������ʾ��item�ڵ�ǰtrans��ΪUnknown

						flag = 1;
						break;
					}
					else
						it3++;
				}
				//--------------------------------------------------------------------------------
			
				if(flag)	//���������־�жϿɷ�ֹ�ظ�ɾ��ͬһtrans��ͬһitem, ��ֹ������ظ�����
				{
					modified_item_num++;
					set_modified_trans_ID.insert(DB_ID);				//���������� ͳ�Ʊ��޸�trans����
					//--------------------------------------------

					cout<<".";
					//cout<<"\nj=="<<j<<"  collection_common_rule_IR.size() == "<<collection_common_rule_IR.size()<<endl;

					set< pair<item,item> >::iterator it4 = collection_common_rule_IR.begin();
					for(; it4 != collection_common_rule_IR.end(); it4++)
					{
						//���ͣ� map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MIN;
						//T_RULE_UNKNOWN a_record;

						pair<item,item> related_rule = *it4;
						T_RULE_UNKNOWN & a_record = bodon_ruleset_all_MIN[*it4];

						int result = is_trans_contain_rule_L_or_R_del(t_copy, *it4, selected_item);

						if(  1 == result  )
						{
							a_record.antecedent_supp -= 1;
							a_record.union_supp -= 1;
							a_record.conf = (double)(a_record.union_supp)/a_record.antecedent_supp; 

							//cout<<"\n֮��:\n";
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
					cout<<"�����ظ�ɾ��item!  Sensitive rules֮����common items!"<<endl;		//getchar();
				}
		}

		//------------------------------------------------------------------------------------------
		//��Block N0 0'��

		cout<<"\n\n ��Block N0 0'��  N0 == "<<N0<<endl;
		cout<<"\n   -- collection_common_rule_IL.size() == "<<collection_common_rule_IL.size()<<endl;

		unsigned N0_iteration;

		if( supporting_trans_vec_NBRS.size() < N0 )  //������õ�partial support trans������ 
		{
			cout<<"   \n��Warning��:The partially supporting trans is not sufficient!!\n";
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

				if( !added_items_vec_encode.empty())										// �ɷ�ֹ�ظ�����
				{

					set_modified_trans_ID.insert(DB_ID);

					modified_item_num = modified_item_num + added_items_vec_encode.size();						//update data side effects
					//-------------------------------------------------------------------------
					cout<<".";
					//cout<<"\nj=="<<j<<"  collection_common_rule_IL.size() == "<<collection_common_rule_IL.size()<<endl;

					set<pair<item,item>>::iterator it5 = collection_common_rule_IL.begin();						//update knowledge side effects
					for( ; it5 != collection_common_rule_IL.end(); it5++  )
					{
						//���ͣ� map<pair<item,item>, T_RULE_UNKNOWN> bodon_ruleset_all_MAX;
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

		supporting_trans_vec.clear();							// ����ա�supporting_trans_vec��ȫ�ֽṹ������
		vector<T_LEN>().swap(supporting_trans_vec);				// �����Ż������ͷ�����ռ�

		supporting_trans_vec_NBRS.clear();
		vector<T_LEN>().swap(supporting_trans_vec_NBRS);		//Free the memory

	}

	//----------------------------------------------------------------------------------------
	
	 write_DB_to_file();										//output transactions from the memory to files

	//---------------------- ͳ��alpha--------------------------------------------------------
	//���ͣ�vector<pair<item,item>> s_set_pair;	
	for(unsigned int j=0;j<s_set_pair.size();j++)		
	{
		item antecedent_encoding, consequent_encoding;

		basket_recode3(s_set_pair[j].first, antecedent_encoding);
		basket_recode3(s_set_pair[j].second, consequent_encoding);				// �����������͡����ۡ�ת��Ϊ�����ʾ����Ϊbodon_rulesetʹ������

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


	/*-------------------- ͳ��beta ----------------------*/

	map<pair<item,item>, T_RULE_UNKNOWN>::iterator it = bodon_ruleset_strong.begin();
	for(; it != bodon_ruleset_strong.end(); it++)
	{
		unsigned min_union_sup = bodon_ruleset_all_MIN[it->first].union_supp;
		double min_conf = bodon_ruleset_all_MIN[it->first].conf;

		if( min_union_sup/(double)db_size < MST || min_conf < MCT )
			beta ++;
	}

	beta = beta -  (s_set_pair.size()-alpha); 


	/*-------------------- ͳ��gamma ----------------------*/

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

	//basket_recode3(antecedent, antecedent_encoding);							//תΪ����
	//basket_recode3(consequent, consequent_encoding);			

	//for( it= a_rule.first.begin(); it!= a_rule.first.end(); it++)
	//{
	//		union_item.insert(*it);
	//}
	//for( it= a_rule.second.begin(); it!= a_rule.second.end(); it++)
	//{
	//		union_item.insert(*it);			//�ϲ�antecedent��consequent
	//}	

	////vector<int> rule_3_attribute = bodon_ruleset[s_set_pair[i]];				// ������д�����������Ҳ���rule����Ϊbodon_rulesetʹ������
	//vector<int> rule_3_attribute = bodon_ruleset[pair<item,item>(antecedent_encoding, consequent_encoding)];
	//	
	//if(rule_3_attribute.empty()) {cout<<"Error--display_rule()! bodon_ruleset���Ҳ��������й���\n"; getchar(); exit(-1);}

	//unsigned long union_sup = rule_3_attribute[2];								// ����rule��֧�ֶ�
	//unsigned long antecedent_sup = rule_3_attribute[0];							// rule���antecedent��֧�ֶ�

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

	outfile<<"�㷨BA 2007  ɾ��items \n\n";

	outfile<<"datafile:                 "<<datafrom<<endl;
	outfile<<"database size:            "<<db_size<<endl;
	outfile<<"Frequent Itemset(1) or Association Rule(2):"<<FIM_or_AR<<endl;;   //��FIM��ARѡ�񿪹ء������������
	outfile<<"Minimum support threshold (MST):      "<<MST<<endl;
	outfile<<"Minimum confidence threshold (MCT):   "<<MCT<<endl;				// �����������
	//outfile<<"Minimum lift threshold (MLT):         "<<MLT<<endl;				// �����������
	outfile<<"MST discount:             "<<MST_discount<<endl;
	outfile<<"MCT discount:             "<<MCT_discount<<endl;				

	outfile<<"Safety Margin (SM) :      "<<SM<<endl;
	outfile<<"The proportion of blocking 0's to 1's (A01):   "<<A01<<endl;


	outfile<<endl;
 
	outfile<<"Number of strong rules calculated by trie tree     =   "<<num_strong_rules_by_trie<<endl; 
	outfile<<"Number of frequent itemset calculated by trie tree =   "<<bodon_itemset_frequent.size()<<endl;
	outfile<<"Number of all itemset calculated by trie tree =   "<<bodon_itemset.size()<<endl;

	outfile<<"\n bodon_itemset.size() == "<<bodon_itemset.size();								//��2014-Aug-4 19:41 ������
	outfile<<"\n bodon_itemset_frequent.size() == "<<bodon_itemset_frequent.size();				//��2014-Sep-2 18:35 ������
	outfile<<"\n num_itemset_by_trie == "<< num_itemset_by_trie << endl;

	outfile<<"\n bodon_ruleset_all_MIN.size() == "<<bodon_ruleset_all_MIN.size()<<endl;		//bodon_ruleset����pre-strong rules;���ͳ�Ʊ�����num_strong_rules_by_trie
	outfile<<"\n bodon_ruleset_all_MAX.size() == "<<bodon_ruleset_all_MAX.size()<<endl;		//bodon_ruleset����pre-strong rules;���ͳ�Ʊ�����num_strong_rules_by_trie
	outfile<<"\n bodon_ruleset_strong.size()  == "<<bodon_ruleset_strong.size()<<endl;		//bodon_ruleset����pre-strong rules;���ͳ�Ʊ�����num_strong_rules_by_trie
	outfile<<"\n bodon_ruleset_NBRS.size()    == "<<bodon_ruleset_NBRS.size()<<endl;		//bodon_ruleset����pre-strong rules;���ͳ�Ʊ�����num_strong_rules_by_trie
				
	outfile<<" num_strong_rules_by_trie == "<<num_strong_rules_by_trie<<endl;

	outfile<<" \nbodon_ruleset_all����pre-strong rules;���ͳ�Ʊ�����num_strong_rules_by_trie\n";
	outfile<<" \n���ڸ�����һ����;��num_itemset_by_trieȡֵ��bodon_itemset_frequent.size()��2�� \n";



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

		basket_recode3(antecedent, antecedent_encoding);	//תΪ����
		basket_recode3(consequent, consequent_encoding);	

		pair<set<int>, set<int>>  a_rule = pair<set<int>, set<int>>(antecedent_encoding, consequent_encoding );
		
		//if(rule_3_attribute.empty()) {cout<<"Error in write_performance()! û����bodon_ruleset���ҵ������й���"<<endl; getchar(); exit(-1);}

		unsigned long union_sup = bodon_ruleset_strong[a_rule].union_supp;					// ����rule��֧�ֶ�
		double conf =  bodon_ruleset_strong[a_rule].conf;

		outfile<<"\n   union_sup = "<<union_sup/(double)db_size ;
		outfile<<"\n   conf      = "<<conf<<endl;


		outfile<<"\n\n   N0 = " << N0_vec[i] << "       N1 = "<< N1_vec[i] <<endl;

	}

	//-----------------------------------------------------------------------------
	end_clock = clock();		//��������ʱ��
	double total_run_time = (double)(end_clock - start_clock)/CLOCKS_PER_SEC;
	double total_run_time_after_apriori = (double)(end_clock - start_clock2)/CLOCKS_PER_SEC;

	outfile<<"total_run_time: \t"<<total_run_time<<endl;
	outfile<<"total_run_time_after_apriori"<<total_run_time_after_apriori<<endl;
	//------------------------------------------------------------------------------
								
	outfile.close();

	cout<<"\n>>> END to write the outcome of App_Intell14 into file!\n";

}


