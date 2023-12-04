
#include "variator.h"
#include "variator_user.h"
#include "variator_internal.h"

#include "old_apriori.h"




void sort_sen_itemset();

void transfer_s_set_pair_2_s_set();

void filter_sup_trans( const item  & goalitem);	//���˳��������ÿ��supporting trans��ID��len, �������ṹ������

void filter_sup_trans2();						//ֻɨ��һ�����ݿ⣻ �ҵ�ÿ��sensitive rules��supporting trans����,�ֱ𱣴�

void filter_sup_trans3();						//ֻɨ��һ�����ݿ⣻ �ҵ�ÿ��sensitive rules��supporting trans����,�ϲ�����

void display_sup_trans() ;

void display_sup_trans_NBRS() ;


void cal_fk_ratio();							//[DSS 2007] ����ÿ��snestive trans��fk ratio; fk=bk/ak


void display_rule(pair<item,item> a_rule);

void display_rule2(pair<item,item> a_rule);

bool is_belong_to_sensitive(set<int> itemset);

int select_item_to_remove(unsigned t_max_fk);


bool is_trans_contain_rule(transaction & goaltrans, pair<item,item>  a_rule);			//Ҫ�����goaltrans�Ѿ������ڲ�����

bool is_trans_only_contain_antecedent(transaction & goaltrans, pair<item,item>  a_rule);	//Ҫ�����goaltrans�Ѿ������ڲ�����

bool is_trans_adding_antecedent_contains_rule(transaction  origin_trans, item  senstive_rule_antecedent, pair<item,item> a_rule );



void hide_rules_process();

void write_DB_to_file();

void printf_file_fk();

void write_performance(int alpha, int beta, int gamma, int modified_trans_num, int num);


void write_stat_TKDE04();




