
#include "variator.h"
#include "variator_user.h"
#include "variator_internal.h"

#include "old_apriori.h"




void sort_sen_itemset();

void transfer_s_set_pair_2_s_set();

void filter_sup_trans( const item  & goalitem);	//过滤出敏感项的每个supporting trans的ID和len, 结果存入结构体向量

void filter_sup_trans2();						//只扫描一遍数据库； 找到每个sensitive rules的supporting trans集合,分别保存

void filter_sup_trans3();						//只扫描一遍数据库； 找到每个sensitive rules的supporting trans集合,合并保存

void display_sup_trans() ;

void display_sup_trans_NBRS() ;


void cal_fk_ratio();							//[DSS 2007] 计算每个snestive trans的fk ratio; fk=bk/ak


void display_rule(pair<item,item> a_rule);

void display_rule2(pair<item,item> a_rule);

bool is_belong_to_sensitive(set<int> itemset);

int select_item_to_remove(unsigned t_max_fk);


bool is_trans_contain_rule(transaction & goaltrans, pair<item,item>  a_rule);			//要求参数goaltrans已经采用内部编码

bool is_trans_only_contain_antecedent(transaction & goaltrans, pair<item,item>  a_rule);	//要求参数goaltrans已经采用内部编码

bool is_trans_adding_antecedent_contains_rule(transaction  origin_trans, item  senstive_rule_antecedent, pair<item,item> a_rule );



void hide_rules_process();

void write_DB_to_file();

void printf_file_fk();

void write_performance(int alpha, int beta, int gamma, int modified_trans_num, int num);


void write_stat_TKDE04();




