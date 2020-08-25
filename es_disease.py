import numpy as np
import pandas as pd
from collections import defaultdict

class ES:
    def __init__(self, gene_set_up, gene_set_down, gene_ranked_list):
        self.gene_set_up = gene_set_up
        self.gene_set_down = gene_set_down
        self.gene_ranked_list = gene_ranked_list

    def get_gene_rank_drug(self, gene_set):
        return np.array([self.gene_ranked_list.index(g) for g in gene_set])

    def normalize_gene_set(self, gene_set):
        n_set = len(gene_set)
        gene_set_rank = np.array(range(n_set)) + 1
        n_gene_set_rank = gene_set_rank/n_set
        return n_gene_set_rank

    def normalize_gene_rank_drug(self, gene_rank_drug):
        return gene_rank_drug/len(self.gene_ranked_list)

    def cal_a_b_ES(self, norm_gene_set, norm_drug_set):
        difference = norm_gene_set - norm_drug_set
        a_es = difference.max()
        difference_2 = -(norm_gene_set - norm_drug_set)+(1/len(norm_gene_set))
        b_es = difference_2.max()
        return a_es, b_es

    def cal_es(self, gene_set):
        gene_rank_drug = self.get_gene_rank_drug(gene_set)
        norm_gene_set = self.normalize_gene_set(gene_set)
        norm_drug_set = self.normalize_gene_rank_drug(gene_rank_drug)
        a_es, b_es = self.cal_a_b_ES(norm_gene_set, norm_drug_set)
        return a_es, b_es

    def cal_es_final(self):
        a_es_up, b_es_up = self.cal_es(self.gene_set_up)
        a_es_down, b_es_down = self.cal_es(self.gene_set_down)

        if a_es_up > b_es_up:
            es_up = a_es_up
        else:
            es_up = -b_es_up

        if a_es_down > b_es_down:
            es_down = a_es_down
        else:
            es_down = -b_es_down

        if es_up*es_down >= 0:
            es = 0
        else:
            es = es_up - es_down
        return es

def get_up_down_genesets(disease_genes_data):
    disease_genes = disease_genes_data[disease_genes_data['padj']<0.01]
    disease_genes = disease_genes.sort_values(by='log2FoldChange', ascending=False)
    gene_up = list(disease_genes[disease_genes['log2FoldChange']>0]['gene'])
    gene_up = [g.split('_')[0] for g in gene_up]
    gene_down = list(disease_genes[disease_genes['log2FoldChange'] < 0]['gene'])
    gene_down = [g.split('_')[0] for g in gene_down][::-1]
    return gene_up, gene_down

def get_drugs_es(gene_drug_left_pd, gene_up_left, gene_down_left ):
    result_es = defaultdict()
    for i in list(range(gene_drug_left_pd.shape[1])):
        rnk = gene_drug_left_pd.iloc[:, i].sort_values(ascending=False)
        gene_list = list(rnk.index)
        drug_name = list(gene_drug_left_pd.columns)[i]
        es_obj = ES(gene_up_left, gene_down_left, gene_list)
        es_score = es_obj.cal_es_final()
        result_es[drug_name] = es_score
    return result_es

def main():
    disease_genes_data = pd.read_csv('covid_genes.csv')  # (17899, 7)
    gene_disease_all = [g.split('_')[0] for g in list(disease_genes_data['gene'])]
    drug_rnk_s_all = pd.read_csv('drug_gene_pd_all.txt', sep='\t') #(12328, 754)
    gene_drug_all = [g for g in list(drug_rnk_s_all.columns)]
    drug_rnk_s_all = drug_rnk_s_all.set_index('drug_name')
    common_gene = [g for g in gene_disease_all if g in gene_drug_all] # 9823
    gene_drug_left_pd = drug_rnk_s_all[common_gene]
    gene_drug_left_pd = gene_drug_left_pd.T
    gene_name = list(gene_drug_left_pd.index)
    gene_up, gene_down = get_up_down_genesets(disease_genes_data) # 626, 626
    gene_up_left = [g for g in gene_up if g in gene_name] #450
    gene_down_left = [g for g in gene_down if g in gene_name] # 463

    # test by one drug
    rnk = gene_drug_left_pd.iloc[:, 3].sort_values(ascending=False)
    gene_list = list(rnk.index)
    es_obj = ES(gene_up_left, gene_down_left, gene_list)
    es_score = es_obj.cal_es_final()

    # get es for all drugs

    es_result = get_drugs_es(gene_drug_left_pd, gene_up_left, gene_down_left)

    es_pd = pd.DataFrame.from_dict(es_result, orient='index')
    es_pd.to_csv('new_es.txt', sep='\t')
