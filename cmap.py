import gseapy as gp
import pandas as pd
import multiprocessing as mp
import os
import numpy as np


# note: multiprocessing may not work on windows
def es_score(rnk, gene_sets_path):
    pre_res = gp.prerank(rnk=rnk, gene_sets=gene_sets_path,
                         processes=4,
                         permutation_num=10,
                         outdir=None, seed=6)
    return pre_res

def es_score_2(i, com_name, gene_sets_path):
    print(i)
    rnk = rnk_s.iloc[:, i].sort_values(ascending=False)
    pre_res = es_score(rnk, com_name[i], gene_sets_path)
    es_list = pre_res.res2d[['es']]
    es_list = es_list.rename(columns={'es': com_name[i]})
    return es_list

def es_score_3(i, gene_sets_path):
    print(i)
    rnk = rnk_s.iloc[:, i].sort_values(ascending=False)
    pre_res = es_score(rnk, gene_sets_path)
    result = pre_res.res2d
    if result.loc['covid_19_up', 'es']*result.loc['covid_19_down', 'es'] < 0:
        es = result.loc['covid_19_up', 'es'] - result.loc['covid_19_down', 'es']
    else:
        es = 0
    return com_name[i], es

def es_permu_n_cmap(i):
    gene_sets_path = 'gene_set_covid.gmt'
    com_name, es_score = es_score_3(i, gene_sets_path)
    es_permu_result = []
    for m in list(range(n)):
        gene_sets_path_r = 'gene_set_covid_{}.gmt'.format(m)
        _, random_es = es_score_3(i, gene_sets_path_r)
        es_permu_result.append(random_es)
    es_mean = np.mean(es_permu_result)
    if es_score > es_mean:
        es_p_value = len([es for es in es_permu_result if es > es_score])/n
    else:
        es_p_value = len([es for es in es_permu_result if es < es_score])/n
    return [com_name, es_score, es_p_value]

# data
def es_score_all(rnk_s, gene_sets_path):
    com_name = list(rnk_s.columns)
    result_es = []
    for i in list(range(rnk_s.shape[1])):
        rnk = rnk_s.iloc[:, i].sort_values(ascending=False)
        pre_res = es_score(rnk, com_name[i], gene_sets_path)
        es_list = pre_res.res2d[['es']]
        es_list = es_list.rename(columns={'es':com_name[i]})
        result_es.append(es_list)

    result_es_all = pd.concat(result_es,
                              axis=1,
                              join='outer',
                              ignore_index=False)
    return result_es_all


def merge_files():
    folder = os.listdir('es_result/')
    folder_left = ['es_result/'+f for f in folder if f.endswith('_ES_score_lung')]
    files = [f +'/'+f_2
             for f in folder_left
             for f_2 in os.listdir(f)
             if f_2.endswith('.report.csv')]

    result = []
    for f in files:
        name = f.split('/')[1].split('_ES_score')[0]
        result_one = pd.read_csv(f)
        result_one['drug_name'] = name
        result_one['term_id'] = result_one.index
        result.append(result_one)

    result_es_all = pd.concat(result,
                              axis=0,
                              ignore_index=True)

    return result_es_all

rnk_s = pd.read_csv("drug_gene_lung_pd_all.txt", sep="\t")
#gene_sets_path = 'gene_set_covid.gmt'
rnk_s = rnk_s.set_index('drug_name')
rnk_s = rnk_s.T
com_name = list(rnk_s.columns)
#result_es_all = es_score_all(rnk_s, gene_sets_path)
n = 100
if __name__ == '__main__':
    pool = mp.Pool(8)
    result_es = []
    funclist = []
    for i in list(range(rnk_s.shape[1])):
        f = pool.apply_async(es_permu_n_cmap, [i])
        funclist.append(f)

    result = []
    for f in funclist:
        result_one = f.get()
        result.append(result_one)

    result_es_all = pd.DataFrame(result)
    result_es_all.to_csv('es_cmap.txt', sep='\t')