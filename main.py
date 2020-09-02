import multiprocessing as mp
from ges import get_expression
import pandas as pd

ges_id = 'GSE70138'
drugs_pd = pd.read_csv('drugtargted_3000_Copy.txt', sep='\t')
gctx_path_GSE92742 = 'GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx'
cmap = pd.read_csv('{}_Broad_LINCS_pert_info.txt'.format(ges_id), sep='\t')
sig_info = pd.read_csv("{}_Broad_LINCS_sig_info.txt".format(ges_id), sep="\t")
gene_info = pd.read_csv("{}_Broad_LINCS_gene_info.txt".format(ges_id), sep="\t")
cmap_left = cmap[cmap['inchi_key'].isin(list(drugs_pd['inchikey']))]
drugs = list(cmap_left["pert_iname"].unique())
print(len(drugs))
per_dict_in = dict(zip(sig_info["sig_id"], sig_info["pert_iname"]))
gene_dict_in = dict(zip(gene_info['pr_gene_id'], gene_info['pr_gene_symbol']))

def get_expression_one(drug_name):
    return get_expression(gctx_path_GSE92742, per_dict_in, gene_dict_in, sig_info, drug_name)



if __name__ == '__main__':
    pool = mp.Pool(8)

    funclist = []
    for d in drugs[1:50]:
        f = pool.apply_async(get_expression_one, [d])
        funclist.append(f)

    result = []
    for f in funclist:
        result_one = f.get()
        result.append(result_one)
    drug_gene_mean_GSE92742 = pd.concat(result)
    drug_gene_mean_GSE92742.to_csv('example.txt', sep='\t')