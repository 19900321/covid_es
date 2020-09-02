from cmapPy.pandasGEXpress.parse_gctx import parse
import multiprocessing as mp

import pandas as pd

# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz

# drug_ids = sig_info["sig_id"][sig_info["pert_iname"].isin(cmap_left["pert_iname"])]
# drug_ids = sig_info["sig_id"][sig_info["pert_iname"] == "vorinostat"]

def get_expression(gctx_path, per_dict, gene_dict, sig_info, cell_need, drug_name):
    drug_ids = list(sig_info[(sig_info["pert_iname"] == drug_name) & (sig_info["cell_id"].isin(cell_need))]["sig_id"])
    data = parse(gctx_path,
                 convert_neg_666=True,
                 cid=drug_ids)
    data_pd = data.data_df
    data_pd.index = [gene_dict[g] if g in gene_dict else None for g in data_pd.index]
    data_pd = data_pd.T
    data_pd['drug_name'] = [per_dict[g] if g in per_dict else None for g in data_pd.index]
    c = list(data_pd.columns)
    c.remove('drug_name')
    data_pd_mean = data_pd.groupby('drug_name')[c].mean()
    return data_pd_mean

def whole_drug_gene_mean(ges_id, drugs_pd, gctx_path):
    cmap = pd.read_csv('{}_Broad_LINCS_pert_info.txt'.format(ges_id), sep='\t')
    gene_info = pd.read_csv("{}_Broad_LINCS_gene_info.txt".format(ges_id), sep="\t", dtype=str)
    sig_info = pd.read_csv("{}_Broad_LINCS_sig_info.txt".format(ges_id), sep="\t")
    cmap_left = cmap[cmap['inchi_key'].isin(list(drugs_pd['inchikey']))]
    drugs = list(cmap_left["pert_iname"].unique())
    cell_info = pd.read_csv('{}_Broad_LINCS_cell_info.txt'.format(ges_id), sep='\t')
    cell_need = list(cell_info[cell_info['primary_site'] == 'lung']['cell_id'].unique())
    print(len(drugs))
    per_dict_in = dict(zip(sig_info["sig_id"], sig_info["pert_iname"]))
    gene_dict_in = dict(zip(gene_info['pr_gene_id'], gene_info['pr_gene_symbol']))
    l_pd = []
    n = 0
    for d_n in drugs:
        a = get_expression(gctx_path, per_dict_in, gene_dict_in, sig_info, cell_need, d_n)
        print(d_n)
        print(n)
        n += 1
        l_pd.append(a)
    drug_gene_mean = pd.concat(l_pd)
    return drug_gene_mean



ges_id = 'GSE70138'
drugs_pd = pd.read_csv('drugtargted_3000_Copy.txt', sep='\t')
gctx_path = 'GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx'
cmap = pd.read_csv('{}_Broad_LINCS_pert_info.txt'.format(ges_id), sep='\t')
sig_info = pd.read_csv("{}_Broad_LINCS_sig_info.txt".format(ges_id), sep="\t")
gene_info = pd.read_csv("{}_Broad_LINCS_gene_info.txt".format(ges_id), sep="\t",  dtype=str)
cmap_left = cmap[cmap['inchi_key'].isin(list(drugs_pd['inchikey']))]
cell_info = pd.read_csv('{}_Broad_LINCS_cell_info.txt'.format(ges_id), sep='\t')
cell_need = list(cell_info[cell_info['primary_site'] == 'lung']['cell_id'].unique())
drugs = list(cmap_left["pert_iname"].unique())
print(len(drugs))
per_dict_in = dict(zip(sig_info["sig_id"], sig_info["pert_iname"]))
gene_dict_in = dict(zip(gene_info['pr_gene_id'], gene_info['pr_gene_symbol']))


def get_expression_one(drug_name):
    return get_expression(gctx_path, per_dict_in, gene_dict_in, sig_info, cell_need, drug_name)

if __name__ == '__main__':
    pool = mp.Pool(8)

    funclist = []
    for d in drugs[5]:
        f = pool.apply_async(get_expression_one, [d])
        funclist.append(f)

    result = []
    for f in funclist:
        result_one = f.get()
        result.append(result_one)
    drug_gene_mean_GSE92742 = pd.concat(result)
    drug_gene_mean_GSE92742.to_csv('lung_drug_test.txt', sep='\t')

#gctx_path_GSE70138 = 'GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx'
#drug_gene_mean_GSE70138 = whole_drug_gene_mean('GSE70138', drugs_pd, gctx_path_GSE70138)

