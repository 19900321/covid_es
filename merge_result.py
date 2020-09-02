import pandas as pd
drug_gene_pd_92742 = pd.read_csv('drug_gene_GSE92742.txt', sep='\t') #(680, 12329
drug_gene_pd_70138 = pd.read_csv('70138_drug_gene_mean.txt', sep='\t') #(372, 12329)
drug_gene_pd_all = pd.concat([drug_gene_pd_92742, drug_gene_pd_70138])
use_col = list(drug_gene_pd_all.columns)
use_col.remove('drug_name')

drug_gene_pd_all = drug_gene_pd_all.groupby('drug_name')[use_col].mean()
drug_gene_pd_all.to_csv('drug_gene_pd_all.txt', sep='\t')


# select the targeted 20 drugs
drug_gene_pd_all = pd.read_csv('drug_gene_pd_all.txt', sep='\t')
drug_name_map = pd.read_csv('drug_new_maped.txt', sep='\t')

es_all = pd.read_csv('es_all.txt', sep = '\t')
drug_name_map = drug_name_map[drug_name_map['pert_iname'].notna()]
drug_targted = pd.read_csv('target_drug.csv')
drug_targted['drugs'] = drug_targted['drugs'].str.lower()
drugs_20 = list(drug_targted['drugs'])
drug_gene_pd_all['drug_name'] = drug_gene_pd_all['drug_name'].str.lower()

drug_name_map['name'] = drug_name_map['name'].str.lower()
drug_name_map['pert_iname'] = drug_name_map['pert_iname'].str.lower()
name_dict = dict(zip(drug_name_map['name'], drug_name_map['pert_iname']))

drug_gene_20 = drug_gene_pd_all[drug_gene_pd_all['drug_name'].isin(drugs_20)]

gene_pd = []
[name_dict[d] for d in drugs_20 if d in name_dict.keys()]
[d for d in drugs_20 if d in list(drug_gene_pd_all['drug_name'])]

len([d for d in drugs_20 if d in list(drug_name_map['name'])])
len([d for d in drugs_20 if d in list(es_all['drug_name'].unique())])

[d for d in drugs_20 if d in list(drug_name_map['name'])]

drug_query = pd.read_csv('cmap_query')
drug_query['name'] = drug_query['name'].str.lower()

drug_in = list(set(drug_name_map[drug_name_map['pert_id'].isin(list(drug_query['pert_id']))]['pert_iname']))

[d for d in drugs_20 if d in list(drug_query['name'])]
drug_gene_20 = drug_gene_pd_all[drug_gene_pd_all['drug_name'].isin(drug_in)]
drug_gene_20.to_csv('drug_15.txt', sep='\t')