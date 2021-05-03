import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

#Eric's code
cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()


clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12

plt.figure(figsize=(8,6))
sns.countplot(data = clinical_data, x = 'PAM50', order = clinical_data['PAM50'].value_counts().index)
plt.title("Distribution of Subtypes in CPTAC-BRCA")
plt.ylabel("Number of Patients")
plt.xlabel("Subtype by PAM50")

clinical_data['Age_in_years'].value_counts(dropna=False)

plt.figure(figsize=(8,6))
sns.countplot(data = clinical_data.loc[clinical_data['Age_in_years'] < 40], x = 'PAM50', 
              order = clinical_data.loc[clinical_data['Age_in_years'] < 40]['PAM50'].value_counts().index,
             color = 'steelblue')
plt.title("Distribution of Subtypes for Young Patients (<40 years) in CPTAC-BRCA")
plt.ylabel("Number of Patients")
plt.xlabel("Subtype by PAM50")

plt.figure(figsize=(8,6))
sns.countplot(data = clinical_data.loc[clinical_data['Age_in_years'] >= 60], x = 'PAM50', 
              order = clinical_data.loc[clinical_data['Age_in_years'] >= 60]['PAM50'].value_counts().index,
             color = 'steelblue')
plt.title("Distribution of Subtypes for Old Patients (>=60 years) in CPTAC-BRCA")
plt.ylabel("Number of Patients")
plt.xlabel("Subtype by PAM50")



clinical_data.shape

young = clinical_data.loc[clinical_data['Age_in_years'] < 40]
old = clinical_data.loc[clinical_data['Age_in_years'] >= 60]

clinical_data.head()

protein_data.head()

protein_young = protein_data.loc[young.index]
protein_old = protein_data.loc[old.index]

protein_young.shape

protein_old.shape

protein_young.mean().sort_values(ascending=False)[0:5]

protein_old.mean().sort_values(ascending=False)[0:5]

protein_data_marked = protein_data.copy()
protein_data_marked['Group'] = np.nan
for i in young.index:
    protein_data_marked.at[i, 'Group'] = "Young"
for j in old.index:
    protein_data_marked.at[j, 'Group'] = "Old"

sns.boxplot(x=protein_data_marked["Group"], y=protein_data_marked["PIK3CA"])
plt.title("PIK3CA Expression in Young and Old Patients in CPTAC-BRCA")

sns.boxplot(x=protein_data_marked["Group"], y=protein_data_marked["TP53"])
plt.title("TP53 Expression in Young and Old Patients in CPTAC-BRCA")

protein_data_marked = protein_data_marked.sort_index()
clinical_data = clinical_data.sort_index()

protein_data_marked['PAM50'] = clinical_data['PAM50']


sns.boxplot(x=protein_data_marked["PAM50"], y=protein_data_marked["PIK3CA"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("PIK3CA Expression by Subtype (PAM50) in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")

sns.boxplot(x=protein_data_marked["PAM50"], y=protein_data_marked["TP53"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("TP53 Expression by Subtype (PAM50) in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")


#protein_data_marked.loc[protein_data_marked['Group'] == "Young"]

sns.boxplot(x=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["PAM50"], y=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["PIK3CA"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("PIK3CA Expression by Subtype (PAM50) for Young Patients in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")

sns.boxplot(x=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["PAM50"], y=protein_data_marked.loc[protein_data_marked['Group'] == "Young"]["TP53"], order=['LumA', 'Basal','LumB','Her2','Normal'])
plt.title("TP53 Expression by Subtype (PAM50) for Young Patients in CPTAC-BRCA")
plt.xlabel("Subtype by PAM50")

clinical_data.shape


from sklearn.manifold import TSNE
import umap

protein_data = protein_data.sort_index()

from sklearn.impute import SimpleImputer
values = protein_data.values
imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
imputed_values = imputer.fit_transform(values)

protein_data_imp = pd.DataFrame(imputed_values, columns = protein_data.columns, index= protein_data.index)

tsne = TSNE(n_components=2, random_state=7, perplexity=35, n_iter=1000, learning_rate=200)
tsne_res = tsne.fit_transform(protein_data_imp)

tsne_df = pd.DataFrame(data = tsne_res, columns = ['Dimension 1', 'Dimension 2'])
tsne_df = tsne_df.set_index(protein_data_imp.index)

final_tsne_df = pd.concat([tsne_df, clinical_data[['PAM50', 'Age_in_years']]], axis = 1)
final_tsne_df = pd.concat([final_tsne_df, protein_data_marked[["Group"]]], axis = 1)

curr_hue = 'PAM50'
plt.figure(figsize=(6,6))
sns.scatterplot(x="Dimension 1", y="Dimension 2", hue = curr_hue, data=final_tsne_df)
plt.title("t-SNE Clustering of Protein Expression Data (Colored by Subtype)")
