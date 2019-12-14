import re
import GEOparse
import pandas
import numpy as np
import pickle

common_genes = [
"DDIT4",    
"GHRL",      
"PER1",      
"EPHX2",     
"GNG2",      
"IL1B",      
"DHRS13",    
"NR1D1",     
"ZNF438",    
"NR1D2",     
"CD38",      
"TIAM2",     
"CD1C",      
"LLGL2",   
"GZMB",      
"CLEC10A",  
"PDK1",      
"GPCPD1",    
"MUM1",      
"STIP1",     
"CHSY1",     
"AK5",       
"CYB561",    
"SLPI",      
"PARP2",     
"PGPEP1",    
"C12orf75",  
"FKBP4",     
"CAMKK1",  
"DTYMK", 
"NPEPL1",  
"MS4A3", 
"IL13RA1", 
"ID3", 
"MEGF6", 
"TCN1",  
"NSUN3", 
"POLH",  
"SYT11",   
"SH2D1B",  
"REM2"]    


gse0 = GEOparse.get_GEO("GSE39445")
#gse0.gpls : GPL15331


def isNaN(num):
    return num != num

id_gene0 = {} # gene symbol : {symbo: ID}
data = gse0.gpls['GPL15331'].table
for i in range(data.shape[0]):
  gene = data.GENE_SYMBOL[i]
  if gene in common_genes:
    id_gene0[data.ID[i]] = gene
  
#  
#id_gene2 = {} # gene symbol : {symbo: ID}
#data = gse2.gpls['GPL10379'].table
#for i in range(data.shape[0]):
#  gene = data.GENE_SYMBOL[i]
#  if gene in common_genes:
#    id_gene2[data.ID[i]] = gene
  

feat_target0 = {}
for i in common_genes:
  feat_target0[i] = {}

feat_target0['target'] = {}

for i, (gsm, value) in enumerate(gse0.phenotype_data[["title", "characteristics_ch1.4.timesampletaken", "characteristics_ch1.5.timesampletaken"]].iterrows()):
  feat_target0['target'][gsm] = 0
  if isNaN(value["characteristics_ch1.4.timesampletaken"]):
    feat_target0['target'][gsm] = value["characteristics_ch1.5.timesampletaken"]
  else:                        
    feat_target0['target'][gsm] = value["characteristics_ch1.4.timesampletaken"]


for name, gsm in gse0.gsms.items():
  print name, len(gse0.gsms)
  for i in range(gsm.table.shape[0]):
    if gsm.table.ID_REF[i] not in id_gene0:
      continue
    gene = id_gene0[gsm.table.ID_REF[i]]
    if gene in common_genes:
       feat_target0[gene][name] = gsm.table.VALUE[i]

f = open("GSE39445.pkl","wb")
pickle.dump(feat_target0,f)
f.close()

pd = pandas.DataFrame(feat_target0)
pd.to_csv('GSE39445.csv')

#f = open("GSE39445.pkl",'rb')
#feat_target0 = pickle.load(f)
#f.close()

print "sample number" , len(feat_target0['target'])
data_augmentation_times = 1001 / len(feat_target0['target']) + 1
print "To > 1000, sample after augmentation " , len(feat_target0['target']) * data_augmentation_times 


batch = {}
for gene in common_genes:
  batch[gene] = feat_target0[gene]

batch['target'] = {}
for gsm, v in feat_target0['target'].items():
  v = v.split(":")
  batch['target'][gsm] = float(v[0]) + float(v[1])/60.0


for i in range(1, data_augmentation_times):
  for gsm in feat_target0['target'].keys():
    batch['target'][gsm+ "_" + str(i)] = batch['target'][gsm]+ 0.01 * i
    for gene in common_genes:    
      v = batch[gene][gsm] + 0.01 * i
      batch[gene][gsm+ "_" + str(i)] = v
  
pd= pandas.DataFrame(batch)
pd.to_csv('GSE39445_automl.csv')
  


#--------------------------- dataset 2--------------------

gse1 = GEOparse.get_GEO("GSE48113")
#gse1.gpls : GPL15331

id_gene1 = {} # gene symbol : {symbo: ID}
data = gse1.gpls['GPL15331'].table
for i in range(data.shape[0]):
  gene = data.GENE_SYMBOL[i]
  if gene in common_genes:
    id_gene1[data.ID[i]] = gene

feat_target1 = {}
for i in common_genes:
  feat_target1[i] = {}

feat_target1['target'] = {}


# GSM1168837  NA ...
for i, (gsm, value) in enumerate(gse1.phenotype_data[["characteristics_ch1.3.time sample taken"]].iterrows()):
  if gsm == 'GSM1168837':
    continue
  feat_target1['target'][gsm] = value["characteristics_ch1.3.time sample taken"]
  for gene in common_genes:
     feat_target1[gene][gsm] = 0


  print gsm, feat_target1['target'][gsm] 

cnt = 0
for name, gsm in gse1.gsms.items():
  if gsm == 'GSM1168837':
    continue
  print name, cnt
  cnt = cnt +1
  for i in range(gsm.table.shape[0]):
    if gsm.table.ID_REF[i] not in id_gene1:
      continue
    gene = id_gene1[gsm.table.ID_REF[i]]
    if gene in common_genes:
       feat_target1[gene][name] = gsm.table.VALUE[i]

f = open("GSE48113.pkl","wb")
pickle.dump(feat_target1,f)
f.close()

pd = pandas.DataFrame(feat_target1)
pd.to_csv('GSE48113.csv')

#f = open("GSE48113.pkl",'rb')
#feat_target1 = pickle.load(f)
#f.close()

batch = {}
batch['target'] = {}
for gene in common_genes:
  batch[gene] = feat_target1[gene]

for gsm, v in feat_target1['target'].items():
  v = v.split(":")
  batch['target'][gsm] = float(v[0]) + float(v[1])/60.0

   
print "sample number" , len(feat_target1['target'])
data_augmentation_times = 1001 / len(feat_target1['target']) + 1
print "To > 1000, sample after augmentation " , len(feat_target1['target']) * data_augmentation_times 

for i in range(1, data_augmentation_times):
  for gsm in feat_target1['target'].keys():
    batch['target'][gsm+ "_" + str(i)] = batch['target'][gsm]+ 0.01 * i
    for gene in common_genes:    
      v = batch[gene][gsm] + 0.01 * i
      batch[gene][gsm+ "_" + str(i)] = v
  

pd= pandas.DataFrame(batch)
pd.to_csv('GSE48113_automl.csv')

#--------------------------- dataset 3--------------------

gse2 = GEOparse.get_GEO("GSE56931")
#gse2.gpls : GPL10379

id_gene2 = {} # gene symbol : {symbo: ID}
data = gse2.gpls['GPL10379'].table
for i in range(data.shape[0]):
  gene = data.GeneSymbol[i]
  if gene in common_genes:
    id_gene2[data.ID[i]] = gene

feat_target2 = {}
for i in common_genes:
  feat_target2[i] = {}

feat_target2['target'] = {}


for i, (gsm, value) in enumerate(gse2.phenotype_data[["characteristics_ch1.4.hour"]].iterrows()):
  feat_target2['target'][gsm] = value["characteristics_ch1.4.hour"]
  print gsm, feat_target2['target'][gsm] 

cnt = 0
for name, gsm in gse2.gsms.items():
  print name, cnt
  cnt = cnt +1
  for i in range(gsm.table.shape[0]):
    if gsm.table.ID_REF[i] not in id_gene2:
      continue
    gene = id_gene2[gsm.table.ID_REF[i]]
    if gene in common_genes:
       feat_target2[gene][name] = gsm.table.VALUE[i]

f = open("GSE56931.pkl","wb")
pickle.dump(feat_target2,f)
f.close()

pd = pandas.DataFrame(feat_target2)
pd.to_csv('GSE56931.csv')

f = open("GSE56931.pkl",'rb')
feat_target2 = pickle.load(f)
f.close()

batch = {}
batch['target'] = {}
for gene in common_genes:
  batch[gene] = feat_target2[gene]

for gsm, v in feat_target2['target'].items():
  print gsm, v
  batch['target'][gsm] = float(v)
   
print "sample number" , len(feat_target2['target'])
data_augmentation_times = 1001 / len(feat_target2['target']) + 1
print "To > 1000, sample after augmentation " , len(feat_target2['target']) * data_augmentation_times 

for i in range(1, data_augmentation_times):
  for gsm in feat_target2['target'].keys():
    batch['target'][gsm+ "_" + str(i)] = batch['target'][gsm]+ 0.01 * i
    for gene in common_genes:    
      v = batch[gene][gsm] + 0.01 * i
      batch[gene][gsm+ "_" + str(i)] = v
  

pd= pandas.DataFrame(batch)
pd.to_csv('GSE56931_automl.csv')



#--------------------------- dataset 4--------------------

gse3 = GEOparse.get_GEO("GSE113883")
#gse3.gpls : GPL18573

id_gene3 = {} # gene symbol : {symbo: ID}
data = gse3.gpls['GPL18573'].table
for i in range(data.shape[0]):
  gene = data.GeneSymbol[i]
  if gene in common_genes:
    id_gene3[data.ID[i]] = gene

feat_target3 = {}
for i in common_genes:
  feat_target3[i] = {}

feat_target3['target'] = {}


for i, (gsm, value) in enumerate(gse3.phenotype_data[["characteristics_ch1.4.hour"]].iterrows()):
  feat_target3['target'][gsm] = value["characteristics_ch1.4.hour"]
  print gsm, feat_target3['target'][gsm] 

cnt = 0
for name, gsm in gse3.gsms.items():
  print name, cnt
  cnt = cnt +1
  for i in range(gsm.table.shape[0]):
    if gsm.table.ID_REF[i] not in id_gene3:
      continue
    gene = id_gene3[gsm.table.ID_REF[i]]
    if gene in common_genes:
       feat_target3[gene][name] = gsm.table.VALUE[i]

f = open("GSE113883.pkl","wb")
pickle.dump(feat_target3,f)
f.close()

pd = pandas.DataFrame(feat_target3)
pd.to_csv('GSE113883.csv')

f = open("GSE113883.pkl",'rb')
feat_target3 = pickle.load(f)
f.close()

batch = {}
batch['target'] = {}
for gene in common_genes:
  batch[gene] = feat_target3[gene]

for gsm, v in feat_target3['target'].items():
  print gsm, v
  batch['target'][gsm] = float(v)
   
print "sample number" , len(feat_target3['target'])
data_augmentation_times = 1001 / len(feat_target3['target']) + 1
print "To > 1000, sample after augmentation " , len(feat_target3['target']) * data_augmentation_times 

for i in range(1, data_augmentation_times):
  for gsm in feat_target3['target'].keys():
    batch['target'][gsm+ "_" + str(i)] = batch['target'][gsm]+ 0.01 * i
    for gene in common_genes:    
      v = batch[gene][gsm] + 0.01 * i
      batch[gene][gsm+ "_" + str(i)] = v
  

pd= pandas.DataFrame(batch)
pd.to_csv('GSE113883_automl.csv')



