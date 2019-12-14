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



data = pandas.read_csv('/home/jingchi/gibraltar/second/deliver/GSE39445_GSE48113_GSE56931.csv').to_dict()

num_sample = len(data['sample'])
# 900 for training and the rest for evaluation
data_train = {}
data_test = {}
for key in data.keys():
  data_train[key] = {}
  data_test[key] = {}
  
for key in data.keys():
  for i in range(900):
    if key == 'target': 
      print data[key][i]
      v = str(data[key][i])
      print v
      v = v.split(":")
      print v
      if len(v) == 2:
        print v
        v = float(v[0]) + float(v[1])/60.0
      else:
        print v
        v = float(v[0])
      data_train[key][i] = v
    else:
      data_train[key][i] = data[key][i]

      
  for i in range(900, num_sample):
    data_test[key][i] = data[key][i]
    if key == 'target': 
      v = str(data[key][i])
      v = v.split(":")
      if len(v) == 2:
        v = float(v[0]) + float(v[1])/60.0
      else:
        v = float(v[0])
      data_test[key][i-900] = v
    else:
      data_test[key][i-900] = data[key][i]


for key in data.keys():
  for i in range(100):
    if key == 'sample': 
      data_train[key][i+900] = data_train[key][i] + "_1"
    else:
      print key, i, data[key][i]
      data_train[key][i+900] = data_train[key][i] + 0.01
    
  

pd= pandas.DataFrame(data_train)
pd.to_csv('GSE39445_GSE48113_GSE56931_train.csv')
pd= pandas.DataFrame(data_test)
pd.to_csv('GSE39445_GSE48113_GSE56931_test.csv')



