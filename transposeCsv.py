import pandas as pd
csv = pd.read_csv("test.csv", header=None)
# use skiprows if you want to skip headers
df_csv = pd.DataFrame(data=csv)
#df_csv.to_csv(index=True, header=True, decimal=',', sep=',',float_format='%.4f')
transposed_csv = df_csv.T
transposed_csv.to_csv('output.csv',decimal=',', sep=',', float_format='.4f')
import csv
with open('output.csv', newline='') as csvfile:
  reader = csv.reader(csvfile)
  i = 0
  for row in reader:
      if i!=0 and i!=7:
          j = 0
          for item in row:
            if j!=0:
                print("{:.4f}".format(float(item)), end=", ")
            j +=1
          print('\n')      
      i+=1
