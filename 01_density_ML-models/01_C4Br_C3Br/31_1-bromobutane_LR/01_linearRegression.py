from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score 
from sklearn.metrics import mean_absolute_percentage_error as skmape

import pandas as pd
import numpy as np
import os
import shutil

import pickle

testmode = True

pwd = os.getcwd()
### create current working directory ###
cwd = os.path.join(pwd, "trainedModels")
#print(f'{cwd}')
if os.path.exists(cwd):
  if testmode:
    shutil.rmtree(cwd)
  else:
    print(f'\nPATH \'{cwd}\' already exists. \n\nExiting without starting or changing anything.\n')
    exit()

os.mkdir(cwd)

rndInts = [678, 147, 561, 237, 588, 951, 490, 395, 877, 297, 721, 711, 985, 171, 75, 16, 669, 530, 999, 794, 936, 111, 816, 968, 48, 986, 829, 996, 272, 759, 390, 930, 633, 928, 854, 554, 562, 78, 222, 294, 725, 582, 731, 249, 791, 35, 180, 510, 593, 634]

testSizes = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]

## grid sampling 1296
data = pd.read_csv('1-bromobutane_trainingsdata.csv')
data = data.drop(data.columns[0], axis=1)
#print(f'{data}')
#data = data.transpose()
#labels = data.iloc[0]
#data = data[1:]
#data.columns = labels
#print(f'{data}')

X = data.drop('density', axis=1)
X = X.drop('density_err', axis=1)
Y = data['density']
#print(f'{data}')
#print(f'{X}')
#print(f'{Y}')

dfStatistics = pd.DataFrame({"ratio": [],
                             "rndint": [],
                             "mape": [],
                             "r2": []})

for thisRatio in testSizes:
  tmp_cwd = os.path.join(cwd, str(thisRatio))
  os.mkdir(tmp_cwd)
  os.chdir(tmp_cwd)
  
  for rndInt in rndInts:
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=thisRatio, random_state=rndInt)
    
    model = LinearRegression()  
    model.fit(X_train, Y_train)
    pickle.dump(model, open(f'trained_model_{thisRatio}_{rndInt}.sav', 'wb'))
    
    Y_prediction = model.predict(X_test)
    mape = skmape(Y_test, Y_prediction)
    r2 = r2_score(Y_test, Y_prediction)
     
    dfEntry = pd.DataFrame({"ratio": [thisRatio],
                            "rndint": [rndInt],
                            "mape": [mape],
                            "r2": [r2]})
    dfStatistics = pd.concat([dfStatistics, dfEntry], ignore_index=True)
    
  os.chdir(cwd)

os.chdir(pwd)
statisticsFileName = 'Stats.csv'
if os.path.exists(statisticsFileName):
  os.remove(statisticsFileName)
  print(f'Removed existing statistics file: \'{statisticsFileName}\'.')
dfStatistics.to_csv(statisticsFileName)
print(f'Wrote statistics to file: \'{statisticsFileName}\'.')
print(f'{dfStatistics}')










