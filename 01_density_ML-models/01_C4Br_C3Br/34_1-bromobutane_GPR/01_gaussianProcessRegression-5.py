### imports ###
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import Matern
from sklearn.gaussian_process.kernels import RationalQuadratic as RQ
from sklearn.gaussian_process.kernels import ExpSineSquared as ESS

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_percentage_error as skmape

import pandas as pd
import numpy as np
import os
import shutil

import pickle

testmode = True

rndInts = [678, 147, 561, 237, 588, 951, 490, 395, 877, 297, 721, 711, 985, 171, 75, 16, 669, 530, 999, 794, 936, 111, 816, 968, 48, 986, 829, 996, 272, 759, 390, 930, 633, 928, 854, 554, 562, 78, 222, 294, 725, 582, 731, 249, 791, 35, 180, 510, 593, 634]

testSizes = [0.30, 0.70]

statisticsFileName = 'Stats-5.csv'



pwd = os.getcwd()
### create current working directory ###
cwd = os.path.join(pwd, "trainedModels")
#print(f'{cwd}')
if os.path.exists(cwd):
  if testmode:
    pass
    #shutil.rmtree(cwd)
  else:
    print(f'\nPATH \'{cwd}\' already exists. \n\nExiting without starting or changing anything.\n')
    exit()
else:
  os.mkdir(cwd)
#os.mkdir(cwd)


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
                             "kernel": [],
                             "length_scale": [],
                             "nu": [],
                             "alpha": [],
                             "periodicity": [],
                             "mape": [],
                             "r2": []})

n_restarts = 9

for thisRatio in testSizes:
  trainSize = 1-thisRatio
  print(f'Training GPR Models with {trainSize}% of the dataset.\n')
  tmp_cwd = os.path.join(cwd, str(thisRatio))
  os.mkdir(tmp_cwd)
  os.chdir(tmp_cwd)
  
  for rndInt in rndInts:
    #print(f'.')
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=thisRatio, random_state=rndInt)
    
    kernel1 = RBF()
    kernel2 = Matern()
    kernel3 = RQ()
    kernel4 = ESS()
    
    model1 = GaussianProcessRegressor(kernel=kernel1, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=n_restarts, random_state=rndInt)
    model2 = GaussianProcessRegressor(kernel=kernel2, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=n_restarts, random_state=rndInt)
    model3 = GaussianProcessRegressor(kernel=kernel3, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=n_restarts, random_state=rndInt)
    model4 = GaussianProcessRegressor(kernel=kernel4, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=n_restarts, random_state=rndInt)
    
    model1.fit(X_train, Y_train)
    pickle.dump(model1, open(f'trained_model_{thisRatio}_{rndInt}_RBF.sav', 'wb'))
    model2.fit(X_train, Y_train)
    pickle.dump(model2, open(f'trained_model_{thisRatio}_{rndInt}_Matern.sav', 'wb'))
    model3.fit(X_train, Y_train)
    pickle.dump(model3, open(f'trained_model_{thisRatio}_{rndInt}_RQ.sav', 'wb'))
    model4.fit(X_train, Y_train)
    pickle.dump(model4, open(f'trained_model_{thisRatio}_{rndInt}_ESS.sav', 'wb'))    

    Y_1prediction = model1.predict(X_test)
    Y_2prediction = model2.predict(X_test)
    Y_3prediction = model3.predict(X_test)
    Y_4prediction = model4.predict(X_test)
    mape1 = skmape(Y_test, Y_1prediction)
    mape2 = skmape(Y_test, Y_2prediction)
    mape3 = skmape(Y_test, Y_3prediction)
    mape4 = skmape(Y_test, Y_4prediction)
    r21 = r2_score(Y_test, Y_1prediction)
    r22 = r2_score(Y_test, Y_2prediction)
    r23 = r2_score(Y_test, Y_3prediction)
    r24 = r2_score(Y_test, Y_4prediction)
    
    dfEntry = pd.DataFrame({"ratio": [thisRatio],
                            "rndint": [rndInt],
                            "kernel": ["RBF"],
                            "length_scale": [model1.kernel.length_scale],
                            "nu": ["-"],
                            "alpha": ["-"],
                            "periodicity": ["-"],
                            "mape": [mape1],
                            "r2": [r21]})
    dfStatistics = pd.concat([dfStatistics, dfEntry], ignore_index=True)
    dfEntry = pd.DataFrame({"ratio": [thisRatio],
                            "rndint": [rndInt],
                            "kernel": ["Matern"],
                            "length_scale": [model2.kernel.length_scale],
                            "nu": [model2.kernel.nu],
                            "alpha": ["-"],
                            "periodicity": ["-"],
                            "mape": [mape2],
                            "r2": [r22]})
    dfStatistics = pd.concat([dfStatistics, dfEntry], ignore_index=True)
    dfEntry = pd.DataFrame({"ratio": [thisRatio],
                            "rndint": [rndInt],
                            "kernel": ["RQ"],
                            "length_scale": [model3.kernel.length_scale],
                            "nu": ["-"],
                            "alpha": [model3.kernel.alpha],
                            "periodicity": ["-"],
                            "mape": [mape3],
                            "r2": [r23]})
    dfStatistics = pd.concat([dfStatistics, dfEntry], ignore_index=True)
    dfEntry = pd.DataFrame({"ratio": [thisRatio],
                            "rndint": [rndInt],
                            "kernel": ["ESS"],
                            "length_scale": [model4.kernel.length_scale],
                            "nu": ["-"],
                            "alpha": ["-"],
                            "periodicity": [model4.kernel.periodicity],
                            "mape": [mape4],
                            "r2": [r24]})
    dfStatistics = pd.concat([dfStatistics, dfEntry], ignore_index=True)

  print(f'\n')
  os.chdir(cwd)

os.chdir(pwd)
if os.path.exists(statisticsFileName):
  os.remove(statisticsFileName)
  print(f'Removed existing statistics file: \'{statisticsFileName}\'.')
dfStatistics.to_csv(statisticsFileName)
print(f'Wrote statistics to file: \'{statisticsFileName}\'.')
print(f'{dfStatistics}')










