import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_percentage_error as skmape

import pandas as pd
import os
import time
import shutil

n_epochs = 10000
learning_rate = 0.001
architecture = f'4-128-64-32-1'
hiddenLayer = 3

trainingDatasetNames = ["1-bromobutane"]
rndInts = [678, 147]
trainRatios = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
testmode = True

class thisModel3(nn.Module):
  def __init__(self, architecture):
    super().__init__()
    try:
      inputdim = int(architecture.split('-')[0])
    except:
      print(f'Can not get input dim from {architecture}')
    try:
      firstlayerdim = int(architecture.split('-')[1])
    except:
      print(f'Can not get firstlayer dim from {architecture}')
    try:
      secondlayerdim = int(architecture.split('-')[2])
    except:
      print(f'Can not get secondlayer dim from {architecture}')
    try:
      thirdlayerdim = int(architecture.split('-')[3])
    except:
      print(f'Can not get thirdlayer dim from {architecture}')
    try:
      outputdim = int(architecture.split('-')[4])
    except:
      print(f'Can not get outputs dim from {architecture}')
      
    self.hidden1 = nn.Linear(inputdim, firstlayerdim)
    self.act1 = nn.LeakyReLU()
    self.hidden2 = nn.Linear(firstlayerdim, secondlayerdim)
    self.act2 = nn.LeakyReLU()
    self.hidden3 = nn.Linear(secondlayerdim, thirdlayerdim)
    self.act3 = nn.LeakyReLU()
    self.output = nn.Linear(thirdlayerdim, outputdim)
    self.act_output = nn.LeakyReLU()

  def forward(self, x):
    x = self.act1(self.hidden1(x))
    x = self.act2(self.hidden2(x))
    x = self.act3(self.hidden3(x))
    x = self.act_output(self.output(x))
    return x

pwd = os.getcwd()
cwd = os.path.join(pwd, "trainedModels")
#print(f'{cwd}')
if os.path.exists(cwd):
  if testmode:
    #shutil.rmtree(cwd)
    pass
  else:
    print(f'\nPATH \'{cwd}\' already exists. \n\nExiting without starting or changing anything.\n')
    exit()
else:
  os.mkdir(cwd)


data = pd.read_csv('1-bromobutane_trainingsdata.csv')
data = data.drop(data.columns[0], axis=1)
#print(f'{data}')
X = data.drop('density', axis=1)
X = X.drop('density_err', axis=1)
Y = data['density']
#print(f'{X}')
#print(f'{Y}')


for thisRatio in trainRatios:
  thisTestRatio = 1-thisRatio
  tmpCWD = os.path.join(cwd, str(thisRatio))
  if os.path.exists(tmpCWD):
    if testmode:
      #shutil.rmtree(tmpCWD)
      pass
    else:
      print(f'\nPATH \'{tmpCWD}\' already exists. \n\nExiting without starting or changing anything.\n')
      exit()
  else:
    os.mkdir(tmpCWD)
  os.chdir(tmpCWD)
    
  for rndInt in rndInts:
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=thisTestRatio, random_state=rndInt)
    X_train = torch.tensor(X_train.to_numpy(), dtype=torch.float32)
    Y_train = torch.tensor(Y_train.to_numpy(), dtype=torch.float32).reshape(-1, 1)
    X_test = torch.tensor(X_test.to_numpy(), dtype=torch.float32)
    Y_test = torch.tensor(Y_test.to_numpy(), dtype=torch.float32).reshape(-1, 1)
    
    dfStatistics = pd.DataFrame({"ratio": [],
                                 "rndint": [],
                                 "learning_rate": [],
                                 "batchsize": [],
                                 "epochs": [],
                                 "loss": [],
                                 "time": [],
                                 "mape": [],
                                 "r2": []})

    for trainingDatasetName in trainingDatasetNames:
      if hiddenLayer == 3:
        model = thisModel3(architecture)
      else:
        print(f'hiddenLayer needs to be 3 (is {hiddenLayer}). Exit.')
        exit()
      
      loss_fn = nn.L1Loss()
      optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        
      batch_size = 5
      while True:
        if len(X_train)%batch_size == 0:
          break
        else:
          batch_size = batch_size + 1
  
      startTime = time.time()
      for epoch in range(n_epochs):
        for i in range(0, len(X_train), batch_size):
          Xbatch = X_train[i:i+batch_size]
          Ypred = model(Xbatch)
          Ybatch = Y_train[i:i+batch_size]
          loss = loss_fn(Ypred, Ybatch)
          optimizer.zero_grad()
          loss.backward()
          optimizer.step()
        
      thisTrainingTime = time.time() - startTime
      torch.save(model.state_dict(), f'trainedModel_{thisRatio}_{rndInt}.pth')
      
      YpredTrain = model(X_train)
      YpredTest = model(X_test)
        
      YpredictionsTest = []
      YexpectationsTest = []
      for i in range(len(YpredTest)):
        YpredictionsTest.append(YpredTest[i].item())
        YexpectationsTest.append(Y_test[i].item())
      thisMAPETest = skmape(YpredictionsTest, YexpectationsTest)
      thisR2Test = r2_score(YpredictionsTest, YexpectationsTest)

      
      dfThisStats = pd.DataFrame({"ratio": [thisRatio],
                                  "rndint": [rndInt],
                                  "learning_rate": [learning_rate],
                                  "batchsize": [batch_size],
                                  "epochs": [n_epochs],
                                  "loss": [loss.item()],
                                  "time": [thisTrainingTime],
                                  "mape": [thisMAPETest],
                                  "r2": [thisR2Test]})
      
      dfStatistics = pd.concat([dfStatistics, dfThisStats], ignore_index=True)

    statisticsFileName = f'StatsPart_{thisRatio}_{rndInt}.csv'
    statsFilePath = os.path.join(pwd, statisticsFileName)
  
    if os.path.exists(statsFilePath):
      os.remove(statsFilePath)
      print(f'Removed existing statistics file: \'{statsFilePath}\'.')
    dfStatistics.to_csv(statsFilePath)
    print(f'Wrote statistics to file: \'{statsFilePath}\'.')
    print(f'{dfStatistics}')



