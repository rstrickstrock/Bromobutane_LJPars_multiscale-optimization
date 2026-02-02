### imports ###
import pandas as pd
import numpy as np
import os
import torch
import torch.nn as nn
import torch.optim as optim

parametersFile = 'this_parameters.csv'
modelConfigFile = 'modelConfig.csv'

#modelFiles = ['trainedModel_1-bromobutane.pth', 'trainedModel_2-bromobutane.pth']
modelFiles = ['trainedModel_2-bromobutane.pth']
predictionFile = 'this_prediction.csv'


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

class thisModel2(nn.Module):
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
      outputdim = int(architecture.split('-')[3])
    except:
      print(f'Can not get outputs dim from {architecture}')

    self.hidden1 = nn.Linear(inputdim, firstlayerdim)
    self.act1 = nn.LeakyReLU()
    self.hidden2 = nn.Linear(firstlayerdim, secondlayerdim)
    self.act2 = nn.LeakyReLU()
    self.output = nn.Linear(secondlayerdim, outputdim)
    self.act_output = nn.LeakyReLU()

  def forward(self, x):
    x = self.act1(self.hidden1(x))
    x = self.act2(self.hidden2(x))
    x = self.act_output(self.output(x))
    return x


thisParameters = pd.read_csv(parametersFile).to_numpy()
thisParameters = torch.tensor(thisParameters, dtype=torch.float32)
#print(f'{thisParameters}')

thisModelConfig = pd.read_csv(modelConfigFile)
#print(f'{thisModelConfig}')
thisNEpochs = thisModelConfig.n_epochs.to_numpy()[0]
thisLearningRate = thisModelConfig.learning_rate.to_numpy()[0]
thisArchitecture = str(thisModelConfig.architecture.to_numpy()[0])
thisHiddenLayer = thisModelConfig.hiddenLayer.to_numpy()[0]
#print(f'{thisNEpochs}')
#print(f'{thisLearningRate}')
#print(f'{thisArchitecture}')
#print(f'{thisHiddenLayer}')

if thisHiddenLayer == 2:
  thisModel = thisModel2(thisArchitecture)
if thisHiddenLayer == 3:
  thisModel = thisModel3(thisArchitecture)
else:
  print(f'thisHiddenLayer needs to be 2 or 3 (is {thisHiddenLayer}). Exit.')
  exit()

if os.path.exists(predictionFile):
  os.remove(predictionFile)
  print(f'Removed existing predictions file: \'{predictionFile}\'.')
with open(predictionFile, 'w') as f:
  pass

for modelFile in modelFiles:
  thisModel.load_state_dict(torch.load(modelFile, weights_only=True))
  thisModel.eval()

  thisPrediction = thisModel(thisParameters).item()
  #print(f'{thisPrediction}')

  with open(predictionFile, 'a') as f:
    f.write(f'{modelFile} {thisPrediction}\n')

