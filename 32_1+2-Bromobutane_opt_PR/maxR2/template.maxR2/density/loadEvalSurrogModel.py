### imports ###
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import pickle
import pandas as pd
import os

parametersFile = 'this_parameters.csv'
predictionFile = 'this_prediction.csv'

modelFiles = ['trainedModel_1-bromobutane.sav', 'trainedModel_2-bromobutane.sav']

thisParameters = pd.read_csv(parametersFile)
#print(f'{thisParameters}')

if os.path.exists(predictionFile):
  os.remove(predictionFile)
  print(f'Removed existing predictions file: \'{predictionFile}\'.')
with open(predictionFile, 'w') as f:
  pass

for modelFile in modelFiles:
  with open(modelFile, "rb") as input_file:
    thisModel = pickle.load(input_file)
    
  thisPrediction = thisModel.predict(thisParameters)
  #print(f'{thisPrediction}')

  with open(predictionFile, 'a') as f:
    f.write(f'{modelFile} {thisPrediction[0]}\n')



  

