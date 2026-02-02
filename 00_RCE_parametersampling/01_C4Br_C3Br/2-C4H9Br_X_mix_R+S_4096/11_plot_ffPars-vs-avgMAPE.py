#import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

fileName = '2-Bromobutane_RCE.csv'
substName = f'{fileName[0:13]}'
titleFontSize = 20
axisLabelFontSize = 16
cBarLabel = f'Avg. MAPE [-]'
figDPI = 300
filterDF = True
maxMAPE = 5.0
parameterCombinations = [[f'SigC', f'SigBr'],
                         [f'SigC', f'EpsC'],
                         [f'SigC', f'EpsBr'],
                         [f'SigBr', f'EpsC'],
                         [f'SigBr', f'EpsBr'],
                         [f'EpsC', f'EpsBr'],]

def getDataAndLabels(x, y):
  if x == f'SigC':
    X = dfResults.SigC800
    xlbl = f'Sigma C [nm]'
  elif x == f'SigBr':
    X = dfResults.SigBr730
    xlbl = f'Sigma Br [nm]'
  elif x == f'EpsC':
    X = dfResults.EpsC800
    xlbl = f'Epsilon C [kJ/mol]'
  elif x == f'EpsBr':
    X = dfResults.EpsBr730
    xlbl = f'Epsilon Br [kJ/mol]'

  if y == f'SigC':
    Y = dfResults.SigC800
    ylbl = f'Sigma C [nm]'
  elif y == f'SigBr':
    Y = dfResults.SigBr730
    ylbl = f'Sigma Br [nm]'
  elif y == f'EpsC':
    Y = dfResults.EpsC800
    ylbl = f'Epsilon C [kJ/mol]'
  elif y == f'EpsBr':
    Y = dfResults.EpsBr730
    ylbl = f'Epsilon Br [kJ/mol]'
  
  return X, xlbl, Y, ylbl



try:
  saveOrShow = sys.argv[1]
except:
  saveOrShow = "show"
if saveOrShow == "save":
  pass
else:
  saveOrShow = "show"
  
if not os.path.isfile(fileName):
  print(f'Can not find and open \'{fileName}\'. Exit.')
  exit()
else:
  dfResults = pd.read_csv(fileName)
  #print(f'{dfResults}')
  try:
    dfResults = dfResults.drop(columns=["Unnamed: 0"])
  except:
    print(f'Something went wrong with\'dfResults = dfResults.drop(columns=["Unnamed: 0"])\'.')
  else:
    #print(f'{dfResults}')
    pass

if filterDF:
  dfResults = dfResults[dfResults["AvgMAPE"] <= maxMAPE]

for combi in parameterCombinations:
  x = combi[0]
  y = combi[1]
  X, xlbl, Y, ylbl = getDataAndLabels(x, y)
  plt.figure()
  plt.scatter(X, Y, c=dfResults.AvgMAPE)
  plt.xlabel(xlbl, fontsize=axisLabelFontSize, weight='bold')
  plt.ylabel(ylbl, fontsize=axisLabelFontSize, weight='bold')
  plt.suptitle(f'{substName}', fontsize=titleFontSize, weight='bold')
  plt.colorbar(label=cBarLabel)
  plt.tight_layout()
  if saveOrShow == "save":
    plt.savefig(f'{substName}_AvgMAPE-RCE_{x}-vs{y}.png', dpi=figDPI)

if saveOrShow == "show":
  plt.show()
