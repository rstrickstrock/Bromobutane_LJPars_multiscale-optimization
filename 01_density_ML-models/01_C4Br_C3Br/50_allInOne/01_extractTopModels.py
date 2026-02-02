import pandas as pd
import os
import shutil

top = 20

def readThisCSV(filename):
  if not os.path.isfile(filename):
    print(f'Can not find and open \'{filename}\'. Exit.')
    exit()
  else:
    dfStatistics = pd.read_csv(filename)
    #print(f'{dfStatistics}')
    try:
      dfStatistics = dfStatistics.drop(columns=["Unnamed: 0"])
    except:
      print(f'Something went wrong with\'dfStatistics = dfStatistics.drop(columns=["Unnamed: 0"])\'.')
    else:
      #print(f'{dfStatistics}')
      pass
      
  return dfStatistics

#statisticsFilesLR = ['Stats_LR_1-bromobutane.csv', 'Stats_LR_2-bromobutane.csv']
#statisticsFilesPR = ['Stats_PR_1-bromobutane.csv', 'Stats_PR_2-bromobutane.csv']
#statisticsFilesRFR = ['Stats_RFR_1-bromobutane.csv', 'Stats_RFR_2-bromobutane.csv']
#statisticsFilesGPR = ['Stats_GPR_1-bromobutane.csv', 'Stats_GPR_2-bromobutane.csv']
#statisticsFilesNNR = ['Stats_NNR_1-bromobutane.csv', 'Stats_NNR_2-bromobutane.csv']
statisticsFiles = ['Stats_LR_1-bromobutane.csv', 'Stats_PR_1-bromobutane.csv', 'Stats_RFR_1-bromobutane.csv', 'Stats_GPR_1-bromobutane.csv', 'Stats_NNR_1-bromobutane.csv', 'Stats_LR_2-bromobutane.csv', 'Stats_PR_2-bromobutane.csv', 'Stats_RFR_2-bromobutane.csv', 'Stats_GPR_2-bromobutane.csv', 'Stats_NNR_2-bromobutane.csv']


for statsFile in statisticsFiles:
  thisOutFile = f'{statsFile.split(".")[0]}_top{top}.csv'
  #print(f'{thisOutFile}')
  thisDF = readThisCSV(statsFile)
  #print(f'{thisDF}')
  
  minMapeIDXs = []
  maxR2IDXs = []

  dfTopModels = pd.DataFrame({"ratio": [],
                              "rndint": [],
                              "mape": [],
                              "r2": []})
                                  
  thisTop = 0
  dfThisStatistics = thisDF
  #print(f'MinMAPE Models:')
  while True:
    idxMinMAPE = dfThisStatistics['mape'].idxmin()
    #print(f'{idxMinMAPE}')
    #minMAPE = dfThisStatistics['mape'].min()
    #print(f'{minMAPE}')
    minMAPERow = dfThisStatistics.loc[idxMinMAPE]
    #print(f'Min MAPE Entry:\n{minMAPERow}\n')
    dfThisStatistics = dfThisStatistics.drop([idxMinMAPE]) 
  
    minMapeIDXs.append(idxMinMAPE)
    dfThisTopModel = pd.DataFrame({"ratio": [minMAPERow.ratio],
                                   "rndint": [minMAPERow.rndint],
                                   "mape": [minMAPERow.mape],
                                   "r2": [minMAPERow.r2]})
                                   
    dfTopModels = pd.concat([dfTopModels, dfThisTopModel], ignore_index=True)
    
    thisTop = thisTop + 1
    if thisTop == top:
      #print(f'')
      break


  thisTop = 0
  dfThisStatistics = thisDF
  #print(f'MaxR2 Models:')
  while True:
    idxMaxR2 = dfThisStatistics['r2'].idxmax()
    #print(f'{idxMaxR2}')
    #maxR2 = dfThisStatistics['r2'].max()
    #print(f'{maxR2}')
    maxR2Row = dfThisStatistics.loc[idxMaxR2]
    #print(f'Max R2 Entry:\n{maxR2Row}\n')
    dfThisStatistics = dfThisStatistics.drop([idxMaxR2])
  
    maxR2IDXs.append(idxMaxR2)
    dfThisTopModel = pd.DataFrame({"ratio": [maxR2Row.ratio],
                                   "rndint": [maxR2Row.rndint],
                                   "mape": [maxR2Row.mape],
                                   "r2": [maxR2Row.r2]})
    dfTopModels = pd.concat([dfTopModels, dfThisTopModel], ignore_index=True)
    

    thisTop = thisTop + 1
    if thisTop == top:
      #print(f'')
      break

  if os.path.exists(thisOutFile):
    os.remove(thisOutFile)
    print(f'Removed existing top{top} statistics file: \'{thisOutFile}\'.')
  dfTopModels.to_csv(thisOutFile)
  print(f'Wrote statistics to file: \'{thisOutFile}\'.')
  #print(f'{dfTopModels}')
  #print(f'')









  

