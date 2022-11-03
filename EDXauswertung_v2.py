# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:33:28 2022

@author: poschmann
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import sys

def getvalue(valueinputtext, valuetype):
    valueinput = ''
    while (type(valueinput) != int) and (type(valueinput) != float):
        print("'e' for exit")
        valueinput = input(valueinputtext)
        try: 
            valueinput = valuetype(valueinput)
        except:
            if valueinput == 'e':
                sys.exit()
            else:
                print('wrong input')
    return valuetype(valueinput)

def getyninput(inputquestion):
    print("'e' for exit")
    ok = input(inputquestion)
    while (ok != 'n') and (ok != 'y'):
        if ok == 'e':
            sys.exit()
        else:
            print('wrong input')
            ok = input(inputquestion)
    return ok


#entry of filepath
os.chdir(input('Enter full path of data folder: '))

#looking for .TXT-file in working directory
#searching quantification files in folders with name containing 'region'
gefundene_dateien = []
path = os.getcwd()
for aktuelles_verzeichnis, unterverzeichnisse, dateien in os.walk(path):
  for datei in dateien:
      #operate on files called quantification.csv
    if datei == 'quantification.csv':  
        file_path = os.path.join(aktuelles_verzeichnis, datei)
        #operate on files in folders with 'region' in the name
        if 'region' in file_path.split('_'):
            gefundene_dateien.append(os.path.join(aktuelles_verzeichnis, datei))

#making a DataFrame of all the data in found quantification files
df = pd.concat((pd.read_csv(f) for f in gefundene_dateien), axis = 0)
#naming of columns
df.columns = ['atomic_number', 'element_symbol', 'element_name', 'atomic_percent', 'weight_percent', 'energy_level']
#calulating the mean values of each elment 
sorted_df_mean = df.groupby(['atomic_number', 'element_symbol'], axis = 0).mean()
#calulating the standard deviation values of each elment 
sorted_df_std = df.groupby(['atomic_number', 'element_symbol'], axis = 0).std()
#making combined dataframe
results = sorted_df_mean.drop(['atomic_percent', 'weight_percent', 'energy_level'], axis = 1)
results['atomic_percent'] = sorted_df_mean['atomic_percent'].round(3)
results['atomic_percent_stdDev']= sorted_df_std['atomic_percent'].round(3)
results['weight_percent']= sorted_df_mean['weight_percent'].round(3)
results['weight_percent_stdDev']= sorted_df_std['weight_percent'].round(3)
#saving combined dataframe
results.to_csv('EDX_results_mean.txt')
gefundene_dateien = []
path = os.getcwd()         
#searching spectrum.emsa files in subfolders with name containing 'region'
gefundene_spektren = []
for aktuelles_verzeichnis, unterverzeichnisse, dateien in os.walk(path):
  for datei in dateien:
      #operate on files called spectrum.emsa
    if datei == 'spectrum.emsa':
        file_path = os.path.join(aktuelles_verzeichnis, datei)
        #operate on files in folders with 'region' in the name
        if 'region' in file_path.split('_'):
            gefundene_spektren.append(os.path.join(aktuelles_verzeichnis, datei))
linetest =0
with open(str(gefundene_spektren[0])) as firstdata:
  for line in firstdata:
      if '#' not in line:
          break
      if '#Spectrum' not in line:
          linetest +=1
      else:
          break

#generate datafram of all EDX spectra
spektren_df = pd.concat((pd.read_csv(f, delimiter = ': ', usecols = [0], skiprows = linetest, skipfooter = 1, header = None, engine = 'python') for f in gefundene_spektren), axis = 1)
#values to generate X-axis (offset might change with recalibration)#values to generate X-axis (offset might change with recalibration)
npoints = 2048
perchan = 9.905
offset = -39.359

#function to normalize spectra
def absolute_maximum_scale(series):
    return series / series.abs().max()
#normalize every spectrum in dataframe
for col in spektren_df.columns:
    spektren_df[col] = absolute_maximum_scale(spektren_df[col])
#assign measurement numbers 0-N
spektren_df.columns = range(1, spektren_df.columns.size +1)
#generate empty list to put in difference values to mean spectrum later
diff_list = []
#generate plot of single spectra
for col in spektren_df.columns:
    spektren_clean = spektren_df[spektren_df[col] >= 0.00005]
    plt.plot(range(len(spektren_clean)), spektren_clean[col] **0.25)
    #calculate difference of each spectrum to mean spectrum and put it in to list
    diff = pd.DataFrame(spektren_clean[col] - spektren_clean.mean(axis=1))   
    #sum up difference spectrum as measure for possible outlier
    diff_list.extend(abs(diff.sum()))
#put out graph
plt.legend(range(1, len(spektren_df.columns) + 1))
plt.xlim([-10, 1200])
plt.ylim([0.1, 1])
plt.show()
#find spectrum that has highest difference to mean spectrum
max_diff_meas = diff_list.index(max(diff_list)) + 1
#give info about spectrum with highest difference to mean spectrum
print('The spectrum ' + str(max_diff_meas) + ' has highest difference to mean spectrum and might be an outlier.')
graph_legend = spektren_df.columns.to_list()
#delete single measurements form dataevaluation
delete_list = []
delete_spectra_q = getyninput('Should a spectrum be removed from evaluation (y/n):') #input('Should a spectrum be removed from evaluation (y/n):')
if delete_spectra_q == 'y' : 
    delete_mult_spectra = delete_spectra_q
    spektren_df.columns = graph_legend
    while delete_mult_spectra == 'y':
        delete_spectra_a = getvalue('Which spectrum should be removed from 1 - ' + str(graph_legend) + ': ', int) #int(input('Which spectrum should be removed from 1 - ' + str(spektren_df.columns.size) + ': '))
        delete_list.append(delete_spectra_a)
        if delete_spectra_a in graph_legend:
            spektren_df = spektren_df.drop(delete_spectra_a, axis = 1) 
            graph_legend = [element for element in graph_legend if element != delete_spectra_a]
            #generate empty list to put in difference values to mean spectrum later
            diff_list = []
            #generate plot of single spectra
            for col in spektren_df.columns:
                spektren_clean = spektren_df[spektren_df[col] >= 0.00005]
                plt.plot(range(len(spektren_clean)), spektren_clean[col] **0.25)
                #calculate difference of each spectrum to mean spectrum and put it in to list
                diff = pd.DataFrame(spektren_clean[col] - spektren_clean.mean(axis=1))   
                #sum up difference spectrum as measure for possible outlier
                diff_list.extend(abs(diff.sum()))
            #put out graph
            plt.legend(graph_legend)
            plt.xlim([-10, 1200])
            plt.ylim([0.1, 1])
            plt.show()
            #find spectrum that has highest difference to mean spectrum
            max_diff_meas = diff_list.index(max(diff_list)) + 1
            #give info about spectrum with highest difference to mean spectrum
            print('The spectrum ' + str(max_diff_meas) + ' has highest difference to mean spectrum and might be an outlier.')
        else :
            print('wrong input')
        delete_mult_spectra = getyninput('Do you want to delete further spectra? (y/n)') #input('Do you want to delete further spectra? (y/n)')
        
            
        
#generate dataframe of mean spectrum with standard deviation
Mean_Spectrum = pd.DataFrame()
Mean_Spectrum['Mean'] = spektren_df.mean(axis=1)
Mean_Spectrum['Std_Dev'] = spektren_df.std(axis=1)
Mean_Spectrum['Energy'] = (spektren_df.index * perchan + offset)/1000

#calculation of y-error range     
y_lower = ((Mean_Spectrum['Mean']) - (Mean_Spectrum['Std_Dev'])) ** 0.5
y_upper = ((Mean_Spectrum['Mean']) + (Mean_Spectrum['Std_Dev'])) ** 0.5

#plot creation
plt.figure(figsize=(10,8))
#shaded standard deviation
plt.fill_between(Mean_Spectrum['Energy'], y_lower, y_upper, alpha=0.3, label = 'Standard Deviation')
#plotting mean EDX spectrum
plt.plot(Mean_Spectrum['Energy'], Mean_Spectrum['Mean']**0.5, label = 'Mean Intensity')
#axis labels
plt.xlabel('Energy /keV')
plt.ylabel('Normalized Square Root Intensity /a.u.')
#legend
plt.legend()
#plot range
plt.xlim([0, 15])
plt.ylim([0, 1])
#plotsaving
plt.savefig('mean_EDX_spectrum.png', dpi = 300)
#plot showing
plt.show()

#searching quantification files in folders with name containing 'region'
for aktuelles_verzeichnis, unterverzeichnisse, dateien in os.walk(path):
  for datei in dateien:
      #operate on files called quantification.csv
    if datei == 'quantification.csv':  
        file_path = os.path.join(aktuelles_verzeichnis, datei)
        #operate on files in folders with 'region' in the name
        if 'region' in file_path.split('_'):
            gefundene_dateien.append(os.path.join(aktuelles_verzeichnis, datei))
        print(gefundene_dateien)
#deleting unwanted measurements from data evaluation with input
if delete_spectra_q == 'y':
    delete_measurement = getyninput('Do you want to delete the spectrum as well from EDX results table? (y/n):')#input('Do you want to delete the spectrum as well from EDX results table? (y/n):')
    if delete_measurement =='y':
        [x-1 for x in delete_list]
        # for loop deletes one after another and uses the index values. because the index changes during deletion process it has to start from the end!
        delete_list.sort(reverse = True)
        for entry in delete_list:
            del gefundene_dateien[entry-1]
#making a DataFrame of all the data in found quantification files
df = pd.concat((pd.read_csv(f) for f in gefundene_dateien), axis = 0)
#naming of columns
df.columns = ['atomic_number', 'element_symbol', 'element_name', 'atomic_percent', 'weight_percent', 'energy_level']
#calulating the mean values of each elment 
sorted_df_mean = df.groupby(['atomic_number', 'element_symbol'], axis = 0).mean()
#calulating the standard deviation values of each elment 
sorted_df_std = df.groupby(['atomic_number', 'element_symbol'], axis = 0).std()
#making combined dataframe
results = sorted_df_mean.drop(['atomic_percent', 'weight_percent', 'energy_level'], axis = 1)
results['atomic_percent'] = sorted_df_mean['atomic_percent'].round(3)
results['atomic_percent_stdDev']= sorted_df_std['atomic_percent'].round(3)
results['weight_percent']= sorted_df_mean['weight_percent'].round(3)
results['weight_percent_stdDev']= sorted_df_std['weight_percent'].round(3)
#saving combined dataframe
results.to_csv('EDX_results_mean.txt')

#generate dataframe with all spectra, mean spectra and standard deviation and output to file
all_spectra = spektren_df
all_spectra['Energy'] = (spektren_df.index * perchan + offset)/1000
all_spectra['Mean'] = Mean_Spectrum['Mean']
all_spectra['Std_Dev']= Mean_Spectrum['Std_Dev']
all_spectra = all_spectra.round(5)
all_spectra.to_csv('all_EDX_spectra.txt')
