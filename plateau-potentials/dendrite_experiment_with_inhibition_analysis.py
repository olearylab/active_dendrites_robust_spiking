import json
import matplotlib.pyplot as plt
from pprint import pprint
import os
import numpy as np
import pandas as pd
import analysis_utils as ana
#from analysis_utils import *
from analysis_utils import tableau
import utils as ut #from utils import *
import seaborn as sns
import time

######################################################
# Analysis for Model2
Major_data = pd.DataFrame(columns = ['Inhib_weight', 'AMPA_num', 'AMPA_locs', 'AMPA_weight',
'NMDA_num', 'NMDA_locs', 'NMDA_weight',
'spike_num','platamp', 'ISI', 'platdur'])

path_to_json = 'PhysFig_with_inhibition/Analysis/'
start_time = time.time()
json_files = [pos_json for pos_json in os.listdir(path_to_json) if pos_json.endswith('.json')]

for index, js in enumerate(json_files):
    with open(os.path.join(path_to_json, js)) as json_file:
        data = json.load(json_file)
        filename = json_files[index]
        AMPA_num = data['SynAMPA']['num']
        AMPA_locs = data['SynAMPA']['locs']
        AMPA_weight = data['SynAMPA']['weight']
        NMDA_num = data['SynNMDA']['num']
        NMDA_locs = data['SynNMDA']['locs']
        NMDA_weight = data['SynNMDA']['weight']
        Inhib_weight = data['Inhibition']['weight']
        # NMDA_Beta = data['SynNMDA']['Beta']
        # NMDA_Cdur = data['SynNMDA']['Cdur']
        spike_num = ana.spike_count(data['recording']['soma']['voltage'])
        ISI, platamp = ana.meas_platamp(data['recording']['soma']['voltage'])
        platdur = ana.meas_platdur(data['recording']['soma']['voltage'])
        # For TTX
        # platamp = TTX_platamp(data['recording']['soma']['voltage'])
        # ISI = 0
        Major_data.loc[index] = [Inhib_weight, AMPA_num, AMPA_locs, AMPA_weight,
        NMDA_num, NMDA_locs, NMDA_weight,
        spike_num, platamp, ISI, platdur]

print("--- %s seconds ---" % (time.time() - start_time))
Major_data = Major_data.sort_values(by = ['Inhib_weight'])
# Add the correct saving path here:
savepath = path_to_json + '/total_results.csv'
Major_data.to_csv(savepath)

# ##### Plotting
File = 'PhysFig_with_inhibition/Analysis/total_results.csv'

df = pd.read_csv(File, index_col = 0)

savepath = 'PhysFig_with_inhibition/'
plt.close()
plt.clf()
plt.figure(figsize = (9,6), dpi = 100)
plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.15)
ax = plt.gca()
plt.style.use("ggplot")
plt.rcParams['axes.edgecolor'] = "black"
plt.rcParams['axes.facecolor'] = "white"
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.grid(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.plot(df['Inhib_weight'], df['platamp'],linestyle = '--', marker = 'o',
markersize = 15, c = tableau(20), linewidth = 2, label = 'Platamp')

plt.xlabel("Inhibitory weight", size = 22, color = 'black')
plt.ylabel("Plateau Amp (mV)", size = 22, color = 'black')
ax.tick_params(labelsize=22, pad = 12, colors = "black")
ax.set_facecolor('white')
plt.xlim([0,10])
plt.ylim([0,25])
ax.set_facecolor('white')
plt.legend(loc = 'best', fontsize = 22)
title = "Platamp_Model_B&W"
ut.save(title, savepath, ext="svg", close=True, verbose=True)

###################################################
plt.close()
plt.clf()
plt.figure(figsize = (9,6), dpi = 100)
plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.15)
ax = plt.gca()
plt.style.use("ggplot")
plt.rcParams['axes.edgecolor'] = "black"
plt.rcParams['axes.facecolor'] = "white"
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.grid(False)
# ax.spines['bottom'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.plot(df['Inhib_weight'] * (0.001 / (0.005 * 0.9)), df['platdur'],linestyle = '--', marker = 'o',
markersize = 15, c = tableau(20), linewidth = 2, label = 'Platdur')

plt.xlabel("g_Inhib / g_NMDA", size = 22, color = 'black')
plt.ylabel("Plateau Duration (ms)", size = 22, color = 'black')
ax.tick_params(labelsize=22, pad = 12, colors = "black")
ax.set_facecolor('white')
plt.xlim([0,1])
plt.ylim([-10,300])
ax.set_facecolor('white')
plt.legend(loc = 'best', fontsize = 22)
title = "Platdur_Model_B&W"
# save(title, savepath, ext="pdf", close=False, verbose=True)
ut.save(title, savepath, ext="png", close=True, verbose=True)