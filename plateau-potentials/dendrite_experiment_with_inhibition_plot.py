import json
import matplotlib.pyplot as plt
from pprint import pprint
import os
import numpy as np
import pandas as pd
#from analysis_utils import *
from analysis_utils import tableau
import utils as ut #from utils import *
import seaborn as sns


###### Load DATA
## Add the correct analysis path here:
path_to_json = 'PhysFig_with_inhibition/Plot/'
json_files = [pos_json for pos_json in os.listdir(path_to_json) if pos_json.endswith('.json')]
nfile = len(json_files)

df = pd.DataFrame(columns = ['soma_trace','dend1_trace','dend2_trace','dend3_trace','labels'])

for index, js in enumerate(json_files):
    with open(os.path.join(path_to_json, js)) as json_file:
        data = json.load(json_file)
        filename = json_files[index]
        time = data['recording']['time'][4000:]
        new_time = [x-100.0 for x in time]
        soma_trace = data['recording']['soma']['voltage'][4000:]
        dend1_trace = data['recording']['basal_34']['voltage_0.8'][4000:]
        dend2_trace = data['recording']['basal_34']['voltage_0.5'][4000:]
        dend3_trace = data['recording']['basal_34']['voltage_0.3'][4000:]
        #labels = data['ExNMDA']['weight']
        labels = data['Inhibition']['weight']
        df.loc[index] = [soma_trace, dend1_trace, dend2_trace, dend3_trace, labels]

df = df.sort_values(by = ['labels'], ascending = True)

############## Shifted Dend0.8 plot
plt.close()
plt.clf()
plt.figure(figsize = (9,6), dpi = 300)
plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.15)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

#new_dend1_trace = []
new_dend2_trace = []
for i in range(nfile-1):
    time1 = [x+330*i for x in new_time]
    #new_dend1_trace.append([x+20*i for x in df['dend1_trace'].iloc[i+1]])
    #plt.plot(time1, new_dend1_trace[i], color = tableau(i+1), label = 'Input w = ' + str('%.2f' % df['labels'].iloc[i+1]), linewidth = 2)
    new_dend2_trace.append([x+20*i for x in df['dend2_trace'].iloc[i+1]])
    plt.plot(time1, new_dend2_trace[i], color = tableau(i+1), label = 'Input w = ' + str('%.2f' % df['labels'].iloc[i+1]), linewidth = 2)
    # plt.legend(loc = 'best', fontsize = 15)
    #plt.title("Voltage Traces at Basal Dendrite (loc:0.8) \n", fontsize = 25)
    plt.title("Voltage Traces at Basal Dendrite (loc:0.5) \n", fontsize = 25)
    plt.ylim([-80, 60])
    # plt.ylabel('mV', fontsize = 23)
    plt.yticks(fontsize=23)
    plt.xlim([50, 2250])
    #plt.xlabel('Time (ms)', fontsize = 15)
    # plt.xticks([])
    # plt.yticks([])
    plt.xticks(fontsize = 23)
    plt.gca().yaxis.grid(False)
    plt.gca().set_facecolor('white')
    #plt.gca().axes.get_xaxis().set_visible(False)
    #plt.gca().axes.get_yaxis().set_visible(False)
    plt.xlabel('Time (ms)', fontsize = 23)

title2 = "Dend0.8_shift"
ut.save(title2, path_to_json, ext="svg", close=True, verbose=True)
# save(title2, path_to_json, ext="ps", close=True, verbose=True)