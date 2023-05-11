"""
Test the effects of asynchronous plateau potentials in the detailed PFC L5 neuron on the soma

The orginal model is from
https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=117207&file=/acker_antic/Model/CA%20229.hoc#tabs-2
Modified by : tsjb2@cam.ac.uk

Run simulation with NMDAmajor.mod file
"""
import CA229 as de # from CA229 import *
import matplotlib.pyplot as plt
from neuron import h
import numpy as np
import utils as ut #from utils import *
import json
import itertools
import time
# import pdb     # For python debugging
# from random import *

h.load_file('stdrun.hoc') # for initialization

##############################
# Function to determine the randomness of activation time of AMPA and NMDA
# Uniform random seen in Fig 2. E1
# Remove the seed to get true randomness
def random_2(low, high, size):
    time_random = np.linspace(low, high, size)
    np.random.seed(10)
    np.random.shuffle(time_random)
    return time_random

# Alpha random seen in Fig 2. E2
def random_beta(low, high, size):
    time_random = np.random.beta(2, 8, size=size)
    np.random.seed(10)
    time_random = (high - low)*time_random + low
    return time_random

################### Test the ratio of different repceptors
def Glu_Stim(TTX = False, Pool1_num = 9, Pool2_num = 9, Syn_w1 = 0.01, Syn_w2 = 0.01, BasalIdx = 34, Loc = [0.2, 0.6]):
    """
    Model the Glumate Stimulation.
    Model the Receptors in 2 pools:
        Pool 1: AMPA + NMDA (same synaptic weight, represent spine conductance)
        Pool 2: NMDA only (represent the extrasyanptic NMDARs)

    Parameters:
    -----------
    TTX: True or False.
        True: setting all the sodium channel conductance to 0 to mimic
              TTX application in experiments
        False: default
    Pool1_num: syanptic AMPA/NMDA numbers
    Pool2_num: extrasyanptic NMDA numbers
    Syn_w1: the syanptic weight of AMPA/NMDA receptors in pool1
    Syn_w2: the syanptic weight of AMPA/NMDA receptors in pool2
    BasalIdx: The indices of the basal dendrites wherein we will place the AMPA and NMDA receptors
    Loc: the input location for AMPA and NMDA receptors in the dendrites
    -----------
    Outputs:
        Figures: recording from soma and 3 different locations from basal dendrites
        json: soma and dendritc voltage recording and parameters info
    """

    Cell = de.CA229()
    # Can adjust channel conductance ratio here:
    # eg. Cell = CA229(KA_ratio = 0.5)
    ###########################################
    timestr = time.strftime("%H%M")
    data = time.strftime("%m_%d")
    # directory = 'Data_' + data +'/'
    directory_root = 'PhysFig/'
    directory = directory_root + 'Major/Plot/'
    # directory = directory_root +'Major/Analysis/'

    if (TTX == True):
        Cell.TTX()
        title = "Major_TTX_Pool1_"+ str(Pool1_num) + "_Pool1_W_" + str(Syn_w1) + \
            "_Pool2_" + str(Pool2_num) + "_Pool2_W_" + str(Syn_w2) + "_" + timestr

    else:
        title = "Major_Pool1_"+ str(Pool1_num) + "_Pool1_W_" + str(Syn_w1) + \
            "_Pool2_" + str(Pool2_num) + "_Pool2_W_" + str(Syn_w2) + "_" + timestr


    ###########################################
    # Adding the stimulus and receptors to the dendrites
    ###########################################

    # Loc and time delay set up
    # For Pool 1:
    loc1 = list(np.linspace(Loc[0], Loc[1], Pool1_num))
    delay1 = random_2(10, 50 + int(Syn_w1*50), Pool1_num)
    # delay1 = random_beta(10, 50 + int(Syn_w1*50), Pool1_num)
    ns = h.NetStim()
    ns.interval = 20
    ns.number = 1
    ns.start = 190
    ns.noise = 0
    # For Pool 2:
    loc2 = list(np.linspace(Loc[0], Loc[1], Pool2_num))
    # delay2 = list(np.linspace(5, 10, Pool2_num))
    delay2 = random_2(15, 55 + int(Syn_w2*60), Pool2_num)
    # delay2 = random_beta(15, 55 + int(Syn_w2*60), Pool2_num)

    # Delays for signals between dendrites:
    delay = 50 # Time between the signals on each dendrite
    dendrite_delays = [(i * delay) for i in range(len(BasalIdx))]

    # Lists for Pool 1
    SynAMPA = []
    nc_AMPA = []
    SynNMDA = []
    # nc_NMDA = []
    # Lists for Pool 2
    ExNMDA = []
    nc_ExNMDA = []
    for idx, basal in enumerate(BasalIdx):
        ###########################################
        # Adding Pool 1
        ###########################################
        for i in range(Pool1_num):
            ###########################
            # Adding AMPA
            SynAMPA.append(h.AMPA(Cell.basal[basal](loc1[i])))
            SynAMPA[-1].gmax = 0.05
            nc_AMPA.append(h.NetCon(ns, SynAMPA[i]))
            nc_AMPA[-1].delay = delay1[i] + dendrite_delays[idx] # delay1[i] #uniform(1,20)
            nc_AMPA[-1].weight[0] = Syn_w1

        ###########################
        # Adding NMDA - Major
        for i in range(Pool1_num):
            tempNMDA = h.nmda(Cell.basal[basal](loc1[i]))
            tempNMDA.gmax = 0.005*Syn_w1
            tempNMDA.onset= delay1[i] + ns.start + dendrite_delays[idx]
            SynNMDA.append(tempNMDA)

        ###########################################
        # Adding Pool 2
        ###########################################
        for i in range(Pool2_num):
            ###########################
            # Adding extrasyanptic NMDA
            tempNMDA2 = h.nmda(Cell.basal[basal](loc2[i]))
            tempNMDA2.gmax = 0.005*Syn_w2
            tempNMDA2.onset= delay2[i] + ns.start + dendrite_delays[idx]
            ExNMDA.append(tempNMDA2)

    ###########################################
    ### Recording
    ###########################################
    t_vec = h.Vector()
    t_vec.record(h._ref_t)
    v_vec_soma = h.Vector()
    v_dend_dict = {}
    for basal in BasalIdx:
        v_dend_dict[basal] = h.Vector()
    cai_soma = h.Vector()

    v_vec_soma.record(Cell.soma[2](0.5)._ref_v)
    for basal in v_dend_dict:
        v_dend_dict[basal].record(Cell.basal[basal](Loc[0])._ref_v)

    cai_soma.record(Cell.soma[2](0.5)._ref_cai)


    ###########################################
    ### Run & Plot
    ###########################################
    h.celsius = 32
    h.v_init =  -73.6927850677
    h.init()
    h.tstop = 1000
    h.run()

#    pdb.set_trace()   #Debugging
    # print v_vec_soma[-1]
    # plt.clf()
    # plt.close()
    # plt.figure(figsize = (16, 6), dpi = 100)
    # plt.plot(t_vec, v_vec_soma, label = 'soma(0.5)', color = 'black')
    # plt.plot(t_vec, v_vec_dend1, label = 'bdend[34](0.8)', color = 'red')
    # plt.plot(t_vec, v_vec_dend2, label = 'Basal[34](0.5)', color = 'blue')
    # plt.plot(t_vec, v_vec_dend3, label = 'Basal[34](0.3)', color = 'green')
    # plt.ylim([-90, 40])
    # plt.xlim([0, 800])
    # plt.legend(loc = 'best')
    # plt.ylabel('mV')
    # plt.xlabel('Time (ms)')
    # plt.title ("Glumate Receptor Activated Plateau Potential")
    #
    # ut.save(title, directory, ext="png", close=True, verbose=True)

    #######################
    # Plot the intracelluar calcium concentration
    # plt.clf()
    # plt.close()
    # plt.figure(figsize = (16, 6), dpi = 100)
    # plt.plot(t_vec, cai_soma, label = 'soma(0.5)', color = 'black')
    # #plt.plot(t_vec, cai_dend, label = 'bdend[34](0.3)', color = 'red')
    #
    # # plt.ylim([-90, 60])
    # plt.xlim([0, 800])
    # plt.legend(loc = 'best')
    # plt.ylabel('mM')
    # plt.xlabel('Time (ms)')
    # plt.title ("Calcium concentration")
    # title1 = "Calcium_" + title
    # ut.save(title1, directory, ext="png", close=True, verbose=True)


    data = ut.Vividict()
    data['SynAMPA']['num'] = Pool1_num
    data['SynAMPA']['locs'] = loc1
    data['SynAMPA']['weight'] = Syn_w1
    data['SynNMDA']['num'] = Pool1_num
    data['SynNMDA']['locs'] = loc1
    data['SynNMDA']['weight'] = Syn_w1
    # data['SynNMDA']['Beta'] = Beta
    # data['SynNMDA']['Cdur'] = Cdur
    data['ExNMDA']['num'] = Pool2_num
    data['ExNMDA']['locs'] = loc2
    data['ExNMDA']['weight'] = Syn_w2
    # data['ExNMDA']['Beta'] = Beta
    # data['ExNMDA']['Cdur'] = Cdur

    data['recording']['time'] = list(t_vec)
    data['recording']['soma']['voltage'] = list(v_vec_soma)
    for basal in v_dend_dict:
        data['recording'][f'basal_{basal}'][f'voltage_{Loc[0]}'] = list(v_dend_dict[basal])
    data['recording']['soma']['ica'] = list(cai_soma)

    ut.savejson(data, title, directory, ext = "json", verbose = False)

######################################################
if __name__ == "__main__":
    print("Running the model")
    start_time = time.time()

    basal_indices = [34, 15, 14]
    loc = [0.9, 0.9]

    weight = 0.5
    Pool_num = 8 + int(20*weight)
    Glu_Stim(False, Pool_num, Pool_num, weight, weight, basal_indices, loc)
    print("Finished.")
    print("--- %s seconds ---" % (time.time() - start_time))