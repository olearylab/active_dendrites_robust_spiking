"""
The analysis functions for EEE stimulation data.
1. Plateau Amplitude
2. Plateau Duration
3. Spike Numbers on plateau
4. Interspike intervals
5. EPSPs amplitude and time stamp

Note: new analysis functions are added on August 1st, 2018
Plateu amplitude and duration measurements in basal dendrites with soma

Author: Peng Penny Gao
penggao.1987@gmail.com

Contributors:
salvadordura@gmail.com
joe.w.graham@gmail.com
"""

import json
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import seaborn as sns

########################################
### Function: to measure plateau duration
########################################
def meas_platdur(data, thresh = 10, dt = 0.025):
    """Measures plateau duration (ms)

    Parameters:
    -----------
        data: voltage trace
        thresh (mV): measure the duration of the plateau which are bigger than baseline + thresh
        dt: the sampling rate (default is 0.025ms or 40kHz)
    Return:
    -----------
        platur (ms): plateau duration in ms
    """
    # Make sure dt = 0.025 and there is no stimulation for the first 50 ms
    # If there is anything before 50ms, change the time points for baseline
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    above = [val for val in stable if val > (baseline + thresh)]
    platdur = dt * len(above)
    return platdur

########################################
### Function: to calculate the number of Spikes
########################################
def spike_count(data, thresh = 0):
    """Measure the number of Spikes

    Parameters:
    -----------
    data: list
        The electrical recording data from soma in the CA229 model.
        It was recorded in h.Vector() originally, then saved as list into json
        file duirng the batch processing of simulation.

    thresh: int
        The threshold to detect spikes. By default: if the membrane potential
        is larger than 0 mV, it will be detected as a spike.

    Return:
    -----------
    count: int
        Number of spikes in the data.
    """
    spike_flag = False
    count = 0
    for idx, val in enumerate(data[:-1]):
        if spike_flag == False and val < thresh and data[idx+1] >= thresh:
            spike_flag = True
        elif spike_flag == True and val >= thresh and data[idx+1] < thresh:
            count += 1
            spike_flag = False
    return count

########################################
### Function: to calculate spike interval
########################################
def IST_spikes (data, dt = 0.025):
    """
    Get the interspike intervals, spike_idx, spike_mvalue and spike_midx
    Return:
    -----------
    spike_mvalue:
        The list of minimum values between each two spikes
    spike_midx:
        The list of timepoints of minimum values
    ISTs:
        The average time of interspike intervals
    """
    spikes = get_EPSPs(data, 60, dt)
    IST = []
    count = len(spikes)
    spike_mvalue = []
    spike_midx = []
    if count <= 1:
        ISTs = 0
    else:
        for i in range(count-1):
            IST.append(spikes[i+1][0]-spikes[i][0])
            idx1 = int(round((spikes[i][0])/dt))
            idx2 = int(round((spikes[i+1][0])/dt))
            temp = np.min(data[idx1:idx2])
            spike_mvalue.append(temp)
            spike_midx.append(np.argmin(data[idx1:idx2]) + idx1)

        ISTs = np.mean(IST)
    return ISTs , spike_mvalue, spike_midx

########################################
### Function: to analyze EPSPs
########################################
def get_EPSPs(data, thresh = 2, dt = 0.025):
    """
    Get the maximum index and value of each EPSP peak
    Return:
    -----------
    spikes:
        A list of turple.
        Each turple has the index of max EPSP, and value of max EPSP.
    """
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    EPSPs = []
    for idx, val in enumerate(stable[:-1]):
        if (val >= (baseline + thresh)) and (val > stable[idx-1]) and (val >= stable[idx+1]):
            time = (idx + 6000) * dt
            EPSPs.append((time, val))
    return EPSPs

########################################
### Function: to analyze spikes
########################################
def single_spike(data, dt = 0.025):
    """
    Get the maximum index and value of each EPSP peak
    Return:
    -----------
    spikes:
        A list of turple.
        Each turple has the index of max EPSP, and value of max EPSP.
    """
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    peak_v = np.amax(stable) - baseline
    peak_t = (np.argmax(stable) + 6000) * dt
    return peak_v, peak_t

######################################
# Get the value index in data when the value is the cloest to target value
def get_closest (data, target):
    """Get the index of value in data which is closest to target.
    """
    t = data.index(min(data, key=lambda x:abs(x-target)))
    return t

########################################
### Function: to measure plateau amplitude
########################################
def meas_platamp(data, dt = 0.025):
    """Measures plateau amplitude (average of voltage - baseline while volt
        trace is above baseline + thresh)"""
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    # above = [val for val in stable if val > (baseline + thresh)]
    spike_num = spike_count(stable)
    if spike_num == 0:
        platamp = max(stable) - baseline
        platdur = 0
        ISI = 0
    elif spike_num == 1: # Maybe filtering would be better?
        spikegap = 5 # ms to skip after spike
        idx = data.index(max(data)) + int(spikegap/dt)
        platamp = data[idx] - baseline
        ISI = 0
    else:
        ISI, spike_mvalue, spike_midx = IST_spikes(data, dt)
        platamp = spike_mvalue[-1] - baseline
    return ISI , platamp

########################################
### Function: to measure plateau amplitude for TTX condition
########################################
def TTX_platamp(data):
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    above = [idx for idx, val in enumerate(stable) if val > (baseline+15)]
    # platamp = max(stable) - baseline
    if len(above) <= 15:
        platamp = max(stable) - baseline
    else:
        idx1 = above[0]
        idx2 = above[-1]
        platamp = stable[int(0.75*idx2 + 0.25*idx1)] - baseline
    return platamp

########################################
### Function: to measure the plateau ampltiude and duration
###           in soma and basal dendrites where the inputs are located.
### Using the time stamp for soma plateau amplitude measurement to determine
### the plateau ampltiude in dendrite.
########################################
def soma_platamp_TTX(data):
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    above = [idx for idx, val in enumerate(stable) if val > (baseline+15)]
    # platamp = max(stable) - baseline
    if len(above) <= 15:
        platamp = np.max(stable) - baseline
        idx = np.argmax(stable)
    else:
        idx1 = above[0]
        idx2 = above[-1]
        idx = int(0.75*idx2 + 0.25*idx1)
        platamp = stable[idx] - baseline
    return idx+6000, platamp

def soma_platdur_TTX(data, dt = 0.025):
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    idx, amp = soma_platamp_TTX(data)
    above = [idx for idx, val in enumerate(stable) if val > (baseline+amp/2.0)]
    platdur = dt * len(above)
    return platdur

def TTX_dend_plat(data, idx, dt = 0.025):
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    platamp = data[idx] - baseline
    threshold = platamp*0.5
    above = [val for val in stable if val > (baseline + threshold)]
    platdur = dt * len(above)
    return platamp, platdur

def soma_plat(data, dt = 0.025):
    """Measures plateau amplitude (average of voltage - baseline while volt
        trace is above baseline + thresh)"""
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    spike_num = spike_count(stable)
    if spike_num == 0:
        platamp = max(stable) - baseline
        platdur = 0
        idx = data.index(max(data))
    elif spike_num == 1: # Maybe filtering would be better?
        spikegap = 5 # ms to skip after spike
        idx = data.index(max(data)) + int(spikegap/dt)
        platamp = data[idx] - baseline
    else:
        ISI, spike_mvalue, spike_midx = IST_spikes(data, dt)
        platamp = spike_mvalue[-1] - baseline
        idx = spike_midx[-1]
    return idx, platamp

def dend_plat(data, idx, dt = 0.025):
    """Measures plateau amplitude (average of voltage - baseline while volt
        trace is above baseline + thresh)"""
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    platamp = data[idx] - baseline
    threshold = platamp*0.5
    above = [val for val in stable if val > (baseline + threshold)]
    platdur = dt * len(above)
    return platamp, platdur

def v_curr_inj(data):
    baseline = np.mean(data[4000:6000])
    stable = data[6000:]
    amp = np.max(stable) - baseline
    return amp

#######################################
# Color
#######################################
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
def tableau (num = 0):
    "The number has to be smaller than 21."
    tableau21 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
                 (0, 0, 0)]
    for i in range(len(tableau21)):
        r, g, b = tableau21[i]
        tableau21[i] = (r / 255., g / 255., b / 255.)
    return tableau21[num]
