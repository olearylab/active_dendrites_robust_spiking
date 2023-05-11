# Active Dendrites Enable Robust Spiking Computations Despite Timing Jitter

Burger, T.S.J., Rule, M.E., O'Leary, T. (2023). Active Dendrites Enable Robust Spiking Computations Despite Timing Jitter. [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.03.22.533815v1), 2023-03. [doi: https://doi.org/10.1101/2023.03.22.533815](https://doi.org/10.1101/2023.03.22.533815).

### Abstract

Dendritic action potentials exhibit long plateaus of many tens of milliseconds, outliving axonal spikes by an order of magnitude. The computational role of these slow events seems at odds with any need to rapidly integrate and relay information throughout large nervous systems. We propose that the timescale of dendritic potentials allows reliable integration of asynchronous inputs. We develop a physiologically grounded model in which the extended duration of dendritic spikes equips each dendrite with a resettable memory of incoming signals. This provides a tractable model for capturing dendritic nonlinearities observed in experiments and in more complex, detailed models. Using this model, we show that long-lived, nonlinear dendritic plateau potentials allow reliable integration of asynchronous spikes. We demonstrate this model supports non-trivial computations in a network solving an arbitrary association/discrimination task using sparse spiking that is subject to timing jitter. This demonstrates a computational role for the specific timecourse of dendritic potentials in situations where decisions occur quickly, reliably, and with a low number of spikes. Our results provide empirically testable hypotheses for the role of dendritic action potentials in cortical function as well as a potential bio-inspired means of realising neuromorphic spiking computations in analog hardware.

### Repository Structure

The code consists of two folders: 

- `plateau-potentials/`: Contains Python code to run the detailed biophysical simulation of Fig. 2 in NEURON.
- `DendriteNetwork/`: Contains Julia code for the abstract model of Fig. 2 and the network simulation of Fig. 4.
