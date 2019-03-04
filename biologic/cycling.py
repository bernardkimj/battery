''' Analysis module to process battery cycling data from Bio-Logic Tester.
    Assumes files have been made by 'Export as txt'

    Author: Bernard Kim
    Principal Investigators: Paul Wright, James Evans
    Institution: University of California, Berkeley '''

from battery.utilities import utilities
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import csv

class GCPL:
    ''' Processes cycling data made with GCPL profile 
        (Galvanostatic Cycling with Potential Limitation)
        Standard CC charge/discharge profile '''

    def __init__(self, file, title=None):
        ''' Loads file and imports raw data into pandas dataframe
            Separates data by cycle into dictionary
            Values per cycle: time, voltage, current, Qdischarge,
            Qcharge, Edischarge, Echarge'''

        self.filename = file
        self.cycles = {}

        rows = list(csv.reader(open(file), delimiter='\t'))

        # All possible headers
        headers = {
            # 'cycle number': 'cycle',
            'time/s': 'time', 
            'Ecell/V': 'voltage', 
            '<I>/mA': 'current', 
            'I/mA': 'current',
            'Q discharge/mA.h': 'Qdischarge',
            'Q charge/mA.h': 'Qcharge', 
            'Energy discharge/W.h': 'Edischarge', 
            'Energy charge/W.h': 'Echarge',
        }

        # Instantiate dict for indices within row (agnostic to export order)
        idxs = {}

        # Get indices for values that were exported
        for header in headers:
            if header in rows[0]:
                idxs[headers[header]] = rows[0].index(header)

        if 'cycle number' in rows[0]:
            cycleidx = rows[0].index('cycle number')
        else:
            cycleidx = None

        self.time = []
        self.voltage = []
        self.current = []
        self.rawcycles = []
        self.cycles = {}

        # Iterate through rows and populate raw lists for time, voltage, current
        # Also fill in dictionary for cycles by creating new key if new cycle 
        # and filling in existing headers for each entry
        for row in rows[1:]:
            self.time.append(float(row[idxs['time']]))
            self.voltage.append(float(row[idxs['voltage']]))
            self.current.append(float(row[idxs['current']]))

            # Check if cycle numbers were exported
            if cycleidx:
                rawcycle = float(row[cycleidx])
                cycle = int(rawcycle)
                self.rawcycles.append(rawcycle)

                # Check if cycle number key has been added to dictionary
                # If not, create it
                if cycle not in self.cycles:
                    self.cycles[cycle] = {}

                # For current row, populate cycle dictionary with existing
                # entries
                for entry in idxs:
                    # Check if entry is in cycle dict, add if not
                    if entry not in self.cycles[cycle]:
                        self.cycles[cycle][entry] = []

                    # Append entry value to appropriate list for given cycle
                    self.cycles[cycle][entry].append(float(row[idxs[entry]]))

        if title:
            self.title = title
        else:
            self.title = file[:-8]
            # titlematch = re.search(r'CELL.*\d{1}_C', self.filename)
            # self.title = titlematch.group(0)[:-2]


    def get_IR_data(self):
        ''' Calculates IR data upon charge, discharge, and rest per cycle 
            Organized in dictionary with keys 'charge', 'discharge', 'rest' 
            Called automatically by 'plot_IR_drop'
            Robust to any number of rest, charge, or discharge steps
            '''

        # Initialize dictionary for raw data
        self.DCIR = {
            cycle: {
                step: {
                    'all': [],
                    'average': None,
                } for step in ['charge', 'discharge', 'rest']
            } for cycle in self.cycles
        }

        # Determine type of step when current changes signs and calculate 
        # instantaneous DC internal resistance
        for idx in range(1, len(self.current)):
            dV = self.voltage[idx] - self.voltage[idx-1]
            dI = self.current[idx] - self.current[idx-1]

            if dI != 0:
                previous = self.current[idx-1]
                present = self.current[idx]

                if np.sign(present) == 1 and np.sign(previous) != 1:
                    key = 'charge'
                elif np.sign(present) == -1 and np.sign(previous) != -1:
                    key = 'discharge'
                elif np.sign(present) == 0 and np.sign(previous) != 0:
                    key = 'rest'

                resistance = np.abs(dV/(dI/1000))

                self.DCIR[int(self.rawcycles[idx])][key]['all'].append(resistance)

        # Fill in average values for all steps for all cycles
        for cycle in self.DCIR:
            for step in self.DCIR[cycle]:
                self.DCIR[cycle][step]['average'] = np.mean(
                    self.DCIR[cycle][step]['all'])


    def plot_setup(self, params):

        # Plotting formatting
        font = {'family': 'Arial', 'size': 20}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if params['xlim']:
            ax.set_xlim(params['xlim'])
        if params['ylim']:
            ax.set_ylim(params['ylim'])

        if params['title']:
            ax.set_title(params['title'])
        else:
            title = self.title + params['titletag']
            ax.set_title(title)

        ax.set_xlabel(params['xlabel'])
        ax.set_ylabel(params['ylabel'])
        ax.grid()

        return fig, ax


    def plot_voltage(self, xlim=None, ylim=None, title=None,
        show=False, save=False):
        ''' Plots voltage vs time '''

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_voltage',
            'xlabel': 'Time [s]',
            'ylabel': 'Voltage [V]'
        }

        fig, ax = self.plot_setup(params=params)

        ax.plot(self.time, self.voltage, color='b', linewidth=2,)

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_efficiency(self, xlim=False, ylim=[90, 110], title=None, 
        show=False, save=False):
        ''' Plots charging efficiency vs cycle '''

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_Efficiency',
            'xlabel': 'Cycle Number',
            'ylabel': 'Coulombic Efficiency [%]'
        }

        cycles = list(self.cycles.keys())[1:]
        efficiency = [np.max(self.cycles[cycle]['Qdischarge'])/np.max(
            self.cycles[cycle]['Qcharge'])*100 for cycle in cycles]

        fig, ax = self.plot_setup(params=params)
        ax.plot(cycles[:-1], efficiency[:-1], color='b', linewidth=2,)

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_capacity(self, discharge=True, charge=True, xlim=False, ylim=[0,1], 
        title=None, show=False, save=False,):
        ''' Plots charge/discharge capacity vs cycle '''

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_Capacity',
            'xlabel': 'Cycle Number',
            'ylabel': 'Capacity ' + r'$[mAh]$'
        }

        fig, ax = self.plot_setup(params=params)

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting
        cycles = list(self.cycles.keys())[1:]

        # Plot both discharge and charge by default
        if discharge:
            Qdischarge = [np.max(self.cycles[cycle]['Qdischarge']) for cycle in cycles]
            ax.plot(cycles[:-1], Qdischarge[:-1], color='b', linewidth=2, label=r'$Q_{discharge}$')
        if charge:
            plotcycles = [cycle-1 for cycle in cycles] # To make charge curve line up with discharge
            Qcharge = [np.max(self.cycles[cycle]['Qcharge']) for cycle in cycles]
            ax.plot(plotcycles, Qcharge, color='r', linewidth=2, label=r'$Q_{charge}$')

        ax.legend()

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_capacity_fade(self, xlim=False, ylim=[0,100], title=None, 
        show=False, save=False,):
        ''' Plots capacity fade vs. cycle
            Similar to plotting capacity vs. cycle but normalized to 
            first cycle capacity '''

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_Fade',
            'xlabel': 'Cycle Number',
            'ylabel': 'Capacity Fade [%]'
        }

        fig, ax = self.plot_setup(params=params)

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting
        cycles = list(self.cycles.keys())[1:]

        baseline = np.max(self.cycles[cycles[0]]['Qdischarge'])
        fade = [np.max(self.cycles[cycle]['Qdischarge'])/baseline*100 for cycle in cycles]

        ax.plot(cycles[:-1], fade[:-1], color='b', linewidth=2)

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_IR_drop(self, steps=['charge', 'discharge', 'rest'], 
        average=True, xlim=None, ylim=None, semilogy=True, 
        title=None, show=False, save=False,):
        ''' Plots extracted IR data vs cycle number
            IRtype accepts kwargs 'charge', 'discharge', 'rest', or 'average' as list
            e.g. ['charge', 'discharge',]
            ''' 

        self.get_IR_data()

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_DCIR',
            'xlabel': 'Cycle Number',
            'ylabel': 'DC Internal Resistance [Ω]'
        }

        fig, ax = self.plot_setup(params=params)

        plotcolors = {'charge':'r', 'discharge':'g', 'rest':'b'}
        plotcycles = list(self.DCIR.keys())

        for step in steps:
            # Plots averaged DCIR per step (if multiples of each step per cycle)
            if average:
                plotIR = [self.DCIR[cycle][step]['average'] 
                    for cycle in self.DCIR]

                if semilogy:
                    ax.semilogy(plotcycles, plotIR,
                        color=plotcolors[step], linewidth=2, label=step)
                else:
                    ax.plot(plotcycles, plotIR,
                        color=plotcolors[step], linewidth=2, label=step)

            else:
                # Plots all raw IR values (multiple steps plotted individually)
                all_IR = np.array([self.DCIR[cycle][step]['all'] 
                    for cycle in self.DCIR]).T

                alphas = np.linspace(1, 0.6, len(all_IR))

                for idp, pulse in enumerate(all_IR):
                    if semilogy:
                        if idp == 0:
                            ax.semilogy(plotcycles, pulse,
                                color=plotcolors[step], alpha=alphas[idp],
                                linewidth=2, label=step)
                        else:
                            ax.semilogy(plotcycles, pulse,
                                color=plotcolors[step], alpha=alphas[idp],
                                linewidth=2)
                    else:
                        if idp == 0:
                            ax.plot(plotcycles, pulse,
                                color=plotcolors[step], alpha=alphas[idp],
                                linewidth=2, label=step)
                        else:
                            ax.plot(plotcycles, pulse,
                                color=plotcolors[step], alpha=alphas[idp],
                                linewidth=2)

        ax.legend()

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_DOD(self, plotcycleinterval=20, charge=True, discharge=True, 
        xlim=[0,1], ylim=None, title=None, show=False, save=False,):
        ''' Plots voltage vs DOD 
            Max Qdischarge is on the y-axis (decreasing to the right)
            plotcycleinterval = cycle interval to include on plot
            '''

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_DOD',
            'xlabel': 'Capacity [mAh/cm' + r'$^2$' + ']',
            'ylabel': 'Voltage [V]',
        }

        fig, ax = self.plot_setup(params=params)

        if charge:
            DOC = {}

            for cycle in list(self.cycles.keys())[1:]:
                if np.remainder(cycle, plotcycleinterval) == 0:
                    DOC[cycle] = {'voltage':[], 'Qcharge':[]}

            # Remove last element of voltage and Qdischarge lists to avoid "tail" in graphs
            # Artifact of splitting up by cycle/tester resolution
            for cycle in DOC:
                for voltage, current, Qcharge in \
                    zip(self.cycles[cycle]['voltage'], \
                        self.cycles[cycle]['current'], \
                        self.cycles[cycle]['Qcharge']):

                    if Qcharge != 0 and current > 0:
                        DOC[cycle]['voltage'].append(voltage)
                        DOC[cycle]['Qcharge'].append(Qcharge)

            coloridx = np.linspace(0.25, 1, len(DOC))

            for idx, cycle in enumerate(DOC):
                ax.plot(
                    DOC[cycle]['Qcharge'], DOC[cycle]['voltage'],
                    linewidth=3, color=plt.cm.Blues(coloridx[idx]), label=cycle)

        if discharge:
            DOD = {}

            for cycle in list(self.cycles.keys())[1:]:
                if np.remainder(cycle, plotcycleinterval) == 0:
                    DOD[cycle] = {'voltage':[], 'Qdischarge':[]}

            # Remove last element of voltage and Qdischarge lists to avoid "tail" in graphs
            # Artifact of splitting up by cycle/tester resolution
            for cycle in DOD:
                for voltage, current, Qdischarge in \
                    zip(self.cycles[cycle]['voltage'], \
                        self.cycles[cycle]['current'], \
                        self.cycles[cycle]['Qdischarge']):

                    if Qdischarge != 0 and current < 0:
                        DOD[cycle]['voltage'].append(voltage)
                        DOD[cycle]['Qdischarge'].append(Qdischarge)

            coloridx = np.linspace(0.25, 1, len(DOD))

            for idx, cycle in enumerate(DOD):
                ax.plot(
                    DOD[cycle]['Qdischarge'][1:-1], DOD[cycle]['voltage'][1:-1],
                    linewidth=2, color=plt.cm.Blues(coloridx[idx]), label=cycle)

        handles, labels = fig.gca().get_legend_handles_labels()
        i=1
        while i < len(labels):
            if labels[i] in labels[:i]:
                del(labels[i])
                del(handles[i])
            else:
                i+=1
        plotlabels = ['Cycle ' +str(label) for label in labels]
        ax.legend(handles, plotlabels)

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


class GCPL5:
    ''' Processes cycling data gathered with GCPL5 profile
        (Galvanostatic Cycling with Potential Limitation)
        Increases sampling frequency immediately after current change
        For pulsed discharging of cells '''


    def __init__(self, file, title=None):
        ''' Loads file and imports raw data into pandas dataframe
            Separates data by cycle into dictionary
            Values per cycle: time, voltage, current, Qdischarge,
            Qcharge, Edischarge, Echarge'''

        self.filename = file
        self.cycles = {}

        rows = list(csv.reader(open(file), delimiter='\t'))

        headers = ['cycle number', 'time/s', 'Ecell/V', 'I/mA', 
            'Q discharge/mA.h','Q charge/mA.h', 'Energy discharge/W.h', 'Energy charge/W.h']
        idxs = []

        for header in headers:
            idxs.append(rows[0].index(header))

        cycleidx = idxs[0]
        timeidx = idxs[1]
        voltageidx = idxs[2]
        currentidx = idxs[3]
        Qdischargeidx = idxs[4]
        Qchargeidx = idxs[5]
        Edischargeidx = idxs[6]
        Echargeidx = idxs[7]

        self.rawcycles = []
        for row in rows[1:]:
            self.rawcycles.append(float(row[cycleidx]))

        cyclenumbers = []
        for entry in set(self.rawcycles):
            cyclenumbers.append(int(entry))
        cyclenumbers.sort()

        for cyclenumber in cyclenumbers:
            self.cycles[cyclenumber] = {
                'time': [],
                'voltage': [],
                'current': [],
                'Qdischarge': [],
                'Qcharge': [],
                'Edischarge': [],
                'Echarge': [],
            }

        self.voltage, self.current = [], []

        for row in rows[1:]:
            self.cycles[int(float(row[cycleidx]))]['time'].append(float(row[timeidx]))
            self.cycles[int(float(row[cycleidx]))]['voltage'].append(float(row[voltageidx]))
            self.cycles[int(float(row[cycleidx]))]['current'].append(float(row[currentidx]))
            self.cycles[int(float(row[cycleidx]))]['Qdischarge'].append(float(row[Qdischargeidx]))
            self.cycles[int(float(row[cycleidx]))]['Qcharge'].append(float(row[Qchargeidx]))
            self.cycles[int(float(row[cycleidx]))]['Edischarge'].append(float(row[Edischargeidx]))
            self.cycles[int(float(row[cycleidx]))]['Echarge'].append(float(row[Echargeidx]))

            self.voltage.append(float(row[voltageidx]))
            self.current.append(float(row[currentidx]))

        self.get_IR_data()

        if title:
            self.title = title
        else:
            titlematch = re.search(r'CELL.*\d{1}_C', self.filename)
            self.title = titlematch.group(0)[:-2]

    def get_IR_data(self):
        ''' Calculates IR data upon charge, discharge, and rest per cycle 
            Organized in dictionary with keys 'charge', 'discharge', 'rest' 
            Called automatically by 'plot_IR_drop'
            '''

        DCIR = {cycle:[] for cycle in self.cycles}

        for idx in range(1, len(self.current)):
            if np.sign(self.current[idx-1]) == 0 and np.sign(self.current[idx]) == -1:
                dV = self.voltage[idx] - self.voltage[idx-1]
                dI = self.current[idx] - self.current[idx-1]
                DCIR[int(self.rawcycles[idx])].append(abs(dV/dI*1000))

        IRarray = [DCIR[cycle] for cycle in DCIR]

        # Omit last cycle if unfinished
        if len(IRarray[0]) != len(IRarray[-1]):
            IRarray.pop()
            DCIR.pop(np.max(list(DCIR.keys())))

        IRarray = np.array(IRarray)

        DCIR_avg = []
        for cycle in DCIR:
            DCIR_avg.append(np.mean(DCIR[cycle]))

        self.DCIR = {
            'cycles': list(DCIR.keys()),
            'array': IRarray,
            'average': DCIR_avg
        }

    def plot_IR_drop(self, average=True, xlim=None, ylim=None, 
        title=None, show=False, save=False, imagetype='png'):
        ''' Plots extracted IR data vs cycle number ''' 

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)
        plotcolors = {'charge':'r', 'discharge':'g', 'rest':'b', 'average':'k'}

        if average:
            ax.plot(self.DCIR['cycles'], self.DCIR['average'], 
                color='b', label='Average')
            plottype = '_average'

        else:
            coloridx = np.linspace(0.4, 1, self.DCIR['array'].shape[1])
            for idx, pulse in enumerate(self.DCIR['array'].T, 1):
                ax.plot(plotcycles, pulse, label='Pulse ' + str(idx), 
                    color=plt.cm.Purples(coloridx[idx-1]), linewidth=2)
            plottype = '_all'

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + plottype)

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DC Internal Resistance [Ω]')
        ax.legend()
        ax.grid()

        if save:
            plt.savefig(self.title + '_DCIR' + plottype + '.' + imagetype)
        if show:
            plt.show()

        plt.close(fig)

    def plot_capacity(self, discharge=True, charge=True, xlim=False, ylim=[0,1], 
        title=None, show=False, save=False, imagetype='png'):
        ''' Plots charge/discharge capacity vs cycle '''

        # Plotting formatting
        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting

        # Plot both discharge and charge by default
        if discharge:
            Qdischarge = [np.max(self.cycles[cycle]['Qdischarge']) for cycle in cycles]
            ax.plot(cycles[:-1], Qdischarge[:-1], color='b', linewidth=2, label=r'$Q_{discharge}$')
        if charge:
            Qcharge = [np.max(self.cycles[cycle]['Qcharge']) for cycle in cycles]
            ax.plot(cycles, Qcharge, color='r', linewidth=2, label=r'$Q_{charge}$')

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_Capacity')

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('Capacity ' + r'$[mAh/cm^2]$')
        ax.legend()
        ax.grid()

        if save:
            plt.savefig(self.title + '_Capacity' + '.'+ imagetype)
        if show:
            plt.show()

        plt.close(fig)


class GCPL_batch:
    ''' Plots batches of samples together
        Uses defined methods from GCPL class '''

    def __init__(self, files):
        ''' Initializes and organizes files 
            Dictionary of GCPL5 objects with overarching title '''

        self.alldata = {}

        for file in files:
            data = GCPL(file)
            samplesearch = re.search(r'_S\d{1,2}', data.filename)
            sample = samplesearch.group(0)[2:]
            self.alldata[int(sample)] = data

        titlesearch = re.search(r'CELL.*\d{8}', \
            self.alldata[list(self.alldata.keys())[0]].title)
        self.title = titlesearch.group(0)


    def plot_capacity(self, discharge=True, charge=True, xlim=False, ylim=[0,1], 
        title=None, show=False, save=False,):
        ''' Plots charge/discharge capacity vs cycle '''

        params = {
            'xlim': xlim,
            'ylim': ylim,
            'title': title,
            'titletag': '_Capacity',
            'xlabel': 'Cycle Number',
            'ylabel': 'Capacity ' + r'$[mAh]$'
        }

        fig, ax = self.plot_setup(params=params)

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting
        cycles = list(self.cycles.keys())[1:]

        # Plot both discharge and charge by default
        if discharge:
            Qdischarge = [np.max(self.cycles[cycle]['Qdischarge']) for cycle in cycles]
            ax.plot(cycles[:-1], Qdischarge[:-1], color='b', linewidth=2, label=r'$Q_{discharge}$')
        if charge:
            plotcycles = [cycle-1 for cycle in cycles] # To make charge curve line up with discharge
            Qcharge = [np.max(self.cycles[cycle]['Qcharge']) for cycle in cycles]
            ax.plot(plotcycles, Qcharge, color='r', linewidth=2, label=r'$Q_{charge}$')

        ax.legend()

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_DCIR_discrete(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, savename=None, imagetype='png'):
        ''' Plots discrete curves for DCIR vs cycle number '''

        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)
        coloridx = np.linspace(0,1,10) # For use with tab10 colormap
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for idx, sample in enumerate(sorted(list(self.alldata.keys()))):
            ax.plot(self.alldata[sample].DCIR['cycles'], self.alldata[sample].DCIR['average'],
                color=plt.cm.tab10(coloridx[idx]), linewidth=3,
                label='S'+str(sample))

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DCIR [Ω]')
        ax.legend()
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        if show:
            plt.show()

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('batch_' + self.title + '_DCIR_discrete' + '.' + imagetype)

        plt.close(fig)

    def plot_DCIR_average(self, confidence=0.95, xlim=None, ylim=None, title=None, 
        show=False, save=False, savename=None, imagetype='png'):
        ''' Plots average DCIR vs cycle number
            Includes confidence interval '''

        # Combines all sample data into list of tuples (x,y) by sample number (len = # of samples)
        data = [(self.alldata[sample].DCIR['cycles'], self.alldata[sample].DCIR['average']) \
            for sample in sorted(list(self.alldata.keys()))]

        cycles, mean, std, lcl, ucl = utilities.batch_average_plot(data, confidence=confidence)

        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for sample in data:
            ax.plot(sample[0], sample[1], color='b', linewidth=3, alpha=0.2)
        ax.plot(cycles, mean, color='b', linewidth=3)

        if confidence:
            ax.plot(cycles, lcl, color='b', linestyle='--', linewidth=2, 
                label=str(confidence)[2:]+'% Confidence Interval (t-test)')
            ax.plot(cycles, ucl, color='b', linestyle='--', linewidth=2)
            ax.fill_between(cycles, ucl, lcl, color='b', alpha=0.2)
        else:
            ax.plot(cycles, mean+std, color='b', linestyle='--', linewidth=2, 
                label='±1 '+r'$\sigma$')
            ax.plot(cycles, mean-std, color='b', linestyle='--', linewidth=2)
            ax.fill_between(cycles, mean+std, mean-std, color='b', alpha=0.2)

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DCIR [Ω]')
        ax.legend()
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('batch_' + self.title + '_DCIR_average' + '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)


class GCPL5_batch:
    ''' Plots batches of samples together
        Uses defined methods from GCPL5 class '''

    def __init__(self, files):
        ''' Initializes and organizes files 
            Dictionary of GCPL5 objects with overarching title '''

        self.alldata = {}

        for file in files:
            data = GCPL5(file)
            samplesearch = re.search(r'_S\d{1,2}', data.filename)
            sample = samplesearch.group(0)[2:]
            self.alldata[int(sample)] = data

        titlesearch = re.search(r'CELL.*\d{8}', \
            self.alldata[list(self.alldata.keys())[0]].title)
        self.title = titlesearch.group(0)

    def plot_DCIR_discrete(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, savename=None, imagetype='png'):
        ''' Plots discrete curves for DCIR vs cycle number '''

        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)
        coloridx = np.linspace(0,1,10) # For use with tab10 colormap
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for idx, sample in enumerate(sorted(list(self.alldata.keys()))):
            ax.plot(self.alldata[sample].DCIR['cycles'], self.alldata[sample].DCIR['average'],
                color=plt.cm.tab10(coloridx[idx]), linewidth=3,
                label='S'+str(sample))

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DCIR [Ω]')
        ax.legend()
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        if show:
            plt.show()

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('batch_' + self.title + '_DCIR_discrete' + '.' + imagetype)

        plt.close(fig)

    def plot_DCIR_average(self, confidence=0.95, xlim=None, ylim=None, title=None, 
        show=False, save=False, savename=None, imagetype='png'):
        ''' Plots average DCIR vs cycle number
            Includes confidence interval '''

        # Combines all sample data into list of tuples (x,y) by sample number (len = # of samples)
        data = [(self.alldata[sample].DCIR['cycles'], self.alldata[sample].DCIR['average']) \
            for sample in sorted(list(self.alldata.keys()))]

        cycles, mean, std, lcl, ucl = utilities.batch_average_plot(data, confidence=confidence)

        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for sample in data:
            ax.plot(sample[0], sample[1], color='b', linewidth=3, alpha=0.2)
        ax.plot(cycles, mean, color='b', linewidth=3)

        if confidence:
            ax.plot(cycles, lcl, color='b', linestyle='--', linewidth=2, 
                label=str(confidence)[2:]+'% Confidence Interval (t-test)')
            ax.plot(cycles, ucl, color='b', linestyle='--', linewidth=2)
            ax.fill_between(cycles, ucl, lcl, color='b', alpha=0.2)
        else:
            ax.plot(cycles, mean+std, color='b', linestyle='--', linewidth=2, 
                label='±1 '+r'$\sigma$')
            ax.plot(cycles, mean-std, color='b', linestyle='--', linewidth=2)
            ax.fill_between(cycles, mean+std, mean-std, color='b', alpha=0.2)

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DCIR [Ω]')
        ax.legend()
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('batch_' + self.title + '_DCIR_average' + '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)



















