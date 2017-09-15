''' Analysis module to process battery cycling data from Bio-Logic Tester.
    Assumes files have been made by 'Export as txt'

    Author: Bernard Kim
    Principal Investigators: Paul Wright, James Evans
    Institution: University of California, Berkeley '''

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

    def __init__(self, file):
        ''' Loads file and imports raw data into pandas dataframe
            Separates data by cycle into dictionary
            Values per cycle: time, voltage, current, Qdischarge,
            Qcharge, Edischarge, Echarge'''

        self.filename = file
        self.cycles = {}

        rows = list(csv.reader(open(file), delimiter='\t'))

        headers = ['cycle number', 'time/s', 'Ecell/V', '<I>/mA', 
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

        titlematch = re.search(r'CYCLE.*\d{1}_C', self.filename)
        self.title = titlematch.group(0)[:-2]

    def get_IR_data(self):
        ''' Calculates IR data upon charge, discharge, and rest per cycle 
            Organized in dictionary with keys 'charge', 'discharge', 'rest' 
            Called automatically by 'plot_IR_drop'
            '''

        self.DCIR = {cycle:{} for cycle in self.cycles}

        for idx in range(1, len(self.current)):
            dV = self.voltage[idx] - self.voltage[idx-1]
            dI = self.current[idx] - self.current[idx-1]

            if np.sign(self.current[idx-1]) == -1 and np.sign(self.current[idx]) == 1:
                key = 'charge'
                self.DCIR[int(self.rawcycles[idx])][key] = abs(dV/dI*1000)
            elif np.sign(self.current[idx-1]) == 0 and np.sign(self.current[idx]) == -1:
                key = 'discharge'
                self.DCIR[int(self.rawcycles[idx])][key] = abs(dV/dI*1000)
            elif np.sign(self.current[idx-1]) == 1 and np.sign(self.current[idx]) == 0:
                key = 'rest'
                self.DCIR[int(self.rawcycles[idx])][key] = abs(dV/dI*1000)

    def plot_efficiency(self, xlim=False, ylim=[90, 110], 
        title=None, show=False, save=False, imagetype='png'):
        ''' Plots charging efficiency vs cycle '''

        cycles = list(self.cycles.keys())[1:]

        efficiency = [np.max(self.cycles[cycle]['Qdischarge'])/np.max(self.cycles[cycle]['Qcharge'])*100 
            for cycle in cycles]

        # Plotting formatting
        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        ax.plot(cycles[:-1], efficiency[:-1], color='b', linewidth=2,)

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_Efficiency')

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('Coulombic Efficiency (%)')
        ax.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_Efficiency' + '.' + imagetype)

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
        cycles = list(self.cycles.keys())[1:]

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

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_Capacity' + '.'+ imagetype)

        plt.close(fig)

    def plot_capacity_fade(self, xlim=False, ylim=[0,100], title=None, 
        show=False, save=False, imagetype='png'):
        ''' Plots capacity fade vs. cycle
            Similar to plotting capacity vs. cycle but normalized to 
            first cycle capacity '''

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting
        cycles = list(self.cycles.keys())[1:]

        baseline = np.max(self.cycles[cycles[0]]['Qdischarge'])
        fade = [np.max(self.cycles[cycle]['Qdischarge'])/baseline*100 for cycle in cycles]

        # Plotting formatting
        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        ax.plot(cycles[:-1], fade[:-1], color='b', linewidth=2)

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_Fade')

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('Capacity Fade (%)')
        ax.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_Fade' + '.' + imagetype)

        plt.close(fig)

    def plot_IR_drop(self, IRtypes=['charge', 'discharge', 'rest', 'average'], 
        xlim=None, ylim=None, semilogy=True, title=None, 
        show=False, save=False, imagetype='png'):
        ''' Plots extracted IR data vs cycle number
            IRtype accepts kwargs 'charge', 'discharge', 'rest', or 'average' as list
            e.g. ['charge', 'discharge',]
            ''' 

        self.get_IR_data()

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)
        plotcolors = {'charge':'r', 'discharge':'g', 'rest':'b', 'average':'k'}

        for IRtype in IRtypes:
            if IRtype == 'average':
                plotcycles = [cycle for cycle in self.DCIR]
                plotIR = [np.mean(list(self.DCIR[cycle].values())) for cycle in plotcycles]

                if semilogy:
                    ax.semilogy(plotcycles[1:-1], plotIR[1:-1], color=plotcolors[IRtype], linewidth=2, label=IRtype)
                else:
                    ax.plot(plotcycles[1:-1], plotIR[1:-1], color=plotcolors[IRtype], linewidth=2, label=IRtype)

            else:
                plotcycles = []

                for cycle in self.DCIR:
                    if IRtype in self.DCIR[cycle]:
                        plotcycles.append(cycle)

                plotIR = [self.DCIR[cycle][IRtype] for cycle in plotcycles]

                if semilogy:
                    ax.semilogy(plotcycles, plotIR, color=plotcolors[IRtype], linewidth=2, label=IRtype)
                else:
                    ax.plot(plotcycles, plotIR, color=plotcolors[IRtype], linewidth=2, label=IRtype)

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_DCIR')

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DC Internal Resistance [Ω]')
        ax.legend()
        ax.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_DCIR' + '.' + imagetype)

        plt.close(fig)

    def plot_DOD(self, plotcycleinterval=20, charge=True, discharge=True, xlim=[0,1], ylim=None, 
        title=None, show=False, save=False, imagetype='png'):
        ''' Plots voltage vs DOD 
            Max Qdischarge is on the y-axis (decreasing to the right)
            plotcycleinterval = cycle interval to include on plot
            '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

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

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_DOD')

        ax.set_xlabel('Capacity [mAh/cm' + r'$^2$' + ']')
        ax.set_ylabel('Voltage [V]')
        ax.set_title(self.title + '_DOD')
        ax.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_DOD' + '.' + imagetype)

        plt.close(fig)

class GCPL5:
    ''' Processes cycling data gathered with GCPL5 profile
        (Galvanostatic Cycling with Potential Limitation)
        Increases sampling frequency immediately after current change
        For pulsed discharging of cells '''


    def __init__(self, file):
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

        titlematch = re.search(r'PULSE.*\d{1}_C', self.filename)
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
        ''' Plots extracted IR data vs cycle number
            IRtype accepts kwargs 'charge', 'discharge', 'rest', or 'average' as list
            e.g. ['charge', 'discharge',]
            ''' 

        self.get_IR_data()

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)
        plotcolors = {'charge':'r', 'discharge':'g', 'rest':'b', 'average':'k'}

        plotcycles = self.DCIR['cycles']

        if average:
            plotIR = self.DCIR['average']

            ax.plot(plotcycles, plotIR, color='b', label='Average')
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

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_DCIR' + plottype + '.' + imagetype)

        plt.close(fig)












