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

class GCPL:
    ''' Processes cycling data made with GCPL profile 
        (Galvanostatic Cycling with Potential Limitation)
        Standard CC charge/discharge profile '''

    def __init__(self, file):
        ''' Loads file and imports raw data into pandas dataframe
            Separates data by cycle into dictionary
            Values per cycle: time, voltage, current, Qdischarge,
            Qcharge, Edischarge, Echarge'''

        df = pd.read_csv(file, sep="\t")

        self.filename = file

        self.time = df['time/s']
        self.rawcyclenumbers = df['cycle number']
        self.voltage = df['Ecell/V']
        self.current = df['<I>/mA']
        self.Qdischarge = df['Q discharge/mA.h']
        self.Qcharge = df['Q charge/mA.h']
        self.Edischarge = df['Energy discharge/W.h']
        self.Echarge = df['Energy charge/W.h']

        self.cyclenumbers = []

        # Create set of cycle numbers to pull out individual cycle numbers
        # Sort list so cycle numbers are in order
        for entry in set(self.rawcyclenumbers):
            self.cyclenumbers.append(int(entry))
        self.cyclenumbers.sort()

        # Populate dictionary 'cycles' to split up raw data by cycle
        self.cycles = {}
        for cyclenumber in self.cyclenumbers:
            idx = np.where(np.array(self.rawcyclenumbers) == cyclenumber)

            self.cycles[cyclenumber] = {
                'time': [self.time[idv] for idv in idx],
                'voltage': [self.voltage[idv] for idv in idx],
                'current': [self.current[idv] for idv in idx],
                'Qdischarge': [self.Qdischarge[idv] for idv in idx],
                'Qcharge': [self.Qcharge[idv] for idv in idx],
                'Edischarge': [self.Edischarge[idv] for idv in idx],
                'Echarge': [self.Echarge[idv] for idv in idx],
            }

    def plot_capacity(self, discharge=True, charge=True, xlim=False, ylim=False, 
        show=False, save=False, imagetype='pdf'):
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
            ax.plot(cycles, Qdischarge, color='b', linewidth=3, label=r'$Q_{discharge}$')
        if charge:
            Qcharge = [np.max(self.cycles[cycle]['Qcharge']) for cycle in cycles]
            ax.plot(cycles, Qcharge, color='r', linewidth=3, label=r'$Q_{charge}$')

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('Capacity ' + r'$[mAh/cm^2]$')
        ax.set_title('Capacity vs. Cycle')
        ax.legend()
        ax.grid()

        if show:
            plt.show()
            plt.close()

        if save:
            plt.savefig(self.filename[:-8] + '_cycle' + '.'+ imagetype)

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
                self.DCIR[int(self.rawcyclenumbers[idx])][key] = abs(dV/dI*1000)
            elif np.sign(self.current[idx-1]) == 0 and np.sign(self.current[idx]) == -1:
                key = 'discharge'
                self.DCIR[int(self.rawcyclenumbers[idx])][key] = abs(dV/dI*1000)
            elif np.sign(self.current[idx-1]) == 1 and np.sign(self.current[idx]) == 0:
                key = 'rest'
                self.DCIR[int(self.rawcyclenumbers[idx])][key] = abs(dV/dI*1000)

    def plot_IR_drop(self, IRtypes=None, xlim=None, ylim=None, semilogy=False, 
        show=False, save=False, imagetype='pdf'):
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
                    ax.semilogy(plotcycles, plotIR, color=plotcolors[IRtype], linewidth=3, label=IRtype)
                else:
                    ax.plot(plotcycles, plotIR, color=plotcolors[IRtype], linewidth=3, label=IRtype)

            else:
                plotcycles = []

                for cycle in self.DCIR:
                    if IRtype in self.DCIR[cycle]:
                        plotcycles.append(cycle)

                plotIR = [self.DCIR[cycle][IRtype] for cycle in plotcycles]

                if semilogy:
                    ax.semilogy(plotcycles, plotIR, color=plotcolors[IRtype], linewidth=3, label=IRtype)
                else:
                    ax.plot(plotcycles, plotIR, color=plotcolors[IRtype], linewidth=3, label=IRtype)

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DC Internal Resistance [Ω]')
        ax.set_title(self.filename[:-8] + '_DCIR')
        ax.legend()
        ax.grid()

        if show:
            plt.show()
            plt.close()

        if save:
            plt.savefig(self.filename[:-8] + '_DCIR_' + '_'.join(IRtypes) + '.' + imagetype)

    def plot_DOD(self, plotcycleinterval=10, ylim=None, show=False, save=False, imagetype='pdf'):
        ''' Plots voltage vs DOD 
            Max Qdischarge is on the y-axis (decreasing to the right)
            plotcycleinterval = cycle interval to include on plot
            '''

        DOD = {}

        for cycle in list(self.cycles.keys())[1:]:
            if np.remainder(cycle, plotcycleinterval) == 0:
                DOD[cycle] = {'voltage':[], 'Qdischarge':[]}

        # Remove last element of voltage and Qdischarge lists to avoid "tail" in graphs
        # Artifact of splitting up by cycle/tester resolution
        for cycle in DOD:
            for voltage, Qdischarge in \
                zip(self.cycles[cycle]['voltage'][0][:-1], self.cycles[cycle]['Qdischarge'][0][:-1]):

                if Qdischarge != 0:
                    DOD[cycle]['voltage'].append(voltage)
                    DOD[cycle]['Qdischarge'].append(Qdischarge)

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        coloridx = np.linspace(0.25, 1, len(DOD))

        for idx, cycle in enumerate(DOD):
            ax.plot(
                DOD[cycle]['Qdischarge'], DOD[cycle]['voltage'],
                linewidth=3, color=plt.cm.Blues(coloridx[idx]), label='Cycle '+str(cycle))

        if ylim:
            ax.set_ylim(ylim)

        ax.set_xlabel('Discharge Capacity [mAh]')
        ax.set_ylabel('Voltage [V]')
        ax.set_title(self.filename[:-8] + '_DOD')
        ax.legend()
        ax.grid()

        if show:
            plt.show()
            plt.close()

        if save:
            plt.savefig(self.filename[:-8] + '_DoD' + '.' + imagetype)

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

        df = pd.read_csv(file, sep="\t")

        self.filename = file

        self.time = df['time/s']
        self.rawcyclenumbers = df['cycle number']
        self.voltage = df['Ecell/V']
        self.current = df['I/mA']
        self.Qdischarge = df['Q discharge/mA.h']
        self.Qcharge = df['Q charge/mA.h']
        self.Edischarge = df['Energy discharge/W.h']
        self.Echarge = df['Energy charge/W.h']

        self.cyclenumbers = []

        # Create set of cycle numbers to pull out individual cycle numbers
        # Sort list so cycle numbers are in order
        for entry in set(self.rawcyclenumbers):
            self.cyclenumbers.append(int(entry))
        self.cyclenumbers.sort()

        # Populate dictionary 'cycles' to split up raw data by cycle
        self.cycles = {}
        for cyclenumber in self.cyclenumbers:
            idx = np.where(np.array(self.rawcyclenumbers) == cyclenumber)

            self.cycles[cyclenumber] = {
                'time': [self.time[idv] for idv in idx],
                'voltage': [self.voltage[idv] for idv in idx],
                'current': [self.current[idv] for idv in idx],
                'Qdischarge': [self.Qdischarge[idv] for idv in idx],
                'Qcharge': [self.Qcharge[idv] for idv in idx],
                'Edischarge': [self.Edischarge[idv] for idv in idx],
                'Echarge': [self.Echarge[idv] for idv in idx],
            }

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
                DCIR[int(self.rawcyclenumbers[idx])].append(abs(dV/dI*1000))

        array = np.array([DCIR[cycle] for cycle in DCIR])

        DCIR_avg = []
        for cycle in DCIR:
            DCIR_avg.append(np.mean(DCIR[cycle]))

        self.DCIR = {
            'cycles': list(DCIR.keys()),
            'array': array,
            'average': DCIR_avg
        }

    def plot_IR_drop(self, average=True, xlim=None, ylim=None, show=False, save=False, imagetype='pdf'):
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
            coloridx = np.linspace(0.25, 1, self.DCIR['array'].shape[1])
            for idx, pulse in enumerate(self.DCIR['array'].T, 1):
                ax.plot(plotcycles, pulse, label='Pulse ' + idx, 
                    color=plt.cm.Purples(coloridx[idx-1]), linewidth=3)
            plottype = '_all'

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        ax.set_xlabel('Cycle Number')
        ax.set_ylabel('DC Internal Resistance [Ω]')
        ax.set_title(self.filename[:-8] + plottype)
        ax.legend()
        ax.grid()

        if show:
            plt.show()
            plt.close()

        if save:
            plt.savefig(self.filename[:-8] + '_DCIR' + plottype + '.' + imagetype)












