''' Analysis module to process battery cycling data from Bio-Logic Tester.
    Assumes files have been made by 'Export as txt'

    Author: Bernard Kim
    Principal Investigators: Paul Wright, James Evans
    Institution: University of California, Berkeley '''

from battery.utilities import utilities
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import re
import xlrd
import datetime

class cycling:
    ''' Processes cycling data from Neware
        Standard CC charge/discharge profile '''

    def __init__(self, files, title=None):
        ''' Loads file and imports raw data into pandas dataframe
            Separates data by cycle into dictionary
            Values per cycle: time, voltage, current, Qdischarge,
            Qcharge, Edischarge, Echarge'''

        self.filename = files[0]
        self.time = []
        self.voltage = []
        self.current = []
        self.rawcycles = []
        self.cycles = {}

        headers = {
            'Cycle ID': 'cycle',
            'Current(mA)': 'current',
            'Voltage(V)': 'voltage',
            'Capacity(mAh)': 'Q',
            'Energy(mWh)': 'W',
            'Cycle': 'cycle',
            'Cur(mA)': 'current',
            'Voltage(V)': 'voltage',
            'CapaCity(mAh)': 'Q',
            'Energy(mWh)': 'W',

            # 'Relative Time(h:min:s.ms)': 'time',
            # 'time/s': 'time', 
            # 'Ecell/V': 'voltage', 
            # 'Ewe/V': 'voltage',
            # '<I>/mA': 'current', 
            # 'I/mA': 'current',
            # 'Q discharge/mA.h': 'Qdischarge',
            # 'Q charge/mA.h': 'Qcharge', 
            # 'Energy discharge/W.h': 'Edischarge', 
            # 'Energy charge/W.h': 'Echarge',
        }

        entries = (
            'voltage',
            'current',
            'Qcharge',
            'Qdischarge',
            'Echarge',
            'Edischarge',
        )

        # Instantiate dict for indices within row 
        # (agnostic to export order)
        idxs = {}

        for file in files:
            wb = xlrd.open_workbook(file)
            sheetnames = [s for s in wb.sheet_names()]

            detailsheets = {}

            for name in sheetnames:
                detailsearch = re.search(r'Detail', name)
                
                if detailsearch is not None:
                    numsearch = re.search(r'__[0-9]{1,2}', name)

                    if numsearch is not None:
                        sheetnum = int(numsearch.group(0)[2:])
                    else:
                        sheetnum = name[-1]

                    detailsheets[sheetnum] = name

            for sheetnum in sorted(detailsheets.keys()):
                sheetname = detailsheets[sheetnum]

                datasheet = wb.sheet_by_name(sheetname)
                num_rows = datasheet.nrows

                if sheetnum == min(detailsheets.keys()):
                    header_row = datasheet.row_values(0)

                    # Get indices for values that were exported
                    for header in headers:
                        if header in header_row:
                            idxs[headers[header]] = header_row.index(header)


                for row_idx in range(1, num_rows):
                    row = datasheet.row_values(row_idx)

                    # time_str = str(row[idxs['time']])
                    # raw_time, dummy = time_str.split('.')
                    # time_s = self.get_sec(raw_time)

                    # self.time.append(float(time_s))
                    self.voltage.append(float(row[idxs['voltage']]))
                    self.current.append(float(row[idxs['current']]))

                    rawcycle = float(row[idxs['cycle']])
                    cycle = int(rawcycle)
                    self.rawcycles.append(rawcycle)

                    # Check if cycle number key has been added to dictionary
                    # If not, create it
                    if cycle not in self.cycles:
                        self.cycles[cycle] = {}

                    # For current row, populate cycle dictionary with existing
                    # entries
                    for entry in entries:
                        # Check if entry is in cycle dict, add if not
                        if entry not in self.cycles[cycle]:
                            self.cycles[cycle][entry] = [0]

                    voltage = float(row[idxs['voltage']])
                    current = float(row[idxs['current']])
                    charge = float(row[idxs['Q']])
                    energy = float(row[idxs['W']])

                    self.cycles[cycle]['voltage'].append(voltage)
                    self.cycles[cycle]['current'].append(current)

                    if np.sign(current) == 1:
                        self.cycles[cycle]['Qcharge'].append(charge)
                        self.cycles[cycle]['Echarge'].append(energy)
                        self.cycles[cycle]['Qdischarge'].append(0)
                        self.cycles[cycle]['Edischarge'].append(0)
                    elif np.sign(current) == -1:
                        self.cycles[cycle]['Qdischarge'].append(charge)
                        self.cycles[cycle]['Edischarge'].append(energy)
                        self.cycles[cycle]['Qcharge'].append(0)
                        self.cycles[cycle]['Echarge'].append(0)

        if title:
            self.title = title
        else:
            endtitle = self.filename.find('.x')
            self.title = self.filename[8:endtitle]
            # titlematch = re.search(r'CELL.*\d{1}_C', self.filename)
            # self.title = titlematch.group(0)[:-2]


    def get_sec(self, time_str):
        """Get Seconds from time."""
        h, m, s = time_str.split(':')
        return int(h) * 3600 + int(m) * 60 + int(s)


    def get_IR_data(self):
        ''' Calculates IR data upon charge, discharge, and rest per cycle 
            Organized in dictionary with keys 'charge', 'discharge', 'rest' 
            Called automatically by 'plot_IR_drop'
            Robust to any number of rest, charge, or discharge steps
            '''

        # Initialize dictionary for raw data
        DCIR = {
            cycle: {
                step: np.nan for step in ['charge', 'discharge', 'rest']
            } for cycle in self.cycles
        }

        # Determine type of step when current changes signs and calculate 
        # instantaneous DC internal resistance
        for idx in range(1, len(self.current)):
            dV = self.voltage[idx] - self.voltage[idx-1]
            dI = self.current[idx] - self.current[idx-1]

            # if dI != 0:
            #     previous = self.current[idx-1]
            #     present = self.current[idx]

            previous = self.current[idx-1]
            present = self.current[idx]

            if np.sign(previous) != np.sign(present):
                if np.sign(present) == 1 and np.sign(previous) != 1:
                    key = 'charge'
                elif np.sign(present) == -1 and np.sign(previous) != -1:
                    key = 'discharge'
                elif np.sign(present) == 0 and np.sign(previous) != 0:
                    key = 'rest'

                resistance = np.abs(dV/(dI/1000)) # Ohms

                DCIR[int(self.rawcycles[idx])][key] = resistance

        # Fill in average values for all steps for all cycles
        self.DCIR = {}

        for cycle in DCIR:
            self.DCIR[cycle] = {}

            for step in DCIR[cycle]:
                self.DCIR[cycle][step] = DCIR[cycle][step]


    def get_dQdV(self, cycles=None):
        ''' Calculates dQ/dV from cycling data
            Length of vector is len(voltage)-1
            Also prepares raw charge/discharge capacity vs. voltage
        '''

        if not cycles:
            cycles = [cycle for cycle in self.cycles]

        # Initialize dictionaries for raw data
        self.dQdV_C = {
            cycle: {
                'dQdV': [],
                'voltage': [],
            } for cycle in cycles
        }

        self.dQdV_D = {
            cycle: {
                'dQdV': [],
                'voltage': [],
            } for cycle in cycles
        }

        self.DoC = {
            cycle: {
                'Qcharge': [],
                'voltage': [],
            } for cycle in cycles
        }

        self.DoD = {
            cycle: {
                'Qdischarge': [],
                'voltage': [],
            } for cycle in cycles
        }

        for cycle in cycles:
            charge, discharge = [], []

            for voltage, current, Qcharge, Qdischarge in zip(
                self.cycles[cycle]['voltage'], \
                self.cycles[cycle]['current'], \
                self.cycles[cycle]['Qcharge'], \
                self.cycles[cycle]['Qdischarge'],
            ):

                if Qcharge != 0 and current > 0:
                    charge.append((voltage, Qcharge))

                if Qdischarge != 0 and current < 0:
                    discharge.append((voltage, Qdischarge))

            for step in (charge, discharge):
                voltage = [point[0] for point in step]
                Q_raw = [point[1] for point in step]

                dQdV, dQdV_voltage = [], []

                for idx in range(1, len(voltage)):
                    dV = voltage[idx] - voltage[idx-1]
                    dQ = Q_raw[idx] - Q_raw[idx-1]

                    if dV != 0:
                        dQdV.append(dQ/dV)
                        dQdV_voltage.append(voltage[idx])

                if step == charge:
                    self.dQdV_C[cycle]['dQdV'] = dQdV
                    self.dQdV_C[cycle]['voltage'] = dQdV_voltage
                    self.DoC[cycle]['voltage'] = voltage
                    self.DoC[cycle]['Qcharge'] = Q_raw
                elif step == discharge:
                    self.dQdV_D[cycle]['dQdV'] = dQdV
                    self.dQdV_D[cycle]['voltage'] = dQdV_voltage
                    self.DoD[cycle]['voltage'] = voltage
                    self.DoD[cycle]['Qdischarge'] = Q_raw


    def get_efficiency(self):
        ''' Calculates efficiency per cycle
            saves to self.cycles[cycle]['efficiency']
        '''

        for cycle in self.cycles:
            self.cycles[cycle]['efficiency'] = np.max(
                self.cycles[cycle]['Qdischarge']) / \
                np.max(self.cycles[cycle]['Qcharge']) * 100


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

        fig, ax = utilities.plot_setup(params=params)

        ax.plot(
            self.time, 
            self.voltage, 
            color='b', 
            linewidth=2,
        )

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
        efficiency = [self.cycles[cycle]['efficiency'] for cycle in cycles]

        fig, ax = utilities.plot_setup(params=params)
        ax.plot(
            cycles[:-1], 
            efficiency[:-1], 
            color='b', 
            linewidth=2,
        )

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def plot_capacity(self, discharge=True, charge=True, xlim=False, ylim=False, 
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

        fig, ax = utilities.plot_setup(item=self, params=params)

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting
        cycles = list(self.cycles.keys())

        # Plot both discharge and charge by default
        if discharge:
            Qdischarge = [np.max(self.cycles[cycle]['Qdischarge']) 
                for cycle in cycles]

            ax.plot(
                cycles, 
                Qdischarge, 
                color='b', 
                linewidth=2, 
                label=r'$Q_{discharge}$'
            )
        if charge:
            Qcharge = [np.max(self.cycles[cycle]['Qcharge']) 
                for cycle in cycles]

            ax.plot(
                cycles, 
                Qcharge, 
                color='r', 
                linewidth=2, 
                label=r'$Q_{charge}$'
            )

        ax.legend()

        if save:
            plt.savefig(self.title + params['titletag']+'.png')
        if show:
            plt.show()

        plt.close(fig)


    def get_capacity_fade(self):
        ''' Calculates capacity fade
            Similar to capacity per cycle but normalized to first cycle capacity
            saves to self.cycles[cycle]['fade']
        '''

        baseline = np.max(self.cycles[cycle[0]]['Qdischarge'])

        for cycle in self.cycles:
            self.cycles[cycle]['fade'] = np.max(
                self.cycles[cycle]['Qdischarge']) / \
                baseline * 100


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

        fig, ax = utilities.plot_setup(params=params)

        # Get list of all cycle numbers
        # Omit cycle 0 from plotting
        cycles = list(self.cycles.keys())[1:]

        # baseline = np.max(self.cycles[cycles[0]]['Qdischarge'])
        fade = [self.cycles[cycle]['fade'] for cycle in cycles]

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

        fig, ax = utilities.plot_setup(params=params)

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

        fig, ax = utilities.plot_setup(params=params)

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
