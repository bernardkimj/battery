''' Analysis module to process electrochemical analysis data from Bio-Logic Tester.
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

class CV:
    ''' Processes cyclic voltammetry data
        Based off code originally developed for analyzing data from Gamry 
        Designed for use with batch processing '''

    def __init__(self, filename):
        ''' Loads file and imports raw data 
            time (s), voltage (V), current (mA) '''

        self.filename = filename
        self.cycles = {}

        rows = list(csv.reader(open(filename), delimiter='\t'))

        headers = ['cycle number', 'time/s', 'Ecell/V', '<I>/mA']
        idxs = []

        for header in headers:
            idxs.append(rows[0].index(header))

        cycleidx = idxs[0]
        timeidx = idxs[1]
        voltageidx = idxs[2]
        currentidx = idxs[3]

        rawcycles = []
        for row in rows[1:]:
            rawcycles.append(float(row[cycleidx]))

        cyclenumbers = []
        for entry in set(rawcycles):
            cyclenumbers.append(int(entry))
        cyclenumbers.sort()

        for cyclenumber in cyclenumbers:
            self.cycles[cyclenumber] = {
                'time': [],
                'voltage': [],
                'current': [],
            }

        for row in rows[1:]:
            self.cycles[int(float(row[cycleidx]))]['time'].append(float(row[timeidx]))
            self.cycles[int(float(row[cycleidx]))]['voltage'].append(float(row[voltageidx]))
            self.cycles[int(float(row[cycleidx]))]['current'].append(float(row[currentidx]))

    def plot_current_voltage(self, cycle_index=list(range(2,11)), imageformat='png', title=None, ylim=None, show=False, save=False):
        ''' Plots current vs voltage for one or all cycles '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.3,1,len(self.cycles)) # for use with Oranges cmap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if cycle_index:
            for index in cycle_index:
                ax.plot(self.cycles[index]['voltage'],
                    self.cycles[index]['current'],
                    marker='.', markersize=3,
                    label = 'Cycle' +str(index), 
                    color=plt.cm.Oranges(coloridx[index-1]))
            # ax.legend(loc='upper left', ncol=3)
            ax.legend()
        else:
            for i in range(1,len(self.cycles)):
                ax.plot(self.cycles[i]['voltage'],
                        self.cycles[i]['current'],
                        marker='.', markersize=3,
                        color=plt.cm.Oranges(coloridx[i-1]),
                        label='Cycle  ' + str(i))
            # ax.legend(loc='upper left', ncol=4)
            ax.legend()

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (mA)')
        plt.grid(b=True, which='major', color='0.9', linestyle='-')

        if ylim:
            plt.ylim(ylim)

        if title:
            plt.title(str(title[:-4]) + '_C' + str(cycle_index))
        else:
            plt.title(str(self.filename[:-4]))

        if show:
            plt.show()

        if save:
            plt.savefig('parts_' + str(self.filename[0:-8]) + '.' + str(imageformat))

        plt.close()

class CV_batch:
    ''' Method for batch processing data from Gamry
    Uses methods defined in CV class

    Author: Bernard Kim
    '''

    def __init__(self, alldata, reduce=False):
        # Accepts lists of class CV
        self.allcycles = {}

        for file in alldata:
            exported = CV(file)

            match = re.search(r'_\d{2}_', file)
            stepidx = int(match.group(0)[1:3])/2
            self.allcycles[stepidx] = exported.cycles

        titlematch = re.search(r'CV_.*_\d{2}_', alldata[0])
        self.title = titlematch.group(0)[:-4]

    def plot_current_voltage(self, cycle_index=10, title=None, ylim=None, imageformat='png', show=False, save=False):
        ''' Plots current vs voltage by cycle with all samples on one graph '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.3,1,len(self.allcycles)) # for Blues colormap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for sample in sorted(self.allcycles):
            if cycle_index:
                ax.plot(self.allcycles[sample][cycle_index]['voltage'],
                        self.allcycles[sample][cycle_index]['current'],
                        marker='.', markersize=8,
                        color=plt.cm.Blues(coloridx[int(sample)-1]),
                        label='Cycle '+str(int((sample-1)*10 + cycle_index)))
            # COMMENTED OUT FOR NOW BECAUSE DOESN'T MAKE SENSE TO INCLUDE FOR NOW
            # else:
            #     for i in range(1,len(self.allcycles[sample])):
            #         ax.plot(self.allcycles[sample][i]['voltage'],
            #                 self.allcycles[sample][i]['current'],
            #                 marker='.', markersize=12)

        ax.legend()
        # ax.legend(loc='upper left', fontsize=10, ncol=3)
        plt.xlabel('Potential (V)')
        plt.ylabel('Current (mA)')
        plt.grid(b=True, which='major', color='0.8', linestyle='-')

        if ylim:
            plt.ylim(ylim)
        # else:
        #     plt.ylim([-25, 25])

        # figtitle = title + '_C' + str(cycle_index) + '_n=' + str(len(self.allcycles))
        if title:
            figtitle = title
        else:
            figtitle = self.title

        plt.title(figtitle)

        if show:
            plt.show()

        if save:
            plt.savefig('batch_' + figtitle + '.' + str(imageformat))

        plt.close()

















































