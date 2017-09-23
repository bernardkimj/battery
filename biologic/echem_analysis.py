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
    ''' Processes cyclic voltammetry data from Biologic tester
        Based off code originally developed for analyzing data from Gamry 
        Designed for use with batch processing '''

    def __init__(self, filename):
        ''' Loads file and imports raw data 
            time (s), voltage (V), current (mA) '''

        self.filename = filename
        self.cycles = {}

        titlesearch = re.search(r'CELL.*_CV', self.filename)
        self.title = titlesearch.group(0)

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

    def plot_current_voltage(self, cycle_index=list(range(2,11)), 
        title=None, ylim=None, show=False, save=False, savename=None, imagetype='png'):
        ''' Plots current vs voltage for one or all cycles '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.3,1,len(self.cycles)) # for use with Oranges cmap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if cycle_index:
            for index in cycle_index:
                ax.plot(self.cycles[index]['voltage'],
                    self.cycles[index]['current'],
                    linewidth=3,
                    label = 'Cycle' +str(index), 
                    color=plt.cm.Oranges(coloridx[index-1]))
            # ax.legend(loc='upper left', ncol=3)
            ax.legend()
        else:
            for i in range(1,len(self.cycles)):
                ax.plot(self.cycles[i]['voltage'],
                        self.cycles[i]['current'],
                        linewidth=3,
                        color=plt.cm.Oranges(coloridx[i-1]),
                        label='Cycle  ' + str(i))
            # ax.legend(loc='upper left', ncol=4)
            ax.legend()

        ax.set_xlabel('Potential (V)')
        ax.set_ylabel('Current (mA)')
        ax.grid(b=True, which='major', color='0.9', linestyle='-')

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
                plt.savefig('parts_' + self.title + '.' + imagetype)

        plt.close(fig)

class CV_batch:
    ''' Method for batch processing data from Biologic
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

        titlematch = re.search(r'CELL_.*_\d{2}_', alldata[0])
        self.title = titlematch.group(0)[:-4]

    def plot_current_voltage(self, cycle_index=10, title=None, ylim=None, 
        imagetype='png', show=False, save=False):
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
            plt.savefig('batch_' + figtitle + '.' + str(imagetype))

        plt.close()

class PEIS:
    ''' Processes potentiostatic EIS data from Biologic tester '''

    def __init__(self, filename=None, thickness=0.001, area=1, zscale='k'):
        ''' Opens file and retrieves data.

        Retrieves frequency, real impedance, imaginary impedance, 
        impedance magnitude, phase angle, time, voltage, and current

        *** NOTE *** Assumsed imaginary impedance is given as positive values

        Unit requirements:
            R_solution [ohm]
            Thickness [cm]
            Area [cm^2]
        '''

        self.filename = filename
        self.thickness = thickness
        self.area = area

        self.zscalestr = zscale
        self.zscaleval = self.get_zscale(zscale)

        titlesearch = re.search(r'CELL_.*PEIS', self.filename)
        self.title = titlesearch.group(0)

        rows = list(csv.reader(open(filename), delimiter='\t'))

        headers = ['freq/Hz', 'Re(Z)/Ohm', '-Im(Z)/Ohm', '|Z|/Ohm', 
            'Phase(Z)/deg', 'time/s', 'Ecell/V', 'I/mA']
        idxs = []

        for header in headers:
            idxs.append(rows[0].index(header))

        freqidx = idxs[0]
        realidx = idxs[1]
        imagidx = idxs[2]
        magnidx = idxs[3]
        phaseidx = idxs[4]
        timeidx = idxs[5]
        voltageidx = idxs[6]
        currentidx = idxs[7]

        self.freq, self.realraw, self.imagraw, self.magnraw, \
            self.phase, self.time, self.voltage, self.current = \
            [], [], [], [], [], [], [], [],

        for row in rows[1:]:
            self.freq.append(float(row[freqidx]))
            self.realraw.append(float(row[realidx]))
            self.imagraw.append(float(row[imagidx]))
            self.magnraw.append(float(row[magnidx]))
            self.phase.append(float(row[phaseidx]))
            self.time.append(float(row[timeidx]))
            self.voltage.append(float(row[voltageidx]))
            self.current.append(float(row[currentidx]))

        self.real = [realraw/self.zscaleval for realraw in self.realraw]
        self.imag = [imagraw/self.zscaleval for imagraw in self.imagraw]
        self.magn = [magnraw/self.zscaleval for magnraw in self.magnraw]

        self.find_r_solution()
        self.conductivity = (self.thickness * 1000) / (self.area * self.r_solution)

    def get_zscale(self, zscale):
        ''' Determines scaling factor for impedance values
            Default value is in Ohms, scales values to kOhms or MOhms '''

        zscaledict = {'k': 1e3, 'M': 1e6}

        return zscaledict[zscale]


    def find_r_solution(self):
        ''' Calculated solution resistance by taking value of real impedance
            closest to real axis intercept '''

        min_imag_val = min(abs(imag) for imag in self.imagraw[:10])
        min_imag_idx = [abs(imag) for imag in self.imagraw].index(min_imag_val)

        self.r_solution = self.realraw[min_imag_idx]

    def plot_nyquist(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, imagetype='png'):
        ''' Plots imaginary vs. real impedance 
            Imaginary and real axes locked to same scale by default '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(12,9), dpi=75)


        ax.plot(self.real, self.imag, color='b', linewidth=2,
            label=r'$R_{solution}$' + ' = ' + '%.2f'%self.r_solution + ' Ω' + '\n' + 
                r'$\sigma$' + ' = ' + '%0.2f'%self.conductivity + ' S')

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        ax.set_aspect('equal', 'datalim')
        ax.set_xlabel(r'$Z_{real}$' + ' [' + self.zscalestr + 'Ω]')
        ax.set_ylabel(r'$Z_{imag}$' + ' [' + self.zscalestr + 'Ω]')
        ax.legend()
        ax.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_nyquist' + '.' + str(imagetype))

        plt.close(fig)

    def plot_bode(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, imagetype='png'):
        ''' Plots imaginary vs. real impedance 
            Imaginary and real axes locked to same scale by default '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, (ax_magn, ax_phase) = plt.subplots(2, sharex=True, figsize=(16,9), dpi=75)

        ax_magn.semilogx(self.freq, self.magn,
            color='b', linewidth=2)
        ax_phase.semilogx(self.freq, self.phase,
            color='b', linewidth=2)

        if xlim:
            ax_phase.set_xlim(xlim)
        if ylim:
            ax_phase.set_ylim(ylim)

        if title:
            ax_magn.set_title(title)
        else:
            ax_magn.set_title(self.title)

        ax_magn.set_ylabel('Magnitude' + ' [' + self.zscalestr + 'Ω]')
        ax_phase.set_ylabel('Phase [°]')
        ax_phase.set_xlabel('Frequency [Hz]')
        ax_magn.grid()
        ax_phase.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_bode' + '.' + str(imagetype))

        plt.close(fig)

class PEIS_batch:
    ''' Batch processes potentiostatic EIS data from Biologic tester
        Plots samples measurements from same batch on same axes 
        For Nyquist and Bode plots 

        Uses defined methods in PEIS class '''

    def __init__(self, alldata, zscale='k'):
        # Accepts lists of class CV
        self.allcycles = {}

        self.zscalestr = zscale
        self.zscaleval = self.get_zscale(zscale=zscale)

        for file in alldata:
            exported = PEIS(file, zscale=self.zscalestr)

            match = re.search(r'_\d{2}_', file)
            stepidx = (int(match.group(0)[1:3])-1)/2
            self.allcycles[stepidx] = exported

        titlematch = re.search(r'CELL_.*_\d{2}_', alldata[0])
        self.title = titlematch.group(0)[:-4]

    def get_zscale(self, zscale):
        ''' Determines scaling factor for impedance values
            Default value is in Ohms, scales values to kOhms or MOhms '''

        zscaledict = {'k': 1e3, 'M': 1e6}

        return zscaledict[zscale]

    def plot_nyquist(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, imagetype='png'):
        ''' Plots Nyquist plot with all samples on one graph '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.3,1,len(self.allcycles)) # for Blues colormap

        fig, ax = plt.subplots(figsize=(12,9), dpi=75)

        for sample in sorted(self.allcycles):
            ax.plot(self.allcycles[sample].real, self.allcycles[sample].imag,
                color=plt.cm.Blues(coloridx[int(sample)]), linewidth=2,
                label='Cycle '+str(int(sample*10))+', '+r'$R_{solution}$'+\
                    ' = '+'%.2f'%self.allcycles[sample].r_solution + ' Ω')

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        ax.set_aspect('equal', 'datalim')
        ax.set_xlabel(r'$Z_{real}$' + ' [' + self.zscalestr + 'Ω]')
        ax.set_ylabel(r'$Z_{imag}$' + ' [' + self.zscalestr + 'Ω]')
        ax.legend()
        ax.grid()

        if show:
            plt.show()

        if save:
            plt.savefig('batch_' + self.title + '_nyquist' + '.' + str(imagetype))

        plt.close(fig)

    def plot_bode(self, xlim=None, ylim=None, title=None,
        show=False, save=False, imagetype='png'):
        ''' Plots Bode plot with all samples on one graph '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.3,1,len(self.allcycles)) # for Blues colormap

        fig, (ax_magn, ax_phase) = plt.subplots(2, sharex=True, figsize=(16,9), dpi=75)

        for sample in sorted(self.allcycles):
            ax_magn.semilogx(self.allcycles[sample].freq, self.allcycles[sample].magn,
                color=plt.cm.Blues(coloridx[int(sample)]), linewidth=2,
                label='Cycle '+str(int(sample*10)))
            ax_phase.semilogx(self.allcycles[sample].freq, self.allcycles[sample].phase,
                color=plt.cm.Blues(coloridx[int(sample)]), linewidth=2,
                label='Cycle '+str(int(sample*10)))

        if xlim:
            ax_phase.set_xlim(xlim)
        if ylim:
            ax_phase.set_ylim(ylim)

        if title:
            ax_magn.set_title(title)
        else:
            ax_magn.set_title(self.title)

        ax_magn.set_ylabel('Magnitude' + ' [' + self.zscalestr + 'Ω]')
        ax_phase.set_ylabel('Phase [°]')
        ax_phase.set_xlabel('Frequency [Hz]')
        ax_magn.grid()
        ax_phase.grid()

        if show:
            plt.show()

        if save:
            plt.savefig('batch_' + self.title + '_bode' + '.' + str(imagetype))

        plt.close(fig)










































