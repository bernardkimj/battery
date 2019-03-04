''' Module for analyzing results exported from Brookfield rheometer
DV3T

Author: Bernard Kim
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley

'''

from battery.utilities import utilities
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
import re

class viscosity:

    ''' Analyzes data from Brookfield rheometer
    Pulls data from .csv file for plotting and analysis.
    '''

    def __init__(self, filename=None):
        ''' Opens file and retrieves data. 
            For viscosity test data from Brookfield DV3T, exported as csv
        '''

        # Step [#], Point [#], Time [s], Viscosity [cP],  Torque [%], Speed [RPM],
        # Shear Stress [dyne/cm^2], Shear Rate [pi/s], Temperature [C], 
        # Density [g/cm^3], Accuracy [+/-cP]


        # All possible headers
        headers = {
            'Step': 'step',
            'Point': 'point',
            'Time': 'time',
            'Viscosity': 'viscosity',
            'Torque': 'torque',
            'Speed': 'speed',
            'Shear Stress': 'shearstress',
            'Shear Rate': 'shearrate',
            'Temperature': 'temperature',
            'Density': 'density',
            'Accuracy': 'accuracy',
        }

        # Instantiate dict for indices within row (agnostic to export order)
        idxs = {}
        self.data = {}

        reader = csv.reader(open(filename, errors='replace'), delimiter=',')
        rows = list(reader)
        headerloc = None
        switch = None

        for index, row in enumerate(rows):
            try:
                if row:
                    # Find the part in the file where the actual data starts
                    if row[0] == 'DATA':
                        headerloc = index + 1
                        switch = index + 2
                # Find the indices of each data column based on header
                # Create a list in self.data for each header that exists
                if index == headerloc:
                    for header in headers:
                        if header in row:
                            idxs[headers[header]] = row.index(header)
                            self.data[headers[header]] = []
                # Write data to appropriate list in self.dict
                # If not a number, write NaN
                if switch is not None and index > switch:
                    for header in self.data:
                        try:
                            self.data[header].append(float(row[idxs[header]]))
                        except ValueError:
                            self.data[header].append(float('nan'))
                            # self.data[header].append(None)

            except Exception:
                raise

        self.plotinfo = {
            'step': {
                'label': 'Step [#]',
                'titlestring': 'Step',
            },
            'point': {
                'label': 'Point [#]',
                'titlestring': 'Point',
            },
            'time': {
                'label': 'Time [s]',
                'titlestring': 'Time',
            },
            'viscosity': {
                'label': 'Viscosity [cP]',
                'titlestring': 'Viscosity',
            },
            'torque': {
                'label': 'Torque [%]',
                'titlestring': 'Torque',
            },
            'speed': {
                'label': 'Speed [RPM]',
                'titlestring': 'Speed',
            },
            'shearstress': {
                'label': 'Shear Stress [dyne/cm'+r'$^2$'+']',
                'titlestring': 'Shear Stress',
            },
            'shearrate': {
                'label': 'Shear Rate [1/s]',
                'titlestring': 'Shear Rate',
            },
            'temperature': {
                'label': 'Temperature [°C]',
                'titlestring': 'Temperature',
            },
            'density': {
                'label': 'Density [g/cm'+r'$^3$'+']',
                'titlestring': 'Density',
            },
            'accuracy': {
                'label': 'Accuracy [+/- cP]',
                'titlestring': 'Accuracy',
            },
        }

        self.title = filename[:-4]

        # self.calculate_yield_stress(show=False)


    def plot_data(self, x=None, y=None, xlim=None, ylim=None, 
        title=None, show=False, save=False, savename=None, imagetype='png'):
        ''' Plots data
            Possible choices for both x and y are
            'step', 'point', 'time', 'viscosity', 'torque', 'speed', 
            'shearstress', 'shearrate', 'temperature', 'density', 'accuracy', 
        '''

        font = {'family':'Arial','size': 32}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        # mask_y = np.isfinite(np.array(self.data[y]).astype(np.double))

        ax.plot(self.data[x], self.data[y],
                marker='o', markersize=8,
                linewidth=3, color='#000000',)

        ax.set_xlabel(self.plotinfo[x]['label'])
        ax.set_ylabel(self.plotinfo[y]['label'])
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + ', ' + self.plotinfo[y]['titlestring'] + 
                ' vs. ' + self.plotinfo[x]['titlestring'])

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig(self.title + '_'+y+'_'+x + '.' + str(imagetype))

        if show:
            plt.show()

        plt.close(fig)


    def calculate_yield_stress(self, show=False, save=False):
        ''' Calculates yield stress based on linear fit of 
            shear stress vs. shear rate data.
            Linearly weights linear fit to give higher priority to higher
            shear rates.
        '''

        shearrate = np.array(self.data['shearrate'])
        shearstress = np.array(self.data['shearstress'])

        coefs = np.polyfit(x=shearrate, y=shearstress, deg=1, w=shearrate)

        self.data['yieldstress'] = coefs[1]

        if show or save:
            line = coefs[0]*shearrate + coefs[1]
            intercept_x = np.array([0, shearrate[0]])
            intercept_y = np.array([coefs[1], line[0]])

            font = {'family':'Arial', 'size':24}
            matplotlib.rc('font', **font)

            fig, ax = plt.subplots(figsize=(16,9), dpi=75)
            ax.plot(shearrate, shearstress, linewidth=3,
                color='k', marker='o', markersize='8',)
            ax.plot(shearrate, line, color='r', linewidth=3)
            ax.plot(intercept_x, intercept_y, linewidth=3, 
                linestyle='--', color='r', 
                label='Yield Stress = '+'%.2f'%coefs[1]+' dyne/cm'+r'$^2$')

            ax.set_xlim([0, np.max(shearrate)])
            ax.grid()
            ax.legend()

            ax.set_xlabel(self.plotinfo['shearrate']['label'])
            ax.set_ylabel(self.plotinfo['shearstress']['label'])
            ax.set_title(self.title)

        if save:
            plt.savefig(self.title+'_yieldstress.png')
        if show:
            plt.show()


class viscosity_batch:
    '''
    Class for working with multiple profiles at once
    Relies on methods defined in 'profile' class
    '''


    def __init__(self, files=None):
        ''' 
        Initializes files according 'profile'
        '''

        # Acceptes lists of filenames
        self.alldata = {}
        titles = []
        plotinfo = []

        first = True
        for file in files:
            exported = viscosity(file)
            exported.calculate_yield_stress(save=True)

            title, sampleidx = utilities.title_search_format(
                filename=file, include=['ballmilling'])

            # title, sampleidx = self.title_search_format(filename=file)

            self.alldata[sampleidx] = exported.data
            titles.append(title)

            if first:
                self.plotinfo = exported.plotinfo
                first = False

        # If titles are identical, groovy
        # Otherwise, print error
        self.title = titles[0]

        if len(set(titles)) is not 1:
            print('Titles do not match!!')


    # def title_search_format(self, filename):
    #     ''' 
    #     Searches filenames for batch 
    #     processing
    #     Returns formatted title with constant sample parameters for 
    #     batch plotting and sample idx
    #     '''

    #     # Searches filename for ball milling parameters, ink name, and 
    #     # casting direction
    #     inkmatch = re.search(r'[A-Z]{3}\d{1}[a-z]{1}\d{2}', filename)
    #     paramsmatch = re.search(r'[0-9]{6}(S)?', filename)
    #     Nmatch = re.search(r'S\d{1,2}', filename)

    #     # Only used for initial qualification of ball mill jars
    #     isoldmatch = re.search(r'old', filename)
    #     isTEGmatch = re.search(r'TEG', filename)

    #     # Determines order of elements in title
    #     titlematches = [inkmatch, paramsmatch]

    #     # Instantiate title
    #     title = ''

    #     # Check and assign sample number
    #     sampleidx = Nmatch.group(0)

    #     # Reconstruct the title in the right order, regardless of 
    #     # filename order
    #     for match in titlematches:
    #         if match is not None:
    #             title = title + match.group(0)

    #     if isoldmatch:
    #         title = title + '_old_jars'
    #     elif isTEGmatch:
    #         title = title + '_TEG_jars'

    #     return title, sampleidx


    def plot_batch(self, x=None, y=None, 
        average=False, confidence=0.95, 
        xlim=None, ylim=None, title=None, 
        show=False, save=False, savename=None, imagetype='png'):
        ''' Plots multiple profiles on same plot
            Possible choices for both x and y are
            'step', 'point', 'time', 'viscosity', 'torque', 'speed', 
            'shearstress', 'shearrate', 'temperature', 'density', 'accuracy', 
        '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0,1,10) # for use with tab10 colormap
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if average:
            data = [(self.alldata[sample][x], self.alldata[sample][y]) \
                for sample in self.alldata]

            indep, mean, std, lcl, ucl = utilities.batch_average_plot(data=data)

            for sample in self.alldata:
                ax.plot(self.alldata[sample][x], self.alldata[sample][y],
                    linewidth=3, color='k', alpha=0.2)

            ax.plot(indep, mean, linewidth=3, color='r', label='mean')
            ax.plot(indep, ucl, color='r', linewidth=3, linestyle='--',)
            ax.plot(indep, lcl, color='r', linewidth=3, linestyle='--', 
                label=str(confidence)[2:]+'% Confidence Interval (t-test)')

        else:
            for idx, sample in enumerate(self.alldata):
                x_plot = np.array(self.alldata[sample][x])
                y_plot = np.array(self.alldata[sample][y])

                if np.isnan(np.min(y_plot)):
                    x_plot = x_plot.astype(np.double)
                    y_plot = y_plot.astype(np.double)
                    mask = np.isfinite(y_plot)

                    ax.plot(x_plot[mask], y_plot[mask],
                        linewidth=3, color=plt.cm.tab10(coloridx[idx]),
                        linestyle='--')

                ax.plot(x_plot, y_plot,
                    linewidth=3, color=plt.cm.tab10(coloridx[idx]), 
                    marker='o', markersize=8, label=str(sample))

        ax.set_xlabel(self.plotinfo[x]['label'])
        ax.set_ylabel(self.plotinfo[y]['label'])
        ax.legend()
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + ', ' + self.plotinfo[y]['titlestring'] + 
                ' vs. ' + self.plotinfo[x]['titlestring'])

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('batch_' + self.title + '_'+y+'_'+x + 
                    '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)












































