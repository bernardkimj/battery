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

        step, point, time, viscosity, torque, speed, shearstress, shearrate, \
            temperature, density, accuracy \
            = [], [], [], [], [], [], [], [], [], [], []

        reader = csv.reader(open(filename, errors='replace'), delimiter=',')
        rows = list(reader)
        switch = None

        for index, row in enumerate(rows):
            try:
                if row:
                    # Find the part in the file where the actual data starts
                    if row[0] == 'DATA':
                        switch = index + 2
                if switch is not None and index > switch:
                    if len(row) > 0:
                        step.append(float(row[0])) # [#]
                        point.append(float(row[1])) # [#]
                        time.append(float(row[2])) # [s]
                        viscosity.append(float(row[3])) # [cP]
                        torque.append(float(row[4])) # [%]
                        speed.append(float(row[5])) # [RPM]
                        shearstress.append(float(row[6])) # [dyne/cm^2]
                        shearrate.append(float(row[7])) # [pi/s]
                        temperature.append(float(row[8])) # [C]
                        density.append(float(row[9])) # [g/cm^3]
                        accuracy.append(float(row[10])) # [+/- cP]
            except Exception:
                raise

        # Save data 
        self.data = {
            'step': step,
            'point': point,
            'time': time,
            'viscosity': viscosity,
            'torque': torque,
            'speed': speed,
            'shearstress': shearstress,
            'shearrate': shearrate,
            'temperature': temperature,
            'density': density,
            'accuracy': accuracy,
        }

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
                'label': 'Shear Rate ['+r'$\pi$'+'/s]',
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


    def plot_data(self, x=None, y=None, xlim=None, ylim=None, 
        title=None, show=False, save=False, savename=None, imagetype='png'):
        ''' Plots data
            Possible choices for both x and y are
            'step', 'point', 'time', 'viscosity', 'torque', 'speed', 
            'shearstress', 'shearrate', 'temperature', 'density', 'accuracy', 
        '''

        font = {'family': 'Arial', 'size': 32}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

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


class viscosity_batch:
    '''
    Class for working with multiple profiles at once
    Relies on methods defined in 'profile' class
    '''


    def __init__(self, alldata=None):
        ''' 
        Initializes files according 'profile'
        '''

        # Acceptes lists of filenames
        self.alldata = {}
        titles = []
        plotinfo = []

        first = True
        for file in alldata:
            exported = viscosity(file)
            title, sampleidx = self.title_search_format(filename=file)
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


    def title_search_format(self, filename):
        ''' 
        Searches filenames for batch 
        processing
        Returns formatted title with constant sample parameters for 
        batch plotting and sample idx
        '''

        # Searches filename for ball milling parameters, ink name, and 
        # casting direction
        inkmatch = re.search(r'[A-Z]{3}\d{1}[a-z]{1}\d{2}', filename)
        paramsmatch = re.search(r'[0-9]{6}', filename)
        Nmatch = re.search(r'S\d{1,2}', filename)

        # Only used for initial qualification of ball mill jars
        isoldmatch = re.search(r'old', filename)
        isTEGpmatch = re.search(r'TEG', filename)

        # Determines order of elements in title
        titlematches = [inkmatch, paramsmatch]

        # Instantiate title
        title = ''

        # Check and assign sample number
        sampleidx = Nmatch.group(0)

        # Reconstruct the title in the right order, regardless of 
        # filename order
        for match in titlematches:
            if match is not None:
                title = title + match.group(0) + '_'

        if isoldmatch:
            title = title + 'old_jars'
        elif isTEGpmatch:
            title = title + 'TEG_jars'

        return title, sampleidx


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
                ax.plot(self.alldata[sample][x], self.alldata[sample][y],
                    linewidth=3, color=plt.cm.tab10(coloridx[idx]), 
                    label=str(sample))

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
                plt.savefig('batch_' + self.title + '_'+y+'_'+x + '.' + 
                    '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)












































