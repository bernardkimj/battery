''' Module for analyzing results retrieved from Dektak profilometer

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

font = {'family': 'Arial', 'size': 16}
matplotlib.rc('font', **font)

class profile:

    ''' Analyzes data from Dektak profilometer

    Pulls data from .csv file for plotting and analysis.

    '''

    def __init__(self, filename=None):
        ''' Opens file and retrieves data. '''

    # x   h
    # um  nm

        x, height = [], []

        rows = list(csv.reader(open(filename), delimiter=','))

        switch = None

        for index, row in enumerate(rows):
            try:
                if row:
                    # Find the part in the file where the actual data starts
                    if row[0][0:2] == 'um' and row[1] == 'A':
                        switch = index
                if switch is not None and index > switch:
                    if len(row) > 0:
                        x.append(float(row[0]) * 0.001)
                        height.append(float(row[1]) * 0.001)

            except Exception:
                raise

        # Save data and convert measurements to mm and um for x and height 
        # respectively
        self.line = {
            'x': x,
            'height': height,
        }

        self.title = filename[:-4]


    def plot_profile(self, xlim=None, ylim=None, 
        title=None, show=False, save=False, savename=None, imagetype='png'):
        ''' Plots profile '''

        font = {'family': 'Arial', 'size': 32}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        ax.plot(self.line['x'], self.line['height'],
                linewidth=3, color='#000000',)

        ax.set_xlabel('Horizontal Position [mm]')
        ax.set_ylabel('Height [µm]')
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
                plt.savefig(self.title + '_plot' + '.' + str(imagetype))

        if show:
            plt.show()

        plt.close(fig)

class profile_batch:
    '''
    Class for working with multiple profiles at once
    Relies on methods defined in 'profile' class
    '''


    def __init__(self, alldata=None):
        ''' 
        Initializes files according 'profile'
        '''

        # Acceptes lists of filenames
        self.lines = {}
        titles = []

        for file in alldata:
            exported = profile(file)
            title, sampleidx = self.title_search_format(filename=file)
            self.lines[sampleidx] = exported.line
            titles.append(title)

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
        isparmatch = re.search(r'pa', filename)
        isperpmatch = re.search(r'pe', filename)

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

        if isparmatch:
            title = title + 'parallel'
        elif isperpmatch:
            title = title + 'perpendicular'

        return title, sampleidx


    def plot_batch_profile(self, average=False, confidence=0.95, 
        xlim=None, ylim=None, title=None, show=False, save=False, 
        savename=None, imagetype='png'):
        '''
        Plots multiple profiles on same plot
        '''

        font = {'family': 'Arial', 'size': 20}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0,1,10) # for use with tab10 colormap
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if average:
            data = [(self.lines[line]['x'], self.lines[line]['height']) \
                for line in self.lines]

            x, mean, std, lcl, ucl = utilities.batch_average_plot(data=data)

            for line in self.lines:
                ax.plot(self.lines[line]['x'], self.lines[line]['height'],
                    linewidth=3, color='k', alpha=0.5)

            ax.plot(x, mean, linewidth=3, color='r', label='mean')
            ax.plot(x, ucl, color='r', linewidth=3, linestyle='--',)
            ax.plot(x, lcl, color='r', linewidth=3, linestyle='--', 
                label=str(confidence)[2:]+'% Confidence Interval (t-test)')

        else:
            for idx, line in enumerate(self.lines):
                ax.plot(self.lines[line]['x'], self.lines['line']['height'],
                    linewidth=3, color=coloridx[idx], label=str(line))

        ax.set_xlabel('Horizontal Position [mm]')
        ax.set_ylabel('Height [µm]')
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
                plt.savefig('batch_' + self.title + '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)












































