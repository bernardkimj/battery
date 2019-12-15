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

        mask = np.isfinite(shearstress)

        if False not in mask:
            pass
        else:
            end = np.where(mask==False)[0][0]
            new = np.zeros(len(mask), dtype=bool)
            new[0:end] = True
            new[end:] = False
            mask = new

        coefs = np.polyfit(
            x=shearrate[mask], 
            y=shearstress[mask], 
            deg=1, 
            w=shearrate[mask],
        )

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


    def __init__(self, files=None, inkrename=None):
        ''' 
        Initializes files according 'profile'
        '''

        # Acceptes lists of filenames
        self.alldata = {}
        self.labels = {}
        titles = []
        inknames = []
        plotinfo = []

        first = True
        for file in files:
            exported = viscosity(file)
            exported.calculate_yield_stress(save=False)

            title, inkname, sampleidx = utilities.title_search_format(
                filename=file, include=['ballmilling'])

            inknames.append(inkname)

            try:
                self.alldata[inkname]
            except KeyError:
                self.alldata[inkname] = {}

            self.alldata[inkname][sampleidx] = exported.data
            titles.append(title)

            if first:
                self.plotinfo = exported.plotinfo
                first = False

        inknames = set(inknames)

        for inkname in inknames:
            if inkrename:
                self.labels[inkname] = inkrename[inkname]
            else:
                self.labels[inkname] = inkname

        # If titles are identical, groovy
        # Otherwise, print error
        if len(inknames) > 1 or len(set(titles)) is not 1:
            self.title = '_'.join(name for name in inknames)
            print('More than one inktype!!')
        else:
            self.title = titles[0]


    def plot_batch(self, x=None, y=None, 
        average=False, confidence=None, 
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
            for idx, inkname in enumerate(self.alldata):
                data = [ (self.alldata[inkname][sample][x], \
                    self.alldata[inkname][sample][y]) \
                    for sample in self.alldata[inkname]
                ]

                indep, mean, std, lcl, ucl = \
                    utilities.batch_average_plot(data=data)

                if not confidence:
                    ax.errorbar(
                        indep,
                        mean,
                        yerr=std,
                        color=plt.cm.tab10(coloridx[idx]),
                        marker='.', 
                        markersize=8, 
                        linewidth=3,
                        capsize=10, 
                        elinewidth=3, 
                        markeredgewidth=3, 
                        label=inkname,
                        )
                else:
                    ax.plot(
                        self.alldata[sample][x], 
                        self.alldata[sample][y],
                        linewidth=3, 
                        color='k', 
                        alpha=0.2
                    )

                    ax.plot(
                        indep, 
                        mean, 
                        linewidth=3, 
                        color='r', 
                        label='mean'
                    )
                    ax.plot(
                        indep, 
                        ucl, 
                        color='r', 
                        linewidth=3, 
                        linestyle='--',
                    )
                    ax.plot(indep, 
                        lcl, 
                        color='r', 
                        linewidth=3, 
                        linestyle='--', 
                        label=str(confidence)[2:]+ 
                            '% Confidence Interval (t-test)'
                    )

        else:
            idx = 0

            for inkname in self.data[inkname]:
                for sample in self.alldata[inkname]:
                    x_plot = np.array(self.alldata[inkname][sample][x])
                    y_plot = np.array(self.alldata[inkname][sample][y])

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

                    idx += 1

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


    def plot_batch_vis_shearstress(self, average=False, confidence=None, 
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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18,10), dpi=100)

        for yval, ax in zip(('viscosity', 'shearstress'), (ax1, ax2)):
            if average:
                for idx, inkname in enumerate(self.alldata):
                    data = [ (self.alldata[inkname][sample]['shearrate'], \
                        self.alldata[inkname][sample][yval]) \
                        for sample in self.alldata[inkname]
                    ]

                    indep, mean, std, lcl, ucl = \
                        utilities.batch_average_plot(data=data)

                    mask = np.isfinite(mean)

                    if not confidence:
                        ax.errorbar(
                            indep[mask],
                            mean[mask],
                            yerr=std[mask],
                            color=plt.cm.tab10(coloridx[idx]),
                            marker='.', 
                            markersize=8, 
                            linewidth=3,
                            capsize=10, 
                            elinewidth=3, 
                            markeredgewidth=3, 
                            label=self.labels[inkname],
                            )
                    else:
                        ax.plot(
                            self.alldata[sample]['shearrate'], 
                            self.alldata[sample][yval],
                            linewidth=3, 
                            color='k', 
                            alpha=0.2
                        )

                        ax.plot(
                            indep, 
                            mean, 
                            linewidth=3, 
                            color='r', 
                            label='mean'
                        )
                        ax.plot(
                            indep, 
                            ucl, 
                            color='r', 
                            linewidth=3, 
                            linestyle='--',
                        )
                        ax.plot(indep, 
                            lcl, 
                            color='r', 
                            linewidth=3, 
                            linestyle='--', 
                            label=str(confidence)[2:]+ 
                                '% Confidence Interval (t-test)'
                        )

            else:
                idx = 0

                for inkname in self.data[inkname]:
                    for sample in self.alldata[inkname]:
                        x_plot = np.array(
                            self.alldata[inkname][sample]['shearrate'])
                        y_plot = np.array(self.alldata[inkname][sample][yval])

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

                        idx += 1

            ax.set_xlabel(self.plotinfo['shearrate']['label'])
            ax.set_ylabel(self.plotinfo[yval]['label'])
            ax.legend()
            # ax.grid()

        if title:
            plt.suptitle(title)
        else:
            plt.suptitle('batch_viscosity_shearstress_' + self.title)

        plt.tight_layout()

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('batch_viscosity_shearstress_' \
                    + self.title + '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)









































