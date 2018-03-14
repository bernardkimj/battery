''' Module for analyzing results from X-ray Diffractometer
Raw data may involve pre-processing

Author: Bernard Kim
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley
'''

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

class xrd:

    ''' Analyzes and plots data from XRD

    Plots data from each XRD run
    *** ADD FUNCTIONALITY FOR BATCH PLOTTING ***
    '''

    def __init__(self, filename=None):
        ''' Processes raw data from Rigaku XRD for furhter use
        XRD as used in MSE 104 instructional lab
        '''

        self.scan = {'angle': [], 'intensity': []}

        with open(filename) as f:
            rows = f.readlines()

            for row in rows:
                if len(row) != 1:
                    line = row.split(',')
                    self.scan['angle'].append(line[0])

                    if line[1][-1] == '\n':
                        self.scan['intensity'].append(line[1][:-1])
                    else:
                        self.scan['intensity'].append(line[1])

        window_size = 15
        raw_intensity = np.asarray(self.scan['intensity'], dtype=np.float64)

        self.smoothed = {}
        self.smoothed['angle'] = self.scan['angle'][int(np.floor(window_size/2-1)):-int(np.ceil(window_size/2))]
        self.smoothed['intensity'] = self.moving_average(raw_intensity, window_size)
        self.title = filename[:-4]

    def plotter(self, ylim=None, title=None, smoothed=False,
        show=False, save=False, savename=None, imagetype='png'):

        ''' Plots XRD data
        Assumes input is two lists as defined by self.scan in rigaku_process
        '''

        # Plotting formatting
        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        ax.plot(self.scan['angle'], self.scan['intensity'],
                 marker='.', markersize=5,
                 color='b', label='raw')
        if smoothed:
            ax.plot(self.smoothed['angle'], self.smoothed['intensity'],
                     color='r', label='smoothed')
            ax.legend()

        ax.set_xlim([0,90])

        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        ax.set_xlabel(r'$2\theta$')
        ax.set_ylabel('Intensity')
        ax.grid()

        if show:
            plt.show()

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig(self.title + '.' + imagetype)

        plt.close(fig)

    def moving_average(self, interval, window_size):
        ''' Calculates moving average to smooth noise from data 
            Based on convolution, so weights must be centered on window
            '''
        n = int(window_size)
        window = np.ones(n)/float(window_size)

        return np.convolve(interval, window, 'valid')



    # def batch_plotter(self, title=None, save=False):

    #     ''' Plots XRD data from multiple runs
    #     Inputs must be processed through self
    #     '''

        
