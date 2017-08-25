''' Module for analyzing results from X-ray Diffractometer
Raw data may involve pre-processing

Author: Bernard Kim
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley
'''

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

    def plotter(self, title=None, save=False, smoothed=False):

        ''' Plots XRD data
        Assumes input is two lists as defined by self.scan in rigaku_process
        '''

        plt.plot(self.scan['angle'], self.scan['intensity'],
                 marker='.', markersize=5,
                 color='b')
        if smoothed:
            plt.plot(self.smoothed['angle'], self.smoothed['intensity'],
                     color='r')

        plt.xlabel(r'$2\theta$')
        plt.ylabel('Intensity')
        plt.title(str(title[:-4]))

        if save:
            plt.savefig('export_' + title[:-4] + '.eps', format='eps')
            plt.clf()
        else:
            plt.show()

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

        
