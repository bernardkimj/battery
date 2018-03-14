''' Module containing common data processing scripts used across
    multiple modules

    Author: Bernard Kim
    Principal Investigators: Prof. Paul Wright, Prof. James Evans
    University: University of California, Berkeley
'''

import numpy as np
import re
import csv
import scipy.stats as stats

class utilities:

    ''' Common functions and data processing scripts '''

    def __init__(self):
        pass

    def batch_average_plot(data, confidence=0.95):
        ''' Shortens x and y vectors to length of shortest vector in group
            Returns point-by-point y-mean, std, lcl, ucl, calculated using t-dist
                instead of z-dist (for small samples)
            Data is list of tuples of lists (x,y)
            Number of x and y vectors must be equal! '''

        # Separates x and y vectors into two different objects
        x = [sample[0] for sample in data]
        y = [sample[1] for sample in data]

        if len(x) != len(y):
            raise ValueError('batch_average_prep: different number of x and y vectors')

        maxlen = 0
        for sample in y:
            if len(sample) > maxlen:
                maxlen = len(sample)

        for sample in x+y:
            if len(sample) < maxlen:
                numnan = maxlen-len(sample)
                for num in range(numnan):
                    sample.append(np.nan)

        indep = np.array([np.mean(row) for row in np.array(x).T])

        dataarray = np.array(y).T
        mean = np.array([np.mean(row) for row in dataarray])
        std = np.array([np.std(row) for row in dataarray])

        lcl, ucl = [], []

        if confidence:
            for m, s in zip(mean, std):
                R = stats.t.interval(alpha=confidence, df=len(y)-1, loc=m, scale=s/np.sqrt(len(y)))
                lcl.append(R[0])
                ucl.append(R[1])

        return indep, mean, std, lcl, ucl

    # def figure_formatting(xlim=xlim, ylim=ylim, title=title, 
    #     show=show, save=save, savename=savename, imagetype=imagetype):

    #     font = {'family': 'Arial', 'size': 28}
    #     matplotlib.rc('font', **font)
    #     fig, ax = plt.subplots(figsize=(16,9), dpi=75)

    #     if xlim:
    #         ax.set_xlim(xlim)
    #     if ylim:
    #         ax.set_ylim(ylim)

    #     if title:
    #         ax.set_title(title)
    #     else:
    #         ax.set_title()



