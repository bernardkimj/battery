''' Module for analyzing results retrieved from Dektak profilometer

Author: Bernard Kim
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley

'''

from battery.utilities import utilities
from operator import itemgetter
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
            'x': np.array(x),
            'height': np.array(height),
        }

        self.title = filename[:-4]


    def find_sample(self, plot=False, save=False, show=False):
        '''
        Finds where the sample actually starts and ends
        Returns index ranges within self.line['x']
        Uses window-gradient with line fit
        '''

        location, height = self.line['x'], self.line['height']

        endlength = 0.25 # length of end of profile to check
        lbsize = 0.025 # window lower bound half size for polyfit on edge
        ubsize = 0.05 # window upper bound half size for polyfit on edge
        windownum = 100 # number of windows to iterate through
        positionnum = 100 # number of positions to iterate through

        front_stop = int(endlength*len(location)) # right index of front range
        back_stop = int(len(location)*(1-endlength)) # left index of back range

        lowerbound = np.round(lbsize*endlength*len(location)) # smallest window size
        upperbound = np.round(ubsize*endlength*len(location)) # largest window size
        windows = np.linspace(lowerbound, upperbound, windownum) # array of window sizes
        # windows = np.arange(lowerbound, upperbound) # array of window sizes

        front = (0, front_stop) # beginning and ending indices of front edge
        back = (back_stop, len(location)-1) # beginning and ending indices of back ege

        allvals = {}

        for edge in [front, back]:
            edge_vals, flat_vals = [], []

            for window in windows:
                values = []

                # array of positions of center of window
                # positions = np.arange(edge[0]+window, edge[1]-window)
                positions = np.linspace(edge[0]+window, edge[1]-window, positionnum)

                for position in positions:
                    sub_loc = location[int(position-window):int(position+window)]
                    sub_height = height[int(position-window):int(position+window)]

                    # apply linear fit to section within window
                    coefs = np.polyfit(sub_loc, sub_height, deg=1)
                    line = np.polyval(coefs, location[edge[0]:edge[1]])
                    residuals = np.sum((line-height[edge[0]:edge[1]])**2)

                    values.append((residuals, coefs[0], coefs[1]))

                values.sort(key=lambda tup:np.abs(tup[1]))

                # for flat line, add tuple with slope closest to zero
                flat_vals.append(values[0])

                if edge == front: # for front, tuple with most positive slope
                    edge_vals.append(max(values, key=itemgetter(1)))
                elif edge == back: # for back, tuple with most negative slope
                    edge_vals.append(min(values, key=itemgetter(1)))

            # get set of coefficients for flat and edge lines with lowest 
            # residual score
            flat_coefs = min(flat_vals, key=itemgetter(0))[1:3]
            edge_coefs = min(edge_vals, key=itemgetter(0))[1:3]

            line_loc = location[edge[0]:edge[1]]
            flat_line = np.polyval(flat_coefs, line_loc)
            edge_line = np.polyval(edge_coefs, line_loc)

            diff = 1
            for loc, flat, steep in zip(line_loc, flat_line, edge_line):
                if np.abs(flat-steep) < diff:
                    diff = np.abs(flat-steep)
                    xval = loc

            idx = np.where(location==xval)[0][0]

            if edge == front:
                key = 'front'
            elif edge == back:
                key = 'back'

            allvals[key] = {
                'idx': idx,
                'line_loc': line_loc,
                'flat_line': flat_line,
                'edge_line': edge_line,
            }

        if plot:
            font = {'family': 'Arial', 'size': 24}
            matplotlib.rc('font', **font)

            fig, ax = plt.subplots(figsize=(16,9), dpi=75)

            ax.plot(location, height,linewidth=2, color='k')
            ymin, ymax = ax.get_ylim()

            for edge in allvals:
                ax.plot(allvals[edge]['line_loc'], 
                    allvals[edge]['edge_line'], linewidth=2, color='r')
                ax.plot(allvals[edge]['line_loc'], 
                    allvals[edge]['flat_line'], linewidth=2, color='g')

            ax.set_ylim([ymin, ymax])

            ax.set_xlabel('Horizontal Position [mm]')
            ax.set_ylabel('Height [µm]')
            ax.set_title(self.title+'_edges')
            ax.grid()

            if save:
                plt.savefig(self.title + '_edges.png')
            if show:
                plt.show()

        return (allvals['front']['idx'], allvals['back']['idx'])


    def find_sample_slope(self, plot=False, save=False, show=False):
        '''
        Finds where the sample actually starts and ends
        Returns index ranges within self.line['x']
        Uses window-gradient with line fit
        '''

        location, height = self.line['x'], self.line['height']

        endlength = 0.25 # length of end of profile to check
        slopelen = 0.05 # size of half window of profile end for slope
        lbsize = 0.025 # window lower bound half size for polyfit on edge
        ubsize = 0.1 # window upper bound half size for polyfit on edge
        windownum = 100 # number of windows to iterate through
        positionnum = 100 # number of positions to iterate through

        front_stop = int(endlength*len(location)) # right index of front range
        back_stop = int(len(location)*(1-endlength)) # left index of back range

        front = (0, front_stop) # indices of front edge
        back = (back_stop, len(location)-1) # indices of back ege

        allvals = {}

        for edge in [front, back]:
            slopewindow = int(slopelen*endlength*len(location)/2)
            slope_pos = np.arange(edge[0]+slopewindow, edge[1]-slopewindow)
            slopes, slope_loc = [], [] # index numbers

            for spos in slope_pos:
                sp1, sp2 = int(spos-slopewindow), int(spos+slopewindow)
                slope = (height[sp2]-height[sp1])/(location[sp2]-location[sp1])
                slopes.append(slope)
                slope_loc.append(location[spos])

            lowerbound = np.round(lbsize*len(slopes)) # smallest window size
            upperbound = np.round(ubsize*len(slopes)) # largest window size
            windows = np.linspace(lowerbound, upperbound, windownum) 
            # array of window sizes
            # windows = np.arange(lowerbound, upperbound) # array of window sizes

            edge_vals, flat_vals = [], []

            for window in windows:
                values = []

                # array of positions of center of window
                # positions = np.arange(edge[0]+window, edge[1]-window)
                positions = np.linspace(int(window), 
                    int(len(slopes)-window), positionnum) # index numbers

                for pos in positions:
                    p1, p2 = int(pos-window), int(pos+window)
                    sub_loc = slope_loc[p1:p2]
                    sub_slope = slopes[p1:p2]

                    # sub_height = height[p1:p2]

                    # apply linear fit to section within window
                    # coefs = np.polyfit(sub_loc, sub_height, deg=1)
                    coefs = np.polyfit(sub_loc, sub_slope, deg=1)
                    line = np.polyval(coefs, slope_loc)
                    residuals = np.sum((line-slopes)**2)

                    values.append((residuals, coefs[0], coefs[1]))

                values.sort(key=lambda tup:np.abs(tup[1]))

                # for flat line, add tuple with slope closest to zero
                flat_vals.append(values[0])

                if edge == front: # for front, tuple with most positive slope
                    edge_vals.append(max(values, key=itemgetter(1)))
                elif edge == back: # for back, tuple with most negative slope
                    edge_vals.append(min(values, key=itemgetter(1)))

            # get set of coefficients for flat and edge lines with lowest 
            # residual score
            flat_coefs = min(flat_vals, key=itemgetter(0))[1:3]
            edge_coefs = min(edge_vals, key=itemgetter(0))[1:3]

            flat_line = np.polyval(flat_coefs, slope_loc)
            edge_line = np.polyval(edge_coefs, slope_loc)

            diff = 1
            for loc, flat, steep in zip(slope_loc, flat_line, edge_line):
                if np.abs(flat-steep) < diff:
                    diff = np.abs(flat-steep)
                    xval = loc

            idx = np.where(location==xval)[0][0]

            if edge == front:
                key = 'front'
            elif edge == back:
                key = 'back'

            allvals[key] = {
                'idx': idx,
                'slopes': slopes,
                'slope_loc': slope_loc,
                'flat_line': flat_line,
                'edge_line': edge_line,
            }

        if plot:
            font = {'family': 'Arial', 'size': 24}
            matplotlib.rc('font', **font)

            fig, ax = plt.subplots(figsize=(16,9), dpi=75)

            # ax.plot(location, slopes,linewidth=2, color='k')

            for edge in allvals:
                ax.plot(allvals[edge]['slope_loc'], 
                    allvals[edge]['slopes'], linewidth=2, color='k')
                ax.plot(allvals[edge]['slope_loc'], 
                    allvals[edge]['edge_line'], linewidth=2, color='r')
                ax.plot(allvals[edge]['slope_loc'], 
                    allvals[edge]['flat_line'], linewidth=2, color='g')

            ax.set_xlabel('Horizontal Position [mm]')
            ax.set_ylabel('Slope [µm/mm]')
            ax.set_title(self.title+'_slopefit')
            ax.grid()

            if save:
                plt.savefig(self.title + '_slopefit.png')
            if show:
                plt.show()

        return (allvals['front']['idx'], allvals['back']['idx'])


    def plot_profile(self, xlim=None, ylim=None, 
        title=None, show=False, save=False, savename=None):
        ''' Plots profile '''

        idx = self.find_sample(plot=True, show=show, save=save)

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        ax.plot(self.line['x'], self.line['height'],
            linewidth=2, color='k', alpha=0.5)
        ax.plot(self.line['x'][idx[0]:idx[1]], 
            self.line['height'][idx[0]:idx[1]],
            linewidth=2, color='k',)

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
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title + '_plot' + '.png')

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
                    linewidth=3, color='k', alpha=0.2)

            ax.plot(x, mean, linewidth=3, color='r', label='mean')
            ax.plot(x, ucl, color='r', linewidth=3, linestyle='--',)
            ax.plot(x, lcl, color='r', linewidth=3, linestyle='--', 
                label=str(confidence)[2:]+'% Confidence Interval (t-test)')

        else:
            for idx, line in enumerate(self.lines):
                ax.plot(self.lines[line]['x'], self.lines[line]['height'],
                    linewidth=3, color=plt.cm.tab10(coloridx[idx]), 
                    label=str(line))

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












































