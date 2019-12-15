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

    def __init__(self, filename=None, isolate=True, 
        showalg=False, savealg=False):
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
                        x.append(float(row[0]))
                        # x.append(float(row[0]) * 0.001)
                        height.append(float(row[1]) * 0.001)

            except Exception:
                raise

        # Save data and convert measurements to mm and um for x and height 
        # respectively
        self.line = {
            'x': np.array(x), # um (for now)
            'height': np.array(height), # um
        }

        self.title = filename[:-4]

        self.idxs = self.find_sample(isolate=isolate, 
            showalg=showalg, savealg=savealg)
        self.get_roughness_waviness()


    def find_sample(self, isolate=True, showalg=False, savealg=False):
        '''
        Finds where the sample actually starts and ends
        Returns index ranges within self.line['x']
        Uses window-gradient with line fit
        '''

        location, height = self.line['x'], self.line['height']

        length = 0.30 # length of end of profile to check
        # window lower/upper bound half size for polyfit
        bounds = (0.025, 0.075)

        allvals = {}

        for edge in ['front', 'back']:
            x_edge, y_edge = self.get_edge(x=location, y=height, 
                length=length, end=edge)

            if edge == 'front':
                slope_sign = 1
            elif edge == 'back':
                slope_sign = -1

            flat_coefs = self.window_gradient(x=x_edge, y=y_edge, 
                bounds=bounds, slope=0)
            steep_coefs = self.window_gradient(x=x_edge, y=y_edge, 
                bounds=bounds, slope=slope_sign)

            flat_line = np.polyval(flat_coefs, x_edge)
            steep_line = np.polyval(steep_coefs, x_edge)

            diff = 1
            for loc, flat, steep in zip(x_edge, flat_line, steep_line):
                if np.abs(flat-steep) < diff:
                    diff = np.abs(flat-steep)
                    xval = loc
            idx = np.where(location==xval)[0][0]

            allvals[edge] = {
                'idx': idx,
                'line_x': x_edge,
                'flat_line': flat_line,
                'steep_coefs': steep_coefs,
                'steep_line': steep_line,
            }

        if isolate:
            # center region of isolated sample to take mean
            meanregion = 0.80
            # length of overlap ends for inner intersection
            overlaplength = 0.30

            sub_idx = (allvals['front']['idx'], allvals['back']['idx'])
            sub_x = location[sub_idx[0]:sub_idx[1]]
            sub_y = height[sub_idx[0]:sub_idx[1]]

            x_mid, y_mid = self.get_edge(x=sub_x, y=sub_y, 
                length=(1-meanregion)/2, end='middle')
            mean_coefs = np.polyfit(x_mid, y_mid, deg=1)

            for edge in ['front', 'back']:
                x_overlap, y_overlap = self.get_edge(x=location, y=height, 
                    length=overlaplength, end=edge)

                mid_line = np.polyval(mean_coefs, x_overlap)
                steep_line = np.polyval(allvals[edge]['steep_coefs'], 
                    x_overlap)

                diff = 1
                for loc, mid, steep in zip(x_overlap, mid_line, steep_line):
                    if np.abs(mid-steep) < diff:
                        diff = np.abs(mid-steep)
                        xval = loc
                mid_idx = np.where(location==xval)[0][0]

                allvals[edge]['mid_x'] = x_overlap
                allvals[edge]['mid_line'] = mid_line
                allvals[edge]['idx'] = mid_idx

        if showalg or savealg:
            font = {'family': 'Arial', 'size': 24}
            matplotlib.rc('font', **font)

            fig, ax = plt.subplots(figsize=(16,9), dpi=75)

            ax.plot(location*0.001, height,linewidth=2, color='k')
            ymin, ymax = ax.get_ylim()

            for edge in allvals:
                ax.plot(allvals[edge]['line_x']*0.001, 
                    allvals[edge]['steep_line'], linewidth=2, color='r')
                ax.plot(allvals[edge]['line_x']*0.001, 
                    allvals[edge]['flat_line'], linewidth=2, color='g')

                if isolate:
                    ax.plot(allvals[edge]['mid_x']*0.001,
                        allvals[edge]['mid_line'], 
                        linewidth=2, color='b')

            ax.set_ylim([ymin, ymax])

            ax.set_xlabel('Horizontal Position [mm]')
            ax.set_ylabel('Height [µm]')
            ax.set_title(self.title+'_edges')
            ax.grid()

            if savealg:
                plt.savefig(self.title + '_edges.png')
            if showalg:
                plt.show()

        return (allvals['front']['idx'], allvals['back']['idx'])


    def get_edge(self, x, y, length, end):
        ''' Isolates ends of raw profile data on which to run window gradient
            algorithm
            Returns x and y vectors corresponding to extracted ends
        '''

        if end == 'front':
            stop = int(length*len(x)) # right index of front range
            idxs = (0, stop) # beginning and ending indices of front range
        elif end == 'back':
            stop = int(len(x)*(1-length)) # left index of back range
            idxs = (stop, len(x)-1) # beginning and ending indices of back edge
        elif end == 'middle':
            leftstop = int(length*len(x)) # left index of middle range
            rightstop = int(len(x)*(1-length)) # right index of middle range
            idxs = (leftstop, rightstop)

        x_end = x[idxs[0]:idxs[1]]
        y_end = y[idxs[0]:idxs[1]]

        return x_end, y_end


    def window_gradient(self, x, y, bounds, slope):
        ''' Applies a window-gradient algorithm to find optimal linear fit
            Provide x, y data, window size boundaries, sign of slope to find
            Returns coefficients of optimal linear fit
        '''

        windownum = 50 # number of windows to iterate through
        positionnum = 50 # number of positions to iterate through

        lbound = np.round(min(bounds)*len(x)) # smallest window size
        ubound = np.round(max(bounds)*len(x)) # largest window size
        # array of window sizes
        windows = np.linspace(lbound, ubound, windownum)

        line_vals = []

        for window in windows:
            values = []

            # array of indices of center of window
            positions = np.linspace(window, len(x)-window, positionnum)

            for position in positions:
                win_x = x[int(position-window):int(position+window)]
                win_y = y[int(position-window):int(position+window)]

                # apply linear fit to section within window
                coefs = np.polyfit(win_x, win_y, deg=1)
                test_line = np.polyval(coefs, x)
                residuals = np.sum((test_line-y)**2)

                # store intermediate values
                values.append((residuals, coefs[0], coefs[1]))

            # depending on desired line slope, save best value for window
            if np.sign(slope) == 0:
                values.sort(key=lambda tup:np.abs(tup[1]))
                line_vals.append(values[0])
            elif np.sign(slope) == 1:
                line_vals.append(max(values, key=itemgetter(1)))
            elif np.sign(slope) == -1:
                line_vals.append(min(values, key=itemgetter(1)))

        # get set of coefficients for lines with lowest residual score
        line_coefs = min(line_vals, key=itemgetter(0))[1:3]

        return(line_coefs)


    def find_sample_slope_LEGACY(self, plot=False, save=False, show=False):
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
            # windows = np.arange(lowerbound, upperbound)

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


    def plot_profile(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, savename=None):
        ''' Plots profile '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        # ax.plot(self.line['x']*0.001, self.line['height'],
        #     linewidth=2, color='k',)

        ax.plot(self.line['x']*0.001, self.line['height'],
            linewidth=2, color='k', alpha=0.3)
        ax.plot(self.line['x'][self.idxs[0]:self.idxs[1]]*0.001, 
            self.line['height'][self.idxs[0]:self.idxs[1]],
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
                plt.savefig(self.title + '_plot' + '.svg')

        if show:
            plt.show()

        plt.close(fig)


    def get_roughness_waviness(self):
        ''' Calculates roughness and waviness of 1-D profile
            Uses FFT
            Adapted from code by Rich Winslow
        '''

        cutoff = 80 # µm

        primary_x = self.line['x'][self.idxs[0]:self.idxs[1]]
        primary = self.line['height'][self.idxs[0]:self.idxs[1]]

        samplewidth = self.line['x'][self.idxs[1]] - \
            self.line['x'][self.idxs[0]]
        samplelength = len(primary)

        flipped = primary[::-1]

        extended = np.concatenate((flipped, primary, flipped))
        f = np.array(np.fft.fft(extended))
        f[1:-1] = f[1:-1]*2

        wavelengths = [
            2*samplewidth/N for N in range(1, samplelength)
        ]

        stop_index = 0
        while wavelengths[stop_index] > cutoff:
            stop_index += 1

        filtered = f
        filtered[stop_index:-1] = 0

        ifft_result = np.real(np.fft.ifft(filtered))
        waviness = ifft_result[samplelength:2*samplelength]

        roughness = primary - waviness

        self.primary = primary
        self.wavelengths = wavelengths
        self.waviness = waviness
        self.roughness = roughness

        self.zero_average_waviness()
        self.calculate_metrics()


    def plot_roughness_waviness(self, xlim=None, ylim=None, title=None,
        show=False, save=False, savename=None):
        ''' Plots roughness and waviness with original trace '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        isolatedtrace = self.line['x'][self.idxs[0]:self.idxs[1]]*0.001

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        ax.plot(self.line['x']*0.001, self.line['height'],
            linewidth=2, color='k', alpha=0.3)
        ax.plot(isolatedtrace, self.line['height'][self.idxs[0]:self.idxs[1]],
            linewidth=2, color='k', label=r'$P_a = $'+'%.3f µm'%self.metrics['P']['Pa'])

        ax.plot(isolatedtrace, self.waviness_zero_avg, linewidth=2, color='b',
            label=r'$W_q = $'+'%.3f µm'%self.metrics['W']['Wq'])
        ax.plot(isolatedtrace, self.roughness, linewidth=2, color='r',
            label=r'$R_q = $'+'%.3f µm'%self.metrics['R']['Rq'])

        ax.legend()
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
            ax.set_title(self.title + '_RW')

        if save:
            if savename:
                plt.savefig(savename + '.svg')
            else:
                plt.savefig(self.title + '_RW' + '.svg')

        if show:
            plt.show()

        plt.close(fig)


    def zero_average_waviness(self):
        ''' Zero-averages the waviness so the center of the waviness
            profile is at zero 

            Done by taking mean of waviness profile and subtracting mean
            at each point. '''

        waviness_avg = np.mean(self.waviness)
        self.waviness_zero_avg = self.waviness - waviness_avg


    def calculate_metrics(self):
        ''' Calculates metrics for roughness and waviness per trace 

        Xa - arithmetic average
        Xq - root mean square
        Xp - maximum height of peaks
        Xv - maximum depth of valleys
        Xz - maximum vertical height difference

        '''

        self.metrics = {}

        for attr in ['P', 'W', 'R']:
            if attr == 'P':
                metric = self.primary
            elif attr == 'W':
                metric = self.waviness_zero_avg
            elif attr == 'R':
                metric = self.roughness

            a = np.mean(np.abs(metric))
            q = np.sqrt(np.sum(np.square(metric))/len(metric))
            p = np.max(metric)
            v = np.min(metric)
            z = np.abs(p) + np.abs(v)

            self.metrics[attr] = {
                attr+'a': a,
                attr+'q': q,
                attr+'p': p,
                attr+'v': v,
                attr+'z': z,
            }


    def level_trace(self, start=None, end=None, show=False, savecsv=False):
        ''' Levels out trace with manually entered start and endpoints 
            Can show new plot, but will not save
            Can save adjusted trace as csv
            '''

        if not start or not end:
            startidx = 0
            endidx = len(self.line['x'])-1
        elif start and end:
            startidx = np.where(self.line['x']==start)[0][0]
            endidx = np.where(self.line['x']==end)[0][0]

        points = [
            (self.line['x'][startidx], self.line['height'][startidx]),
            (self.line['x'][endidx], self.line['height'][endidx]),
        ]

        c1, c0 = np.polyfit(
            np.array([points[0][0], points[1][0]]), 
            np.array([points[0][1], points[1][1]]), 
            deg=1
        )

        level_line = np.array([c1*x+c0 for x in self.line['x']])
        new_trace = np.array( [h-h_new for h, h_new in \
            zip(self.line['height'], level_line)] )

        self.line['height'] = new_trace

        if show:
            fig, ax = plt.subplots()
            ax.plot(
                self.line['x'],
                self.line['height'],
            )

            plt.show()

        if savecsv:
            with open(self.title+'_leveled.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',',)

                for x, y in zip(self.line['x'], self.line['height']):
                    writer.writerow([x, y*1000])


class profile_batch:
    '''
    Class for working with multiple profiles at once
    Relies on methods defined in 'profile' class
    '''

    def __init__(self, alldata=None, include=['casting']):
        ''' 
        Initializes files according 'profile'
        '''

        # Acceptes lists of filenames
        self.lines = {}
        self.idxs = {}
        self.metrics = {}
        titles = []

        for file in alldata:
            exported = profile(file)
            title, inkname, samplenum = utilities.title_search_format(filename=file, 
                include=include)
            self.lines[samplenum] = exported.line
            self.idxs[samplenum] = exported.idxs
            self.metrics[samplenum] = exported.metrics
            titles.append(title)

        # If titles are identical, groovy
        # Otherwise, print error
        self.title = titles[0]

        if len(set(titles)) is not 1:
            print('Titles do not match!!')

        utilities.save_json(data=self.metrics, filename=self.title+'.json')


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
            rous, wavs, pa = [], [], []

            for idx in self.lines:
                rous.append(self.metrics[idx]['R']['Rq'])
                wavs.append(self.metrics[idx]['W']['Wq'])
                pa.append(self.metrics[idx]['P']['Pa'])

                start = self.idxs[idx][0]
                end = self.idxs[idx][1]

                isolatedtrace = self.lines[idx]['x'][start:end]*0.001 - \
                    self.lines[idx]['x'][start]*0.001

                ax.plot(
                    isolatedtrace, 
                    self.lines[idx]['height'][start:end],
                    linewidth=2, 
                    color=plt.cm.tab10(coloridx[int(idx[1:])-1]), 
                    label=str(idx)
                )

            mean_r = np.mean(rous)
            mean_w = np.mean(wavs)
            mean_a = np.mean(pa)

            stats = {
                'r': r'$R_q = $'+'%.3f µm'%mean_r,
                'w': r'$W_q = $'+'%.3f µm'%mean_w,
                'a': r'$P_a = $'+'%.3f µm'%mean_a,
            }

            stats_text = stats['r']+'\n'+stats['w']+'\n'+stats['a']

            x_tlim, y_tlim = ax.get_xlim(), ax.get_ylim()

            x_tpos = 0.05*(x_tlim[1]-x_tlim[0]) + x_tlim[0]
            y_tpos = 0.80*(y_tlim[1]-y_tlim[0]) + y_tlim[0]

            ax.text(x_tpos, y_tpos, stats_text, bbox=dict(facecolor='w', 
                edgecolor='k'))

        ax.set_xlabel('Horizontal Position [mm]')
        ax.set_ylabel('Height [µm]')
        ax.legend(loc='upper right')
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


class profile_combined:
    '''
    Class for plotting both scan directions for one sample on same figure
    Relies on methods defined in 'profile' class
    '''

    def __init__(self, files=None, include=['casting']):
        ''' 
        Initializes files according 'profile'
        '''

        # Acceptes lists of filenames
        self.lines = {}
        self.idxs = {}
        self.metrics = {}
        self.titles = {}

        for filename in files:
            exported = profile(filename)

            dir_search = re.search(r'pe|pa', filename)
            if dir_search.group(0) == 'pe':
                direction = 'perp'
            elif dir_search.group(0) == 'pa':
                direction = 'par'

            title, inkname, samplenum = utilities.title_search_format(
                filename=filename, include=include)

            self.lines[direction] = exported.line
            self.idxs[direction] = exported.idxs
            self.metrics[direction] = exported.metrics
            self.titles[direction] = title


    def plot_profile_combined(self, ylim=None, title=None, 
        show=False, save=False, savename=None, imagetype='png'):
        '''
        Plots multiple profiles on same plot
        '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
            figsize=(18,10), dpi=100)

        for direction, ax in zip(('par', 'perp'), (ax1, ax2)):
            start = self.idxs[direction][0]
            end = self.idxs[direction][1]

            midpoint = int((end-start)/2 + start)
            x_mid = self.lines[direction]['x'][midpoint]

            x_plot = (self.lines[direction]['x'] - x_mid)*0.001

            # full trace, 20% opacity
            ax.plot(
                x_plot,
                self.lines[direction]['height'],
                color='#000000',
                alpha=0.3,
                linewidth=3,
            )

            # isolated trace, full opacity
            ax.plot(
                x_plot[start:end],
                self.lines[direction]['height'][start:end],
                color='#000000',
                linewidth=3,
            )

            ax.set_ylabel('Height [µm]')

        ax1.set_xlim([-6, 6])

        if ylim:
            ax.set_ylim(ylim)
        else:
            y_lim = ax1.get_ylim()
            ax1.set_ylim([-10, y_lim[1]*1.1])
            ax2.set_ylim([-10, y_lim[1]*1.1])

        for direction, ax in zip(('par', 'perp'), (ax1, ax2)):
            stats = {
                'r': r'$R_q = $'+'%.2f µm'%self.metrics[direction]['R']['Rq'],
                'w': r'$W_q = $'+'%.2f µm'%self.metrics[direction]['W']['Wq'],
                'a': r'$P_a = $'+'%.2f µm'%self.metrics[direction]['P']['Pa'],
            }

            stats_text = stats['r']+'\n'+stats['w']+'\n'+stats['a']

            x_tlim, y_tlim = ax.get_xlim(), ax.get_ylim()

            x_tpos = 0.02*(x_tlim[1]-x_tlim[0]) + x_tlim[0]
            y_tpos = 0.95*(y_tlim[1]-y_tlim[0]) + y_tlim[0]

            ax.text(
                x_tpos, 
                y_tpos, 
                stats_text, 
                bbox=dict(
                    facecolor='w', 
                    edgecolor='k',
                ),
                verticalalignment='top',
                horizontalalignment='left',
                fontsize=20,
            )

        ax1.set_title('Parallel')
        ax2.set_title('Perpendicular')
        ax2.set_xlabel('Horizontal Position [mm]')

        if title:
            fig.suptitle(title)

        # plt.tight_layout()


        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('combined_' + self.titles['par'] + '.' + imagetype)

        if show:
            plt.show()

        plt.close(fig)









































