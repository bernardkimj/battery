''' Module for analyzing results retrieved from Gamry

Author: Rich Winslow
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley

Further edits: Bernard Kim
'''

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import re

font = {'family': 'Arial', 'size': 16}
matplotlib.rc('font', **font)

class EIS:

    ''' Analyzes data from Gamry EIS

    Pulls data from .dta file for plotting and analysis. Can look at the
    Nyquist plot and determine the DC resistance for ionic conductivity.

    '''

    def __init__(self, filename=None, thickness=0.001, area=1):
        ''' Opens file and retrieves data.

        Retrieves time, frequency, real impedance, imaginary impedance,
        magnitude, and phase. Assumes that the first part of the file is an
        OCV test and that the header for the table consists of two lines.

        Unit requirements:
            R_solution [ohm]
            Thickness [cm]
            Area [cm^2]
        '''

        self.time = []
        self.freq = []
        self.real = []
        self.imag = []
        self.phaz = []
        self.magn = []

        self.r_solution = None
        self.conductivity = None

        self.filename = filename
        self.thickness = thickness
        self.area = area

        with open(filename, errors='replace') as f:
            rows = f.readlines()

            switch = False
            for index, row in enumerate(rows):
                row = row.split()
                try:
                    if row:
                        if row[0] == 'ZCURVE':
                            switch = index + 2

                        if (self.is_num(row[0]) and switch and index > switch):
                            self.time.append(float(row[1]))
                            self.freq.append(float(row[2]))
                            self.real.append(float(row[3]))
                            self.imag.append(float(row[4]))
                            self.magn.append(float(row[6]))
                            self.phaz.append(float(row[7]))
                except Exception:
                    raise

            try:
                self.calculate_conductivity()
            except Exception:
                raise

    def calculate_conductivity(self):
        try:
            max_imag_index = self.imag.index(max(self.imag))
            self.r_solution = self.real[max_imag_index]
            self.conductivity = (
                (self.thickness * 1000) / (self.area * self.r_solution))
        except Exception:
            raise

    def list_metrics(self):
        print('File: ' + self.filename)
        print('R_solution: ' + str(self.r_solution) + ' ohm')
        print('Thickness: ' + str(self.thickness) + ' cm')
        print('Area: ' + str(self.area) + ' cm^2')
        print('Conductivity: ' + str(self.conductivity) + ' mS/cm')
        print('--')

    def plot_nyquist(self, log_plot=False, ylim=None, save=False):
        ''' Plots real impedance vs negative imaginary impedance '''

        self.calculate_conductivity()
        conductivity_string = "{0:.3f}".format(self.conductivity)

        if log_plot:
            plt.loglog(self.real, [-1 * v for v in self.imag],
                       marker='.', markersize=15)
        else:
            plt.plot(self.real, [-1 * v for v in self.imag],
                     marker='.', markersize=15)

        if ylim:
            plt.ylim(ylim[0], ylim[1])

        plt.xlabel('Z_real (ohm)')
        plt.ylabel('(-) Z_imag (ohm)')

        if save:
            plt.savefig(self.filename + '_nyquist.pdf')
        else:
            plt.title('Nyquist Plot - ' +
                      conductivity_string + ' mS/cm - ' + self.filename)
            plt.show()

    def plot_bode(self, save=False):

        fig, (ax_magn, ax_phaz) = plt.subplots(2, sharex=True)
        ax_magn.semilogx(self.freq, self.magn)
        ax_phaz.semilogx(self.freq, self.phaz)

        ax_magn.set_ylabel('Magnitude [Ω]')
        ax_phaz.set_ylabel('Phase [°]')
        ax_phaz.set_xlabel('Frequency')

        if save:
            plt.savefig(self.filename + '_bode.pdf')

        plt.show()




    def is_num(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

class CV:

    ''' Analyzes data from Gamry cyclic voltammetry data

    Pulls data from .dta file for plotting and analysis.

    '''

    def __init__(self, filename=None, reduce=False):
        ''' Opens file and retrieves data. '''

    # Pt  T     Vf    Im  Vu  Sig Ach IERange Over
    # s   V vs. Ref.  A   V   V   V   #       bits

        self.cycles = {}

        with open(filename, errors='replace') as f:
            rows = f.readlines()

            switch = False
            current_cycle_time = []
            current_cycle_voltage = []
            current_cycle_current = []

            for index, row in enumerate(rows):
                row = row.split()
                try:
                    if row:
                        if row[0][0:5] == 'CURVE':
                            curve_number = int(row[0][5::])
                            switch = index + 2
                            if current_cycle_time:
                                self.cycles[curve_number-1] = {
                                    'time': current_cycle_time,
                                    'voltage': current_cycle_voltage,
                                    'current': current_cycle_current,
                                }
                                current_cycle_time = []
                                current_cycle_voltage = []
                                current_cycle_current = []

                        if (self.is_num(row[0]) and switch and index > switch):
                            # Save data and convert current to mA
                            current_cycle_time.append(float(row[1]))
                            current_cycle_voltage.append(float(row[2]))
                            current_cycle_current.append(float(row[3]) * 1000)

                except Exception:
                    raise

            # Save data and convert current to mA
            self.cycles[curve_number] = {
                'time': current_cycle_time,
                'voltage': current_cycle_voltage,
                'current': current_cycle_current * 1000,
            }

        self.title = filename[:-4]

        if reduce:
            CV.reduce_file(self)

    def is_num(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def reduce_file(self):
        ''' Reduces file size of oversampled data '''

        stepsize = 10 # ideal stepsize, mV/s

        for curve_number in self.cycles:
            if curve_number == min(list(self.cycles)) and len(self.cycles[curve_number]['time']) > 6000/stepsize + 1:
                samplerate = (len(self.cycles[curve_number]['time'])-1)/(6000/stepsize)
                CV.undersample(self, samplerate, curve_number)
            elif curve_number == max(list(self.cycles)) and len(self.cycles[curve_number]['time']) > 2000/stepsize:
                samplerate = len(self.cycles[curve_number]['time'])/(2000/stepsize)
                CV.undersample(self, samplerate, curve_number)
            elif len(self.cycles[curve_number]['time']) > 8000/stepsize:
                samplerate = len(self.cycles[curve_number]['time'])/(8000/stepsize)
                CV.undersample(self, samplerate, curve_number)

    def undersample(self, samplerate, curve_number):
        ''' Undersampling algorithm '''

        reducedidx = []
        for idx in list(range(len(self.cycles[curve_number]['time']))):
            if idx%samplerate == 0:
                reducedidx.append(idx)

        reducedtime = []; reducedvoltage = []; reducedcurrent = []
        for idx in reducedidx:
            reducedtime.append(self.cycles[curve_number]['time'][idx])
            reducedvoltage.append(self.cycles[curve_number]['voltage'][idx])
            reducedcurrent.append(self.cycles[curve_number]['current'][idx])

        self.cycles[curve_number]['time'] = reducedtime
        self.cycles[curve_number]['voltage'] = reducedvoltage
        self.cycles[curve_number]['current'] = reducedcurrent

    def find_min_max(self, list, param, extreme):
        ''' For list with dictionary entries [{}, {},], finds maximum or 
        minimum value based on specified dictionary key'''
        values = [item[param] for item in list]
        if extreme == 'max':
            max_value = max(values)
            max_idx = values.index(max_value)
            return list[max_idx]
        elif extreme == 'min':
            min_value = min(values)
            min_idx = values.index(min_value)
            return list[min_idx]

    def moving_average(self, interval, window_size, weight=None):
        ''' Calculates moving average to smooth noise from data 
            Choose weighting option based on kwarg 'weight'
            Based on convolution, so weights must be centered on window
            '''
        n = int(window_size)

        if weight == None or weight == 'unweighted':
            window = np.ones(n)/float(window_size)
        # Weighting for linear and exp weights is not properly centered, 
        # sample average is biased towards future numbers
        elif weight == 'linear':
            window = np.arange(1,n+1)/(n*(n+1)/2)
        elif weight == 'exp':
            alpha = 2/(n+1)
            coeffs = np.zeros(n)

            for idx, coeff in enumerate(coeffs):
                coeffs[idx] = (1 - alpha)**idx

            window = coeffs/np.sum(coeffs)

        return np.convolve(interval, window, 'valid')

    def plot_current_voltage(self, cycle_index=0, xlim=None, ylim=None, 
        title=None, show=False, save=False, imagetype='png'):
        ''' Plots current vs voltage for one or all cycles '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.4,1,len(self.cycles)) # for use with Blues cmap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if cycle_index:
            ax.plot(self.cycles[cycle_index]['voltage'],
                    self.cycles[cycle_index]['current'],
                    linewidth=2
                    color=plt.cm.Blues(coloridx[cycle_index-1]))
        else:
            for i in range(1,len(self.cycles)):
                ax.plot(self.cycles[i]['voltage'],
                        self.cycles[i]['current'],
                        linewidth=2,
                        color=plt.cm.Blues(coloridx[i-1]),
                        label='Cycle '+str(i))
            ax.legend()

        ax.set_xlabel('Potential [V]')
        ax.set_ylabel('Current [mA]')
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_C' + str(cycle_index))

        if show:
            plt.show()

        if save:
            plt.savefig('single_' + self.title + '_C' + str(cycle_index) + '.' + str(imagetype))

        plt.close(fig)

    def extrema_smoothed_wd(self, title=None, save_csv=False, showplot=False, saveplot=False):
        ''' Extracts maximum and minumum points on each CV curve for each cycle
            Uses gradient-window (g-w) method on smoothed curves
            Written by Bernard Kim '''

        del self.cycles[max(self.cycles)] # removes last element of self.cycles (incomplete cycles)
        # !!!!! MAY WANT TO CHANGE THIS LATER FOR LINEAR EXTRAPOLATION PURPOSES !!!!!

        voltagethresh = 2 # positive and negative hard voltage bounds, no effect when == 2
        smoothsize = 0.10 # size of smoothing window, as percent of individual scan length
        lbsize = 0.10 # window lower bound size (for g-w), as percent of scan length
        ubsize = 0.25 # window upper bound size (for g-w), as percent of scan length

        allPairs = {}

        for cycle in self.cycles:
            voltages = self.cycles[cycle]['voltage']
            currents = self.cycles[cycle]['current']
            times = self.cycles[cycle]['time']

            # Split voltage, current, and time lists into cathodic and anodic scans
            scan_switch = voltages.index(max(voltages))

            anodic = {
                'V_raw': np.array(voltages[:scan_switch]),
                'I_raw': np.array(currents[:scan_switch]),
            }

            cathodic = {
                'V_raw': np.array(voltages[scan_switch:]),
                'I_raw': np.array(currents[scan_switch:]),
            }

            allPairs[cycle] = {}

            # Smooth and run gradient-window algorithm on each scan separately
            for scan in [anodic, cathodic]:

                # Smooth current data using moving average
                smoothing_window = round(smoothsize*len(scan['V_raw']))
                scan['I_smoothed'] = self.moving_average(scan['I_raw'], smoothing_window)

                # Extract corresponding smoothed voltage points to smoothed current points
                overhang = len(scan['I_raw']) - len(scan['I_smoothed'])
                lower_pad = np.ceil(overhang/2)
                upper_pad = np.floor(overhang/2)
                scan['V_smoothed'] = scan['V_raw'][int(lower_pad):int(-upper_pad)]

                # Set window size as a function of smoothed scan data length
                lowerbound = round(lbsize*len(scan['V_smoothed']))
                upperbound = round(ubsize*len(scan['V_smoothed']))
                # 1/2 window size, symmetric positive and negative offsets
                windows = list(range(lowerbound,upperbound)) # index of center of window

                totalExtrema = [[] for window in windows]

                for window in windows:
                    # Create list of valid positions based on window size
                    positions = list(range(window,len(scan['V_smoothed'])-window))
                    points = [[] for position in positions]

                    for position in positions:
                        # Extract voltage, current at left, middle, and right window positions
                        backvoltage, backcurrent = scan['V_smoothed'][position-window], \
                                                   scan['I_smoothed'][position-window],
                        frontvoltage, frontcurrent = scan['V_smoothed'][position+window], \
                                                     scan['I_smoothed'][position+window],
                        middlevoltage, middlecurrent = scan['V_smoothed'][position], \
                                                       scan['I_smoothed'][position],

                        backslope = (middlecurrent-backcurrent)/(middlevoltage-backvoltage)
                        frontslope = (frontcurrent-middlecurrent)/(frontvoltage-middlevoltage)

                        # Create dictionary per point with slope products, middle points, and endpoint voltages
                        # ** Of SMOOTHED voltage and current, not raw
                        points[position-window] = {
                            'slopeprod': backslope*frontslope,
                            'pair': {
                                'voltage': middlevoltage, 
                                'current': middlecurrent,
                            },
                            'raw_pair':{
                                'voltage': scan['V_raw'][int(position+lower_pad)],
                                'current': scan['I_raw'][int(position+lower_pad)],
                            },
                            'endpoints': {
                                'back': backvoltage,
                                'front': frontvoltage,
                            },
                        }

                    # Fill with dummy values to prevent empty sequence errors
                    extremaPairs = [{'voltage':0, 'current':0}]

                    for point in points:
                        # Negative slope product means slopes have opposite signs and surround local extrema
                        if point['slopeprod'] < 0 and abs(point['pair']['voltage']) < voltagethresh:
                            # Local maxima/minima
                            extremaPairs.append(point['raw_pair'])

                    # Pull out point for maximum/minumum current for given window size
                    if scan is anodic:
                        windowExtreme = self.find_min_max(extremaPairs, 'current', 'max')
                    elif scan is cathodic:
                        windowExtreme = self.find_min_max(extremaPairs, 'current', 'min')

                    # Populate master list of extreme currents for all window sizes
                    totalExtrema[window - windows[0]] = windowExtreme

                # Find most optimal extrema for all window and gradient possibilities and
                # populate master list of extrema per cycle
                if scan is anodic:
                    realExtreme = self.find_min_max(totalExtrema, 'current', 'max')
                    allPairs[cycle]['anodic'] = realExtreme
                elif scan is cathodic:
                    realExtreme = self.find_min_max(totalExtrema, 'current', 'min')
                    allPairs[cycle]['cathodic'] = realExtreme

                # Compute potential difference (E_p_ox - E_p_red) and current ratio (I_p_ox/I_p_red)
                # and add to master list

                # allPairs[cycle]['']

            if showplot:
                fig, ax = plt.subplots(figsize=(16,9), dpi=60)
                ax.plot(cathodic['V_raw'], cathodic['I_raw'], color='#800000', linewidth=4, label='raw')
                ax.plot(anodic['V_raw'], anodic['I_raw'], color='#800000', linewidth=4)
                ax.plot(cathodic['V_smoothed'], cathodic['I_smoothed'], color='#005DBB', linewidth=2, label='smoothed')
                ax.plot(anodic['V_smoothed'], anodic['I_smoothed'], color='#005DBB', linewidth=2)
                ax.plot(allPairs[cycle]['cathodic']['voltage'], allPairs[cycle]['cathodic']['current'],
                    color='#00D48D', marker="o", markersize=20, markeredgewidth=2, markerfacecolor='None')
                ax.plot(allPairs[cycle]['anodic']['voltage'], allPairs[cycle]['anodic']['current'],
                    color='#00D48D', marker="o", markersize=20, markeredgewidth=2, markerfacecolor='None')
                ax.set_xlabel('Voltage [V]')
                ax.set_ylabel('Current [mA]')
                ax.legend()
                plt.title(str(title[0:-4]) + '_C' + str(cycle))
                ax.grid()

                plt.show()
                plt.close()

                if saveplot:
                    plt.savefig(title[0:-4] + '_C' + str(cycle) + 'extrema' + '.pdf', format='pdf')

        if save_csv:
            # write data to csv file
            filename = str(title)

            with open('extrema_' + filename[:-4] + '.csv', 'w', newline='\n') as f:
                fwrite = csv.writer(f, delimiter=',',quotechar='"')
                fwrite.writerow(['Cycle','E_p_ox (V)','I_p_ox (mA)','E_p_red (V)', 'I_p_red (mA)'])

                for allPair in allPairs:
                    E_p_ox = allPairs[allPair]['anodic']['voltage']
                    I_p_ox = allPairs[allPair]['anodic']['current']
                    E_p_red = allPairs[allPair]['cathodic']['voltage']
                    I_p_red = allPairs[allPair]['cathodic']['current']
                    cyclenum = allPair
                    row = [cyclenum, E_p_ox, I_p_ox, E_p_red, I_p_red]
                    fwrite.writerow(row)

                f.close

    def extrema_smoothed_radius_LEGACY(self, title=None):
        ''' Extracts maximum and minimum points on each CV curve for each cycle
            Uses gradient/window method 
            Written by Bernard Kim '''

        del self.cycles[7] # removes last element of self.cycles (incomplete cycle)
        # !!!!!! MAY WANT TO CHANGE THIS LATER FOR LINEAR EXTRAPOLATION PURPOSES !!!!!!

        allPairs = []

        for cyclenum,cycle in enumerate(self.cycles,1):
            voltages = self.cycles[cycle]['voltage']
            currents = self.cycles[cycle]['current']
            times = self.cycles[cycle]['time']

            # Split voltage, current, and time lists into cathodic and anodic scans
            scan_switch = voltages.index(max(voltages))

            cathodic = {
                'V_raw': np.array(voltages[:scan_switch]),
                'I_raw': np.array(currents[:scan_switch]),
            }

            anodic = {
                'V_raw': np.array(voltages[scan_switch:]),
                'I_raw': np.array(currents[scan_switch:]),
            }

            # Smooth current data using unweighted moving average
            smoothing_window = round(0.10*len(times)/2) # 10 percent of each scan length

            for electrode in [cathodic, anodic]:
                electrode['I_smoothed'] = self.moving_average(electrode['I_raw'], smoothing_window, 'exp')

                # Calculate curvature, kappa, at each point using first and second derivatives
                electrode['d_one'] = np.gradient(electrode['I_smoothed'])
                electrode['d_two'] = np.gradient(electrode['d_one'])
                electrode['kappa'] = (np.abs(electrode['d_two']))/((1 + electrode['d_one']**2)**(3/2))

                # Find index of maximum curvature to find points of local extrema for I_smoothed
                idx = np.argmax(electrode['kappa'])

                # Compensate for truncated vector to pull corresponding data point from raw data vectors
                overhang = len(electrode['I_raw'])-len(electrode['I_smoothed'])
                lower_pad = np.ceil(overhang/2)
                upper_pad = np.floor(overhang/2)
                electrode['V_smoothed'] = electrode['V_raw'][int(lower_pad):int(-upper_pad)]

                electrode['peak'] = {
                    'V_peak': electrode['V_raw'][int(idx+lower_pad)],
                    'I_peak': electrode['I_raw'][int(idx+lower_pad)],
                }

            allPairs.append({
                'cathodic': cathodic['peak'],
                'anodic': anodic['peak']
            })

        print(cathodic['kappa'])
        plt.plot(voltages,currents,color='b')
        plt.plot(cathodic['V_smoothed'],cathodic['I_smoothed'],color='r')
        plt.plot(anodic['V_smoothed'],anodic['I_smoothed'],color='r')
        plt.show()

        # write data to csv file
        filename = str(title)

        with open('extrema_' + filename[:-4] + '.csv', 'w', newline='\n') as f:
            fwrite = csv.writer(f, delimiter=',',quotechar='"')
            fwrite.writerow(['Cycle','E_p_red (V)', 'I_p_red (mA)','E_p_ox (V)','I_p_ox (mA)'])

            for idx,allPair in enumerate(allPairs,1):
                E_p_red = allPair['cathodic']['V_peak']
                I_p_red = allPair['cathodic']['I_peak']
                E_p_ox = allPair['anodic']['V_peak']
                I_p_ox = allPair['anodic']['I_peak']
                row = [idx, E_p_red, I_p_red, E_p_ox, I_p_ox]
                fwrite.writerow(row)

            f.close

    def extrema_raw_wd_LEGACY(self, title=None):
        ''' Extracts maximum and minimum points on each CV curve for each cycle
            Uses gradient/window method 
            Written by Bernard Kim '''

        # LEGACY CODE - uses window/gradient algorithm to find current extrema

        del self.cycles[7] # removes last element of self.cycles (incomplete cycle)
        # !!!!!! MAY WANT TO CHANGE THIS LATER FOR LINEAR EXTRAPOLATION PURPOSES !!!!!!

        voltagethresh = 1.5 # positive and negative hard voltage bounds

        allPairs = []

        for cyclenum,cycle in enumerate(self.cycles,1):
            voltages = self.cycles[cycle]['voltage']
            currents = self.cycles[cycle]['current']
            times = self.cycles[cycle]['time']

            # Set window size as a function of cycle data length
            lowerbound = round(0.10*len(times)) # set to 1% of total cycle length
            upperbound = round(0.25*len(times)) # set to 10% of total cycle length
            windows = list(range(lowerbound,upperbound)); # index of center of window
            # 1/2 window size, symmetric positive and negative offsets

            totalMax = [[] for window in windows]
            totalMin = [[] for window in windows]

            for window in windows:
                # create list of valid positions based on window size
                positions = list(range(window,len(times)-window))
                points = [[] for position in positions]

                for position in positions:
                    # pull out voltage, current at left, middle, and right window positions
                    backvoltage, backcurrent = voltages[position-window], currents[position-window]
                    frontvoltage, frontcurrent = voltages[position+window], currents[position+window]
                    middlevoltage, middlecurrent = voltages[position], currents[position]

                    if frontvoltage == middlevoltage or backvoltage == middlevoltage: # otherwise divide by 0
                        backslope = 0
                        frontslope = 0
                    else:
                        # calculate secant 1st derivative based on lower and upper combinations of points
                        backslope = (middlecurrent-backcurrent)/(middlevoltage-backvoltage)
                        frontslope = (frontcurrent-middlecurrent)/(frontvoltage-middlevoltage)

                    # for each point, create dictionary containing slope products, middle points, and 
                    # endpoint voltages
                    points[position-window] = {
                        'slopeprod': backslope*frontslope,
                        'pair': {
                            'voltage': middlevoltage, 
                            'current': middlecurrent,
                        },
                        'endpoints': {
                            'back': backvoltage,
                            'front': frontvoltage,
                        },
                    }

                # fill with dummy values to prevent empty sequence errors
                maxPairs = [{'voltage':0, 'current':0}]
                minPairs = [{'voltage':0, 'current':0}]

                for point in points:
                    # negative slope product means slopes have opposite signs and surround local extrema
                    if point['slopeprod'] < 0:
                        # local maxima, monotonically increasing voltage
                        if (
                                point['endpoints']['back'] < point['pair']['voltage'] and 
                                point['pair']['voltage'] < point['endpoints']['front'] 
                                and point['pair']['voltage'] < voltagethresh
                            ):
                            maxPairs.append(point['pair']) # adds valid maxPair to maxPairs list
                        # local minima, monotonically decreasing voltage
                        elif (
                                point['endpoints']['back'] > point['pair']['voltage'] and 
                                point['pair']['voltage'] > point['endpoints']['front'] 
                                and point['pair']['voltage'] > -voltagethresh
                            ):
                            minPairs.append(point['pair']) # adds valid minPair to minPairs list

                # pull out point for max/min current for given window size
                windowMax = max(maxPairs, key=lambda x:x['current'])
                windowMin = min(minPairs, key=lambda x:x['current'])

                # populate master list of max/min currents for all window sizes
                totalMax[window - windows[0]] = windowMax
                totalMin[window - windows[0]] = windowMin

            # find most optimal max/min for all window and gradient possibilities
            realMax = max(totalMax, key=lambda x:x['current'])
            realMin = min(totalMin, key=lambda x:x['current'])

            # populate master list of max/min per cycle (list of dictionaries with subdictionaries)
            allPairs.append(
                {
                'Max': realMax,
                'Min': realMin,
                }
            )

        # write data to csv file
        filename = str(title)

        with open('extrema_' + filename[:-4] + '.csv', 'w', newline='\n') as f:
            fwrite = csv.writer(f, delimiter=',',quotechar='"')
            fwrite.writerow(['Cycle','Max Voltage (V)', 'Max Current (mA)','Min Voltage (V)','Min Current (mA)'])

            for idx,allPair in enumerate(allPairs,1):
                maxVolt = allPair['Max']['voltage']
                maxCurr = allPair['Max']['current']
                minVolt = allPair['Min']['voltage']
                minCurr = allPair['Min']['current']
                row = [idx, maxVolt, maxCurr, minVolt, minCurr]
                fwrite.writerow(row)

            f.close

class CV_batch:
    ''' Method for batch processing data from Gamry
    Uses methods defined in CV class

    Author: Bernard Kim
    '''

    def __init__(self, alldata, reduce=False):
        # Accepts lists of class CV
        self.allcycles = {}

        for file in alldata:
            exported = CV(file, reduce=reduce)

            match = re.search(r'S\d{1,2}', file)
            sampleidx = int(match.group(0)[1:])
            self.allcycles[sampleidx] = exported.cycles

        titlesearch = re.search(alldata[0])

    def plot_current_voltage(self, cycle_index=0, xlim=None, ylim=None,
        title=None, show=False, save=False):
        ''' Plots current vs voltage by cycle with all samples on one graph '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0.4,1,len(self.allcycles)) # for use with Blues colormap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for idx, sample in enumerated(sorted(self.allcycles)):
            if cycle_index:
                ax.plot(self.allcycles[sample][cycle_index]['voltage'],
                        self.allcycles[sample][cycle_index]['current'],
                        linewidth=2,
                        color=plt.cm.Blues(coloridx[int(idx)-1]),
                        label='S'+str(sample))

        ax.legend()
        ax.set_xlabel('Potential [V]')
        ax.set_ylabel('Current [mA/cm' +r'$^2$' +']')
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title()


        figtitle = title + '_C' + str(cycle_index) + '_n=' + str(len(self.allcycles))
        plt.title(figtitle)

        if show:
            plt.show()

        if save:
            plt.savefig('export_' + figtitle + '.pdf', format='pdf')

        plt.close()

    def plot_extrema(files, title=None, samplen=None, averaged=False, errorbar=None, save=False):
        ''' Plots exported extrema data into two plots of average max and min 
            voltage and currents with stdev error bars '''

        # instantiate lists for relevant data
        # these need to be lists of lists, where sublists are for different cycles?
        cycles = list(range(1,7)) # hardcoded as 6 cycles for now, possibly make variable later
        V_red, I_red = [[] for cycle in cycles], [[] for cycle in cycles]
        V_ox, I_ox = [[] for cycle in cycles], [[] for cycle in cycles]

        V_red_avg, V_red_std = [], []
        I_red_avg, I_red_std = [], []
        V_ox_avg, V_ox_std = [], []
        I_ox_avg, I_ox_std = [], []

        # construct dictionary 'data' compiling data from csv files
        for file in files:
            with open(file) as f:
                reader = csv.DictReader(f)
                for idx,row in enumerate(reader):
                    V_ox[idx].append(float(row['E_p_ox (V)']))
                    I_ox[idx].append(float(row['I_p_ox (mA)']))
                    V_red[idx].append(float(row['E_p_red (V)']))
                    I_red[idx].append(float(row['I_p_red (mA)']))

        V_ox_arr = np.array(V_ox)
        V_ox_arr[V_ox_arr == 0] = 'nan'
        V_ox_rows, V_ox_cols = V_ox_arr.shape

        I_ox_arr = np.array(I_ox)
        I_ox_arr[I_ox_arr == 0] = 'nan'
        I_ox_rows, I_ox_cols = I_ox_arr.shape

        V_red_arr = np.array(V_red)
        V_red_arr[V_red_arr == 0] = 'nan'
        V_red_rows, V_red_cols = V_red_arr.shape

        I_red_arr = np.array(I_red)
        I_red_arr[I_red_arr == 0] = 'nan'
        I_red_rows, I_red_cols = I_red_arr.shape

        markers = ['.','o','v','^','s','p','d','h','*','x','8','D']

        ## Begin plotting section
        font = {'family': 'Arial', 'size': 13}
        matplotlib.rc('font', **font)

        fig, (ax_V_ox, ax_V_red) = plt.subplots(2, 1, sharex=True)

        ax_V_ox.set_title('Anodic Scan')
        ax_I_ox = ax_V_ox.twinx()
        ax_V_ox.set_ylabel('Voltage (V)', color='b')
        ax_V_ox.tick_params('y', colors='b')
        ax_I_ox.set_ylabel('Current (mA)', color='r')
        ax_I_ox.tick_params('y', colors='r')

        ax_V_red.set_title('Cathodic Scan')
        ax_I_red = ax_V_red.twinx()
        ax_V_red.set_xlabel('Cycle')
        ax_V_red.set_ylabel('Voltage (V)', color='b')
        ax_V_red.tick_params('y', colors='b')
        ax_I_red.set_ylabel('Current (mA)', color='r')
        ax_I_red.tick_params('y', colors='r')

        if averaged: # For averaged batch processing of all samples
            for an_V,an_I,ca_V,ca_I in zip(V_ox_arr,I_ox_arr,V_red_arr,I_red_arr):
                V_ox_avg.append(np.nanmean(an_V))
                V_ox_std.append(np.nanstd(an_V, ddof=1))
                I_ox_avg.append(np.nanmean(an_I))
                I_ox_std.append(np.nanstd(an_I, ddof=1))
                V_red_avg.append(np.nanmean(ca_V))
                V_red_std.append(np.nanstd(ca_V, ddof=1))
                I_red_avg.append(np.nanmean(ca_I))
                I_red_std.append(np.nanstd(ca_I, ddof=1))

            if errorbar: # With errorbars shown (y axes may be blown out of scale)
                ax_V_ox.errorbar(
                    cycles, V_ox_avg, yerr=V_ox_std, 
                    color='b', marker='.', markersize=10,
                    capsize=10, elinewidth=2)
                ax_I_ox.errorbar(
                    cycles, I_ox_avg, yerr=I_ox_std, 
                    color='r', marker='.', markersize=10,
                    capsize=10, elinewidth=2)
                ax_V_red.errorbar(
                    cycles, V_red_avg, yerr=V_red_std, 
                    color='b', marker='.', markersize=10,
                    capsize=10, elinewidth=2)
                ax_I_red.errorbar(
                    cycles, I_red_avg, yerr=I_red_std, 
                    color='r', marker='.', markersize=10,
                    capsize=10, elinewidth=2)
                figtitle = 'batch_avg_eb_' + title + '_n='+ str(samplen)

            else: # Without errorbars shown
                ax_V_ox.plot(cycles, V_ox_avg, color='b')
                ax_V_ox.set_ylim([0, 1.5])
                ax_I_ox.plot(cycles, I_ox_avg, color='r')
                ax_I_ox.set_ylim([0, 20])
                ax_V_red.plot(cycles, V_red_avg, color='b')
                ax_V_red.set_ylim([-1.5, 0])
                ax_I_red.plot(cycles, I_red_avg, color='r')
                ax_I_red.set_ylim([-20, 0])
                figtitle = 'batch_avg_' + title + '_n='+ str(samplen)

        else: # For all samples on the same plot
            for V_ox_col in range(V_ox_cols):
                ax_V_ox.plot(cycles, V_ox_arr[:,V_ox_col], color='b',marker=markers[V_ox_col])
            for I_ox_col in range(I_ox_cols):
                ax_I_ox.plot(cycles, I_ox_arr[:,I_ox_col], color='r',marker=markers[I_ox_col])
            for V_red_col in range(V_red_cols):
                ax_V_red.plot(cycles, V_red_arr[:,V_red_col], color='b',marker=markers[V_red_col])
            for I_red_col in range(I_red_cols):
                ax_I_red.plot(cycles, I_red_arr[:,I_red_col], color='r',marker=markers[I_red_col])
            figtitle = 'batch_' + title + '_n='+ str(samplen)

        plt.suptitle(figtitle)

        if save:
            plt.savefig(figtitle + '.eps', format='eps')

        plt.show()