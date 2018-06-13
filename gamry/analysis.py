''' Module for analyzing results retrieved from Gamry

Author: Bernard Kim
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley

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

    def __init__(self, filename=None, thickness=0.0365, area=1, zscale='k'):
        ''' Opens file and retrieves data.

        Retrieves time, frequency, real impedance, imaginary impedance,
        magnitude, and phase. Assumes that the first part of the file is an
        OCV test and that the header for the table consists of two lines.

        Unit requirements:
            R_solution [ohm]
            Thickness [cm]
            Area [cm^2]
        '''

        self.filename = filename
        self.thickness = thickness
        self.area = area

        self.zscalestr, self.zscaleval = self.get_zscale(zscale)

        titlesearch = re.search(r'EXDTA_.*_S\d{1,2}', self.filename)
        self.title = titlesearch.group(0)[6:]

        self.time = []
        self.freq = []
        self.realraw = []
        self.imagraw = []
        self.phase = []
        self.magnraw = []

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
                            self.realraw.append(float(row[3]))
                            self.imagraw.append(float(row[4]))
                            self.magnraw.append(float(row[6]))
                            self.phase.append(float(row[7]))
                except Exception:
                    raise

        self.real = [realraw/self.zscaleval for realraw in self.realraw]
        self.imag = [-1*imagraw/self.zscaleval for imagraw in self.imagraw]
        self.magn = [magnraw/self.zscaleval for magnraw in self.magnraw]

        self.find_r_solution()
        self.conductivity = (self.thickness * 1000) / (self.area * self.r_solution)

    def is_num(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def get_zscale(self, zscale):
        ''' Determines scaling factor for impedance values
            Default value is in Ohms, scales values to kOhms or MOhms '''

        zscaledict = {
            None: {
                'zscalestr': 'Ω', 
                'zscaleval': 1e0,
            },
            'k': {
                'zscalestr': 'kΩ', 
                'zscaleval': 1e3,
            },
            'M': {
                'zscalestr': 'MΩ', 
                'zscaleval': 1e6,
            },
        }

        return zscaledict[zscale]['zscalestr'], zscaledict[zscale]['zscaleval']

    def find_r_solution(self):
        ''' Calculated solution resistance by taking value of real impedance
            closest to real axis intercept '''

        min_imag_val = min(abs(imag) for imag in self.imagraw[:10])
        min_imag_idx = [abs(imag) for imag in self.imagraw].index(min_imag_val)

        self.r_solution = self.realraw[min_imag_idx]

    def plot_nyquist(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, imagetype='png', savename=None):
        ''' Plots imaginary vs. real impedance 
            Imaginary and real axes locked to same scale by default '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(12,9), dpi=75)

        ax.plot(self.real, self.imag, color='b', linewidth=2,
            label=r'$R_{solution}$' + ' = ' + '%.2f'%self.r_solution + ' Ω' + '\n' + 
                r'$\sigma$' + ' = ' + '%0.2f'%self.conductivity + ' S')

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        ax.set_aspect('equal', 'datalim')
        ax.set_xlabel(r'$Z_{real}$' + ' [' + self.zscalestr + ']')
        ax.set_ylabel(r'$-Z_{imag}$' + ' [' + self.zscalestr + ']')
        ax.legend()
        ax.grid()

        if show:
            plt.show()

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig(self.title + '_nyquist' + '.' + imagetype)

        plt.close(fig)

    def plot_bode(self, xlim=None, ylim=None, title=None, 
        show=False, save=False, imagetype='png'):
        ''' Plots imaginary vs. real impedance 
            Imaginary and real axes locked to same scale by default '''

        font = {'family': 'Arial', 'size': 16}
        matplotlib.rc('font', **font)

        fig, (ax_magn, ax_phase) = plt.subplots(2, sharex=True, figsize=(16,9), dpi=75)

        ax_magn.loglog(self.freq, self.magn,
            color='b', linewidth=2)
        ax_phase.semilogx(self.freq, self.phase,
            color='b', linewidth=2)

        if xlim:
            ax_phase.set_xlim(xlim)
        if ylim:
            ax_phase.set_ylim(ylim)

        if title:
            ax_magn.set_title(title)
        else:
            ax_magn.set_title(self.title)

        ax_magn.set_ylabel('Magnitude' + ' [' + self.zscalestr + ']')
        ax_phase.set_ylabel('Phase [°]')
        ax_phase.set_xlabel('Frequency [Hz]')
        ax_magn.grid()
        ax_phase.grid()

        if show:
            plt.show()

        if save:
            plt.savefig(self.title + '_bode' + '.' + str(imagetype))

        plt.close(fig)

class CV:

    ''' Analyzes data from Gamry cyclic voltammetry data

    Pulls data from .dta file for plotting and analysis.

    '''

    def __init__(self, filename=None):
        ''' Opens file and retrieves data. '''

    # Pt  T     Vf    Im  Vu  Sig Ach IERange Over
    # s   V vs. Ref.  A   V   V   V   #       bits

        self.cycles = {}
        cyclenumbers = []

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
                        if row[0][0:4] == 'AREA':
                            self.area = float(row[2])
                        elif row[0][0:5] == 'CURVE' and len(row[0]) > 5:
                            curve_number = int(row[0][5::])
                            cyclenumbers.append(curve_number)
                            switch = index + 2

                            # Save previous cycle's information
                            if current_cycle_time:
                                self.cycles[curve_number-1] = {
                                    'time': current_cycle_time,
                                    'voltage': current_cycle_voltage,
                                    'current': current_cycle_current,
                                }

                                # Erase previous cycle's temporary data
                                current_cycle_time = []
                                current_cycle_voltage = []
                                current_cycle_current = []

                        # Save current cycle's data into temporary lists
                        if (self.is_num(row[0]) and switch and index > switch):
                            current_cycle_time.append(float(row[1]))
                            current_cycle_voltage.append(float(row[2]))
                            current_cycle_current.append(float(row[3]) 
                                * 1000/self.area)

                except Exception:
                    raise

        self.title = filename[:-4]

    def is_num(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False


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
        title=None, show=False, save=False, savename=None):
        ''' Plots current vs voltage for one or all cycles '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        # coloridx = np.linspace(0.4,1,len(self.cycles)) # for use with Blues cmap
        coloridx = np.linspace(0,0.9,len(self.cycles)) # for use with viridis cmap

        fig, ax = plt.subplots(figsize=(16,10), dpi=75)

        if cycle_index:
            for cycle in cycle_index:
                ax.plot(self.cycles[cycle]['voltage'],
                        self.cycles[cycle]['current'],
                        linewidth=3,
                        color=plt.cm.inferno(coloridx[cycle-1]),
                        label = 'Cycle '+str(cycle))
            ax.legend()
        else:
            for cycle in self.cycles:
                ax.plot(self.cycles[cycle]['voltage'],
                        self.cycles[cycle]['current'],
                        linewidth=3,
                        color=plt.cm.inferno(coloridx[cycle-1]),
                        label='Cycle '+str(cycle))
            ax.legend()

        ax.set_xlabel('Potential [V]')
        ax.set_ylabel('Current [mA/cm$^2$]')
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_C' + str(cycle_index))

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig('single_' + self.title + '_C' + str(cycle_index) + '.png')

        if show:
            plt.show()

        plt.close(fig)


    def extrema_smoothed_wd(self, save_csv=False, showplot=False, saveplot=False):
        ''' Extracts maximum and minumum points on each CV curve for each cycle
            Uses gradient-window (g-w) method on smoothed curves
            Written by Bernard Kim '''

        del self.cycles[max(self.cycles)] # removes last element of self.cycles (incomplete cycles)
        # !!!!! MAY WANT TO CHANGE THIS LATER FOR LINEAR EXTRAPOLATION PURPOSES !!!!!

        voltagethresh = 2 # positive and negative hard voltage bounds, no effect when == 2
        smoothsize = 0.05 # size of smoothing window, as fraction of individual scan length
        lbsize = 0.10 # window lower bound size (for g-w), as percent of scan length
        ubsize = 0.25 # window upper bound size (for g-w), as percent of scan length

        allPairs = {}

        for cycle in self.cycles:
            voltages = self.cycles[cycle]['voltage']
            currents = self.cycles[cycle]['current']
            times = self.cycles[cycle]['time']

            # Split voltage, current, and time lists into cathodic and anodic scans
            scan_switch = voltages.index(min(voltages))

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
                        windowExtreme = self.find_min_max(extremaPairs, 'current', 'min')
                    elif scan is cathodic:
                        windowExtreme = self.find_min_max(extremaPairs, 'current', 'max')

                    # Populate master list of extreme currents for all window sizes
                    totalExtrema[window - windows[0]] = windowExtreme

                # Find most optimal extrema for all window and gradient possibilities and
                # populate master list of extrema per cycle
                if scan is anodic:
                    realExtreme = self.find_min_max(totalExtrema, 'current', 'min')
                    allPairs[cycle]['anodic'] = realExtreme
                elif scan is cathodic:
                    realExtreme = self.find_min_max(totalExtrema, 'current', 'max')
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
                ax.set_title('wd_' + self.title + '_C' + str(cycle))
                ax.grid()

                plt.show()
                # plt.close()

                if saveplot:
                    plt.savefig('wd_' + self.title + '_C' + str(cycle) + '.png')

        if save_csv:
            # write data to csv file
            filename = str(self.title)

            with open('extrema_' + filename + '.csv', 'w', newline='\n') as f:
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


    def get_charges(self):
        ''' Determines total amount of charge passed for both cathodic and 
            anodic scans per cycle. Integrates current wrt time using midpoint 
            (rectangular) rule. '''

        # Instantiate dictionary to export
        allcharge = {
            cycle: {
                'oxidation': 0,
                'reduction': 0,
            } for cycle in self.cycles
        }

        for cycle in self.cycles:
            times = self.cycles[cycle]['time']
            currents = self.cycles[cycle]['current']

            # Get time step value, dt
            t_start = min(times)
            t_end = max(times)
            dt = (t_end-t_start)/(len(times)-1)

            oxidation = 0 # positive currents
            reduction = 0 # negative currents

            for current in currents:
                charge = current*dt * (1/3600) # [mAh], mA*s*(1h/3600 s)

                if np.sign(current) == 1:
                    oxidation += np.abs(charge)
                elif np.sign(current) == -1:
                    reduction += np.abs(charge)

            allcharge[cycle]['oxidation'] = oxidation
            allcharge[cycle]['reduction'] = reduction

        self.charges = allcharge


    def plot_charge(self, show=False, save=False, title=None, savename=None):
        ''' Plots total oxidation and reduction charge per cycle of CV '''

        self.get_charges()

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
            figsize=(16,10), dpi=75)

        plotcycles = [cycle for cycle in self.charges]
        oxidation = [self.charges[cycle]['oxidation'] for cycle in self.charges]
        reduction = [self.charges[cycle]['reduction'] for cycle in self.charges]

        efficiency = []
        for ox, red in zip(oxidation, reduction):
            efficiency.append(ox/red * 100)

        ax1.plot(plotcycles, oxidation, 
            linewidth=3, marker='o', markersize=10, color='r', 
            label='Oxidation')
        ax1.plot(plotcycles, reduction, 
            linewidth=3, marker='o', markersize=10, color='b', 
            label='Reduction')

        ax2.plot(plotcycles, efficiency,
            linewidth=3, marker='o', markersize=10, color='k',
            label='oxidation/reduction')

        ax1.set_ylabel('Charge [mAh/cm'+r'$^2$'+']')
        ax2.set_ylabel('Cycle Efficiency [%]')
        ax2.set_xlabel('Cycle Number')

        ax2.set_xlim(xmin=1, xmax=max(plotcycles))

        ax1.legend()
        ax2.legend()
        ax1.grid()
        ax2.grid()

        if title:
            fig.suptitle(title)
        else:
            fig.suptitle(self.title + '_charge')

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title+'_charge'+'.png')

        if show:
            plt.show()

        plt.close(fig)




    def reduce_file_LEGACY(self):
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

    def undersample_LEGACY(self, samplerate, curve_number):
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

class CV_batch:
    ''' Method for batch processing data from Gamry
    Uses methods defined in CV class

    Author: Bernard Kim
    '''

    def __init__(self, alldata, var):
        # Accepts list of filenames
        self.allcycles = {}
        titles = []

        for file in alldata:
            exported = CV(file)
            title, sampleidx = self.title_search_format(filename=file, var=var)
            self.allcycles[sampleidx] = exported.cycles
            titles.append(title)

        # If all titles are identical, great
        # Otherwise, throw soft error
        self.title = titles[0]

        if len(set(titles)) is not 1:
            print('Titles do not match!!')

    def title_search_format(self, filename, var=None):
        ''' 
        Searches filenames depending on experiment data type for batch 
        processing
        Returns formatted title with constant sample parameters for 
        batch plotting and sample idx
        var can be "IL", "mol", "nu", or "N"
        '''

        # Searches filename for GPE vs. IL sample type, IL identity, 
        # concentration, scan rate, and sample number
        isGPEmatch = re.search(r'GPE', filename)
        isILmatch = re.search(r'IL', filename)
        ILmatch = re.search(r'[A-Z]{1}MIM[a-zA-Z\d]{3,4}', filename)
        molmatch = re.search(r'0,\d{1}m', filename)
        numatch = re.search(r'\d{2,3}mVs', filename)
        Nmatch = re.search(r'v\d{1,2}', filename)

        # Determines order of elements in title
        titlematches = [isGPEmatch, isILmatch, ILmatch, molmatch, numatch]

        # Instantiate title, samplenum, and samplevar as placeholders
        title = ''
        samplenum, samplevar = None, None

        # Check if sample number exists and assign if it does
        if Nmatch:
            samplenum = Nmatch.group(0)

        # Determine sampleidx whether or not samplenum exists
        if var == 'N':
            sampleidx = samplenum
        else:
            if var == 'IL':
                samplevar = ILmatch.group(0)
                if Nmatch:
                    sampleidx = samplevar + '_' + samplenum
                else:
                    sampleidx = samplevar
            elif var == 'mol':
                samplevar = molmatch.group(0)
                if Nmatch:
                    sampleidx = samplevar + '_' + samplenum
                else:
                    sampleidx = samplevar
            elif var == 'nu':
                samplevar = numatch.group(0)
                if Nmatch:
                    sampleidx = samplevar[:-3] + ' ' + 'mV/s ' + samplenum
                else:
                    sampleidx = samplevar[:-3] + ' ' + 'mV/s'
            elif var is None:
                samplevar = None
                sampleidx = filename

        # Extract samplevar from title and don't add it back into the title
        # Otherwise, reconstruct the title in the right order, regardless of 
        # filename order
        for match in titlematches:
            if match is not None:
                if match.group(0) == samplevar:
                    continue
                else:
                    title = title + match.group(0) + '_'

        title = title[:-1]

        return title, sampleidx

    def plot_current_voltage(self, cycle_index=0, xlim=None, ylim=None,
        title=None, show=False, save=False, savename=None, imagetype='png'):
        ''' Plots current vs voltage by cycle with all samples on one graph '''

        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)

        # coloridx = np.linspace(0.5,1,len(self.allcycles)) 
        # for use with Blues colormap
        coloridx = np.linspace(0,1,10) # for use with tab10 colormap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for idx, sample in enumerate(self.allcycles):
            if cycle_index:
                ax.plot(self.allcycles[sample][cycle_index]['voltage'],
                        self.allcycles[sample][cycle_index]['current'],
                        linewidth=3,
                        color=plt.cm.tab10(coloridx[int(idx)]),
                        label=str(sample))

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
            ax.set_title('batch_' + self.title + '_C' + str(cycle_index))

        if save:
            if savename:
                fig.savefig(savename + '.' + imagetype)
            else:
                fig.savefig('batch_' + self.title + '_C' + str(cycle_index) + 
                    '.' + str(imagetype))

        if show:
            plt.show()

        plt.close(fig)

class CV_extrema_plot:
    ''' Method for plotting exported extrema data 
        Plots individual peak coordinates on same I-V curve '''

    def __init__(self, files):
        ''' Process raw csv files into useable data structures
            Assumes csv files have been generated from CV class '''

        self.alldata = {}

        titlesearch = re.search(r'(GPE|IL).*_S', files[0])
        self.title = titlesearch.group(0)[:-2]

        # scansearch = re.search(r'\d{1,3}mVs',files[0])
        # self.scan = scansearch.group(0)

        # Based off headers from generated csv files from CV class
        # headers = ['Cycle', 'Epa [V]', 'Ipa [mA]', 'Epc [V]', 'Ipc [A]']
        headers = ['Cycle', 'E_p_ox (V)', 'I_p_ox (mA)', 'E_p_red (V)', 'I_p_red (mA)']

        for file in files:
            # Get sample number
            samplesearch = re.search(r'_S\d{1,2}', file)
            sample = samplesearch.group(0)[2:]

            rows = list(csv.reader(open(file), delimiter=','))
            idxs = []

            for header in headers:
                idxs.append(rows[0].index(header))

            # Robustness to be column order-agnostic within csv
            cycleidx = idxs[0]
            Epaidx, Ipaidx = idxs[1], idxs[2]
            Epcidx, Ipcidx = idxs[3], idxs[4]

            # # Ensuring cycle numbers are sorted
            # cyclenumbers = [int(row[cycleidx]) for row in rows[1:]]
            # cyclenumbers.sort()

            cycles, Epa, Ipa, Epc, Ipc = [], [], [], [], []

            ### Iterate over rows to get each column as separate list
            for row in rows[1:]:
                cycles.append(row[cycleidx])
                Epa.append(row[Epaidx])
                Ipa.append(row[Ipaidx])
                Epc.append(row[Epcidx])
                Ipc.append(row[Ipcidx])

            self.alldata[int(sample)] = {
                'cycles': cycles,
                'Epa': Epa,
                'Ipa': Ipa,
                'Epc': Epc,
                'Ipc': Ipc,
            }

    def plot_extrema(self, xlim=None, ylim=None, title=None,
        show=False, save=False, savename=None, imagetype='png'):
        ''' Plots extrema points from CV on i-v plot
            Samples are plotted as lines for each cycle '''

        colors = [plt.cm.viridis(idx) for idx in np.linspace(0.1,1,8)]
        # colors = [plt.cm.plasma(idx) for idx in np.linspace(0,1,len(self.alldata))]

        font = {'family': 'Arial', 'size': 28}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        marker = '.'
        markersize = 25
        markeredgecolor = 'k'
        markercolorfirst = 'w'
        markercolorlast = 'k'

        for sample in sorted(list(self.alldata.keys())):
            ax.plot(self.alldata[sample]['Epa'], self.alldata[sample]['Ipa'],
                color=colors[sample-1], linewidth=4, label='S'+str(sample),)
            ax.plot(self.alldata[sample]['Epc'], self.alldata[sample]['Ipc'],
                color=colors[sample-1], linewidth=4,)

            # ax.plot(self.alldata[sample]['Epa'], self.alldata[sample]['Ipa'],
            #     color=colors[coloridx], linewidth=4, label='S'+str(sample),)
            # ax.plot(self.alldata[sample]['Epc'], self.alldata[sample]['Ipc'],
            #     color=colors[coloridx], linewidth=4,)

            ax.plot([self.alldata[sample]['Epa'][0]], [self.alldata[sample]['Ipa'][0]],
                color=markercolorfirst, 
                marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)
            ax.plot([self.alldata[sample]['Epa'][-1]], [self.alldata[sample]['Ipa'][-1]],
                color=markercolorlast, 
                marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)
            ax.plot([self.alldata[sample]['Epc'][0]], [self.alldata[sample]['Ipc'][0]],
                color=markercolorfirst,
                marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)
            ax.plot([self.alldata[sample]['Epc'][-1]], [self.alldata[sample]['Ipc'][-1]],
                color=markercolorlast, 
                marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)

        ax.plot([4], [0], color=markercolorfirst, linewidth=0, 
                marker=marker, markersize=markersize, markeredgecolor=markeredgecolor,
                label='First Cycle')
        ax.plot([4], [0], color=markercolorlast, linewidth=0, 
                marker=marker, markersize=markersize, markeredgecolor=markeredgecolor,
                label='Last Cycle')

        ax.set_xlabel('Potential [V]')
        ax.set_ylabel('Current [mA]')
        ax.legend(fontsize=24)
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim([-2,2])

        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title)

        if show:
            plt.show()

        if save:
            if savename:
                plt.savefig(savename + '.' + imagetype)
            else:
                plt.savefig('extrema_batch_' + self.title + '.' + imagetype)

        plt.close(fig)


