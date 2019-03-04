''' Module for analyzing results retrieved from Gamry

Author: Bernard Kim
Principal Investigators: Prof. Paul Wright, Prof. James Evans
University: University of California, Berkeley

'''

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
import re
import bisect
from decimal import Decimal
from battery.utilities import utilities

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
        self.conductivity = (self.thickness*1000)/(self.area*self.r_solution)

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
            label=r'$R_{solution}$' + ' = ' + '%.2f'%self.r_solution + ' Ω' + 
            '\n' + r'$\sigma$' + ' = ' + '%0.2f'%self.conductivity + ' S')

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

        fig, (ax_magn, ax_phase) = plt.subplots(
            2, sharex=True, figsize=(16,9), dpi=75)

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
            cycle_time = []
            cycle_voltage = []
            cycle_current = []

            for index, row in enumerate(rows):
                row = row.split()
                try:
                    if row:
                        if row[0][0:4] == 'AREA':
                            self.area = float(row[2])
                        elif row[0][0:3] == 'EOC':
                            self.Eoc = float(row[2])
                        elif row[0][0:5] == 'CURVE' and len(row[0]) > 5:
                            curve_number = int(row[0][5::])
                            cyclenumbers.append(curve_number)
                            switch = index + 2

                            # Save previous cycle's information
                            if cycle_time:
                                self.cycles[curve_number-1] = {
                                    'time': cycle_time,
                                    'voltage': cycle_voltage,
                                    'current': cycle_current,
                                }

                                # Erase previous cycle's temporary data
                                cycle_time = []
                                cycle_voltage = []
                                cycle_current = []

                        # Save current cycle's data into temporary lists
                        if (self.is_num(row[0]) and switch and index > switch):
                            cycle_time.append(float(row[1]))
                            cycle_voltage.append(float(row[2]))
                            cycle_current.append(float(row[3]) 
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
            minimum value based on specified dictionary key 
        '''

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


    def plot_current_voltage(self, cycle_index=0, 
        xlim=None, ylim=None, corrected=False,
        title=None, show=False, save=False, savename=None):
        ''' Plots current vs voltage for one or all cycles '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        # coloridx = np.linspace(0.4,1,len(self.cycles)) # use with Blues cmap
        coloridx = np.linspace(0,0.9,len(self.cycles)) # use with viridis cmap

        fig, ax = plt.subplots(figsize=(16,10), dpi=75)

        if corrected:
            # Get reference shift if not calculated already
            try:
                if self.E_corrected:
                    pass
            except AttributeError:
                self.calculate_reference_shift()

            voltage = self.E_corrected
            correctedstr = '_corrected'
        else:
            voltage = {cycle: self.cycles[cycle]['voltage'] 
                for cycle in self.cycles}
            correctedstr = ''

        if cycle_index:
            for cycle in cycle_index:
                ax.plot(voltage[cycle],
                        self.cycles[cycle]['current'],
                        linewidth=3,
                        color=plt.cm.inferno(coloridx[cycle-1]),
                        label = 'Cycle '+str(cycle))
            ax.legend()
        else:
            for cycle in self.cycles:
                if cycle == min(self.cycles) or cycle == max(self.cycles) \
                    or cycle%50 == 0:
                    ax.plot(voltage[cycle],
                            self.cycles[cycle]['current'],
                            linewidth=3,
                            color=plt.cm.inferno(coloridx[cycle-1]),
                            label='Cycle '+str(cycle))
                elif cycle%25 == 0:
                    ax.plot(voltage[cycle],
                        self.cycles[cycle]['current'],
                        linewidth=3,
                        color=plt.cm.inferno(coloridx[cycle-1]),
                        )
            ax.legend()

        ax.set_xlabel('Potential, vs Zn/Zn$^{2+}$ [V]')
        ax.set_ylabel('Current [mA/cm$^2$]')
        ax.grid()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_C' + str(cycle_index) + correctedstr)

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title + '_C' + str(cycle_index) + \
                    correctedstr + '.png')

        if show:
            plt.show()

        plt.close(fig)


    def extrema_smoothed_wd(self, save_json=False, showplot=False, 
        saveplot=False):
        ''' Extracts maximum and minumum points on each CV curve for each cycle
            Uses gradient-window (g-w) method on smoothed curves
            Written by Bernard Kim '''

        # # removes last element of self.cycles (incomplete cycles)
        # del self.cycles[max(self.cycles)]

        voltagethresh = 2 # positive and negative hard voltage bounds
        # as percent of scan length
        smoothsize = 0.05 # size of smoothing window
        lbsize = 0.10 # window lower bound size (for w-g)
        ubsize = 0.25 # window upper bound size (for w-g)

        allPairs = {}

        for cycle in self.cycles:
            voltages = self.cycles[cycle]['voltage']
            currents = self.cycles[cycle]['current']
            times = self.cycles[cycle]['time']

            # Split voltage, current, and time lists into 
            # cathodic and anodic scans
            scan_switch = voltages.index(min(voltages))

            cathodic = {
                'V_raw': np.array(voltages[:scan_switch]),
                'I_raw': np.array(currents[:scan_switch]),
            }

            anodic = {
                'V_raw': np.array(voltages[scan_switch:]),
                'I_raw': np.array(currents[scan_switch:]),
            }

            allPairs[cycle] = {}

            # Smooth and run gradient-window algorithm on each scan separately
            for scan in [cathodic, anodic]:

                # Smooth current data using moving average
                smoothing_window = round(smoothsize*len(scan['V_raw']))
                scan['I_smoothed'] = self.moving_average(
                    scan['I_raw'], smoothing_window)

                # Extract corresponding smoothed voltage points to 
                # smoothed current points
                overhang = len(scan['I_raw']) - len(scan['I_smoothed'])
                lower_pad = np.ceil(overhang/2)
                upper_pad = np.floor(overhang/2)
                scan['V_smoothed'] = scan['V_raw'][
                    int(lower_pad):int(-upper_pad)]

                # Set window size as a function of smoothed scan data length
                lowerbound = round(lbsize*len(scan['V_smoothed']))
                upperbound = round(ubsize*len(scan['V_smoothed']))
                # 1/2 window size, symmetric positive and negative offsets
                # index of center of window
                windows = list(range(lowerbound,upperbound))

                totalExtrema = [[] for window in windows]

                for window in windows:
                    # Create list of valid positions based on window size
                    positions = list(
                        range(window,len(scan['V_smoothed'])-window))
                    points = [[] for position in positions]

                    for position in positions:
                        # Extract voltage, current at left, middle, and right 
                        # window positions
                        backvoltage, backcurrent = \
                            scan['V_smoothed'][position-window], \
                            scan['I_smoothed'][position-window],
                        frontvoltage, frontcurrent = \
                            scan['V_smoothed'][position+window], \
                            scan['I_smoothed'][position+window],
                        middlevoltage, middlecurrent = \
                            scan['V_smoothed'][position], \
                           scan['I_smoothed'][position],

                        backslope = (middlecurrent-backcurrent)/(
                            middlevoltage-backvoltage)
                        frontslope = (frontcurrent-middlecurrent)/(
                            frontvoltage-middlevoltage)

                        # Create dictionary per point with slope products, 
                        # middle points, and endpoint voltages
                        # ** Of SMOOTHED voltage and current, not raw
                        points[position-window] = {
                            'slopeprod': backslope*frontslope,
                            'pair': {
                                'voltage': middlevoltage, 
                                'current': middlecurrent,
                            },
                            'raw_pair':{
                                'voltage': scan['V_raw'][
                                    int(position+lower_pad)],
                                'current': scan['I_raw'][
                                    int(position+lower_pad)],
                            },
                            'endpoints': {
                                'back': backvoltage,
                                'front': frontvoltage,
                            },
                        }

                    # Fill with dummy values to prevent empty sequence errors
                    extremaPairs = [{'voltage':0, 'current':0}]

                    for point in points:
                        # Negative slope product means slopes have opposite 
                        # signs and surround local extrema
                        if point['slopeprod'] < 0 and abs(
                            point['pair']['voltage']) < voltagethresh:
                            # Local maxima/minima
                            extremaPairs.append(point['raw_pair'])

                    # Pull out point for maximum/minumum current for given 
                    #  window size
                    if scan is cathodic:
                        windowExtreme = self.find_min_max(
                            extremaPairs, 'current', 'min')
                    elif scan is anodic:
                        windowExtreme = self.find_min_max(
                            extremaPairs, 'current', 'max')

                    # Populate master list of extreme currents for all 
                    # window sizes
                    totalExtrema[window - windows[0]] = windowExtreme

                # Find most optimal extrema for all window and gradient 
                # possibilities and populate master list of extrema per cycle
                if scan is cathodic:
                    realExtreme = self.find_min_max(
                        totalExtrema, 'current', 'min')
                    allPairs[cycle]['cathodic'] = realExtreme
                elif scan is anodic:
                    realExtreme = self.find_min_max(
                        totalExtrema, 'current', 'max')

                    # Find linear fit of beginning linear portion of 
                    # anodic scan to determine adjusted peak current value
                    coefs = np.polyfit(
                        anodic['V_raw'][0:int(0.5*len(anodic['V_raw']))], 
                        anodic['I_raw'][0:int(0.5*len(anodic['I_raw']))], 
                        deg=1,)

                    I_unadjusted = realExtreme['current']
                    I_diff = coefs[1] + coefs[0]*realExtreme['voltage']
                    I_adjusted = I_unadjusted - I_diff

                    realExtreme['current'] = I_adjusted
                    allPairs[cycle]['anodic'] = realExtreme


            # Compute potential difference (E_p_ox - E_p_red) and 
            # current ratio (I_p_ox/I_p_red) and add to master list

            # allPairs[cycle]['']

            if showplot and cycle == 1:
                fig, ax = plt.subplots(figsize=(16,9), dpi=60)
                ax.plot(
                    cathodic['V_raw'], cathodic['I_raw'], 
                    color='#800000', linewidth=4, label='raw')
                ax.plot(
                    anodic['V_raw'], anodic['I_raw'], 
                    color='#800000', linewidth=4)
                # ax.plot(cathodic['V_smoothed'], cathodic['I_smoothed'], 
                #    color='#005DBB', linewidth=2, label='smoothed')
                # ax.plot(anodic['V_smoothed'], anodic['I_smoothed'], 
                #    color='#005DBB', linewidth=2)
                ax.plot(
                    anodic['V_raw'], coefs[1]+coefs[0]*anodic['V_raw'], 
                    color='k')

                ax.plot(
                    allPairs[cycle]['cathodic']['voltage'], 
                    allPairs[cycle]['cathodic']['current'],
                    color='#00D48D', marker="o", markersize=20, 
                    markeredgewidth=2, markerfacecolor='None')
                ax.plot(
                    allPairs[cycle]['anodic']['voltage'], 
                    allPairs[cycle]['anodic']['current'],
                    color='#00D48D', marker="o", markersize=20, 
                    markeredgewidth=2, markerfacecolor='None')
                ax.set_xlabel('Voltage [V]')
                ax.set_ylabel('Current [mA]')
                ax.legend()
                ax.set_title('wd_' + self.title + '_C' + str(cycle))
                ax.grid()

                # plt.show()
                # plt.close()

                if saveplot:
                    plt.savefig('wd_' +self.title+ '_C' + str(cycle) + '.png')

        if save_json:
            # write data to csv file
            filename = str(self.title)

            peaks = {}

            for cycle in allPairs:
                peaks[cycle] = {
                    'E_p_ox': allPairs[cycle]['anodic']['voltage'], 
                    'I_p_ox': allPairs[cycle]['anodic']['current'], 
                    'E_p_red': allPairs[cycle]['cathodic']['voltage'], 
                    'I_p_red': allPairs[cycle]['cathodic']['current'], 
                }

            utilities.save_json(data=peaks, filename=filename+'.json')


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
            dt = (t_end-t_start)/(len(times))

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

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, 
            figsize=(16, 12), dpi=75)

        plotcycles = [cycle for cycle in self.charges]
        oxidation = [self.charges[cycle]['oxidation'] for cycle in self.charges]
        reduction = [self.charges[cycle]['reduction'] for cycle in self.charges]

        # Initialize arrays b/c reduction_all and oxidation_all need initial
        # value
        efficiency = [oxidation[0]/reduction[0] * 100]
        reduction_all = [reduction[0]]
        oxidation_all = [oxidation[0]]

        for ox, red in zip(oxidation[1:], reduction[1:]):
            efficiency.append(ox/red * 100)
            reduction_all.append(reduction_all[-1] + red)
            oxidation_all.append(oxidation_all[-1] + ox)

        ax1.plot(plotcycles, reduction, 
            linewidth=3, marker='o', markersize=10, color='b', 
            label='Reduction')
        ax1.plot(plotcycles, oxidation, 
            linewidth=3, marker='o', markersize=10, color='r', 
            label='Oxidation')

        ax2.plot(plotcycles, reduction_all, 
            linewidth=3, marker='o', markersize=10, color='b', 
            label='Reduction')
        ax2.plot(plotcycles, oxidation_all, 
            linewidth=3, marker='o', markersize=10, color='r', 
            label='Oxidation')

        ax3.plot(plotcycles, efficiency,
            linewidth=3, marker='o', markersize=10, color='k',
            label='oxidation/reduction')

        # ax1.set_ylabel('Charge [mAh/cm'+r'$^2$'+']')
        ax1.set_ylabel('Charge \n [mAh]')
        ax2.set_ylabel('Cumulative Charge \n [mAh]')
        ax3.set_ylabel('Coulombic Efficiency \n [%]')
        ax3.set_xlabel('Cycle Number')

        ax2.set_xlim(xmin=1, xmax=max(plotcycles))

        # ax1.set_ylim([0, 0.055])
        ax3.set_ylim([0, 100])

        ax1.legend()
        ax2.legend()
        ax1.grid()
        ax2.grid()
        ax3.grid()

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


    def calculate_reference_shift(self, plot=False):
        ''' Calculates possible shift of reference electrode potential due to 
            unequal amounts of oxidation and reduction occuring at working 
            electrode.

            Unequal amounts of oxidation and reduction per cycle suggests that 
            the solution composition is shifting during cycling, which should
            then shift the reference potential of the reference electrode due 
            to Faraday's law.
        '''

        R = 8.314 # Universal gas constant, [J/mol*K]
        T = 25+273.15 # Temperature, [K]
        F = 96485 # Faraday's constant, [coulombs]
        MW_Zn = 65.38 # MW of Zn, [g/mol]
        MW_ZnOtf = 363.51 # MW of ZnOtf, [g/mol]
        n = 2 # electrons involved in redox rxn
        s = 1 # stoichiometric constant of Zn/Zn2+

        m_total_init = 1 # [g] initial mass, assume about 1 g

        rho_IL = {
            'EMIMOTF': 1.387, # [g/mL]
            'BMIMOTF': 1.292, # [g/mL]
        }

        MW_IL = {
            'EMIMOTF': 260.23, # [g/mol]
            'BMIMOTF': 288.29, # [g/mol]
        }

        # Get charges per cycle if not gathered already
        try:
            if self.charges:
                pass
        except AttributeError:
            self.get_charges()

        titlematches = utilities.title_search_echem(self.title, test='CV')
        molsearch = re.search(r'0,\d{1}', titlematches['mol'])

        # Get concentration and IL from filename
        molal = float('0.'+molsearch.group(0)[2]) # mol/kg
        IL = titlematches['IL']

        # mass of salt [g]
        m_ZnOtf_init = m_total_init/(1 + (1/(molal*0.001*MW_ZnOtf))) 
        m_Zn_init = MW_Zn/MW_ZnOtf * m_ZnOtf_init # [g]
        m_IL_init = m_total_init - m_ZnOtf_init # [g]

        m_Zn = {0: m_Zn_init}

        for cycle in self.charges:
            delta_Q = self.charges[cycle]['reduction'] - \
                self.charges[cycle]['oxidation'] # [mAh]

            # If delta_Q is positive, more red than ox, Zn2+ conc decreases
            # If delta_Q is negative, more ox than red, Zn2+ conc increases

            # (mol * g/mol * mAh)/(Coul / 3.6Coul/mAh)
            # Measure of how much Zn is going into solution
            # Positive value means Zn is entering solution
            # Negative value means Zn is leaving solution
            delta_m = -s*MW_Zn*delta_Q/(n*F/3.6) # [g]
            m_Zn[cycle] = m_Zn[cycle-1] + delta_m

        mol_Zn, molal_Zn = {}, {}

        for cycle in m_Zn:
            mol_Zn[cycle] = \
                (m_Zn[cycle]/MW_Zn)/((m_IL_init/rho_IL[IL])*0.001) # [mol/L]
            molal_Zn[cycle] = (m_Zn[cycle]/MW_Zn)/(m_IL_init*0.001)

        E_shift, E_corrected = {}, {}

        for cycle in molal_Zn:
            E_shift[cycle] = -(R*T)/(n*F)*np.log(molal_Zn[cycle])
            E_corrected[cycle] = self.Eoc + E_shift[cycle]

        self.mol_Zn = mol_Zn
        self.molal_Zn = molal_Zn
        self.E_shift = E_shift
        self.E_corrected = E_corrected

        if plot:
            font = {'family': 'Arial', 'size': 24}
            matplotlib.rc('font', **font)

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16,10), dpi=75)

            cycles = [cycle for cycle in self.E_corrected]
            molals = [self.molal_Zn[cycle] for cycle in self.molal_Zn]
            potentials = [self.E_corrected[cycle] for cycle in self.E_corrected]

            ax1.plot(
                cycles, molals,
                linewidth=3, color='k')
            ax2.plot(
                cycles, potentials,
                linewidth=3, color='k')

            ax2.set_xlabel('Cycle Number')
            ax1.set_ylabel('Zn$^{2+}$ Concentration [mol/kg]')
            ax2.set_ylabel('Potential [V]')
            ax1.grid()
            ax2.grid()

            fig.suptitle('Reference Electrode Potential Shift')
            plt.savefig(self.title+'mass_change.png')

            plt.show()
            plt.close(fig)


    def extrema_smoothed_radius_LEGACY(self, title=None):
        ''' Extracts maximum and minimum points on each CV curve for each cycle
            Uses gradient/window method 
            Written by Bernard Kim '''

        # removes last element of self.cycles (incomplete cycle)
        del self.cycles[7]

        allPairs = []

        for cyclenum,cycle in enumerate(self.cycles,1):
            voltages = self.cycles[cycle]['voltage']
            currents = self.cycles[cycle]['current']
            times = self.cycles[cycle]['time']

            # Split voltage, current, and time lists into 
            # cathodic and anodic scans
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
            # 10 percent of each scan length
            smoothing_window = round(0.10*len(times)/2) 

            for electrode in [cathodic, anodic]:
                electrode['I_smoothed'] = self.moving_average(
                    electrode['I_raw'], smoothing_window, 'exp')

                # Calculate curvature, kappa, at each point using first 
                # and second derivatives
                electrode['d_one'] = np.gradient(electrode['I_smoothed'])
                electrode['d_two'] = np.gradient(electrode['d_one'])
                electrode['kappa'] = (np.abs(electrode['d_two']))/(
                    (1 + electrode['d_one']**2)**(3/2))

                # Find index of maximum curvature to find points of 
                # local extrema for I_smoothed
                idx = np.argmax(electrode['kappa'])

                # Compensate for truncated vector to pull corresponding 
                # data point from raw data vectors
                overhang = len(electrode['I_raw'])-len(electrode['I_smoothed'])
                lower_pad = np.ceil(overhang/2)
                upper_pad = np.floor(overhang/2)
                electrode['V_smoothed'] = electrode['V_raw'][
                    int(lower_pad):int(-upper_pad)]

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
            fwrite.writerow(['Cycle', 'E_p_red (V)', 
                'I_p_red (mA)','E_p_ox (V)','I_p_ox (mA)'])

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

    def __init__(self, files):
        # Accepts list of filenames
        self.alldata = {}
        self.Eoc = {}

        params = {}

        for file in files:
            params_match = utilities.title_search_echem(filename=file, 
                test='CV')

            for key in params_match:
                if key not in params:
                    params[key] = []
                else:
                    params[key].append(params_match[key])

        for key in params:
            if len(set(params[key])) is not 1:
                samplekey = key

        del params[samplekey]
        titlekeys = params.keys()

        self.title = ''
        for key in titlekeys:
            self.title = self.title + str(params[key][0]) + '_'

        self.title = self.title[:-1]

        for file in files:
            exported = CV(file)
            params_match = utilities.title_search_echem(filename=file, 
                test='CV')

            self.alldata[params_match[samplekey]] = exported.cycles
            self.Eoc[params_match[samplekey]] = exported.Eoc


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


    def plot_current_voltage(self, cycle_index=0, normalized=False, 
        xlim=None, ylim=None, title=None, show=False, save=False, 
        savename=None, imagetype='png'):
        ''' Plots current vs voltage by cycle with all samples on one graph '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        # coloridx = np.linspace(0.5,1,len(self.allcycles)) 
        # for use with Blues colormap
        coloridx = np.linspace(0,1,10) # for use with tab10 colormap

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for idx, sample in enumerate(self.alldata):
            if cycle_index:
                if normalized:
                    potential = [E - self.Eoc[sample] \
                        for E in self.alldata[sample][cycle_index]['voltage']]
                else:
                    potential = self.alldata[sample][cycle_index]['voltage']
                ax.plot(
                    potential, self.alldata[sample][cycle_index]['current'],
                    linewidth=3, color=plt.cm.tab10(coloridx[int(idx)]),
                    label=str(sample))
            else:
                print('Can\'t let you do that, Star Fox!')
                print('Must specify cycle index!!')
                raise

        ax.legend()
        ax.set_xlabel('Potential, vs. Zn/Zn2+ [V]')
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
        '''

    def __init__(self, allfiles):
        ''' Process raw csv files into useable data structures
            Assumes csv files have been generated from CV class '''


        self.alldata = {}
        self.nus = None

        for filegroup in allfiles:
            nus = []

            for file in filegroup:
                params_match = utilities.title_search_echem(filename=file, 
                    test='CV')

                molsearch = re.search(r'0,\d{1}', params_match['mol'])
                mol = float('0.'+molsearch.group(0)[-1])

                if mol not in self.alldata:
                    self.alldata[mol] = {}

                nusearch = re.search(r'\d{2,3}', params_match['nu']) # mV/s
                nu = int(nusearch.group(0))
                nus.append(nu)

                self.alldata[mol][nu] = utilities.read_json(filename=file)

            try:
                if not self.nus:
                    self.nus = nus
                elif nus == self.nus:
                    pass
            except nus != self.nus:
                print('Scan rates don\'t match!!')
                raise 


    def plot_diffusivity(self, cycle_index=0, ILtype=None, title=None, 
        show=False, save=False, savename=None,):
        ''' Plots extrema points from CV on i-v plot
            Samples are plotted as lines for each cycle '''

        try:
            if cycle_index:
                pass
        except Exception:
            print('Must specify cycle number!!')
            raise

        font = {'family': 'Arial', 'size': 20}
        matplotlib.rc('font', **font)

        # coloridx = np.linspace(0.5,1,len(self.allcycles)) 
        # for use with Blues colormap
        coloridx = np.linspace(0,1,10) # for use with tab10 colormap
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,9), dpi=75)

        if ILtype is 'EMIMOTF':
            rho = 1.387 # g/mL
        elif ILtype is 'BMIMOTF':
            rho = 1.292 # g/mL

        I_ox, I_red = {}, {}

        for mol in self.alldata:
            I_ox[mol], I_red[mol] = [], []

            for nu in self.nus:
                I_ox[mol].append(
                    self.alldata[mol][nu][str(cycle_index)]['I_p_ox'])
                I_red[mol].append(
                    self.alldata[mol][nu][str(cycle_index)]['I_p_red'])

        nu_sqrt = [(nu*0.001)**0.5 for nu in self.nus] # V/s^0.5

        for idx, mol in enumerate(I_ox):
            coefs = np.polyfit(nu_sqrt, I_ox[mol], deg=1)
            diffusivity = (coefs[0]/(2.69e5 * 2**1.5 * mol/rho))**2

            ax1.plot(
                nu_sqrt, I_ox[mol], 
                label=str(mol)+' m, D_o = %.3e cm^2/s'%diffusivity,
                linewidth=3, color=plt.cm.tab10(coloridx[int(idx)]))

        for idx, mol in enumerate(I_red):
            coefs = np.polyfit(nu_sqrt, I_red[mol], deg=1)
            diffusivity = (coefs[0]/(2.69e5 * 2**1.5 * mol/rho))**2

            ax2.plot(
                nu_sqrt, I_red[mol], 
                label=str(mol)+' m, D_o = %.3e cm^2/s'%diffusivity,
                linewidth=3, color=plt.cm.tab10(coloridx[int(idx)]))

        ax1.legend()
        ax1.set_xlabel(r'$\nu^{1/2} [V/s]^{1/2}$')
        ax1.set_ylabel('Peak Oxidation Current [mA/cm' +r'$^2$' +']')
        ax1.grid()

        ax2.legend()
        ax2.set_xlabel(r'$\nu^{1/2} [V/s]^{1/2}$')
        ax2.set_ylabel('Peak Reduction Current [mA/cm' +r'$^2$' +']')
        ax2.grid()


        if title:
            fig.suptitle(title)
        else:
            fig.suptitle(ILtype + '_C' + str(cycle_index))

        if save:
            if savename:
                fig.savefig(savename + '.' + imagetype)
            else:
                fig.savefig(ILtype + '_C' + str(cycle_index) + '.png')

        if show:
            plt.show()

        plt.close(fig)





        # self.alldata = {}

        # titlesearch = re.search(r'(GPE|IL).*_S', files[0])
        # self.title = titlesearch.group(0)[:-2]

        # # scansearch = re.search(r'\d{1,3}mVs',files[0])
        # # self.scan = scansearch.group(0)

        # # Based off headers from generated csv files from CV class
        # # headers = ['Cycle', 'Epa [V]', 'Ipa [mA]', 'Epc [V]', 'Ipc [A]']
        # headers = ['Cycle', 'E_p_ox (V)', 'I_p_ox (mA)', 'E_p_red (V)', 'I_p_red (mA)']

        # for file in files:
        #     # Get sample number
        #     samplesearch = re.search(r'_S\d{1,2}', file)
        #     sample = samplesearch.group(0)[2:]

        #     rows = list(csv.reader(open(file), delimiter=','))
        #     idxs = []

        #     for header in headers:
        #         idxs.append(rows[0].index(header))

        #     # Robustness to be column order-agnostic within csv
        #     cycleidx = idxs[0]
        #     Epaidx, Ipaidx = idxs[1], idxs[2]
        #     Epcidx, Ipcidx = idxs[3], idxs[4]

        #     # # Ensuring cycle numbers are sorted
        #     # cyclenumbers = [int(row[cycleidx]) for row in rows[1:]]
        #     # cyclenumbers.sort()

        #     cycles, Epa, Ipa, Epc, Ipc = [], [], [], [], []

        #     ### Iterate over rows to get each column as separate list
        #     for row in rows[1:]:
        #         cycles.append(row[cycleidx])
        #         Epa.append(row[Epaidx])
        #         Ipa.append(row[Ipaidx])
        #         Epc.append(row[Epcidx])
        #         Ipc.append(row[Ipcidx])

        #     self.alldata[int(sample)] = {
        #         'cycles': cycles,
        #         'Epa': Epa,
        #         'Ipa': Ipa,
        #         'Epc': Epc,
        #         'Ipc': Ipc,
        #     }

    # def plot_extrema(self, xlim=None, ylim=None, title=None,
    #     show=False, save=False, savename=None, imagetype='png'):
    #     ''' Plots extrema points from CV on i-v plot
    #         Samples are plotted as lines for each cycle '''

    #     colors = [plt.cm.viridis(idx) for idx in np.linspace(0.1,1,8)]
    #     # colors = [plt.cm.plasma(idx) for idx in np.linspace(0,1,len(self.alldata))]

    #     font = {'family': 'Arial', 'size': 28}
    #     matplotlib.rc('font', **font)
    #     fig, ax = plt.subplots(figsize=(16,9), dpi=75)

    #     marker = '.'
    #     markersize = 25
    #     markeredgecolor = 'k'
    #     markercolorfirst = 'w'
    #     markercolorlast = 'k'

    #     for sample in sorted(list(self.alldata.keys())):
    #         ax.plot(self.alldata[sample]['Epa'], self.alldata[sample]['Ipa'],
    #             color=colors[sample-1], linewidth=4, label='S'+str(sample),)
    #         ax.plot(self.alldata[sample]['Epc'], self.alldata[sample]['Ipc'],
    #             color=colors[sample-1], linewidth=4,)

    #         # ax.plot(self.alldata[sample]['Epa'], self.alldata[sample]['Ipa'],
    #         #     color=colors[coloridx], linewidth=4, label='S'+str(sample),)
    #         # ax.plot(self.alldata[sample]['Epc'], self.alldata[sample]['Ipc'],
    #         #     color=colors[coloridx], linewidth=4,)

    #         ax.plot([self.alldata[sample]['Epa'][0]], [self.alldata[sample]['Ipa'][0]],
    #             color=markercolorfirst, 
    #             marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)
    #         ax.plot([self.alldata[sample]['Epa'][-1]], [self.alldata[sample]['Ipa'][-1]],
    #             color=markercolorlast, 
    #             marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)
    #         ax.plot([self.alldata[sample]['Epc'][0]], [self.alldata[sample]['Ipc'][0]],
    #             color=markercolorfirst,
    #             marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)
    #         ax.plot([self.alldata[sample]['Epc'][-1]], [self.alldata[sample]['Ipc'][-1]],
    #             color=markercolorlast, 
    #             marker=marker, markersize=markersize, markeredgecolor=markeredgecolor)

    #     ax.plot([4], [0], color=markercolorfirst, linewidth=0, 
    #             marker=marker, markersize=markersize, markeredgecolor=markeredgecolor,
    #             label='First Cycle')
    #     ax.plot([4], [0], color=markercolorlast, linewidth=0, 
    #             marker=marker, markersize=markersize, markeredgecolor=markeredgecolor,
    #             label='Last Cycle')

    #     ax.set_xlabel('Potential [V]')
    #     ax.set_ylabel('Current [mA]')
    #     ax.legend(fontsize=24)
    #     ax.grid()

    #     if xlim:
    #         ax.set_xlim(xlim)
    #     else:
    #         ax.set_xlim([-2,2])

    #     if ylim:
    #         ax.set_ylim(ylim)

    #     if title:
    #         ax.set_title(title)
    #     else:
    #         ax.set_title(self.title)

    #     if show:
    #         plt.show()

    #     if save:
    #         if savename:
    #             plt.savefig(savename + '.' + imagetype)
    #         else:
    #             plt.savefig('extrema_batch_' + self.title + '.' + imagetype)

    #     plt.close(fig)


class CA:

    ''' Analyzes data from Gamry chronoamperometry data

    Pulls data from .dta file for plotting and analysis.

    '''

    def __init__(self, filename=None):
        ''' Opens file and retrieves data. '''

        # Pt  T     Vf    Im  Vu  Sig Ach IERange Over
        # s   V vs. Ref.  A   V   V   V   #       bits

        v_tolerance = 0.005 # [V], tolerance for determining voltage steps

        with open(filename, errors='replace') as f:
            rows = f.readlines()

            switch = False
            time = []
            voltage = []
            current = []

            for index, row in enumerate(rows):
                row = row.split()
                try:
                    if row:
                        # if row[0] == 'VPRESTEP':
                        #     self.Vstep0 = float(row[2])
                        if row[0] == 'NOTES':
                            try:
                                note = rows[index+2][1]
                                mass_search = re.search(r'\d{1}.\d{3}g', note)
                                self.mass = float(mass_search.group(0)[0:-1])
                            except Exception:
                                self.mass = 1

                        
                        elif row[0] == 'VSTEP1':
                            self.Vstep1 = float(row[2])
                        elif row[0] == 'VSTEP2':
                            self.Vstep2 = float(row[2])
                        elif row[0] == 'AREA':
                            self.area = float(row[2])
                        elif row[0] == 'EOC':
                            self.Eoc = float(row[2])
                        elif row[0] == 'CURVE':
                            switch = index + 2

                        # Save current cycle's data into temporary lists
                        if (self.is_num(row[0]) and switch and index > switch):
                            time.append(float(row[1]))
                            voltage.append(float(row[2]))
                            current.append(float(row[3]) 
                                * 1000/self.area)

                except Exception:
                    raise

            self.time = time
            self.voltage = voltage
            self.current = current

        Vstep1_mask, Vstep2_mask = [], []

        for time, voltage in zip(self.time, self.voltage):
            if time > 0 and np.abs(voltage-self.Vstep1) < v_tolerance:
                Vstep1_mask.append(1)
                Vstep2_mask.append(0)
            elif time > 0 and np.abs(voltage-self.Vstep2) < v_tolerance:
                Vstep1_mask.append(0)
                Vstep2_mask.append(1)
            else:
                Vstep1_mask.append(0)
                Vstep2_mask.append(0)

        Vstep1_pos = np.where(np.array(Vstep1_mask)==1)[0]
        Vstep2_pos = np.where(np.array(Vstep2_mask)==1)[0]

        self.Vstep1_idx = [Vstep1_pos[0], Vstep1_pos[-1]]
        self.Vstep2_idx = [Vstep2_pos[0], Vstep2_pos[-1]]

        self.ss_current1 = self.current[self.Vstep1_idx[1]]
        self.ss_current2 = self.current[self.Vstep2_idx[1]]

        self.mol_Zn, self.molal_Zn = utilities.get_init_concentrations(
            title=filename, mass=self.mass, test='CA')

        self.title = filename[:-4]


    def is_num(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False


    def plot_current_voltage(self, title=None, 
        show=False, save=False, savename=None):
        ''' Plots current vs voltage for one or all cycles '''

        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
            figsize=(12,9), dpi=75)

        ax1.plot(self.time, self.voltage,
                 linewidth=3, marker='o', markersize=5, color='#000000')
        ax2.plot(self.time, self.current,
                 linewidth=3, marker='o', markersize=5, color='#000000',)

        i_text = 'i$^1_{ss}$ = %.3f mA/cm$^2$'%self.ss_current1 + '\n' +\
            'i$^2_{ss}$ = %.3f mA/cm$^2$'%self.ss_current2
        x_tlim, y_tlim = ax2.get_xlim(), ax2.get_ylim()

        x_tpos = 0.67*(x_tlim[1]-x_tlim[0]) + x_tlim[0]
        y_tpos = 0.68*(y_tlim[1]-y_tlim[0]) + y_tlim[0]

        ax2.text(x_tpos, y_tpos, i_text, bbox=dict(facecolor='w', 
            edgecolor='k'))

        ax1.set_ylabel('Potential vs. Ref [V]')
        ax2.set_ylabel('Current [mA/cm$^2$]')
        ax2.set_xlabel('Time [s]')

        if title:
            fig.suptitle(title)
        else:
            fig.suptitle(self.title)

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title + '.png')

        if show:
            plt.show()

        plt.close(fig)


    def calculate_diffusivity(self):
        ''' Calculates diffusivity using Cottrell equation '''

        n = 2 # electrons involved in redox rxn
        F = 96485 # Faraday's constant, [coulombs]
        A = self.area # [cm^2]
        c0 = self.mol_Zn*1000 # [mol/cm^3]

        # bounds for window gradient to find coefs on Cottrell plot
        bounds = (0.05, 0.45)


        self.step1 = {'idx': self.Vstep1_idx, 'slope': -1}
        self.step2 = {'idx': self.Vstep2_idx, 'slope': 1}
        # self.step1 = {'idx': self.Vstep1_idx, 'slope': np.sign(self.Vstep1)}
        # self.step2 = {'idx': self.Vstep2_idx, 'slope': np.sign(self.Vstep2)}

        for step in [self.step1, self.step2]:
            # make step values at Vstep = 0V read as positive (oxidation)
            if step['slope'] == 0:
                step['slope'] == 1

            idx = step['idx']
            slope = step['slope']

            time_raw = self.time[idx[0]:idx[1]]
            time = [(t-time_raw[0])**(-0.5) for t in time_raw[1:]]
            current = self.current[idx[0]+1:idx[1]]

            coefs = utilities.window_gradient(x=time, y=current,
                bounds=bounds, slope=slope, lsq=True)
            line = [coefs[1] + t*coefs[0] for t in time]
            D = slope * coefs[0] * np.pi/((n*F*A*c0)**2)

            step['time'] = time
            step['current'] = current
            step['line'] = line
            step['D'] = D


    def plot_cottrell(self, show=False, save=False, title=None, 
        savename=None):
        ''' Plots Cottrell plots for both steps '''

        # Get diffusivity if not gathered already
        try:
            if self.step1 and self.step2:
                pass
        except AttributeError:
            self.calculate_diffusivity()

        font = {'family': 'Arial', 'size':20}
        matplotlib.rc('font', **font)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,9), dpi=75)

        ax1.plot(self.step1['time'], self.step1['current'], 
            linewidth=3, marker='o', markersize=5, color='#000000')
        ax1.plot(self.step1['time'], self.step1['line'], 
            linewidth=3, color='#FF0000',
            label='D$_0$ = %.3E cm$^2$/s'% Decimal(self.step1['D'])+\
            '\n'+'V = '+str(self.Vstep1)+'V')
        ax2.plot(self.step2['time'], self.step2['current'],  
            linewidth=3, marker='o', markersize=5, color='#000000')
        ax2.plot(self.step2['time'], self.step2['line'],
            linewidth=3, color='#FF0000',
            label='D$_0$ = %.3E cm$^2$/s'% Decimal(self.step2['D'])+\
            '\n'+'V = '+str(self.Vstep2)+'V')

        ax1.set_xlabel('t$^{-0.5}$, [s$^{-0.5}$]')
        ax1.set_ylabel('Current, I [mA/cm$^2$]')
        ax1.legend()

        ax2.set_xlabel('t$^{-0.5}$, [s$^{-0.5}$]')
        ax2.set_ylabel('Current, I [mA/cm$^2$]')
        ax2.legend()

        if title:
            fig.suptitle(title)
        else:
            fig.suptitle(self.title + '_diffusivity')

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title + '_diffusivity' + '.png')

        if show:
            plt.show()


class CA_batch:
    ''' Method for batch processing data from Gamry
        Uses methods defined in CA class

        Author: Bernard Kim
    '''

    def __init__(self, files, redox=None, date=None):
        # redox accepts 'red' or 'ox'
        # version is used for title only, accepts string, 
        # e.g. 'new', 'cycled', etc.

        # Accepts list of filenames
        self.alldata = {}
        self.Eoc = {}
        self.ss_currents = {}

        self.redox = redox

        params = {}

        for file in files:
            # Determining file name
            params_match = utilities.title_search_echem(filename=file, 
                test='CA')

            for key in params_match:
                if key not in params:
                    params[key] = []
                else:
                    params[key].append(params_match[key])

            # Data exporting and formatting for use
            exported = CA(file)

            if self.redox is 'red':
                V_step = exported.Vstep1
                ss_current = exported.ss_current1
            elif self.redox is 'ox':
                V_step = exported.Vstep2
                ss_current = exported.ss_current2

            self.alldata[V_step] = exported
            self.ss_currents[V_step] = ss_current

        self.title = ''
        for key in params:
            self.title = self.title + str(params[key][0]) + '_'

        self.title = self.title + self.redox + '_' + date


    def plot_voltage(self, title=None, 
        show=False, save=False, savename=None):
        ''' Plots current vs voltage '''

        font = {'family': 'Arial', 'size': 22}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0, 1, 10) # for use with tab10 colormap

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
            figsize=(12,9), dpi=75)

        for idx, V_step in enumerate(self.alldata):
            label = 'i$_{ss}$ = %.3f mA/cm$^2$'%self.ss_currents[V_step]

            ax1.plot(self.alldata[V_step].time, self.alldata[V_step].voltage,
                linewidth=3, marker='o', markersize=5, 
                color=plt.cm.tab10(coloridx[int(idx)]))
            ax2.plot(self.alldata[V_step].time, self.alldata[V_step].current,
                linewidth=3, marker='o', markersize=5, 
                color=plt.cm.tab10(coloridx[int(idx)]),
                label=label)

        ax1.set_ylabel('Potential vs. Ref [V]')
        ax2.set_ylabel('Current [mA/cm$^2$]')
        ax2.set_xlabel('Time [s]')

        ax2.legend()

        if title:
            fig.suptitle(title)
        else:
            fig.suptitle(self.title)

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title + '.png')

        if show:
            plt.show()

        plt.close(fig)


    def plot_cottrell(self, title=None,show=False, save=False, savename=None):
        ''' Plots Cottrell plots of CA curve for multiple samples'''

        for data in self.alldata:
            try:
                if self.alldata[data].step1 and self.alldata[data].step2:
                    pass
            except AttributeError:
                self.alldata[data].calculate_diffusivity()


        font = {'family': 'Arial', 'size': 22}
        matplotlib.rc('font', **font)

        coloridx = np.linspace(0, 1, 10) # for use with tab10 colormap

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,9), dpi=75)
        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        for idx, data in enumerate(self.alldata):
            if self.redox is 'red':
                ax.plot(self.alldata[data].step1['time'], 
                    self.alldata[data].step1['current'], 
                    linewidth=3, marker='o', markersize=5, color='#000000')
                ax.plot(self.alldata[data].step1['time'], 
                    self.alldata[data].step1['line'], 
                    linewidth=3, 
                    color=plt.cm.tab10(coloridx[int(idx)]),
                    label='D$_0$ = %.3E cm$^2$/s'% Decimal(
                        self.alldata[data].step1['D'])+\
                        '\n'+'V = '+str(self.alldata[data].Vstep1)+'V')
            elif self.redox is 'ox':
                ax.plot(self.alldata[data].step2['time'], 
                    self.alldata[data].step2['current'], 
                    linewidth=3, marker='o', markersize=5, color='#000000')
                ax.plot(self.alldata[data].step2['time'], 
                    self.alldata[data].step2['line'], 
                    linewidth=3, 
                    color=plt.cm.tab10(coloridx[int(idx)]),
                    label='D$_0$ = %.3E cm$^2$/s'% Decimal(
                        self.alldata[data].step2['D'])+\
                        '\n'+'V = '+str(self.alldata[data].Vstep2)+'V')

        ax.set_xlabel('t$^{-0.5}$, [s$^{-0.5}$]')
        ax.set_ylabel('Current, I [mA/cm$^2$]')
        ax.legend()

        if title:
            ax.set_title(title)
        else:
            ax.set_title(self.title + '_cottrell_batch_'+ self.redox)

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(self.title + '_cottrell_batch_' + \
                    self.redox + '.png')

        if show:
            plt.show()





































