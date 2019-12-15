''' Module containing common data processing scripts used across
    multiple modules

    Author: Bernard Kim
    Principal Investigators: Prof. Paul Wright, Prof. James Evans
    University: University of California, Berkeley
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re
import csv
from operator import itemgetter
import scipy.stats as stats
import json

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
            raise ValueError('batch_average_prep: different number of \
                x and y vectors')

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
                R = stats.t.interval(alpha=confidence, df=len(y)-1, loc=m, 
                    scale=s/np.sqrt(len(y)))
                lcl.append(R[0])
                ucl.append(R[1])

        return indep, mean, std, lcl, ucl

    def title_search_format(filename, include=None):
        ''' 
        Searches filenames for batch processing for multiple samples
        Returns formatted title with constant sample parameters for 
        batch plotting and sample idx

        Must specify what to search for within each title string, will depend 
        on type of data (e.g. Dektak, Brookfield, Biologic, etc.)

        include must take the form of a list of strings
        Possibilities for strings are ballmilling', 'casting'

        '''

        allmatches = {}
        orderedparams = {}

        # General title strings
        ink = re.search(r'[A-Z]{3}\d{1}[a-z]{1}\d{2}', filename) # ink name
        samplenum = re.search(r'S\d{1,2}', filename, re.IGNORECASE) # sample number
        date = re.search(r'[0-9]{8}', filename) # date

        allmatches['inktype'] = {
            'ink': ink, 
            'date': date,  
        }
        orderedparams['inktype'] = ['ink', 'date']

        if samplenum:
            sampleidx = samplenum.group(0)
        else:
            raise ValueError('No sample number!')

        if include:
            # ball milling parameters, in order of ball diameter, frequency, 
            # time
            if 'ballmilling' in include:
                BMparams = re.search(r'[_,-][0-9]{6}(S)?', filename) 
                sieving = re.search(r'M[0-9]{2,3}', filename)
                isold = re.search(r'old', filename) # if old ball milling jars
                isTEG = re.search(r'TEG', filename) # if TEG jars

                allmatches['ballmilling'] = {
                    'BMparams': BMparams, 
                    'sieving': sieving, 
                    'isold': isold, 
                    'isTEG': isTEG, 
                }
                orderedparams['ballmilling'] = [
                    'BMparams', 'sieving', 'isold', 'isTEG']

            # blade coater
            if 'casting' in include:
                # casting direction
                pa = re.search(r'pa', filename)
                pe = re.search(r'pe', filename)

                # doctor blade height, stencil thickness, casting speed
                DB = re.search(r'\d{2,3}umDB', filename, re.IGNORECASE)
                ST = re.search(r'\d{2,3}umST', filename, re.IGNORECASE)
                speed = re.search(r'\d{2,3}mms', filename, re.IGNORECASE)

                allmatches['casting'] = {
                    'DB': DB, 
                    'ST': ST, 
                    'speed': speed, 
                    'pa': pa, 
                    'pe': pe,
                }
                orderedparams['casting'] = ['DB', 'ST', 'speed', 'pa', 'pe']

        grouporder = ['inktype'] + [group for group in include]
        title = ''

        for group in grouporder:
            for item in orderedparams[group]:
                if allmatches[group][item] is not None:
                    append = allmatches[group][item].group(0)

                    front, back = 0, len(append)
                    if append[0] == '_' or append[0] == '-':
                        front = 1
                    if append[-1] == '_' or append[-1] == '-':
                        back = -1

                    append = append[front:back]

                    title = title + append + '-'
                # with allmatches[group][item] as match:
                #     if match is not None:
                #         title = title + match.group(0) + '_'

        title = title[:-1] # cut off last underscore
        inkname = ink.group(0)

        return title, inkname, sampleidx

    def title_search_echem(filename, test=None):
        '''
        Searches batches of filenames for batch processing
        Returns formatted title depending on variable parameter
        Only for use with gamry module

        Must specify type of electrochemical test: 'CV' or 'CA'
        '''

        allmatches = {}

        if test is 'CV':
            IL = re.search(r'[A-Z]{1}MIM[a-zA-Z\d]{3,4}', filename).group(0)
            mol = re.search(r'(0,\d{1}m|neat)', filename).group(0)
            nu = re.search(r'\d{2,3}mVs', filename).group(0)

            allmatches = {
                'IL': IL,
                'mol': mol,
                'nu': nu,
            }

        elif test is 'CA':
            IL = re.search(r'[A-Z]{1}MIM[a-zA-Z\d]{3,4}', filename).group(0)
            mol = re.search(r'(0,\d{1}m|neat)', filename).group(0)

            allmatches = {
                'IL': IL,
                'mol': mol,
            }

        return allmatches

    def plot_setup(item, params):

        # Plotting formatting
        font = {'family': 'Arial', 'size': 20}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(16,9), dpi=75)

        if params['xlim']:
            ax.set_xlim(params['xlim'])
        if params['ylim']:
            ax.set_ylim(params['ylim'])

        if params['title']:
            ax.set_title(params['title'])
        else:
            title = item.title + params['titletag']
            ax.set_title(title)

        ax.set_xlabel(params['xlabel'])
        ax.set_ylabel(params['ylabel'])
        ax.grid()

        return fig, ax

    def save_json(data, filename=None):
        ''' Creates .json file to save desireable metrics
            Accepts dictionary format of data to save and filename

            Checks first if filename exists and adds to it
            Otherwise, creates new file

        '''

        output = json.dumps(data, separators=(',', ': '), indent=4)

        with open(filename, 'w') as f:
            f.write(output)

    def read_json(filename=None):
        ''' Reads and saves .json file
            Converts back into dictionary (if saved as dictionary)
        '''

        with open(filename) as f:
            data = json.load(f)

        return data

    def get_init_concentrations(title=None, mass=1, test=None):

        R = 8.314 # Universal gas constant, [J/mol*K]
        T = 25+273.15 # Temperature, [K]
        F = 96485 # Faraday's constant, [coulombs]
        MW_Zn = 65.38 # MW of Zn, [g/mol]
        MW_ZnOtf = 363.51 # MW of ZnOtf, [g/mol]
        n = 2 # electrons involved in redox rxn
        s = 1 # stoichiometric constant of Zn/Zn2+

        m_total_init = mass # [g] initial mass, assume about 1 g

        rho_IL = {
            'EMIMOTF': 1.387, # [g/mL]
            'BMIMOTF': 1.292, # [g/mL]
        }

        MW_IL = {
            'EMIMOTF': 260.23, # [g/mol]
            'BMIMOTF': 288.29, # [g/mol]
        }

        titlematches = utilities.title_search_echem(title, test=test)
        molsearch = re.search(r'0,\d{1}', titlematches['mol'])

        # Get concentration and IL from filename
        molal = float('0.'+molsearch.group(0)[2]) # mol/kg
        IL = titlematches['IL']

        # mass of salt [g]
        m_ZnOtf_init = m_total_init/(1 + (1/(molal*0.001*MW_ZnOtf))) 
        m_Zn_init = MW_Zn/MW_ZnOtf * m_ZnOtf_init # [g]
        m_IL_init = m_total_init - m_ZnOtf_init # [g]

        mol_Zn = (m_Zn_init/MW_Zn)/((m_IL_init/rho_IL[IL])*0.001) # [mol/L]
        molal_Zn = (m_Zn_init/MW_Zn)/(m_IL_init*0.001) # [mol/kg]

        return mol_Zn, molal_Zn

    def window_gradient(x, y, bounds, slope, windownum=50, positionnum=50, 
        lsq=False):
        ''' Applies a window-gradient algorithm to find optimal linear fit
            Provide x, y data, window size boundaries, sign of slope to find
            bounds is a tuple (small, large) that specifies window size
            Returns coefficients of optimal linear fit
        '''

        linear_tol = 0.05 # tolerance value for linear residual fit

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
                
                if lsq:
                    residuals = np.sum((test_line-y)**2)
                else:
                    residuals = 0

                    for raw, fit in zip(y, test_line):
                        if np.abs(raw/fit) > linear_tol:
                            residuals += 1

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

    def moving_average(interval, size, weight=None):
        ''' Calculates moving average to smooth noise from data 
            Choose weighting option based on kwarg 'weight'
            Based on convolution, so weights must be centered on window
            'size' is desired fraction of interval to be considered for smoothing
            '''

        n = int(np.ceil(len(interval)*size))

        if weight == None or weight == 'unweighted':
            window = np.ones(n)/n
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

        smoothed = np.convolve(interval, window, 'valid')

        overhang = len(interval) - len(smoothed)
        frontpad = int(np.ceil(overhang/2))
        endpad = int(np.floor(overhang/2))
        idxs = (frontpad, len(interval)-endpad)

        return smoothed, idxs

    def check_forward_monotonicity(series, type, length):
        ''' Checks if list-type object is monotonically increasing 
            or decreasing over a specified length
            type accepts 'increasing' or decreasing'
            series must be 1D list-type object
            length accepts integer value
            returns boolean truth value
        '''

        checklen = length
        checklist = []

        for idc in range(1, checklen+1):
            if type == 'increasing':
                if (series[idc] > series[idc-1]):
                    checklist.append(True)
                else:
                    checklist.append(False)
            elif type == 'decreasing':
                if (series[idc] < series[idc-1]):
                    checklist.append(True)
                else:
                    checklist.append(False)

        if np.sum(checklist) == checklen:
            return True
        else:
            return False


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


















