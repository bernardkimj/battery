''' Analysis module to calculate sheet resistance and related conductivity
    of printed electrodes.

    Author: Bernard Kim
    Principal Investigators: Paul Wright, James Evans
    Institution: University of California, Berkeley '''

from battery.utilities import utilities
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import re


class resistivity:

    ''' Analyzes sheet resistance data
        Measured and recorded by hand into custom Excel spreadsheets '''

    def __init__(self, filenames):
        ''' Loads datasheets and organizes in dictionary '''

        self.data = {}

        for filename in filenames:
            # Get experiment date from filename
            date = re.search(r'\d{8}', filename).group(0)

            # Load Excel file
            xls_file = pd.ExcelFile(filename)

            # Load data sheets
            df = xls_file.parse('Samples')
            df = df.dropna()
            inks = xls_file.parse('Ink Recipes')
            inks = inks.dropna()

            # Identify distinct ink names using a set
            inknames = list(set(inks['Ink Name']))

            # Gather all samples by ink name
            for inkname in inknames:
                # Instantiate data entry for new ink names
                if inkname not in self.data:
                    self.data[inkname] = {'dates': {}}

                # Store sample data and ink recipes for each datasheet by date and ink name
                self.data[inkname]['dates'][date] = {
                    'data': df.loc[df['Ink Name'] == inkname],
                    'recipedf': inks.loc[inks['Ink Name'] == inkname],
                }

            # print('Completed importing ' + filename)

        resistivity.get_statistics(self)
        resistivity.get_recipe(self)

    def get_statistics(self):
        ''' Gathers all resistivity and thickness measurements by ink name 
            and calculates statistics 
            Stores values in dictionary by inkname'''

        for inkname in self.data:
            alldf = pd.DataFrame()

            for date in self.data[inkname]['dates']:
                # Compile all data for each inkname for all dates 
                # and add date column to df
                df = self.data[inkname]['dates'][date]['data']
                df.insert(len(list(df)), 'date', date)
                alldf = alldf.append(df)

                self.data[inkname]['dates'][date]['statistics'] = {
                    'n_samples': len(df),
                    'thickness': {
                        'alldata': df['Measured Thickness (mm)'],
                        'mean': np.mean(df['Measured Thickness (mm)']),
                        'std': np.std(df['Measured Thickness (mm)']),
                    },
                    'resistivity': {
                        'parallel': {
                            'alldata': df['Resistivity, Parallel (Ωm)'],
                            'mean': np.mean(df['Resistivity, Parallel (Ωm)']),
                            'std': np.std(df['Resistivity, Parallel (Ωm)']),
                        },
                        'perpendicular': {
                            'alldata': df['Resistivity, Perpendicular (Ωm)'],
                            'mean': np.mean(df['Resistivity, Perpendicular (Ωm)']),
                            'std': np.std(df['Resistivity, Perpendicular (Ωm)']),
                        },
                    }
                }

            thickness = alldf['Measured Thickness (mm)']
            res_parallel = alldf['Resistivity, Parallel (Ωm)']
            res_perpendicular = alldf['Resistivity, Perpendicular (Ωm)']

            self.data[inkname]['statistics'] = {
                'raw_data': alldf,
                'n_samples': len(alldf),
                'thickness': {
                    'alldata': thickness,
                    'mean': np.mean(thickness),
                    'std': np.std(thickness),
                },
                'resistivity': {
                    'parallel': {
                        'alldata': res_parallel,
                        'mean': np.mean(res_parallel),
                        'std': np.std(res_parallel),
                    },
                    'perpendicular': {
                        'alldata': res_perpendicular,
                        'mean': np.mean(res_perpendicular),
                        'std': np.std(res_perpendicular),
                    },
                },
            }

            # print('Finished importing ' + inkname)

    def electrode_type(self, inkname=None):
        ''' Returns string for active material based on electrode type '''

        if re.search(r'\ACFA', inkname) is not None:
            activematerial = 'MnO2'
        elif re.search(r'\AAFA', inkname) is not None:
            activematerial = 'Zn'

        return activematerial

    def get_recipe(self):
        ''' Organizes ink recipe to organize data as a percentage of desired
            material '''

        AB = 'Acetylene Black'

        for inkname in self.data:
            activematerial = self.electrode_type(inkname=inkname)

            for date in self.data[inkname]['dates']:
                rawrecipe = self.data[inkname]['dates'][date]['recipedf'].set_index('Species')

                self.data[inkname]['dates'][date]['recipe'] = {
                    'ideal': {
                        activematerial: rawrecipe.loc[activematerial, 'Ideal Amount (g)'],
                        AB: rawrecipe.loc[AB, 'Ideal Amount (g)'],
                        'PVDF-HFP': rawrecipe.loc['PVDF-HFP', 'Ideal Amount (g)'],
                        'NMP': rawrecipe.loc['NMP (in gel)', 'Ideal Amount (g)'] + 
                            rawrecipe.loc['NMP (additional)', 'Ideal Amount (g)'],
                    },
                    'actual': {
                        activematerial: rawrecipe.loc[activematerial, 'Actual Amount (g)'],
                        AB: rawrecipe.loc[AB, 'Actual Amount (g)'],
                        'PVDF-HFP': rawrecipe.loc['PVDF-HFP', 'Actual Amount (g)'],
                        'NMP': rawrecipe.loc['NMP (in gel)', 'Actual Amount (g)'] + 
                            rawrecipe.loc['NMP (additional)', 'Actual Amount (g)'],
                    },
                    'temperature': rawrecipe.loc[activematerial, 'Room Temperature (°C)'],
                    'humidity': rawrecipe.loc[activematerial, 'Room Humidity (%RH)'],
                }

                ideal = self.data[inkname]['dates'][date]['recipe']['ideal']
                actual = self.data[inkname]['dates'][date]['recipe']['actual']
                
                if 'fractions' not in self.data[inkname]:
                    sumfraction = np.sum([
                        ideal[activematerial], 
                        ideal[AB], 
                        ideal['PVDF-HFP'],
                        ])

                    fractions = {
                        activematerial: ideal[activematerial]/sumfraction,
                        AB: ideal[AB]/sumfraction,
                        'PVDF-HFP': ideal['PVDF-HFP']/sumfraction,
                        'solidphase': sumfraction/(sumfraction + ideal['NMP']),
                    }

                    self.data[inkname]['fractions'] = fractions


                error = {
                    activematerial: (actual[activematerial]-ideal[activematerial])/ideal[activematerial],
                    AB: (actual[AB]-ideal[AB])/ideal[AB],
                    'PVDF-HFP': (actual['PVDF-HFP']-ideal['PVDF-HFP'])/ideal['PVDF-HFP'],
                    'NMP': (actual['NMP']-ideal['NMP'])/ideal['NMP'],
                }

                self.data[inkname]['dates'][date]['recipe']['error'] = error

    def compile_table(self, inknames=None):
        ''' Compiles all relevant data based on specified inknames for easy plotting '''

        # Checks first inkname to determine electrode type (cathode or anode)
        activematerial = self.electrode_type(inkname=inknames[0])

        headers = [
            'inkname',
            'n_samples',
            activematerial,
            'Acetylene Black',
            'PVDF-HFP',
            'solidphase',
            'thickness_mean',
            'thickness_std',
            'res_par_mean',
            'res_par_std',
            'res_perp_mean',
            'res_perp_std',
        ]

        rawtable = []

        for inkname in inknames:
            stats = self.data[inkname]['statistics']
            fractions = self.data[inkname]['fractions']

            row = (
                inkname, 
                stats['n_samples'], 
                fractions[activematerial]*100, 
                fractions['Acetylene Black']*100, 
                fractions['PVDF-HFP']*100, 
                fractions['solidphase']*100, 
                stats['thickness']['mean'], 
                stats['thickness']['std'], 
                stats['resistivity']['parallel']['mean'], 
                stats['resistivity']['parallel']['std'], 
                stats['resistivity']['perpendicular']['mean'], 
                stats['resistivity']['perpendicular']['std'], 
            )

            rawtable.append(row)

        table = pd.DataFrame.from_records(rawtable, columns=headers)
        table.to_csv(activematerial+'_statistics.csv')


        return table


    def plot_resistivity(self, inknames=None, xaxis=None, 
        title=None, show=False, save=False, savename=None):
        ''' Plots resistivity as a function of specified x-axis
            'AB', 'PVDF-HFP', 'activematerial', 'solidphase', 'thickness' 
            Must specify inknames to include in plot via 'inknames' kwarg '''

        # Import compiled table based on desired inknames to plot
        table = self.compile_table(inknames)
        
        # Sets which electrode type specified inks belong to
        activematerial = self.electrode_type(inkname=inknames[0])

        if activematerial == 'Zn':
            electrode = 'Anode'
            par_color = '#6699ff'
            perp_color = '#0000ff'
        elif activematerial == 'MnO2':
            electrode = 'Cathode'
            par_color = '#ef2e2e'
            perp_color = '#420d0d'

        # Sorts imported table based on desired xaxis param
        sortedtable = table.sort_values(xaxis)

        # Set strings for plot title and labels
        plotstrings = {
            'Acetylene Black': {
                'xlabel': 'wt% Acetylene Black',
                'titlelabel': 'Acetylene Black fraction',
            },
            'PVDF-HFP': {
                'xlabel': 'wt% PVDF-HFP',
                'titlelabel': 'PVDF-HFP fraction',
            },
            'activematerial': {
                'xlabel': 'wt% ' + activematerial,
                'titlelabel': activematerial + ' fraction',
            },
            'solidphase': {
                'xlabel': 'wt% Solid Phase',
                'titlelabel': 'Solid Phase',
            },
            'thickness': {
                'xlabel': 'Sample Thickness (mm)',
                'titlelabel': 'Sample Thickness',
            },
        }

        res_par_mean = sortedtable['res_par_mean']
        res_par_std = sortedtable['res_par_std']
        res_perp_mean = sortedtable['res_perp_mean']
        res_perp_std = sortedtable['res_perp_std']

        labels = sortedtable['n_samples']

        # Plotting mumbo-jumbo
        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(18,10), dpi=100)

        ax.errorbar(
            sortedtable[xaxis], 
            res_par_mean, 
            yerr=res_par_std,
            color=par_color, 
            marker='.', 
            markersize=8, 
            linewidth=3,
            capsize=10, 
            elinewidth=3, 
            markeredgewidth=3, 
            label='Parallel'
        )

        ax.errorbar(
            sortedtable[xaxis], 
            res_perp_mean, 
            yerr=res_perp_std,
            color=perp_color, 
            marker='.', 
            markersize=8, 
            linewidth=3,
            capsize=10, 
            elinewidth=3, 
            markeredgewidth=3, 
            label='Perpendicular'
        )

        ax.legend()
        ax.set_xlabel(plotstrings[xaxis]['xlabel'])
        ax.set_ylabel('Resistivity, ' + r'$\Omega-m$')

        ymin, ymax = ax.get_ylim()
        ydelta = (ymax-ymin)
        ax.set_ylim([0, ydelta*1.1])

        heights = [max(par_mean+par_std, perp_mean+perp_std) 
            for par_mean, par_std, perp_mean, perp_std 
            in zip(res_par_mean, res_par_std, res_perp_mean, res_perp_std)]

        for xval, height, label in zip(sortedtable[xaxis], heights, labels):
            ax.text(xval, height + 0.05*ydelta, 'n=' + str(label), 
                ha='center', va='bottom', fontsize=20)

        # if title:
        #     ax.set_title(title)
        # else:
        #     ax.set_title(electrode + ' resistivity by ' + plotstrings[xaxis]['titlelabel'])

        plt.tight_layout()

        if save:
            if savename:
                plt.savefig(savename + '.png')
            else:
                plt.savefig(electrode + xaxis + '.png')

        if show:
            plt.show()


    def plot_resistivity_bar(self, inknames=None, xaxis=None, show=False, save=False):
        ''' Plots resistivity as a function of specified x-axis
            'AB', 'PVDF-HFP', 'activematerial', 'solidphase', 'thickness' 
            Must specify inknames to include in plot via 'inknames' kwarg '''

        # Import compiled table based on desired inknames to plot
        table = self.compile_table(inknames)
        
        # Sets which electrode type specified inks belong to
        activematerial = self.electrode_type(inkname=inknames[0])

        if activematerial == 'Zn':
            electrode = 'Anode'
        elif activematerial == 'MnO2':
            electrode = 'Cathode'

        # Sorts imported table based on desired xaxis param
        sortedtable = table.sort_values(xaxis)

        # Set strings for plot title and labels
        plotstrings = {
            'Acetylene Black': {
                'xlabel': 'wt% Acetylene Black',
                'titlelabel': 'Acetylene Black fraction',
            },
            'PVDF-HFP': {
                'xlabel': 'wt% PVDF-HFP',
                'titlelabel': 'PVDF-HFP fraction',
            },
            'activematerial': {
                'xlabel': 'wt% ' + activematerial,
                'titlelabel': activematerial + ' fraction',
            },
            'solidphase': {
                'xlabel': 'wt% Solid Phase',
                'titlelabel': 'Solid Phase',
            },
            'thickness': {
                'xlabel': 'Sample Thickness (mm)',
                'titlelabel': 'Sample Thickness',
            },
        }

        res_par_mean = sortedtable['res_par_mean']
        res_par_std = sortedtable['res_par_std']
        res_perp_mean = sortedtable['res_perp_mean']
        res_perp_std = sortedtable['res_perp_std']

        labels = sortedtable['n_samples']

        # Plotting mumbo-jumbo
        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots(figsize=(18,10), dpi=75)
        width = 0.15

        switch = 0
        # labels = []

        rect_par = ax.bar(
            sortedtable[xaxis] - width/2, 
            res_par_mean,
            width, 
            yerr=res_par_std,
            color='#6699ff',
            label='Parallel'
        )

        for xpos, mean, std in zip(sortedtable[xaxis], res_par_mean, res_par_std):
            ax.text(
                xpos - width/2,
                mean + std,
                # '%d \n'%par_avg + r'$\pm$' +'%d'%par_std,
                '%d '%mean + r'$\pm$' +' %d'%std,
                ha='center',
                va='bottom',
                rotation=80,
                fontsize=18,
            )

        rect_perp = ax.bar(
            sortedtable[xaxis] + width/2, 
            res_perp_mean,
            width, 
            yerr=res_perp_std,
            color='#0000ff',
            label='Perpendicular'
        )

        for xpos, mean, std in zip(sortedtable[xaxis], res_perp_mean, res_perp_std):
            ax.text(
                xpos + width/2,
                mean + std,
                # '%d \n'%par_avg + r'$\pm$' +'%d'%par_std,
                '%d '%mean + r'$\pm$' +' %d'%std,
                ha='center',
                va='bottom',
                rotation=80,
                fontsize=18,
            )

        if not switch:
            ax.legend()
            switch +=1

        xticks = sortedtable[xaxis]
        xlabels = sortedtable[xaxis]

        plt.xticks(xticks, xlabels)

        # ax.set_ylim([0, 17500])
        # ax.set_ylim([0, 21500])

        ax.set_xlabel(plotstrings[xaxis]['xlabel'])
        ax.set_ylabel(r'$Resistivity, \Omega-m$')

        if save:
            plt.savefig(electrode + xaxis + '.png')
        if show:
            plt.show()


    def plot_thickness(self, inkname=None, date=None, save=False):
        ''' Plots resistivity as a function of thickness for batch of samples
            Must specify inknames to include in plot via 'inknames' kwarg '''

        thickness = self.data[inkname]['dates'][date]['statistics']\
            ['thickness']['alldata']
        res_par = self.data[inkname]['dates'][date]['statistics']\
            ['resistivity']['parallel']['alldata']
        res_perp = self.data[inkname]['dates'][date]['statistics']\
            ['resistivity']['perpendicular']['alldata']

        # Plotting mumbo-jumbo
        font = {'family': 'Arial', 'size': 14}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots()

        ax.scatter(thickness*1000, res_par, label='parallel')
        ax.scatter(thickness*1000, res_perp, label='perpendicular')

        ax.legend()

        figtitle = inkname + ' Thickness vs. Resistivity'

        ax.set_ylim([0, np.max([res_par, res_perp])*1.10])
        ax.set_xlabel('Sample Thickness (µm)')
        ax.set_ylabel(r'$Resistivity, \Omega-m$')
        ax.set_title(figtitle)
        ax.grid()

        if save:
            plt.savefig('thickness_' + inkname + '_' + date + '.pdf', format='pdf')

        plt.show()


class cc_resistivity:

    ''' Analyzes sheet resistance and conductivity data
        Measured and recorded by hand into custom Excel spreadsheet '''

    def __init__(self, filename):
        ''' Loads datasheets and organizes in dictionary '''

        # Assumes csv with column order of 
        # 'Ink Name', 'Sample Number', 'Thickness (mm)', 
        # 'R_perp (Ωm)', 'R_par (Ωm)'

        # Instantiate dict for indices within row (agnostic to export order)
        idxs = {}
        self.data = {}

        reader = csv.reader(open(filename, errors='replace'), delimiter=',')
        rows = list(reader)

        headers = {
            'Thickness (mm)': 'thickness',
            'R_perp (Ωm)': 'r_perp',
            'R_par (Ωm)': 'r_par',
            'C_perp (S/m)': 'c_perp',
            'C_par (S/m)': 'c_par',
        }

        for index, row in enumerate(rows):
            if index == 0:
                for header in headers:
                    if header in row:
                        idxs[headers[header]] = row.index(header)

            if index > 0:
                inkname = row[0]
                # samplenum = row[1]
                # thickness = row[2]
                # r_perp = row[3]
                # r_par = row[4]

                freq = int(inkname[0:2])
                ball = int(inkname[2:4])
                time = int(inkname[5:])

                if freq not in self.data:
                    self.data[freq] = {}

                if ball not in self.data[freq]:
                    self.data[freq][ball] = {}

                if time not in self.data[freq][ball]:
                    self.data[freq][ball][time] = {
                        headers[header]: [] for header in headers
                    }

                # sample = {headers[header]: [] for header in headers}

                for header in headers:
                    try:
                        self.data[freq][ball][time][headers[header]].append(
                            float(row[idxs[headers[header]]]))
                    except ValueError:
                        self.data[freq][ball][time][headers[header]].append(
                            float('nan'))

        for freq in self.data:
            for ball in self.data[freq]:
                for time in self.data[freq][ball]:
                    # with self.data[freq][ball][time] as sample:
                    self.data[freq][ball][time]['thickness_avg'] = \
                        np.nanmean(self.data[freq][ball][time]['thickness'])
                    self.data[freq][ball][time]['thickness_std'] = \
                        np.nanstd(self.data[freq][ball][time]['thickness'])
                    self.data[freq][ball][time]['r_perp_avg'] = \
                        np.nanmean(self.data[freq][ball][time]['r_perp'])
                    self.data[freq][ball][time]['r_perp_std'] = \
                        np.nanstd(self.data[freq][ball][time]['r_perp'])
                    self.data[freq][ball][time]['r_par_avg'] = \
                        np.nanmean(self.data[freq][ball][time]['r_par'])
                    self.data[freq][ball][time]['r_par_std'] = \
                        np.nanstd(self.data[freq][ball][time]['r_par'])
                    self.data[freq][ball][time]['c_perp_avg'] = \
                        np.nanmean(self.data[freq][ball][time]['c_perp'])
                    self.data[freq][ball][time]['c_perp_std'] = \
                        np.nanstd(self.data[freq][ball][time]['c_perp'])
                    self.data[freq][ball][time]['c_par_avg'] = \
                        np.nanmean(self.data[freq][ball][time]['c_par'])
                    self.data[freq][ball][time]['c_par_std'] = \
                        np.nanstd(self.data[freq][ball][time]['c_par'])


    def plot_conductivity_bar(self, show=False, save=False):

        ''' Plots bar graphs of conductivities
            Intended for three separate figures organized in heirarchy of 
            Frequency >> ball size >> milling time '''

        # Plotting formatting
        font = {'family': 'Arial', 'size': 24}
        matplotlib.rc('font', **font)

        timenum = 3
        ballnum = 4

        for freq in sorted(self.data.keys()):
            fig, ax = plt.subplots(figsize=(18,10), dpi=75)
            width = 0.5

            # # each direction gets 1/2 a unit, space between each ball size
            # # gets one fill unit
            # ind = np.arange((3+1)*(4))

            switch = 0
            labels = []

            for idb, ball in enumerate(sorted(self.data[freq].keys()), 1):
                for idt, time in enumerate(sorted(self.data[freq][ball].keys()), 1):
                    x_loc = 1 + (ballnum*(idb-1) + idt)

                    x_loc_par = x_loc - width/2
                    x_loc_perp = x_loc + width/2

                    par_avg = self.data[freq][ball][time]['c_par_avg']
                    par_std = self.data[freq][ball][time]['c_par_std']
                    perp_avg = self.data[freq][ball][time]['c_perp_avg']
                    perp_std = self.data[freq][ball][time]['c_perp_std']


                    rect_par = ax.bar(
                        x_loc_par, 
                        par_avg,
                        width, 
                        yerr=par_std,
                        color='#6699ff',
                        label='Parallel'
                    )

                    rect_perp = ax.bar(
                        x_loc_perp, 
                        perp_avg,
                        width, 
                        yerr=perp_std,
                        color='#0000ff',
                        label='Perpendicular'
                    )

                    ax.text(
                        x_loc_par,
                        par_avg + par_std + 500,
                        # '%d \n'%par_avg + r'$\pm$' +'%d'%par_std,
                        '%d '%par_avg + r'$\pm$' +' %d'%par_std,
                        ha='center',
                        va='bottom',
                        rotation=80,
                        fontsize=18,
                    )

                    ax.text(
                        x_loc_perp,
                        perp_avg + perp_std + 500,
                        # '%d \n'%perp_avg + r'$\pm$' +'%d'%perp_std,
                        '%d '%perp_avg + r'$\pm$' +' %d'%perp_std,
                        ha='center',
                        va='bottom',
                        rotation=80,
                        fontsize=18,
                    )

                    if not switch:
                        ax.legend()
                        switch +=1

                    if idt == 2:
                        labels.append((x_loc, '%sh \n %smm' %(time, ball)))
                    else:
                        labels.append((x_loc, '%sh' %time))

            # for xval, height, label in zip(sortedtable[xaxis], par_height, labels):
            #     ax.text(xval, height, 'n=' + str(label), ha='center', va='bottom', fontsize=24)

            xticks = [label[0] for label in labels]
            xlabels = [label[1] for label in labels]

            plt.xticks(xticks, xlabels)

            # ax.set_ylim([0, 17500])
            ax.set_ylim([0, 21500])

            # ax.set_xlabel('Balls and Times')
            ax.set_ylabel('Conductivity [S/m]')
            ax.set_title('%s Hz'%freq)

            if save:
                plt.savefig('Niballmilling_%sHz.png'%freq)
            if show:
                plt.show()























