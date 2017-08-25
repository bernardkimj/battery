''' Analysis module to calculate sheet resistance and related conductivity
    of printed electrodes.

    Author: Bernard Kim
    Principal Investigators: Paul Wright, James Evans
    Institution: University of California, Berkeley '''

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

        return table

    def plot_resistivity(self, inknames=None, xaxis=None, save=False):
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

        # Plotting mumbo-jumbo
        font = {'family': 'Arial', 'size': 14}
        matplotlib.rc('font', **font)

        fig, ax = plt.subplots()

        ax.errorbar(
            sortedtable[xaxis], res_par_mean, yerr=res_par_std,
            color='b', marker='.', markersize=8,
            capsize=10, elinewidth=2, label='parallel')
        ax.errorbar(
            sortedtable[xaxis], res_perp_mean, yerr=res_perp_std,
            color='g', marker='.', markersize=8,
            capsize=10, elinewidth=2, label='perpendicular')
        ax.legend()

        figtitle = electrode + ' resistivity by ' + plotstrings[xaxis]['titlelabel']

        ax.set_xlabel(plotstrings[xaxis]['xlabel'])
        ax.set_ylabel(r'$Resistivity, \Omega-m$')
        ax.set_title(figtitle)
        ax.grid()

        ymin, ymax = ax.get_ylim()
        ydelta = (ymax-ymin)
        ax.set_ylim([ymin, ydelta*1.1])

        par_height, perp_height = [], []
        for mean, std in zip(res_par_mean, res_par_std):
            par_height.append(mean+std)
        for mean, std in zip(res_perp_mean, res_perp_std):
            perp_height.append(mean+std)

        textheights = []
        for par, perp in zip(par_height, perp_height):
            textheights.append(max(par, perp) + ydelta*0.025)

        labels = sortedtable['n_samples']

        for xval, height, label in zip(sortedtable[xaxis], textheights, labels):
            ax.text(xval, height, 'n=' + str(label), ha='center', va='bottom', fontsize=12)

        if save:
            plt.savefig(electrode + xaxis + '.pdf', format='pdf')

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
















