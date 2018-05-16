import os
import glob
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy
# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

# favorite Seaborn settings for notebooks
rc={'lines.linewidth': 2,
    'axes.labelsize' : 16,
    'axes.titlesize' : 18,
    'axes.facecolor' : 'F4F3F6',
    'axes.edgecolor' : '000000',
    'axes.linewidth' : 1.2,
    'xtick.labelsize' : 13,
    'ytick.labelsize' : 13,
    'grid.linestyle' : ':',
    'grid.color' : 'a6a6a6'}
sns.set_context('notebook', rc=rc)
sns.set_style('darkgrid', rc=rc)
sns.set_palette("deep", color_codes=True)

#===============================================================================
# define variables to use over the script
date = 20170627
username = 'sbarnes'
# run = 'r1'
#operator = np.array(['01new007', '01new010', '01new011', '01new013','01new014'])

# read the CSV file with the mean fold change
df = pd.read_csv('output/' + str(date) + '_' + 'O1new' + \
        '_lacI_titration_MACSQuant.csv')
rbs = df.rbs.unique()
operator = df.operator.unique()

def set_color_cycle(self, clist=None):
    if clist is None:
        clist = rcParams['axes.color_cycle']
    self.color_cycle = itertools.cycle(clist)

#===============================================================================

# plot the curve for the 0 IPTG cultures
plt.figure(figsize=(7, 7))
binding_energy = df.binding_energy.unique()
for i, op in enumerate(operator):
    repressor_array = np.logspace(0, 3, 200)
    fc_theory = 1 / (1 + 2 * repressor_array / 5E6 * np.exp(- binding_energy[i]))

    plt.plot(repressor_array, fc_theory)

#reset color cycle to match colors
plt.gca().set_color_cycle(None)

#no_iptg = df.groupby('operator').get_group(0)
for op in operator:
    plt.plot(df[df.operator == op].repressors, df[df.operator == op].fold_change_A, \
     marker='o', linewidth=0, label = op)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('repressor copy number')
plt.ylabel('fold-change')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('output/'  + 'lacI_titration_ctrl.png')

# Produce a second plot showing the cultures on a semi-log scale.

plt.figure(figsize=(7, 7))
binding_energy = df.binding_energy.unique()
for i, op in enumerate(operator):
    repressor_array = np.logspace(0, 3, 200)
    fc_theory = 1 / (1 + 2 * repressor_array / 5E6 * np.exp(- binding_energy[i]))

    plt.plot(repressor_array, fc_theory)

#reset color cycle to match colors
plt.gca().set_color_cycle(None)

#no_iptg = df.groupby('operator').get_group(0)
for op in operator:
    plt.plot(df[df.operator == op].repressors, df[df.operator == op].fold_change_A, \
     marker='o', linewidth=0, label = op)
plt.xscale('log')
plt.xlabel('repressor copy number')
plt.ylabel('fold-change')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('output/'  + 'lacI_titration_semilog.png')
