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

# read the CSV file with the mean fold change
df = pd.read_csv('output/20160829_lacI_titration_MACSQuant.csv')
rbs = df.rbs.unique()
#=============================================================================== 

# compute the theoretical repression level
repressor_array = np.logspace(1, 3, 100)
epsilon_array = np.array([-15.3, -13.9, -9.7, -17])
operators = np.array(['O1new010', 'O1new011', 'O1new013', 'O1new014', 'O1new007'])

colors = sns.hls_palette(len(operators), l=.3, s=.8)
# plot theoretical curve
plt.figure(figsize=(7, 7))
for i, o in enumerate(operators):
    fold_change_theor = 1 / (1 + 2 * repressor_array / 5E6 \
            * np.exp(-epsilon_array[i]))
    plt.plot(repressor_array, fold_change_theor, label=o,
            color=colors[i])
    plt.plot(df[(df.operator == o) & (df.rbs != 'auto') & \
            (df.rbs != 'delta')].repressors, 
            df[(df.operator == o) & (df.rbs != 'auto') & \
            (df.rbs != 'delta')].fold_change,
            marker='o', linewidth=0, color=colors[i], 
            label=o + ' flow cytometer')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('repressor copy number')
plt.ylabel('fold-change')
plt.xlim(left=10)
plt.legend(loc='lower left')
plt.tight_layout()
plt.savefig('output/lacI_titration.png')
