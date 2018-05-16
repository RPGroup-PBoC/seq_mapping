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

# list the directory with the data
datadir = '../../../data/flow/csv/'
files = np.array(os.listdir(datadir))
csv_bool = np.array(['csv' in f for f in files])
files = files[np.array(csv_bool)]

# define the patterns in the file names to read them
operators = np.array(['O1new010', 'O1new011', 'O1new013', 'O1new014', 'O1new007'])
energies = np.array([-15.3, -13.9, -9.7, -17, -10])
rbs = np.array(['auto', 'delta', 'RBS1L', 'RBS1', 'RBS1027', 'RBS446', 'RBS1147', 'HG104'])
repressors = np.array([0, 0, 870, 610, 130, 62, 30, 11])

#=============================================================================== 
fsc_range = [5E3, 3E4]
ssc_range = [2E4, 1E5]

# scatter-plot of the side scattering vs the frontal scattering
df_example = pd.read_csv(datadir + files[0])
df_sample = df_example.sample(frac=0.05)

plt.figure()
plt.scatter(df_sample['FSC-A'], df_sample['SSC-A'], alpha=0.2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.axvline(x=fsc_range[0], color='darkred')
plt.axvline(x=fsc_range[1], color='darkred')
plt.axhline(y=ssc_range[0], color='darkred')
plt.axhline(y=ssc_range[1], color='darkred')
plt.tight_layout()
plt.savefig('output/FSC_SSC.png')

#=============================================================================== 
# initialize the DataFrame to save the mean expression levels
df = pd.DataFrame()
# read the files and compute the mean YFP value
for i, operator in enumerate(operators):
    for j, strain in enumerate(rbs):
        # find the file
        r_file = glob.glob(datadir + '*' + operator + '_' + strain + '_' + '*csv')
        print(r_file)
        for repeat, filename in enumerate(r_file):
            # read the csv file
            data = pd.read_csv(filename)
            # filter by side and front scattering
            data = data[(data['SSC-A'] > ssc_range[0]) & \
                    (data['FSC-A'] > fsc_range[0]) & \
                    (data['SSC-A'] < ssc_range[1]) & \
                    (data['FSC-A'] < fsc_range[1])]
            # filter by excluding the bottom 2.5 and top 97.5 percentiles
            percentile = np.percentile(data['FITC-A'], [2.5, 97.5])
            data = data[(data['FITC-A'] > percentile[0]) & \
                        (data['FITC-A'] < percentile[1])]

           # compute the mean and append it to the data frame along the
            # operator and strain
            df = df.append([[20160829, 'nbellive', operator, energies[i], 
                        strain, repressors[j], 0, data['FITC-A'].mean()]],
                        ignore_index=True)

# rename the columns of the data_frame
df.columns = ['date', 'username', 'operator', 'binding_energy', \
        'rbs', 'repressors', 'IPTG', 'mean_YFP']

# initialize pandas series to save the corrected YFP value
mean_bgcorr = pd.Series()
# correct for the autofluorescence background
for operator in df.operator.unique():
# substract the background grouping by strain
    mean_bgcorr = \
    mean_bgcorr.append(df.groupby('operator').get_group(operator).mean_YFP - \
    float(df[(df.rbs == 'auto') & \
    (df.operator == operator)].mean_YFP))

        
mean_bgcorr.name = 'mean_YFP_bgcorr'
df = pd.concat([df, mean_bgcorr], join_axes=[df.index], axis=1, join='inner')

# initialize pandas series to save the fold-change 
mean_fc = pd.Series()
# compute the fold-change
for operator in df.operator.unique():
# substract the background grouping by strain
    mean_fc = \
mean_fc.append(df.groupby('operator').get_group(operator).mean_YFP_bgcorr / \
        float(df[(df.rbs == 'delta') & \
    (df.operator == operator)].mean_YFP_bgcorr))

mean_fc.name = 'fold_change'
df = pd.concat([df, mean_fc], join_axes=[df.index], axis=1, join='inner')

# write
df.to_csv('output/20160829_lacI_titration_MACSQuant.csv', index=False)
#=============================================================================== 
# Add the comments to the header of the data file
filenames = ['./comments.txt', 'output/20160829_lacI_titration_MACSQuant.csv']
with open('../../../data/20160712_lacI_titration_MACSQuant.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())
