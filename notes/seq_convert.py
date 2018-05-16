import numpy as np
import pandas as pd

seqs = pd.read_csv('temp_seqs.csv')
seqs = np.array(seqs['seqs'])
seqs = [seq.upper() for seq in seqs]
np.savetxt('caps_seqs.csv', seqs, delimiter=',', fmt="%s")
