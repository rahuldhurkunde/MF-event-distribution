import h5py
import sys
import numpy as np

no = sys.argv[1]

filename = 'noise_%s.hdf' %no
f = h5py.File(filename, 'r')

#Get the data

filen = 'noisefiles/noise_%s' %no

data = np.array(f['data'])
np.savetxt(filen, data)
