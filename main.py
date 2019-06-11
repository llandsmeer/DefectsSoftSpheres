import numpy as np
raw = np.loadtxt('./sim/log')
nparts = (raw[:,0]==0).sum()
first = raw[raw[:,0]==0]

run = raw[:,0].astype(int)
p = raw[:,1:4]
c = raw[:,4:7]
n = raw[:,7:11]
