
import os
import numpy as np
import numpy.ma as ma
import pylab as plt


'''

'''

# Defire some parameters
nrow = 770
ncol = 405
nlay = 3
nts = 101
ncells = nrow*ncol*nlay

hed = np.empty([nrow*ncol*nlay*2, nts])

ifile = 'f:\projects\capture\codes\input\dataset_head.dat' 
ofile = 'f:\projects\capture\codes\output\dataset_head_mod.npy' 


# Define some option
opt_read_and_convert_dataset = True # Delete test lines
opt_get_dry_cell_id = False

if opt_read_and_convert_dataset: 
    count = 0
    count2 = 0
    ts = 0
    with open(ifile) as f:
        for line in f:
            val = line.split()            
            #print(f'{val} \n')
            #print(f'{val[0].isdigit()} \n')
            if val[0].isdigit():
                hed[count,ts] = float(val[0])
                count+=1
                count2+=1
                if count == ncol*nrow*nlay*2:
                    count=0
                    ts+=1
                
                #ofile.write(line)
    np.save(ofile, hed)



if opt_get_dry_cell_id:
    for ts in range(0,2,1):
        hed2 = hed[ncells:ncells*2,ts]
        hed3 = np.reshape(hed2,[nrow*ncol,nlay])
        hed_lay1 = hed3[:,0]
        hed_lay1_rs = np.reshape(hed_lay1, [nrow,ncol])
        plt.imshow(hed_lay1_rs)
        
        plt.show()
        #z2 = np.reshape(np.flipud(z), [nrow*ncol, 1])  # flipup for GMS dataset
