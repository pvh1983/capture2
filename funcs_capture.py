import numpy as np
import os
import matplotlib.pyplot as plt
from funcs_capture import *
from shutil import copyfile


def plot_hist(data1, time_step):
    # print abs(data)
    plt.hist(abs(data), 10, density=True, facecolor='g', alpha=0.75)
    plt.xlabel('Withdrawal Rate (ft^3/day)')
    plt.ylabel('Probability')
    plt.title('Time step = %2d. Nsamples = %d \n.' %
              (time_step+1, data1.shape[0]))
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([8, 10, 0, 2])
    plt.grid(True)
    ofile2 = 'Qact_TS_' + str(time_step) + '_.png'
    plt.savefig(ofile2)
    # plt.show()


def plot_hist_err(data2):
    # print abs(data)
    plt.hist(abs(data2), 10, density=True, facecolor='g', alpha=0.75)
    plt.xlabel('Water Balance Errors (ft^3/day)')
    plt.ylabel('Probability')
    plt.title('Nsamples = %d \n.' % (data2.shape[0]))
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([8, 10, 0, 2])
    # plt.grid(True)
    plt.savefig('Err.png')
    # plt.show()

def print_2D_data():
    # Print x y and CF to a file
    for k2 in range(nlay):
        fid3 = open('CF.tmp', 'w')
        fid3.write('X, Y, Z \n')
        ofile_cf = 'CF2D_fID_' + \
            str(feature_id) + '_ts_' + str(ts + 1) + \
            '_lay_' + str(k2 + 1) + '.txt'
        for i2 in range(nrow):
            for j2 in range(ncol):
                # print k,i,j
                fid3.write('%9.1f, %9.1f, %4.3f \n' %
                           (x_coor[j2], y_coor[i2], c[k2, i2, j2]))
        fid3.close()
        copyfile('CF.tmp', ofile_cf)

    # Interpolation


def interpolation_cell(c, nlay, nrow, ncol):
    c2 = c
    d9 = np.empty([9, 1])
    for k2 in range(nlay):
        # print 'k2 = ', k2
        for i2 in range(1, nrow-1, 1):
            for j2 in range(1, ncol-1, 1):
                d9 = [c[k2, i2-1, j2-1], c[k2, i2-1, j2], c[k2, i2-1, j2+1],
                      c[k2, i2, j2-1], c[k2, i2, j2], c[k2, i2, j2+1],
                      c[k2, i2+1, j2-1], c[k2, i2+1, j2], c[k2, i2+1, j2+1]]
                # print 'd9 = ',d9
                d9n = [i for i in d9 if i >= 0.0]  # Find value -999
                if d9n:  # If d9 is NOT empty
                    # print feature_id,ts+1,k2,i2,j2,len(d9n),sum(d9n)/len(d9n)
                    c2[k2, i2, j2] = sum(d9n)/len(d9n)