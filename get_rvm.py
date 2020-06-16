#!/usr/bin/env python3




import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
import loaders



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    args = parser.parse_args()

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    phase = []

    for filename in args.filenames:
        domain_radius = h5py.File(filename, 'r')['run_config']['domain_radius'][()]
        time          = h5py.File(filename, 'r')['time'][()]

        r     = loaders.get_dataset(filename, 'radius')
        dA    = loaders.get_dataset(filename, 'cell_area')
        vr    = loaders.get_dataset(filename, 'radial_velocity')
        vp    = loaders.get_dataset(filename, 'phi_velocity')
        vx    = loaders.get_dataset(filename, 'x_velocity')
        phi   = loaders.get_dataset(filename, 'phi')
        sigma = loaders.get_dataset(filename, 'sigma')



        #ECCENTRICITY CALCULATION MM08

        radial_cut = (r > 1.5) * (r < 5.0)

        m = 1
        top = np.sum((dA * sigma * vr * np.exp(1.j * m * phi)) * radial_cut)
        bottom = np.sum((dA * sigma * vp) * radial_cut)
        e_m = np.absolute(top/bottom)
        print e_m
        phase.append(np.angle(top/bottom))
       

        # ax1.hist(
        #     vx.flat,
        #     weights=(dA * (r < 3) * (r > 2)).flat,
        #     bins=100,
        #     density=True,
        #     histtype='step',
        #     label=r'$\rm{{orbit}} = {:.01f}$'.format(time / 2 / np.pi))
    
    dphi = np.diff(phase)
    dphi[dphi < -0.9 * 2 * np.pi] += 2 * np.pi
    dphi[dphi > +0.9 * 2 * np.pi] -= 2 * np.pi
    cumphase = np.cumsum(dphi)

    #for i in range(len(cumphase)):
        #print cumphase[i]

    

    #ax1.set_xlim(-0.2, 1.2)
    #ax1.set_xlabel(r'$e^2$')
    #ax1.set_ylabel(r'$P(e^2)$')
    #ax1.set_yscale('log')
    #ax1.legend()
    #plt.show()
