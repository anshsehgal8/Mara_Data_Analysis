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
    e_m = []

    for filename in args.filenames:
        domain_radius = h5py.File(filename, 'r')['run_config']['domain_radius'][()]
        time          = h5py.File(filename, 'r')['time'][()]
        mass_ratio    = h5py.File(filename, 'r')['run_config']['mass_ratio'][()]
        run           = mass_ratio * 1000

        r     = loaders.get_dataset(filename, 'radius')
        dA    = loaders.get_dataset(filename, 'cell_area')
        vr    = loaders.get_dataset(filename, 'radial_velocity')
        vp    = loaders.get_dataset(filename, 'phi_velocity')
        vx    = loaders.get_dataset(filename, 'x_velocity')
        phi   = loaders.get_dataset(filename, 'phi')
        sigma = loaders.get_dataset(filename, 'sigma')

        v_kep = np.sqrt(1./r)



        #ECCENTRICITY CALCULATION MM08

        radial_cut = (r > 1.5) * (r < 5.0)

        m = 1
        top = 2. * np.sum((dA * sigma * vr * np.exp(1.j * m * phi)) * radial_cut)
        bottom = np.sum((dA * sigma * v_kep) * radial_cut)
        e_m.append(np.absolute(top/bottom))
        #print e_m
        phase.append(np.angle(top/bottom))
    


    
    dphi = np.diff(phase)
    dphi[dphi < -0.9 * 2 * np.pi] += 2 * np.pi
    dphi[dphi > +0.9 * 2 * np.pi] -= 2 * np.pi
    cumphase = np.cumsum(dphi)
    np.savetxt('phasedata%iM_J.txt' %run,cumphase)
