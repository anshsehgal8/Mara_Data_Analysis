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
        dM    = sigma * dA 


        GM    = 1.0                 # binary total mass
        omega    = (GM / r**3)**0.5    # Keplerian orbital frequency of the gas parcel
        L     = r * vp              # specific angular momentum and energy of the gas parcel
        E     = 0.5 * (vp**2 + vr**2) - GM / r
        e_squared = 1.0 - (0.5 * omega * L / E)**2
        e = np.sqrt(e_squared)
        i = np.where(~np.isnan(e))

        


        #radial_cut = (r > 1.5) * (r < 5.0)
        #avg_e = (e[i] * dA[i] * radial_cut[i]).sum() / (dA[i] * radial_cut[i]).sum()
        #print avg_e
        #print radial_cut.shape

        #plt.figure()
        #plt.hist(e_squared.flat,weights=radial_cut.flat,bins=100)
        #plt.show()




