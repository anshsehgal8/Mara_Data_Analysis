#!/usr/bin/env python3

import argparse
import os
import pickle
import numpy as np
import h5py
import loaders




def gas_eccentricity_diagnostics(filename):
    r     = loaders.get_dataset(filename, 'radius')
    dA    = loaders.get_dataset(filename, 'cell_area')
    vr    = loaders.get_dataset(filename, 'radial_velocity')
    vp    = loaders.get_dataset(filename, 'phi_velocity')
    vx    = loaders.get_dataset(filename, 'x_velocity')
    phi   = loaders.get_dataset(filename, 'phi')
    sigma = loaders.get_dataset(filename, 'sigma')
    dM    = dA * sigma

    domain_radius = h5py.File(filename, 'r')['run_config']['domain_radius'][()]
    rbins = np.linspace(0.0, domain_radius, 256)

    GM     = 1.0
    E      = +0.5 * (vp**2 + vr**2) - GM / r
    a      = -0.5 * GM / E
    omega2 = GM / a**3
    L      = r * vp
    e2     = 1.0 - (0.5 * L / E)**2 * omega2

    e2dM_binned, _ = np.histogram(r, weights=dM * e2, bins=rbins)
    dM_binned,   _ = np.histogram(r, weights=dM,      bins=rbins)
    dA_binned,   _ = np.histogram(r, weights=dA,      bins=rbins)

    vrm1dM, _ = np.histogram(r, weights=dM * vr * np.exp(1.j * phi), bins=rbins)
    vpdM,   _ = np.histogram(r, weights=dM * vp,                     bins=rbins)

    # For each checkpoint, we need to generate a file containing these measurements:
    # 1. Surface density radial profile
    # 2. MM08 as a function of r
    # 3. KD06 as a function of r
    # 4. MM08 globally
    # 5. KD06 globally
    # 6. Phase angle as a function of r
    # 7. Phase angle globally
    # 8. Filename (diagnostics file that was analyzed, and the time)

    radial_cut = (rbins[1:] > 1.5) * (rbins[1:] < 5.0)
    sigma_vs_r = dM_binned / dA_binned

    # MM08 measurements:
    mm08_vs_r         = vrm1dM / vpdM
    mm08_global       = vrm1dM[radial_cut].sum() / vpdM[radial_cut].sum()
    mm08_e_vs_r       = np.absolute(mm08_vs_r)
    mm08_phase_vs_r   = np.angle(mm08_vs_r)
    mm08_e_global     = np.absolute(mm08_global)
    mm08_phase_global = np.angle(mm08_global)

    # KD06 measurements:
    kd06_e_vs_r   = (e2dM_binned / dM_binned)**0.5
    kd06_e_global = (e2dM_binned[radial_cut].sum() / dM_binned[radial_cut].sum())**0.5

    
    path = os.getcwd()
    stringname = '/eccentricity_measures.%04d.h5'%(counter)
    newfilename = path + stringname
    #print newfilename
    newfilename = str(newfilename).zfill(4)
    h5f = h5py.File(newfilename, 'w')
    h5f['sigma_vs_r']        = sigma_vs_r
    h5f['mm08_vs_r']         = mm08_vs_r
    h5f['mm08_e_vs_r']       = mm08_e_vs_r
    h5f['mm08_phase_vs_r']   = mm08_phase_vs_r
    h5f['mm08_e_global']     = mm08_e_global
    h5f['mm08_phase_global'] = mm08_phase_global
    h5f['kd06_e_vs_r']       = kd06_e_vs_r
    h5f['kd06_e_global']     = kd06_e_global
    h5f['time']              = h5py.File(filename, 'r')['time'][()]
    h5f['filename']          = filename




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    args = parser.parse_args()
    counter = 0

    for filename in args.filenames:
        gas_eccentricity_diagnostics(filename)
        counter +=1
