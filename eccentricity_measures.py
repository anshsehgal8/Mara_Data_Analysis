#!/usr/bin/env python3

import argparse
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

    GM    = 1.0
    E     = +0.5 * (vp**2 + vr**2) - GM / r
    a     = -0.5 * GM / E
    omega = (GM / a**3)**0.5
    L     = r * vp
    e2    = 1.0 - (0.5 * omega * L / E)**2

    e2dM_binned, _ = np.histogram(r, weights=dM * e2, bins=rbins)
    dM_binned,   _ = np.histogram(r, weights=dM,      bins=rbins)

    vrm1dM, _ = np.histogram(r, weights=dM * vr * np.exp(1.j * phi), bins=rbins)
    vpdM,   _ = np.histogram(r, weights=dM * vp,                     bins=rbins)

    kd06_vs_r = (e2dM_binned / dM_binned)**0.5
    mm08_vs_r = vrm1dM / vpdM




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    args = parser.parse_args()

    for filename in args.filenames:
        gas_eccentricity_diagnostics(filename)
