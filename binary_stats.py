#!/usr/bin/env python3




import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt




class HelperFunctions(object):
    """
    Contains functions to convert datasets in an HDF5 diagnostics file into more
    useful data.
    """
    def cell_center_x_array(self, vertices):
        x = vertices[:,:,0]
        xc = 0.25 * (x[:-1,:-1] + x[1:,:-1] + x[:-1,1:] + x[1:,1:])
        return xc

    def cell_center_y_array(self, vertices):
        y = vertices[:,:,1]
        yc = 0.25 * (y[:-1,:-1] + y[1:,:-1] + y[:-1,1:] + y[1:,1:])
        return yc

    def area_array(self, vertices):
        x = vertices[:,:,0]
        y = vertices[:,:,1]
        dx = np.diff(x, axis=0)
        dy = np.diff(y, axis=1)
        dx_m = 0.5 * (dx[:,:-1] + dx[:,1:])
        dy_m = 0.5 * (dy[:-1,:] + dy[1:,:])
        return dx_m * dy_m

    def radius_array(self, vertices):
        xc = self.cell_center_x_array(vertices)
        yc = self.cell_center_y_array(vertices)
        return (xc**2 + yc**2)**0.5

    def phi_array(self, vertices):
        xc = self.cell_center_x_array(vertices)
        yc = self.cell_center_y_array(vertices)
        return np.arctan2(yc, xc)




def get_dataset(fname, key):
    """
    Generates a dataset from an HDF5 diagnostics file, either by loading it
    directly, or performing a simple operation on the raw data.
    
    :type       fname:     str
    :param      fname:     The name of the HDF5 file to load (or generate) from
    :type       key:       str
    :param      key:       The name of the dataset to load (or generate). Can be
                           one of [sigma, radial_velocity, phi_velocity, radius,
                           phi, cell_area, mass].

    :returns:   The dataset
    :rtype:     a 3D numpy array with shape (num_blocks, block_size, block_size)
    
    :raises     KeyError:  If the key was not known
    """
    h5f = h5py.File(fname, 'r')
    helpers = HelperFunctions()

    if key in ['sigma', 'radial_velocity', 'phi_velocity']:
        return np.array([g[...] for g in h5f[key].values()])

    if key == 'radius':
        return np.array([helpers.radius_array(g[...]) for g in h5f['vertices'].values()])

    if key == 'phi':
        return np.array([helpers.phi_array(g[...]) for g in h5f['vertices'].values()])

    if key == 'cell_area':
        return np.array([helpers.area_array(g[...]) for g in h5f['vertices'].values()])

    if key == 'mass':
        return get_dataset(fname, 'sigma') * get_dataset(fname, 'cell_area')

    if key == 'x_velocity':
        vr = get_dataset(fname, 'radial_velocity')
        vp = get_dataset(fname, 'phi_velocity')
        phi = get_dataset(fname,'phi')
        return vr * np.cos(phi) - vp * np.sin(phi)

    if key == 'y_velocity':
        vr = get_dataset(fname, 'radial_velocity')
        vp = get_dataset(fname, 'phi_velocity')
        phi = get_dataset(fname,'phi')
        return vr * np.sin(phi) + vp * np.cos(phi)

    raise KeyError('unknown dataset: ' + key)




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    args = parser.parse_args()

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)

    for filename in args.filenames:
        domain_radius = h5py.File(filename, 'r')['run_config']['domain_radius'][()]
        vx = get_dataset(filename, 'x_velocity')
        dA = get_dataset(filename, 'cell_area')
        r  = get_dataset(filename, 'radius')
        cut = (r < 1.5) * (r > 0.5) # this cut should isolate the gap
        ax1.hist(vx.flat, weights=(dA * cut).flat, bins=200, histtype='step')

    plt.show()
