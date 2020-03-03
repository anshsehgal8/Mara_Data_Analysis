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
        xc = 0.5 * (x[:-1,:-1] + x[1:,:-1] + x[:-1,1:] + x[1:,1:])
        return xc

    def cell_center_y_array(self, vertices):
        y = vertices[:,:,1]
        yc = 0.5 * (y[:-1,:-1] + y[1:,:-1] + y[:-1,1:] + y[1:,1:])
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
    	return get_dataset(fname, 'radial_velocity') * np.cos(get_dataset(fname,'phi')) - np.sin(get_dataset(fname,'phi')) * get_dataset(fname,'phi_velocity')
    if key == 'y_velocity':
    	return get_dataset(fname, 'radial_velocity') * np.sin(get_dataset(fname,'phi')) + np.cos(get_dataset(fname,'phi')) * get_dataset(fname,'phi_velocity')

    raise KeyError('unknown dataset: ' + key)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument("inclination", type=float)
    parser.add_argument("varpi", type=float)
    args = parser.parse_args()

    #print(get_dataset(args.filenames[0], 'sigma').shape)
    #print(get_dataset(args.filenames[0], 'radius').shape)
    #print(get_dataset(args.filenames[0], 'phi').shape)
    radius = get_dataset(args.filenames[0],'radius')
    mass   = get_dataset(args.filenames[0],'mass')
    vphi   = get_dataset(args.filenames[0], 'phi_velocity')
    vr     = get_dataset(args.filenames[0], 'radial_velocity')
    area   = get_dataset(args.filenames[0], 'cell_area')

    def line_of_sight_vector():
    	ux = np.sin(args.inclination) * np.cos(args.varpi)
    	uy = np.sin(args.inclination) * np.sin(args.varpi)
    	uz = np.cos(args.inclination)
    	return [ux,uy,uz]

    

    
    

    Vx = get_dataset(args.filenames[0],'x_velocity')
    Vy = get_dataset(args.filenames[0],'y_velocity')
    Vz = 0

    #V = np.sqrt(Vx**2 + Vy**2)


    vel_array = np.array([Vx,Vy,Vz])
    #LOS = V * np.sin(1.05)
    V = np.dot(line_of_sight_vector(),vel_array)
    count,bins = np.histogram(V,weights=area,bins=1000)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(bins[:-1],count)

    #ax1.set_xlim(-0.4,0.4)
    plt.show()

    

