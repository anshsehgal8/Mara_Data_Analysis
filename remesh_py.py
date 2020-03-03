#!/usr/bin/env python3

from itertools import chain
import h5py
import numpy as np
import matplotlib.pyplot as plt

def cell_center_arrays(vertices):
    x = vertices[:,:,0]
    y = vertices[:,:,1]
    xc = 0.5 * (x[:-1,:-1] + x[1:,:-1] + x[:-1,1:] + x[1:,1:])
    yc = 0.5 * (y[:-1,:-1] + y[1:,:-1] + y[:-1,1:] + y[1:,1:])
    return xc, yc

def area_array(vertices):
    x = vertices[:,:,0]
    y = vertices[:,:,1]
    dx = np.diff(x, axis=0)
    dy = np.diff(y, axis=1)
    dx_m = 0.5 * (dx[:,:-1] + dx[:,1:])
    dy_m = 0.5 * (dy[:-1,:] + dy[1:,:])
    return dx_m * dy_m


def radius_array(vertices):
    xc, yc = cell_center_arrays(vertices)
    return (xc**2 + yc**2)**0.5


def theta_array(vertices):
    xc, yc = cell_center_arrays(vertices)
    return np.arctan2(yc, xc)


def mass_array(sigma, vertices):
    return sigma * area_array(vertices)

fname = 'diagnostics.0300.h5'
h5f = h5py.File(fname, 'r')

sigma_list    = [g[...] for g in h5f['sigma'].values()]
vertices_list = [g[...] for g in h5f['vertices'].values()]
mass     = lambda: [mass_array(*i).flat  for i in zip(sigma_list, vertices_list)]
radius   = lambda: [radius_array(i).flat for i in vertices_list]
theta    = lambda: [theta_array (i).flat for i in vertices_list]


def accept(q):
    return q < 0.5 * np.pi  or q > 1.5 * np.pi 


filtered_radius = [m for m, q in zip(chain(*radius()), chain(*theta())) if accept(q)]
filtered_mass   = [m for m, q in zip(chain(*mass  ()), chain(*theta())) if accept(q)]
count, bins = np.histogram(filtered_radius, weights=filtered_mass, bins=250)
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.step(bins[:-1], count, where='post')
ax1.set_yscale('log')
plt.show()