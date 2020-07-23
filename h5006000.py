from netCDF4 import Dataset

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LightSource
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle, Polygon
from matplotlib.colors import LogNorm
#import matplotlib.colors as colors

#from scipy.interpolate import interp1d

from mpl_toolkits.basemap import Basemap

import numpy as np

import cartopy.crs as crs
import cartopy._crs as _ccrs
from cartopy.feature import NaturalEarthFeature

import wrf
from wrf import (to_np, getvar, smooth2d, get_cartopy,get_basemap, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

import warnings



def altura(array , altura = 500, it = 0):
    print('ll')
    i500 = 0 
    i6000 = 0
    h500 = np.zeros_like(array[it,0,...])
    h6000 = h500  
    for i in range(305):
        for j in range(305):
            hgt = array[it,:,j,i]

            hgtm = hgt/9.8

            z500diff = np.abs(hgtm - altura)
            mindiff = np.min(z500diff)
            lz500diff = list(z500diff)
            iz500 = lz500diff.index(mindiff)
            
            h500 [j,i] = iz500
            
        
    plt.imshow(h500)
    plt.colorbar()
    plt.show()
    
    return h500