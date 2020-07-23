from netCDF4 import Dataset

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LightSource
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle, Polygon
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.colors as colors

#from scipy.interpolate import interp1d

from mpl_toolkits.basemap import Basemap

import numpy as np

import pandas as pd

import cartopy.crs as crs
import cartopy._crs as _ccrs
from cartopy.feature import NaturalEarthFeature

import wrf
from wrf import (to_np, getvar, smooth2d, get_cartopy,get_basemap, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

import time

import os

import glob


#from metpy.units import units
#from metpy.calc import temperature_from_potential_temperature

import warnings


warnings.filterwarnings("ignore")