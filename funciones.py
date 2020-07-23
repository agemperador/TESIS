import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.patches as patches
import matplotlib.path as mpath
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import math as mt



def emagrama(time1,lat1,lon1,t1,td1,p1,u1,v1,guardar,nombre):

    latix,lonix=lat_idx(lat1),lon_idx(lon1)
    l=time1%24
    d=int(time1/24)+8
    tope=25

    vientox=u[time1,0:tope,latix,lonix].round(decimals=0)
    vientoy=v[time1,0:tope,latix,lonix].round(decimals=0)

    intensidad=mpcalc.wind_speed(vientox,vientoy)

    dirección=mpcalc.wind_direction(u[time1,:,latix,lonix],v[time1,:,latix,lonix])

    p1=p1/100
    pr = p1[0:tope] * units.hPa
    T = t1[time1,0:tope,latix,lonix]-273 
    temp=T.round(decimals=0) * units.degC
    Td = td1[time1,0:tope,latix,lonix] -273
    tempd=Td.round(decimals=0) * units.degC

    lcl_pressure, lcl_temperature = mpcalc.lcl(pr[0], temp[0], tempd[0])

    parcel_prof = mpcalc.parcel_profile(pr, temp[0], tempd[0]).to('degC')

    fig = plt.figure(figsize=(14,13))

    gs = gridspec.GridSpec(3, 3)

    skew = SkewT(fig,rotation=30)

    skew.plot(pr, temp, 'r', linewidth=2)
    skew.plot(pr, tempd, 'g', linewidth=2)
    skew.plot_barbs(pr,vientox, vientoy)
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
    skew.plot(pr, parcel_prof.round(decimals=1), 'k', linewidth=2)
    skew.shade_cin(pr, temp, parcel_prof.round(decimals=1))
    skew.shade_cape(pr, temp, parcel_prof.round(decimals=1))

    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
    skew.ax.set_ylim(1000, 30)
    skew.ax.set_xlim(-120, 80)
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    print(parcel_prof[parcel_prof>temp].shape,temp[parcel_prof>temp].shape,pr[parcel_prof>temp].shape)
    CAPE=-np.trapz(9.81*(parcel_prof[parcel_prof>temp ]-temp[parcel_prof>temp])/temp[parcel_prof>temp],pr[parcel_prof>temp])
    CIN=-np.trapz(9.81*(parcel_prof[parcel_prof<temp]-temp[parcel_prof<temp])/temp[parcel_prof<temp],pr[parcel_prof<temp])

    plt.text(-1,70,"CAPE = %s"%CAPE.round(decimals=1))
    plt.text(-3,65,"CIN = %s"%CIN.round(decimals=1))

    plt.title(nombre+"-lat=%s lon=%s - %s-%s utc" %(lat1,lon1,d,l),fontsize=20)

        
    ax_hod = inset_axes(skew.ax, '45%', '45%', loc=1)
    h = Hodograph(ax_hod, component_range=100.)
    h.add_grid(increment=20)
    h.plot_colormapped(vientox*2, vientoy*2,intensidad) 
    if guardar==True:
        plt.savefig("img/"+nombre+"-lat %s lon %s - %s-%s utc.png" %(lat1,lon1,d,l),format='png')
        plt.close()
    else: plt.show()

    

def dewpoint_approximation(T,RH):
 
    a = 17.271
    b = 237.7 # degC
    Td = (b * gamma(T,RH)) / (a - gamma(T,RH))
 
    return Td
 
def gamma(T,RH):
    a = 17.271
    b = 237.7 # degC
    g = (a * T / (b + T)) + np.log(RH/100.0)
 
    return g
def lon_idx(lng):
    idxlon=abs(lng+360-245)*4
    return(int(idxlon))

def lat_idx(ltd):
    idxlat=abs(20-ltd)*4
    return(int(idxlat))


