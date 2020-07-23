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


def somb_cont(var_somb,var_cont,it = 0,milon = mlon,Malon = Mlon,Malat = Mlat,milat = mlat,
              lev_somb = 0,lev_cont = 0,cmap_somb = 'coolwarm',cmap_cont = "rainbow",nombre = ' ',
              viento = False, lev_viento = 0, sep_barb = 1, cmap_viento = 'nada',
              estados = True,guardar = False):
        
    
    ###REVISAR LA TD PORQUE A VECES DA MAYOR QUE T

    """
    #### --------------------------------------#######

    # SOMBREADO Y CONTORNO
    # somb_cont() grafica una variable sombreada y una contorno
    # En orden el imput es:
    # variable sombreada, meter directamente la variable entera
    # variable contorno, lo mismo
    # tiempo, en indices, si es "-1" toma todos los tiempos
    # longitud mínima, tomando como 0 el meridiano de greenwich
    # longitud máxima
    # latitud mínima tomando como 0 el eucador
    # latitud máxima 
    # nivel de altura del sombreado
    # nivel de altura del contorno
    # cmap del sombreado
    # cmap del contorno
    # nombre
    # True si se quiere el basemap de las provincias, False sino
    # True si se quiere guardar las imágenes
    # EJEMPLO:  somb_cont(t,q,-1,-70,-55,-40,-20,1,1,"rainbow","Blues","Prueba Q T",True,True)
    """
    
    
    print(type(it))
    if milon==999:
        mixlon=0
    else: mixlon=abs(milon+360-245)*4
    if  Malon==999:
        Mixlon=len(lon)
    else: Mixlon= abs(Malon+360-245)*4
    if milat==999:
        mixlat=0
    else: mixlat=abs(20-milat)*4
    if Malat==999:
        Mixlat=len(lat)
    else: Mixlat=abs(20-Malat)*4 

    if it==-1:
        it=ti    
    elif type(it)==int:
        it=[it]
        
    mxlon,Mxlon, mxlat, Mxlat = lon_idx(milon), lon_idx(Malon), lat_idx(milat), lat_idx(Malat)
    if len(var_somb.shape) == 4: var_somb =  var_somb[:,lev_somb,mxlat:Mxlat,mxlon:Mxlon]
    else: var_somb = var_somb[:,mxlat:Mxlat,mxlon:Mxlon]
    if len(var_cont.shape) == 4: var_cont = var_cont[:,lev_somb,mxlat:Mxlat,mxlon:Mxlon]
    else: var_cont =  var_cont[:,mxlat:Mxlat,mxlon:Mxlon]

    lons, lats = np.meshgrid(lon[mxlon:Mxlon], lat[mxlat:Mxlat])
    mlon,Mlon=np.min(lons),np.max(lons)
    mlat,Mlat=np.min(lats),np.max(lats)

    m = Basemap(projection='cyl', resolution="l",llcrnrlon=mlon,llcrnrlat=mlat, urcrnrlon=Mlon, urcrnrlat=Mlat)
    xi, yi = m(lons, lats)

    
    somb_min=np.min(var_somb[:,:,:])
    somb_max=np.max(var_somb[:,:,:])
    ncont=25
    clevs=np.linspace(somb_min,somb_max,ncont)

    cont_min=np.min(var_cont[:,:,:])
    cont_max=np.max(var_cont[:,:,:])
    clevc=np.linspace(cont_min,cont_max,ncont)
    

    for i in it:
        l=i%24
        d=int(i/24)+8

        fig=plt.figure(figsize=(10,10))
        
        s = m.contourf(xi,yi,var_somb[i,:,:],clevs,cmap=plt.get_cmap(cmap_somb),extend='both')
        cbar=m.colorbar(s,location='right',pad=0.05,size=0.2)
        cbar.ax.tick_params(labelsize=12,direction="out",length = 4,pad=3)

        c= m.contour(xi,yi,var_cont[i,:,:],clevc,cmap=plt.get_cmap(cmap_cont))
        clab = plt.clabel(c,clevc,fontsize=14,fmt='%.1f')

        if viento == True:
            a=sep_barb
            v_x = u[i,lev_viento,mxlat:Mxlat:a,mxlon:Mxlon:a]
            v_y = v[i,lev_viento,mxlat:Mxlat:a,mxlon:Mxlon:a]
            x_v = xi[::a,::a]
            y_v = yi[::a,::a]
            print(v_x.shape,v_y.shape,x_v.shape,y_v.shape)
            #vi=m.barbs(x_v, y_v, v_x, v_y, np.sqrt((v_x*2) ** 2 + (v_y*2) ** 2),
            #           cmap = cmap_viento if cmap_viento != 'nada' else None,
            #           length=6)
            vi=m.streamplot(x_v, y_v, v_x, v_y, cmap = plt.get_cmap('coolwarm'))
        
        m.drawparallels(np.arange(mlat, Mlat,5), labels=[0.3,0,0,0], fontsize=10,linewidth=0.4)
        m.drawmeridians(np.arange(mlon, Mlon,10), labels=[0,0,0,0.3], fontsize=10,linewidth=0.4)
        m.drawcoastlines(linewidth=0.4)
        m.drawcountries(linewidth=0.7)
        

        
        if estados==True:
            m.drawstates(linewidth=1)
        plt.title(nombre +" -2018/11/%s-%s UTC" %(d,l),fontsize=20, y=1.02,loc="center")
        
        dl=2.
        [ptlat,ptlon] = [-31.6-dl/2.,-64.0-dl/2.]
        
        m.plot([ptlon,ptlon],[ptlat,ptlat+dl],'-k',linewidth=2)
        m.plot([ptlon,ptlon+dl],[ptlat+dl,ptlat+dl],'-k',linewidth=2)
        m.plot([ptlon+dl,ptlon+dl],[ptlat+dl,ptlat],'-k',linewidth=2)
        m.plot([ptlon+dl,ptlon],[ptlat,ptlat],'-k',linewidth=2)
        if guardar==True:
            #plt.savefig('img/' +nombre+" -%s-%s utc" %(d,l))
            plt.savefig("../Latex/img/"+nombre.replace(' ','').replace('(','').replace(')','').replace('/','') +
                        "-%s-%sutc" %(d,l))
        else:  plt.show()




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


