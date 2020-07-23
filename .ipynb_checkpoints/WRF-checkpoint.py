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

import time

import warnings


warnings.filterwarnings("ignore")

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


def hovmoller(archivo, path,var_s = 'W',i= 150,j= 150, guardar = False, show = True):
    """
    Solo sirve para muchos tiempos
    """
    print(path+archivo)
    data = Dataset(path+archivo,'r')

    opress = data.variables['P']
    if opress.shape[0]==1: return 0
    opressb = data.variables['PB']

    ot = data.variables['T']
    otb = data.variables['T00']

    ow = data.variables[var_s]
    oh = data.variables['HGT']

    lats = data.variables['XLAT']
    lons = data.variables['XLONG']

    topo = np.asarray(data.variables['HGT'])
    lozada = [-31.651943,-64.07947] #Lat Lon de Lozada, Córdoba
    
    print('i = %f - j = %f '%(i,j))
    print(lats[0,i,j],lons[0,i,j])
    dz = opress.shape[1]
    press = opress[0,:,j,i] + opressb[0,:,j,i]

    t = ot[:,:,j,i] + 290.

    w = ow[:,1:,j,i]
    wu = ow.units
    dzt = t[:,1:dz]  - t[:,0:dz-1]


    time = np.arange(t.shape[0])
    time2D, press2D = np.meshgrid(time, press)

    pmin = press.min()
    pmax = press.max()

    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)
    ww = ax.pcolor(time2D, press2D, w.transpose(), cmap=plt.get_cmap('seismic'),zorder = -1,
                  vmin = -5,vmax = 5)
    cbar = plt.colorbar(ww)
    cbar.set_label(wu)

    var_c = dzt.transpose()

    ncont = 20

    c_min, c_max =np.min(var_c),np.max(var_c)

    if c_min == c_max:
        return 0

    clev=np.linspace(c_min,c_max,ncont)

    c = ax.contour(time2D[1:dz,:], press2D[1:dz,:], 
                   var_c ,clev, color='k' ,linestyles=np.where(clev >= 0, "-", "--"),
                  )
    clab = plt.clabel(c,clev,fontsize=12,fmt='%.0f',colors='k')


    if arch[32:-30]== '/test01/': carpeta = '/control/'
    else: carpeta = arch[32:-30]

    file = arch[-30:]
    plt.ylim(80000.,10000.)
    plt.xlabel('time')
    plt.ylabel('press')
    plt.title(arch[-30:]+'%s - %s'%(str(lats[0,i,j]),str(lons[0,i,j])))

    
    if show == True:
        plt.show()

    if guardar == True:
        plt.savefig('./img%sHovmoller/Hovmoller%s%s'%(carpeta,carpeta[1:-1], archivo[-30:]), format = 'png')
        plt.close(fig)

        #plt.show()
    data.close()

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def corte_vert(archivo,tiempos = [0], n_var_s = 'no', n_var_c='no',y_corte = 200, n_var_x = 'XLONG',
               ini = 0, end=305,
               set_topo = True,
               guardar = False, show = True, mask_s = 1, minimo = -15,maximo = 15, cmap ='coolwarm'):
    

    data = Dataset(archivo,'r')
    #opress = data.variables['P_PL']
    #ot = data.variables['T_PL']
    #ow = data.variables['']
    
    opress = data.variables['P']
    opressb = data.variables['PB']

    #ot = data.variables['T']
    #otb = data.variables['T00']


    if set_topo == True:
        oh = data.variables['HGT']
        topo = np.asarray(data.variables['HGT'])
        topografia = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,topo)))

    
    if n_var_s != 'no': ovars = data.variables[n_var_s]

    if n_var_c != 'no': ovarc = data.variables[n_var_c]

    lats = data.variables['XLAT']
    lons = data.variables['XLONG']

    x = data.variables[n_var_x]
    print (archivo)
    
    if type(tiempos)!=list: tiempos = [tiempos]
        
    
    
    for t in tiempos:
        
        
        
        time = data.variables['Times']
                
        if time.shape[0]<=t: 
            
            continue
    
        time = ''.join(str(time[t,...]).replace('b','').replace("'",'').replace('[','').replace(']','').replace('\n','').replace(' ',''))
    
        if time == '--------------------------------------':
            print('Tiempo invalido')
            return 0 
        
        if t >= lats.shape[0]:
            print('Se paso de tiempos')
            continue

        lozada = [-31.651943,-64.07947] #Lat Lon de Lozada, Córdoba

        i = y_corte

        press = np.array(opress[t,:,i,i]) + np.array(opressb[t,:,i,i])

        #w = ot[0,:,i,:] +290
    
        #np.max(ow)
        xtrm = np.max([minimo, maximo])
        
        if n_var_s=='W':
            var_s = ovars[t,:-1,i,:]
        else: var_s = ovars[t,:,i,:]
            
        if n_var_c == 'W':
            var_c = ovarc [t,:-1,i,:]
        else: var_c = ovarc [t,:,i,:]
            
            

        x2D, press2D = np.meshgrid(x[t,i,:], press)


        fig = plt.figure(figsize=(15,7))

        ax = fig.add_subplot(111)

        colmiss='#AAAAAA'
        wcmap = plt.get_cmap(cmap)
        #wcmap.set_bad(color=colmiss)
        #wnorm = Normalize(-xtrm,xtrm)
        for j in range(305):
            
            if set_topo ==True:
                top = topografia[t,i,j]
            else: top = 100000
                

            ##### HAY UN PROBLEMA CON  LA MASCARA
            var_s[:,j] = np.ma.array(var_s[:,j], mask = press2D[:,j] > top)
            var_c[:,j] = np.ma.array(var_c[:,j], mask = press2D[:,j] > top)

        var_s = np.ma.array(var_s , mask = np.abs(var_s) < mask_s)
        #v = np.ma.array(w , mask = np.abs(w) < 2)
        s = ax.pcolormesh(x2D[:,ini:end],press2D[:,ini:end],var_s[:,ini:end], cmap = wcmap,vmax =maximo,vmin =minimo)
        
        ncont = 5
        c_min, c_max =np.round(np.min(var_c)),np.round(np.max(var_c))
        if c_min != c_max:
            clev=np.linspace(c_min,c_max,ncont)
            clev= np.delete(clev,np.where(clev==0))
            pc = plt.contour(x2D[:,ini:end],press2D[:,ini:end],var_c[:,ini:end],clev, colors = 'k', )

            clab = plt.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')
        
        if set_topo == True: plt.plot(x2D[0,ini:end],topografia[0,i,ini:end],linewidth = 5, color = 'k' )

        
        #plt.vlines(lozada[1],ymin=97500.,ymax=1000)

        plt.colorbar(s,extend ='both')

        plt.ylim (96600.,1000)
        plt.ylabel('Presión (hPa)')

        plt.xlabel('Lon °')
        plt.title('Corte vertical longitudinal - lat = %f'%lats[0,i,i])
        plt.suptitle(time)
        carpeta = archivo[32:-30]
        if carpeta == '/test01/': carpeta='/control/'

        if show == True: plt.show()

        if guardar == True:
            plt.savefig('./img%sZsec/CorteVert_%s_%s_%s_%i' %(n_var_s,carpeta,carpeta[1:-1],time,t//4),format='png')
            print('Imagen guardada %i'%t)
            plt.close(fig)
        
    data.close()

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def carga_var(data,name_var,time,lev_var,ini = 0,end = 305):
    """
    SOLO SIRVE CON WRF
    
    Función para cargar las variables comunes.
    
    No carga viento ni lluvia
    
    Se ingresa el archivo netcdf 'data', el nombre de la variable 'name_var',
    el tiempo 'time', nivel 'lev_var', lat lon inicial y final 'ini' 'end'
    """
    #### Si hay nombre de la variable
    if name_var != '0':
        
        if (type (name_var)== str) and (type(time) == int):
            #cargo la variable del wrf
            var = getvar(data,name_var, timeidx = time)
                
            #si es una variable 4d me quedo solo con el nivel que quiero
            if len (var.shape)>2: var = var[lev_var,...]
            smooth_var = smooth2d(var, 3, cenweight=4)
        else:
            print("""\033[1;32m La variable tiene que estar escrita como string y el tiempo como int \n 
                  Si no escribis ts por default es 0 \033[0m""")
            return 0
        
    else:
        smooth_var = 0
    
    return smooth_var
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def carga_v_rain(data, i,lev_viento = 0, lev_rich = 15, modo = 'normal',
          viento = False,sep_barb = 2, lluvia = False, topografia = False,
          ini =0,end = 305 , var_lluvia = 'RAINNC'): 
    """
    SOLO SIRVE CON WRF
    
    Esta función carga las variables especiales como viento, lluvia y los indices:
    
    necesita las cosas tipicas: data, tiempo, nivel.
    'modo' me dice si es en nivel eta(normal), presión o indices.
    viento, lluvia y topografia son variables bool
    'var_lluvia' se puede cambiar si se desea otro tipo de lluvia
    
    Carga las variables segun corresponde:
    
    Al viento lo carga segun lo que se pida
    
    """
    
    
    """""""""""""""""""""""""""""""""VIENTO"""""""""""""""""""""""""""""""""
    # En modo presion cargo los vientos en niveles de presion
    if modo == 'presion': var_viento = ['U_PL', 'V_PL']
    #si es nivel 0 tomo viento a 10 m
    elif lev_viento ==0: var_viento = ['U10', 'V10']
    #sino tomo la variable en nivel eta
    else: var_viento = ['U', 'V']
    if viento == True:
        
 
        v_x = getvar(data, var_viento[0], timeidx= i ) 
        v_y = getvar(data, var_viento[1], timeidx= i) 
        
        #me quedo solo con un nivel
        if len(np.asarray(v_x).shape)==2: 
            v_x = v_x[ini:end,ini:end]
            v_y = v_y[ini:end,ini:end]
        else: 
            v_x = v_x[lev_viento, ini:end,ini:end]
            v_y = v_y[lev_viento, ini:end,ini:end]
        
        v = [v_x,v_y]
        
    else: v = 0
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""  
    """""""""""""""""""""""""""""LLUVIA"""""""""""""""""""""""""""""
    if lluvia == True:
        raini = getvar(data,var_lluvia,timeidx = i)[ini:end,ini:end]
        rain_ = getvar(data,var_lluvia,timeidx = i-1)[ini:end,ini:end]

        rain = raini-rain_
        
    else: rain = 0

    if modo == 'indices':
        #indices es cape y richardson
        var_viento = ['U', 'V']
        
        v_top = 17 #Altura de aprox 300hpa donde h = 6km
        v_bot = 1 # Altura de aprox 950hpa donde h = 500m
        
        #Cargo los cuatro vientos para calcular richardson
        v_x_top = getvar(data, var_viento[0], timeidx= i ) [v_top, ini:end, ini:end]
        v_y_top = getvar(data, var_viento[1], timeidx= i ) [v_top, ini:end, ini:end]
        v_x_bot = getvar(data, var_viento[0], timeidx= i ) [v_bot, ini:end, ini:end]
        v_y_bot = getvar(data, var_viento[1], timeidx= i ) [v_bot, ini:end, ini:end]
        
        v_x_bot,v_y_bot, v_x_top,v_y_top = to_np(v_x_bot), to_np(v_y_bot),to_np(v_x_top), to_np(v_y_top)
        
        #Calculo la cortante
        shear = np.asarray(list(map( lambda xtop,xbot,ytop,ybot: ((xtop-xbot)**2+(ytop-ybot)**2)/2., v_x_top,v_x_bot, v_y_top, v_y_bot )))
        #Cargo el cape
        cape = np.array(getvar(data,'AFWA_CAPE', timeidx = i))[ini:end,ini:end]
        #cape = np.expand_dims(cape,axis = 0)
        #cape = np.repeat(cape,59,axis = 0) sirve para repetir la matriz en un tercer eje
        #v_media = v_x**2+v_y**2
        
        #Calculo el richardson
        ri = cape/shear

        print('Richardson activado')
    else: ri,cape, shear =0,0,0
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    return to_np(v),rain, ri, cape, shear

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def set_size(w,h, ax=None):
    """ 
    Seteo el tamalo del gráfico para que quede bien
    w, h: width, height in inches 
    """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w+0.2)/(r-l)
    figh = float(h+0.2)/(t-b)
    ax.figure.set_size_inches(figw, figh)
    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def sombreado (x,y, var_somb,ax, ini = 0, end = 305, cmap_s = 'coolwarm', zorder = 0, k=1, vmin = 0, vmax = 0):
    """
    Gráfico de la variable sombreada
    """    
    #var_s = set_var(var_somb)
    
    if vmin == vmax: 
        vmin = np.min(var_somb)
        vmax = np.max(var_somb)
    
    ps = ax.pcolormesh(x, y, var_somb[ini:end,ini:end]*k, cmap = plt.get_cmap(cmap_s,10),zorder = zorder, vmin = vmin, vmax = vmax )

    return ps

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def contorno (x,y, var_cont,ax, ini = 0, end = 305, cmap_c = 'hot_r', zorder = 1,ncont = 12, n_var_c = 'no'):
    """
    Gráfico de la variable de contorno
    """
    #var_c = set_var(var_cont)
    var_c = var_cont
    c_min, c_max =np.min(var_c),np.max(var_c)
    clev=np.linspace(c_min,c_max,ncont)
    if n_var_c == 'UP_HELI_MAX': 
        clev= np.delete(clev,np.where(clev==0))

    if c_min >= c_max: clev = [0.0]
    pc = ax.contour(x,y, var_c[ini:end,ini:end],clev, cmap = plt.get_cmap(cmap_c),zorder = zorder)
    clab = plt.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')

    return pc
 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_viento(x,y,v,bm,ini=0,end=305,cmap_viento = 'no',zorder = 2, density = 1):
    """
    Gráfico de viento 
    """ 
    ### Cmap personalizado para el viento con jet a partir de los 13m/s
    """""""""""""""CMAP VIENTO"""""""""""""""
    if cmap_viento == 'no':
        cdict = {
        'blue': [(0,1,0.5),(0.4,0.5,0.4) , (1, 1, 1)],
        'green': [(0,0,0.2),(0.4,0.2,0), (1, 0, 0)],
        'red': [(0,0,0) ,(0.4,0,0.4), (1, 1, 1)]}
        cmap_viento = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,10)
        k = 3

    else: k = 1
    """"""""""""""""""""""""""""""""""""""""""
    x_v = x[0,ini:end]
    y_v = y[ini:end,0]
    

    v_x = to_np(v[0])[ini:end,ini:end]
    v_y = to_np(v[1])[ini:end,ini:end]
    
    #Calculo la velocidad
    speed = np.sqrt((v_x*2) ** 2 + (v_y*2) ** 2)#[ini:end,ini:end]
    #Esto es para que varie el ancho de las lineas de corriente
    
    lw = k*speed / speed.max()
    
    # Esto verifica que ninguna variable tenga dimension 0
    if len(x_v)<=0 or len(y_v)<=0 or len(v_x)<=0 or len(v_y)<=0 or len(speed)<=0:
        print('Esta todo mal con esto', x_v.shape,
          y_v.shape, v_x.shape, v_y.shape, speed.shape, lw.shape)
        return 0
    #Si la velocidad máxima es menor a este numero no grafico las lineas de corriente
    if speed.max()>0.001:
        vi= bm.streamplot(x_v,y_v,  v_x, v_y, color = speed , density=density
                     ,cmap = cmap_viento, linewidth= lw, arrowsize = 1.5,
                         arrowstyle = '->', zorder = zorder)
    else: 
        print('No hay viento')
        vi = None
    return vi

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_precip(x,y, rain, bm, ini=0,end = 305, cmap_lluvia='rainbow', zorder = 3, vminr = 0, vmaxr = 0):
    """
    Gráfico de precipitación
    """
    #print(vminr, vmaxr)
    if vminr == vmaxr:
        vminr = np.min(rain)
        vmaxr = np.max(rain)
    #Enmascaro la lluvia menor a 0.5
    r = np.ma.array(rain, mask= rain < 0.5)
    ra = bm.pcolormesh(x,y,r[ini:end,ini:end], cmap = plt.get_cmap(cmap_lluvia,10),zorder = zorder, vmin = vminr, vmax = vmaxr)
    return ra    

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_recuadro(x,ax = False, ancho= 0.1,alto=0.1):
    """
    EN CONSTRUCCIÓN
    
    Esto va a graficar un recuadro en una región
    """
    lats = [-34,-32,-32,-34]
    lons = [-65,-65,-63,-63]
    x,y = ax(lats,lons)
    xy = zip(x,y)
    currentAxis = plt.gca()
    poly = Polygon( list(xy), facecolor='red', alpha=0.4 )
    currentAxis.add_patch(poly)
    #plt.gca().add_patch(poly)
    #else:   currentAxis = ax
    #currentAxis.add_patch(Rectangle((  0.5 - ancho/2, 0.5 - alto/2), ancho, alto, fill = False, alpha=1,linewidth=3))
    
    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def info(var = '0',lev = '0', viento = False, lluvia = False, topografia=False):
    """
    Esta función genera un texto con todo lo que se va a graficar
    """
    string = ''
    if var != '0':
        string = string + var + str(lev) 
    
    elif (viento and lluvia and topografia): string = 'No'
    else:
        if viento: 
            string = string + 'Viento en ' + str (lev)
        if lluvia: 
            string = string + ', Lluvia activada '
        if topografia:
            string = string + ', Topografía activada'

    return string

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def mapa_somb_cont2(path = '../../../../media/agustin/Linux/salidas_wrf/test01/',
               archivo = 'wrfout_d02_2018-11-10_12:00:00', modo = 'normal',
               ini = 0, end = 100,rows = 2,cols = 2,
               var_s = '0',  ts = [0],   lev_s = 0  , us = 0 , cmap_s = 'coolwarm', vmin = 0,vmax = 0,mask = 'no',
               var_c1 = '0', lev_c1 = 0 , uc1 = 0, cmap_c1 = 'hot_r', 
               var_c2 = '0', lev_c2 = 0 , uc2 = 0, cmap_2 = 0, cmap_c2 = 'winter',
               viento = True, lev_viento = 0,sep_barb = 1,density = 1, cmap_viento = 'no',
               titulo = 'Sin titulo definido', guardar = False, savepath = './img/', show = True,
               topografia = True, topo_mode = 'soft',
               lluvia = True, cmap_lluvia = 'gist_ncar', vminr = 0, vmaxr = 0,
               lev_rich = 15, lluvia_inicial = -1
                   ):   
    """
    FUNCIÓN PRINCIPAL DE GRAFICOS SOMBREADOS Y 2 CONTORNOS CON LLUVIA Y VIENTO
    
    Se le pueden dar parametros de lo que se te ocurra
    Se pueden hacer gráficos en niveles eta y de presión.
    Se pueden hacer gráfico de variables comunes o de variables especiales como los indices
    Se puede graficar viento y lluvia    
    
    ##############################
    #-----EJEMPLO DE MAPA--------#
    ##############################
    
    mapa_somb_cont2(path= PATH_DEL_ARCHIVO,archivo=STRING_DEL_ARCHIVO, var_s='Q2',vmin = 0.008,vmax = 0.022,
            var_c1='TH2',
            guardar=False,savepath = STRING_DEL_PATH_DE_SALIDA, show = True,ini = 0, end = 304,rows = 1, cols = 1,
            ts = [t], lluvia= True, topografia=False, modo = 'normal', #lev_rich=k,
            titulo="Q 2m (s) - $\Theta$ 2m (c) -  V sup (slines) - Lluvia (mm/hr)",
            cmap_s = 'BrBG',cmap_c1='hot_r', viento = True,
            cmap_lluvia = 'rainbow')
    """
    
    # Si le dan un solo tiempo lo vuelve lista
    if type(ts) != list: ts = [ts]
    """CARGA DE DATOS"""
    #Archivo ncdf
    data = Dataset(path+archivo,'r')
    
    #Lat y lon
    lats, lons = getvar(data,'XLAT'), getvar(data,'XLONG')
    if topografia==True: 
        #Cargo la topografia
        topo = np.asarray(getvar(data,'HGT'))
        
    #Estos son los tiempos escritos en palabras
    times =  data.variables['Times']
    ###Un poco de control###

    #Si el ultimo de los tiempos supera el máximo tiempo, esta todo mal y corta
    if ts[-1]>np.asarray(times).shape[0]: 
        print('Hay un tiempo que esta fuera de rango')
        return 0
    
    if np.asarray(lats).shape[0] <= end or np.asarray(lons).shape[0] <= end:
        end = min(np.asarray(lats).shape[0],np.asarray(lons).shape[0])-1
        
    if np.any(np.asarray(ts)>=len(times)):
        for i in range (len(ts)):
            ts[i] = i*len(times)//len(ts) -1
    #########################
    tiempos = []
    #Guardo los tiempos en texto
    for j,i in enumerate(ts):    
        tiempos.append( ''.join(str(times[i]).replace('b','').replace("'",'').replace('[','').replace(']','').replace('\n','').replace(' ','')))
        #print(tiempos[-1])
    #Pasan estas cosas
    if tiempos ==['--------------------------------------']:
        print('Problemas con los tiempos')
        return 0

    # Descripción corta de lo que se está graficando
    if modo == 'normal':
        print ('Mapa con límites lat: %f - %f, lon: %f - %f. \n Variables: %s somb, %s cont 1, %s cont 2,  %s. \n Archivo: %s. \n Tiempos: '
          %(lats[ini,ini], lats[end,end], lons[ini,ini], lons[end,end], 
            info(var_s,lev_s), info(var_c1,lev_c1),info(var_c2,lev_c2),
            info(viento = viento, lev = lev_viento, lluvia = lluvia, topografia = topografia),archivo), tiempos)
    else: 
        print('Modo %s'%modo)
    
    # Cargo el basemap de las latitudes
    bm = get_basemap(lats[ini:end,ini:end])
        
    x, y = bm(to_np(lons)[ini:end,ini:end], to_np(lats)[ini:end,ini:end])
    

    #Por si le pifia a los tiempos en función de las columnas
    if len(ts)!= rows*cols:
        if np.sqrt(len(ts)) == round(np.sqrt(len(ts))):
            rows = int(np.sqrt(len(ts)))
            cols = rows
        else: 
            rows = int(np.sqrt(len(ts)))+1
            cols = rows
            
    zorderR = 3
                
    #Guardo los nombres de las variables
    n_var_s, n_var_c1, n_var_c2 = var_s, var_c1, var_c2
    
    if modo == 'indices': 
        n_var_s = 'AFWA_CAPE'
        n_var_c2 = 'Richarson'
        zorderR = 0
        if vmax<vmin:
            vmin = np.min(to_np(data.variables[n_var_s]))
            vmax = np.max(to_np(data.variables[n_var_s]))
            print('Valores máximos de CAPE %f -- Valores mínimos de CAPE %f'%(vmax, vmin))
    
    """
    if modo == 'normal' or modo == 'presion':
        if vmax<vmin:
            _ = to_np(data.variables[n_var_s])
            if len(_) == 4: _ = _[:,lev_s,:,:]
            vmax = np.max(_)
            vmin = np.min(_)
    """
    lozada = [-31.651943,-64.07947] #Lat Lon de Lozada, Córdoba
    
    lats, lons = np.array(lats),np.array(lons)
    
    #print('Lozada: ',lats[lats == lozada[0]], lons[lons == lozada[1]])

    ###############################
    ######  GENERO LA IMAGEN ######
    ###############################
    
    fig = plt.figure()#,constrained_layout=True)
    axx = []
    for i, j in zip(ts,range(len(ts))):
        
        ## Un subplot por cada mapa uso las rows cols y la enumeración j
        ax = fig.add_subplot(int(rows*100+cols*10 + j+1))
        
        #Aca guardo los subplots
        axx.append(ax)
        
        """""""""""""""FUNCION DE CARGA DE DATOS"""""""""""""""
        
        ## Esta carga las cosas raras
        v,rain, ri, cape,shear = carga_v_rain(data = data, i = i,lev_viento = lev_viento, lev_rich=lev_rich, modo = modo,
                              viento = viento, sep_barb = sep_barb, lluvia = lluvia, ini =ini,end = end)
        
        if type(lluvia_inicial) != int:
            rain = lluvia_inicial
            
            
        ## Esta carga las cosas comunes
        if modo == 'normal' or modo =='presion': 
            smooth_var_s  = carga_var(data,n_var_s,i ,lev_s,ini,end)
            smooth_var_c2 = carga_var(data,n_var_c2,i, lev_c2,ini,end)
        
        elif modo == 'indices':
            smooth_var_s = cape
            smooth_var_c2 = ri
            
        if mask != 'no': smooth_var_s = np.ma.array(smooth_var_s, mask = smooth_var_s < mask)

            
        smooth_var_c1 = carga_var(data,n_var_c1,i,lev_c1,ini,end)
        

        ## Me genero los bools para despues ver que grafico y que no
        if type(smooth_var_s) == int: somb = False
        else: somb = True
        if type(smooth_var_c1) == int: cont1 = False
        else: cont1 = True
        if type(smooth_var_c2) == int: cont2 = False
        else: cont2 = True
        
        

        """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        """""""""ACA EMPIEZA EL MAPA"""""""""
        # DEFINO DIA Y HORA
        l=i%24
        d=int(i/24)+12 


        """""""""""""""GRAFICO SOMBREADO"""""""""""""""
        if somb == True or modo == 'indices':
           
            ps = sombreado(x,y, smooth_var_s,ax,ini =ini, end = end,cmap_s= cmap_s, vmin = vmin, vmax = vmax)

            ###TOPOGRAFÍA EN CONTORNO###
        if topografia == True: 
            
            ax.contour(x,y, topo[:-1,:-1],cmap = 'gray_r', linewidths = 2, zorder = 0)
                
        """""""""""""""""""""GRAFICO CONTORNO 1"""""""""""""""""""""
        if cont1 == True:
            
            pc1 = contorno(x,y, smooth_var_c1,ax, ini =ini, end = end, cmap_c = cmap_c1, zorder = 1,ncont = 12, n_var_c = n_var_c1)
            
        """""""""""""""""""""GRAFICO CONTORNO 2"""""""""""""""""""""
        if cont2 == True or modo =='indices':

            pc2 = contorno(x,y, smooth_var_c2,ax, ini =ini, end = end, cmap_c = cmap_c2, zorder = 1,ncont = 5)

        """""""""""""""""""""""""""GRAFICO VIENTO"""""""""""""""""""""""""""
        if viento == True:
            
            vi = plot_viento(x,y,v,bm,ini,end,cmap_viento=cmap_viento,zorder= 4, density = density)
            
            
        """""""""""""""""""""""""""""LLUVIA"""""""""""""""""""""""""""""
        if lluvia == True:
            
            ra = plot_precip(x,y, rain, ax, ini,end, cmap_lluvia, zorder = zorderR, vminr = vminr, vmaxr = vmaxr)


        """""""""""""""BASEMAP"""""""""""""""""""""""""""""""""""""""""""""
        bm.drawcoastlines(linewidth=0.25)
        bm.drawstates(linewidth=1)
        bm.drawcountries(linewidth=0.25)
        """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        """
        #plot_recuadro(lozada,ax =bm)

        #if recuadro != [-1,-1]:
            #plot_recuadro(recuadro,)

        ax.add_patch(Rectangle((lats[100,100],lons[100,100] ), 1, 1, fill = False, alpha=1))
       """
            
        ## Muestro los tiempos de cada plot para ver donde se corta por cualquier error
        #print(tiempos[j] + ' OK')
        
        ## Titulo de cada gráfico
        ax.set_title("\n %s" %tiempos[j], fontsize = 12)

        titulo_ = titulo.replace(' ','').replace('(','').replace(')','').replace('/','')
        

    ##################################
    #Esto es para que quede acomodado#
    ##################################
    #fig.tight_layout()
    set_size(12,8)
    ##################################
    
    ##################################
    ###"""ACA VAN LAS COLORBARS """###
    ##################################
    if somb == True: 
        caxv = plt.axes([0.8, 0.15, 0.015, 0.7])

        cbars = fig.colorbar(ps, ax=axx,cax=caxv, orientation = 'vertical')
        #Defino los limites del gráfico
        #ps.set_clim(vmin, vmax)
        cbars.ax.set_ylabel(n_var_s + ' (%s)' %data.variables[n_var_s].units,fontsize = 12)
        #cbars.draw_all()
        
    if lluvia == True:
        caxh2 = plt.axes([0.2, 0.15, 0.01, 0.7])
        cbar = fig.colorbar(ra, ax = axx, cax = caxh2 ,orientation = 'vertical')
        #cbar.ax.set_yticklabels([vminr,10,20,30,40,vmaxr])
        #Defino los limites del gráfico
        #ra.set_clim(0, 40)
        #cbar.draw_all()
        cbar.ax.set_ylabel('Lluvia (%s)' %data.variables['RAINNC'].units, labelpad = -70 , fontsize = 12)
        
    if viento == True and vi is not None:
        caxh = plt.axes([0.2,0.08, 0.6, 0.02])
        cbarv = fig.colorbar(vi.lines,ax=axx,cax=caxh, orientation = 'horizontal')    
        minvalue = 0
        maxvalue = 30
        if lev_viento != 0:
            minvalue = 0
            maxvalue = 80
        #Defino los limites del gráfico
        cbarv.set_clim(minvalue, maxvalue)
        cbarv.ax.set_xlabel('Velocidad del viento (%s)' %data.variables['U10'].units, fontsize = 12)
        cbarv.draw_all()
        

    ###################################
    
    ## El titulo de toda la imagen
    fig.suptitle(titulo+ '\n',fontsize = 18, y = 0.95)
    
    if guardar==True: 
        plt.savefig("%s%s-%s-utc.png" %(savepath,titulo_,tiempos[j]))
        print('Imagen guardada en %s'%savepath)

    if show == True: plt.show()
    data.close()
    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
class mapa ():
    
    def __init__(self, archivo,archivo_prev,wb_v,wb_c,wb_r,wb_i,wb_s):
        super().__init__()
        
        self.data = Dataset(archivo,'r')
        self.data_wb_v = Dataset(wb_v,'r')
        self.data_wb_c = Dataset(wb_c,'r')
        self.data_wb_r = Dataset(wb_r,'r')
        self.data_wb_i = Dataset(wb_i,'r')
        self.data_wb_s = Dataset(wb_s,'r')
        self.archivo_prev = archivo_prev

        

        
        
    def cargar_var(self):
        lat,lon = getvar(self.data,'XLAT'), getvar(self.data,'XLONG')
        self.bm = get_basemap(lat)
        self.lat,self.lon = to_np(lat),to_np(lon)
        topo = np.array(getvar(self.data,'HGT'))
        self.topo = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,topo)))
        
        self.times = self.data.variables['Times'][:]


        self.x, self.y  = self.bm(self.lon,self.lat)
        
        rain =  np.array(self.data.variables['RAINNC'])
        
        data_prev = Dataset(self.archivo_prev,'r')
        rain[1:] = rain[1:,...]-rain[:-1,...]
        rain[0] = rain[0,...]-np.array(data_prev.variables['RAINNC'][-1,...])
        
        data_prev.close()

        self.t = np.array(self.data.variables['T'])

        self.rain = np.ma.masked_array(rain, rain <2)
        
        self.u,self.v,self.w = np.array(self.data.variables['U']),np.array(self.data.variables['V']), np.array(self.data.variables['W'])
        
        pres = np.array(self.data.variables['P']) 
        presb = np.array(self.data.variables['PB'])
        self.press = (pres + presb)/100
        self.pmax, self.pmin = 966,100
        
        self.window = 20
        self.stride = 11
        self.stride_lon = 16
        
        
        self.wcmap = plt.get_cmap('seismic')
        


        return
    

    


    def cargar_var_especies(self, i,k,lat, lon):


        self.qice = np.array(self.data.variables['QICE'][i,:,lat,lon])
        self.qs =  np.array(self.data.variables['QSNOW'][i,:,lat,lon])
        self.qc =  np.array(self.data.variables['QCLOUD'][i,:,lat,lon])
        self.qr =  np.array(self.data.variables['QRAIN'][i,:,lat,lon])
        self.qv =  np.array(self.data.variables['QVAPOR'][i,:,lat,lon])

        self.qvtend = np.array(self.data_wb_v.variables['QVTTEND'][k,:,lat,lon])
        self.qctend = np.array(self.data_wb_c.variables['QCTTEND'][k,:,lat,lon])
        self.qitend = np.array(self.data_wb_i.variables['QITTEND'][k,:,lat,lon])
        self.qrtend = np.array(self.data_wb_r.variables['QRTTEND'][k,:,lat,lon])
        self.qstend = np.array(self.data_wb_s.variables['QSTTEND'][k,:,lat,lon])

        self.qvh_adv = np.array(self.data_wb_v.variables['QVH_ADV'][k,:,lat,lon])
        self.qch_adv = np.array(self.data_wb_c.variables['QCH_ADV'][k,:,lat,lon])
        self.qih_adv = np.array(self.data_wb_i.variables['QIH_ADV'][k,:,lat,lon])
        self.qrh_adv = np.array(self.data_wb_r.variables['QRH_ADV'][k,:,lat,lon])
        self.qsh_adv = np.array(self.data_wb_s.variables['QSH_ADV'][k,:,lat,lon])

        self.qvz_adv = np.array(self.data_wb_v.variables['QV_ZADV'][k,:,lat,lon])
        self.qcz_adv = np.array(self.data_wb_c.variables['QC_ZADV'][k,:,lat,lon])
        self.qiz_adv = np.array(self.data_wb_i.variables['QI_ZADV'][k,:,lat,lon])
        self.qrz_adv = np.array(self.data_wb_r.variables['QR_ZADV'][k,:,lat,lon])
        self.qsz_adv = np.array(self.data_wb_s.variables['QS_ZADV'][k,:,lat,lon])


    def carga_sing_var(self,data,nombreVariable,i,lev,lat,lon):

        return np.array(data.variables[nombreVariable][i,lev,lat,lon])







    def crear_mapa_especies(self,i,k,gs=[1,2], save = False, carpeta = './'):

        self.i = i
        self.k = k

        fig = plt.figure(figsize=(20,8))

        self.gs = fig.add_gridspec(gs[0],gs[1])

        ax_lat = fig.add_subplot(self.gs[0,0])

        self.eje_especies_lat(ax_lat, nombre = 'Especies Completo')

        ax_lon = fig.add_subplot(self.gs[0,1])

        self.eje_especies_lon(ax_lon, nombre = 'Especies Completo')


        #plt.colorbar(self.q_plot_somb)

        plt.tight_layout()

        
        if save == True:
            plt.savefig(carpeta+'/WB_Especies%s.png'%self.get_time(self.i))


        plt.show()    
          
        
                

    def crear_mapa_WB (self,i,k,gs = [5,4], save = False, carpeta = './'):

        self.k = k    
        self.i = i

        variables = {
            'QRAIN':'Greens',
            'QVAPOR':'Purples',
            'QCLOUD':'Blues',
            'QICE':'Greys',
            'QSNOW':'Reds'
            }

        fig = plt.figure(figsize = (20,15))

        self.gs = fig.add_gridspec(gs[0],gs[1])
        

        for i,var in enumerate(variables.keys()):
            ax_lat = fig.add_subplot(self.gs[i,0:2])

            self.eje_especies_ind_lat(ax_lat,var,somb_cmap= variables[var],nombre = 'Especies %s Lat'%var[:2])

            ax_lon = fig.add_subplot(self.gs[i,2:])

            self.eje_especies_ind_lon(ax_lon,var,somb_cmap= variables[var],nombre = 'Especies %s Lon'%var[:2])


        plt.tight_layout()

        plt.text(-0.2,4,self.get_time(self.i))
        
        if save == True:
            plt.savefig(carpeta+'/WB_%s.png'%self.get_time(self.i))


        plt.show()    
          
    def crear_mapa_SC_Madura(self,i,gs = [5,6], save = False, carpeta = './'):
     
       
        self.i = i

        fig = plt.figure(figsize = (20,15))

        self.gs = fig.add_gridspec(gs[0],gs[1])

        ax_conv_lat = fig.add_subplot(self.gs[0,:2])

        self.eje_conv_maduraLat (ax_conv_lat,i)

        ax_conv_lon = fig.add_subplot(self.gs[0,2:4])

        self.eje_conv_maduraLon (ax_conv_lon,i)


        #QR ADV Z
        cax_1 = fig.add_axes([0.665, 0.83 , 0.01, 0.15])
        cbar_1 = plt.colorbar(self.plot_1, ax= [ax_conv_lat,ax_conv_lon],cax = cax_qvh1, orientation = 'vertical')  
        cbar_1.set_label('QV ADV Z')
        
        #QV TEND  
        cax_2 = fig.add_axes([0.715, 0.83 , 0.01, 0.15])
        cbar_2 = plt.colorbar(self.plot_2, ax= [ax_conv_lat,ax_conv_lon],cax = cax_qvz1, orientation = 'vertical')    
        cbar_2.set_label('QV TEND')
        
        
        
                
        ax2 = fig.add_subplot(self.gs[1:3,:4])
        
        self.eje_xz(ax2,i)

        ax4 = fig.add_subplot(self.gs[3:5,:4])
        
        self.eje_yz(ax4,i)

        
        cax_ascensos = fig.add_axes([0.765, 0.83 , 0.01, 0.15])
        cbar_ascensos = plt.colorbar(self.ascensos, ax= [ax_conv_lat,ax_conv_lon],cax = cax_ascensos, orientation = 'vertical')    
        cbar_ascensos.set_label('W (m/s)')
        
        
                
        ax3 = fig.add_subplot(self.gs[0,5:6])
        
        self.hodografa(ax3,i)
        
        cax_hod = fig.add_axes([1., 0.83 , 0.01, 0.15])
        cbar_hod = plt.colorbar(self.hod, ax= ax3,cax = cax_hod, orientation = 'vertical',)
        cbar_hod.ax.set_yticklabels(np.round(np.linspace(self.press[self.i,8,self.max_lat,self.max_lon],self.press[self.i,28,self.max_lat,self.max_lon],10 ),decimals= 0))
        cbar_hod.set_label('Nivel de Presión (HPa)')
        
        
        ax5 =fig.add_subplot(self.gs[3:5,4:6])
        
        self.eje_centrado_xy(ax5,i)
        
        ax6 = fig.add_subplot(self.gs[1:3,4:6])
        
        self.eje_mapa_xy (ax6,i)
      
    
        cax_lluvia = fig.add_axes([0.815, 0.83 , 0.01, 0.15])
        cbar_lluvia = plt.colorbar(self.ps, ax= [ax_conv_lat,ax_conv_lon],cax = cax_lluvia, orientation = 'vertical')    
        cbar_lluvia.set_label('Lluvia (mm)')
        

        plt.tight_layout()

        plt.show()
    
    def crear_mapa_SC_CI (self,i,k,gs = [5,6], save = False, carpeta = './'):
        
        self.k = k    
        self.i = i
        
        fig = plt.figure(figsize = (20,15))

        self.gs = fig.add_gridspec(gs[0],gs[1])

        ax_esp_lat = fig.add_subplot(self.gs[0,:2])

        self.eje_xz(ax_esp_lat,i)

        ax_esp_lon = fig.add_subplot(self.gs[0,2:4])

        self.eje_yz(ax_esp_lon,i)




        #cax_qvh1 = fig.add_axes([0.665, 0.83 , 0.01, 0.15])
        ##cbar_qvh1 = plt.colorbar(self.qvh_plot, ax= [ax_conv_lat,ax_conv_lon],cax = cax_qvh1, orientation = 'vertical')  
        #cbar_qvh1.set_label('QV ADV H')
        
        
        #cax_qvz1 = fig.add_axes([0.715, 0.83 , 0.01, 0.15])
        #cbar_qvz1 = plt.colorbar(self.qvz_plot, ax= [ax_conv_lat,ax_conv_lon],cax = cax_qvz1, orientation = 'vertical')    
        #cbar_qvz1.set_label('QV ADV Z')
        
                
        ax2 = fig.add_subplot(self.gs[1:3,:4])
        
        
        self.eje_especies_lat(ax2, nombre = 'Especies Completo Lat')


        ax4 = fig.add_subplot(self.gs[3:5,:4])
        
        self.eje_especies_lon(ax4, nombre = 'Especies Completo Lat')


        
        cax_ascensos = fig.add_axes([0.765, 0.83 , 0.01, 0.15])
        cbar_ascensos = plt.colorbar(self.ascensos, ax= [ax_esp_lat,ax_esp_lon],cax = cax_ascensos, orientation = 'vertical')    
        cbar_ascensos.set_label('W (m/s)')
        
        
                
        ax3 = fig.add_subplot(self.gs[0,5:6])
        
        self.hodografa(ax3,i)
        
        cax_hod = fig.add_axes([1., 0.83 , 0.01, 0.15])
        cbar_hod = plt.colorbar(self.hod, ax= ax3,cax = cax_hod, orientation = 'vertical',)
        cbar_hod.ax.set_yticklabels(np.round(np.linspace(self.press[self.i,8,self.max_lat,self.max_lon],self.press[self.i,28,self.max_lat,self.max_lon],10 ),decimals= 0))
        cbar_hod.set_label('Nivel de Presión (HPa)')
        
        
        ax5 =fig.add_subplot(self.gs[3:5,4:6])
        
        self.eje_centrado_xy(ax5,i)
        
        ax6 = fig.add_subplot(self.gs[1:3,4:6])
        
        self.eje_mapa_xy (ax6,i)
      
    
        cax_lluvia = fig.add_axes([0.815, 0.83 , 0.01, 0.15])
        cbar_lluvia = plt.colorbar(self.ps, ax= [ax_esp_lat,ax_esp_lon],cax = cax_lluvia, orientation = 'vertical')    
        cbar_lluvia.set_label('Lluvia (mm)')
        

        plt.tight_layout()

        #plt.suptitle('Inicio de la convección profunda, tiempo:')
        
        if save == True:
            plt.savefig(carpeta+'/SC_CI_WB_Mod%s.png'%self.get_time(i))
            

        plt.show()

    def crear_mapa_SC_formacion (self,i,gs = [5,6], save = False, carpeta = './'):
        
        self.i = i
        
        fig = plt.figure(figsize = (20,15))
        
        self.gs = fig.add_gridspec(gs[0],gs[1])
        
        ax_formacion_lat = fig.add_subplot(self.gs[0,:2])
        
        self.eje_formacionLat (ax_formacion_lat,i)
        
        ax_formacion_lon = fig.add_subplot(self.gs[0,2:4])
        
        self.eje_formacionLon (ax_formacion_lon,i)
        
        
        cax_qvh1 = fig.add_axes([0.7, 0.81, 0.015, 0.18])
        cbar_qvh1 = plt.colorbar(self.qvh_plot_1, ax= [ax_formacion_lat,ax_formacion_lon],cax = cax_qvh1, orientation = 'vertical')    
                
        
        ax2 = fig.add_subplot(self.gs[1:3,:4])
        
        self.eje_xz(ax2,i)
        
        
        ax3 = fig.add_axes([0.45,0.3,0.155,0.12])
        
        self.hodografa(ax3,i)

        ax4 = fig.add_subplot(self.gs[3:5,:4])
        
        self.eje_yz(ax4,i)
        
        
        ax5 =fig.add_subplot(self.gs[3:5,4:6])
        
        self.eje_centrado_xy(ax5,i)
        
        ax6 = fig.add_subplot(self.gs[1:3,4:6])
        
        self.eje_mapa_xy (ax6,i)
      
        
        
        plt.tight_layout()

        #plt.suptitle ('Formación de la celda, tiempo :')
        
        
        plt.show()
    
    def crear_mapa_cv (self,i, gs = [6,2], save = False, carpeta = './'):
        
        
        fig = plt.figure(figsize=(10,12))
        self.gs = fig.add_gridspec(gs[0],gs[1])

        
        ax1 = fig.add_subplot(self.gs[0:2,1])
        
        #ax_bar_pp = fig.add_axes([0.83, 0.66, 0.02, 0.23])
        
        self.eje_mapa_xy(ax1,i)

        ax0 = fig.add_subplot(self.gs[0:2,0])
        
        self.eje_centrado_xy(ax0,i,0)
        
        ax2 = fig.add_subplot(self.gs[2:4,:])
        
        self.eje_xz(ax2,i)
        
        
        ax3 = fig.add_axes([0.135,0.46,0.155,0.12])
        
        self.hodografa(ax3,i)

        ax4 = fig.add_subplot(self.gs[4:6,:])
        
        self.eje_yz(ax4,i)
        
        
        caxr = fig.add_axes([1., 0.66, 0.02, 0.28])
        
        cbar_pp = plt.colorbar(self.p_rain, ax = [ax0,ax1],cax = caxr, orientation = 'vertical', pad = 0.03)
        cbar_pp.ax.set_ylabel('Precipitación hr')

        caxq = fig.add_axes([1.0, 0.11, 0.02, 0.45])

        cbar_qv = plt.colorbar(self.q_plot, ax = [ax2,ax4], cax = caxq, orientation = 'vertical' )
        cbar_qv.ax.set_ylabel('qvapor kg/kg')
        
        caxw = fig.add_axes([0.1,0.06,0.65,0.02])
        
        cbar_w = plt.colorbar(self.ascensos, ax = ax4, cax = caxw, orientation='horizontal')
        cbar_w.ax.set_xlabel('Ascensos y descensos (m/s)')
        
        
        
        plt.tight_layout()
        
        plt.savefig(carpeta+'collage_track.png')
        
        plt.show()

######################## Este lo voy a terminar sacando  ###############
    def crear_mapa_wb(self,i,gs = [2,2], save = False, carpeta = './'):
        

        
        fig = plt.figure(figsize=(10,12))
        self.gs = fig.add_gridspec(gs[0],gs[1])

        
        ax1 = fig.add_subplot(self.gs[0,0:1])
        
        self.eje_wb_1(ax1,i)
        
        ax2 = fig.add_subplot(self.gs[0,1:2])
        
        self.eje_wb_2(ax2,i)
        
        ax3 = fig.add_subplot(self.gs[1,:])
        
        self.eje_wb_3(ax3,i)
        
        
        cax_qten = fig.add_axes([0.93, 0.28, 0.02, 0.15])
        cbar_qten = plt.colorbar(self.qten_plot, ax= ax3,cax = cax_qten, orientation = 'vertical')


        cax_qvh = fig.add_axes([0.93, 0.12, 0.02, 0.15])
        cbar_qvh = plt.colorbar(self.qvh_plot, ax = ax3, cax = cax_qvh, orientation = 'vertical')
        
        plt.tight_layout()
        
        plt.savefig(carpeta+'collage_wb_%i.png'%i)
        plt.show()     
########################################################################

    def eje_especies_ind_lat(self,ax,q_name ,i='a', k = 'a', somb_cmap= 'Greys', nombre =''):

        data_q = {
            'QRAIN': self.data_wb_r,
            'QICE': self.data_wb_i,
            'QSNOW': self.data_wb_s,
            'QCLOUD': self.data_wb_c,
            'QVAPOR': self.data_wb_v,        
        }

        print(type(data_q['QRAIN']))

        i,k = self.set_ik(i,k)

        x2D, press2D = self.xp2D()

        lev = list(range(0,len(press2D[:,0])))
        lat = list(range(   self.max_lat-self.window,self.max_lat+self.window, 1 ))

        lon = self.max_lon

        qtend_str = q_name[:2]+'TTEND'

        qh_str = q_name[:2]+'H_ADV'

        qv_str = q_name[:2]+'_ZADV'

        qvar = self.carga_sing_var(self.data,q_name,i,lev,lat,lon   )

        qtend = self.carga_sing_var(data_q[q_name],qtend_str,k,lev,lat,lon   )
        qhadv = self.carga_sing_var(data_q[q_name],qh_str,k,lev,lat,lon   )
        qvadv = self.carga_sing_var(data_q[q_name],qv_str,k,lev,lat,lon   )


        qvar = self.setVarVert (qvar,press2D,)
        qtend = self.setVarVert (qvar,press2D)
        qhadv = self.setVarVert (qvar,press2D)
        qvadv = self.setVarVert (qvar,press2D)

        ax.pcolormesh(x2D,press2D, qvar, cmap=plt.get_cmap(somb_cmap))

        mapa.plot_contour_especies(self,ax,x2D,press2D,qtend,color='k')

        ax.quiver(x2D,press2D,qhadv,qvadv,scale= 5, color ='k')

        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        ax.set_facecolor('white')



    def eje_especies_ind_lon(self,ax,q_name ,i='a', k = 'a',somb_cmap='Greys', nombre =''):

        data_q = {
            'QRAIN': self.data_wb_r,
            'QICE': self.data_wb_i,
            'QSNOW': self.data_wb_s,
            'QCLOUD': self.data_wb_c,
            'QVAPOR': self.data_wb_v,        
        }

        i,k = self.set_ik(i,k)

        y2D, press2D = self.yp2D()

        lev = list(range(0,len(press2D[:,0])))
        lon = list(range(   self.max_lon-self.window,self.max_lon+self.window, 1 ))

        lat = self.max_lat

        qtend_str = q_name[:2]+'TTEND'

        qh_str = q_name[:2]+'H_ADV'

        qv_str = q_name[:2]+'_ZADV'

        qvar = self.carga_sing_var(self.data,q_name,i,lev,lat,lon   )

        qtend = self.carga_sing_var(data_q[q_name],qtend_str,k,lev,lat,lon   )
        qhadv = self.carga_sing_var(data_q[q_name],qh_str,k,lev,lat,lon   )
        qvadv = self.carga_sing_var(data_q[q_name],qv_str,k,lev,lat,lon   )

        qvar = self.setVarVert (qvar,press2D)
        qtend = self.setVarVert (qvar,press2D)
        qhadv = self.setVarVert (qvar,press2D)
        qvadv = self.setVarVert (qvar,press2D)

        ax.pcolormesh(y2D,press2D, qvar, cmap=plt.get_cmap(somb_cmap))

        mapa.plot_contour_especies(self,ax,y2D,press2D,qtend,color='k')

        ax.quiver(y2D,press2D,qhadv,qvadv,scale= 5, color ='k')

        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        ax.set_facecolor('white')

        
   


    def plot_contour_especies(self,ax,x2D,press2D, var, color):



        def mult(x):
            return x*100000
        var = mult(var)
        
        m = np.max([abs(np.min(var)),abs(np.max(var))])

        levels =  np.linspace(-m,m,10)
        print(levels)
        cs = ax.contour(x2D,press2D,var,levels = levels,colors = color) 
        cs.collections[5].remove()
        cs.collections[4].remove()

        #plt.setp( zc, linewidth=5)
        plt.clabel(cs,inline=1, fmt='%1.2f')

        l,_ = cs.legend_elements()
        
        return l

    def plot_barbs(self,ax,x2D,press2D):


        a=4 #Separacion de flechas

        end_v = 30
        ax.barbs(x2D[:end_v,3:-3:a],press2D[:end_v,3:-3:a], u[:end_v,3:-3:a],w[:end_v,3:-3:a],
                           color ='#393939',  pivot = 'middle',length = 7, linewidth =1 , zorder = 4 ) 



    def eje_especies_lat(self,ax, i='a', k = 'a', nombre =''):


        i,k = self.set_ik(i,k)

        x2D, press2D = self.xp2D()

        mapa.cargar_var_especies(self,i,k,self.max_lat,list(range(self.max_lon-self.window,self.max_lon+self.window,1)))
        
        mapa.eje_especies(self,ax,x2D,press2D)
        
        plt.title("Corte vertical XZ: %s (Contorno tend .10**-5)"%nombre)  


    def eje_especies_lon(self,ax, i='a', k = 'a', nombre =''):

        i,k = self.set_ik(i,k)

        y2D, press2D = self.yp2D()

        mapa.cargar_var_especies(self,i,k,list(range(self.max_lat-self.window,self.max_lat+self.window,1)), self.max_lon)
       
        mapa.eje_especies(self,ax,y2D,press2D)
        
        plt.title("Corte vertical YZ: (Contorno tend .10**-5) %s"%nombre)  

    def eje_especies (self,ax,x2D,press2D):

        alpha = 0.6

        mask = 0.000001

        ice = self.setVarVert (self.qice,press2D, mask = mask)
        cloud = self.setVarVert (self.qc,press2D, mask = mask)
        rain = self.setVarVert (self.qr,press2D, mask = mask)
        vapor = self.setVarVert (self.qv,press2D, mask = mask)
        snow = self.setVarVert (self.qs,press2D, mask = mask)

        mask = 0   

        ice_ten = self.setVarVert (self.qitend,press2D,mask = mask)
        cloud_ten = self.setVarVert (self.qctend,press2D, mask = mask)
        rain_ten = self.setVarVert (self.qrtend,press2D, mask = mask)
        vapor_ten = self.setVarVert (self.qvtend,press2D, mask = mask)
        snow_ten = self.setVarVert (self.qstend,press2D, mask = mask)
        
        i_hadv = self.setVarVert (self.qih_adv,press2D,mask = mask)
        c_hadv = self.setVarVert (self.qch_adv,press2D, mask = mask)
        r_hadv = self.setVarVert (self.qrh_adv,press2D, mask = mask)
        v_hadv = self.setVarVert (self.qvh_adv,press2D, mask = mask)
        s_hadv = self.setVarVert (self.qsh_adv,press2D, mask = mask)

        i_vadv = self.setVarVert (self.qiz_adv,press2D,mask = mask)
        c_vadv = self.setVarVert (self.qcz_adv,press2D, mask = mask)
        r_vadv = self.setVarVert (self.qrz_adv,press2D, mask = mask)
        v_vadv = self.setVarVert (self.qvz_adv,press2D, mask = mask)
        s_vadv = self.setVarVert (self.qsz_adv,press2D, mask = mask)

        self.q_plot_somb = ax.pcolormesh(x2D,press2D, cloud, cmap=plt.get_cmap('Blues'), alpha = alpha)
        self.q_plot_somb = ax.pcolormesh(x2D,press2D, vapor, cmap=plt.get_cmap('Purples'), alpha = alpha)
        self.q_plot_somb = ax.pcolormesh(x2D,press2D, snow, cmap=plt.get_cmap('Reds'), alpha = alpha)   
        self.q_plot_somb = ax.pcolormesh(x2D,press2D, ice, cmap=plt.get_cmap('Greys'), alpha = 0.3)
        self.q_plot_somb = ax.pcolormesh(x2D,press2D, rain, cmap=plt.get_cmap('Greens'), alpha = 0.7)

    

        l1 = mapa.plot_contour_especies(self,ax, x2D, press2D,ice_ten,color = 'grey')

        l2 = mapa.plot_contour_especies(self,ax, x2D, press2D,cloud_ten,color = 'b')

        l3 = mapa.plot_contour_especies(self,ax,x2D,press2D,rain_ten,color='g')

        l4 = mapa.plot_contour_especies(self,ax,x2D,press2D,vapor_ten,color='indigo')

        l5 = mapa.plot_contour_especies(self,ax,x2D,press2D,snow_ten,color='r')

        #q1 = ax.quiver(y2D,press2D,r_hadv,v_vadv,alpha = 0.8,scale= 5, color = 'g')
        #q2 = ax.quiver(y2D,press2D,v_hadv,v_vadv,alpha = 0.1,scale= 5, color='indigo')
        #q3 = ax.quiver(y2D,press2D,i_hadv,i_vadv,alpha = 0.8,scale= 5, color='k')




        ax.legend([l1[-1],l2[-1],l3[-1],l4[-1],l5[-1]], 
                ['Hielo','Nube','Lluvia','Vapor', 'Nieve'])


        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        ax.set_facecolor('white')


    def eje_q_lat(self, ax,vartend,varh,varz,i='a',k='a', nombre='', hor = True, vert = True, somb_scale = 100000 ):
        
        i,k = self.set_ik(i,k)

        x2D, press2D = self.xp2D()
        
        varh_adv = self.setVarVert( varh [k,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window]*100,press2D,mask = 0.1)     
        varz_adv = self.setVarVert( varz [k,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window]*100,press2D,mask = 0.1)    
        var_tend = self.setVarVert( vartend  [k,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window]*somb_scale,press2D,mask = 0.1)
        

        self.q_plot_somb = ax.pcolormesh(x2D,press2D, var_tend, cmap=plt.get_cmap('RdBu_r'), vmin = -5, vmax = 5)        

        
        ## Esto esta provisorio
        
        clev = [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2]

        if vert == True:
            pc = ax.contour(x2D,press2D, varz_adv,levels = 8,linewidth = 3,fmt='%0f' ,cmap = plt.get_cmap('PuOr'), vmin = -5, vmax = 5, zorder = 10 )
            plt.clabel(pc, inline=True, fontsize=8,fmt='%1.0f', colors = 'k')

        if hor == True: 
            pc = ax.contour(x2D,press2D,varh_adv, levels=8, linewidth = 3, cmap =  plt.get_cmap('PuOr'), vmin = -5, vmax = 5, zorder = 10 )
            plt.clabel(pc, inline=True, fontsize=8,fmt='%1.0f',colors = 'k')
                        
        ax.plot(x2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical XZ: %s"%nombre)

    def eje_q_lon(self,ax,vartend,varh,varz,i='a',k='a', nombre='', hor = True, vert = True, somb_scale = 100000):
        
        i,k = self.set_ik(i,k)

        y2D,press2D  = self.yp2D()
        
        varh_adv = self.setVarVert( varh[k,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon]*100,press2D,mask = 0.1)          
        varz_adv = self.setVarVert( varz[k,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon]*100,press2D,mask = 0.1)        
        var_tend = self.setVarVert( vartend[k,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon]*somb_scale,press2D,mask = 0.1)

        ax.pcolormesh(y2D,press2D, var_tend, cmap=plt.get_cmap('RdBu_r'), vmin = -5, vmax = 5)  
        

        ## Esto esta provisorio
        if hor == True:
            pc = ax.contour(y2D,press2D, varz_adv,levels =8, linewidth = 3, cmap = plt.get_cmap('PuOr'),vmin = -5, vmax = 5, zorder = 10 )
            plt.clabel(pc, inline=True, fontsize=8,fmt='%1.0f', colors = 'k')

        if vert == True:
            pc = ax.contour(y2D,press2D,varh_adv, linewidth = 3, cmap =  plt.get_cmap('PuOr'), vmin = -5, vmax = 5, zorder = 10 )
            plt.clabel(pc, inline=True, fontsize=8, fmt='%1.0f', colors = 'k')

            
        ax.plot(y2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical YZ: %s"%nombre)

    def eje_conv_maduraLon (self,ax,i):
    
        y2D,press2D  = self.yp2D()
        
        qv_tend = self.setVarVert( self.qvtend[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D, mask =0.01)             
        qr_tend = self.setVarVert( self.qrtend[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D, mask=0.01)
        qrz_adv = self.setVarVert(self.qrz_adv[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D,mask = 0.005)
        
        
        qv_tend = np.ma.masked_array(qv_tend,mask =  press2D<600)
        qr_tend = np.ma.masked_array(qr_tend, mask = press2D<800)
        
        ax.pcolormesh(y2D,press2D, qrz_adv, cmap = plt.get_cmap('BrBG'))        
        ax.pcolormesh(y2D,press2D,qv_tend, cmap= plt.get_cmap('RdBu_r'))
        ax.pcolormesh(y2D,press2D,qr_tend,cmap = plt.get_cmap('seismic'))
        
        
        ax.hlines(800, y2D[0,0],y2D[0,-1])

        ax.hlines(600, y2D[0,0],y2D[0,-1])

        ax.plot(y2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical YZ: Adv H Qv (somb inferior), \n Adv Z Qv (somb medio), Tend. Qi (somb superior)  ")

    def eje_conv_maduraLat(self,ax,i):
        
        x2D,press2D = self.xp2D()
        
        
        qv_tend = self.setVarVert( self.qvtend [i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window],press2D, mask =0.01)     
        qrz_adv = self.setVarVert( self.qrz_adv[i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window],press2D, mask =0.01)     
        qr_tend = self.setVarVert( self.qrtend [i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window],press2D, mask =0.01)    
        
        
        
        qrz_adv = np.ma.masked_array(qrz_adv, mask = press2D<600)
        qv_tend = np.ma.masked_array(qv_tend, mask = press2D<800)
        
        
        ax.pcolormesh(x2D,press2D,qrz_adv,cmap = plt.get_cmap('seismic'))
        ax.pcolormesh(x2D,press2D, qv_tend, cmap=plt.get_cmap('BrBG'))
        ax.pcolormesh(x2D,press2D,qr_tend, cmap= plt.get_cmap('RdBu_r'))
        

        
        ax.hlines(800, x2D[0,0],x2D[0,-1])

        ax.hlines(600, x2D[0,0],x2D[0,-1])
                        
        ax.plot(x2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical XZ: Adv H Qv (somb inferior),\n Adv Z Qv (somb medio), Tend. Qi (somb superior)  ")

    def eje_convIniLon(self,ax,i):
        
        y2D,press2D  = self.yp2D()
        
        qvh_adv = self.setVarVert( self.qvh_adv[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D, mask =0.01)
        
        qvh_adv = np.ma.masked_array(qvh_adv,mask =  press2D<800)
        
        
        qvz_adv = self.setVarVert( self.qvz_adv[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D, mask=0.01)
        
        qvz_adv = np.ma.masked_array(qvz_adv, mask = press2D>800)
        qvz_adv = np.ma.masked_array(qvz_adv, mask = press2D<400)
        
        qc_tend = self.setVarVert(self.qctend[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D,mask = 0.005)
        qc_tend = np.ma.masked_array(qc_tend,mask =  press2D>400)
        
        ax.pcolormesh(y2D,press2D,qc_tend, cmap= plt.get_cmap('RdBu_r'))
        
        self.qvz_plot = ax.pcolormesh(y2D,press2D, qvz_adv, cmap = plt.get_cmap('BrBG'), vmin = -0.15, vmax = 0.15)
        
        self.qvh_plot = ax.pcolormesh(y2D,press2D,qvh_adv,cmap = plt.get_cmap('RdBu_r'), vmin = -0.15, vmax = 0.15)
        
        ax.hlines(800, y2D[0,0],y2D[0,-1])

        ax.hlines(400, y2D[0,0],y2D[0,-1])

        ax.plot(y2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical YZ: Adv H Qv (somb inferior),\n Adv Z Qv (somb medio), Tend. Qi (somb superior)  ")

    def eje_convIniLat(self,ax,i):
        
        x2D,press2D = self.xp2D()
        
        qvh_adv = self.setVarVert( self.qvh_adv[i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window], press2D, mask =0.005)
        
        qvz_adv = self.setVarVert( self.qvz_adv[i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window], press2D, mask=0.005)
        
        qvz_adv = np.ma.masked_array(qvz_adv, mask = press2D>800)
        qvz_adv = np.ma.masked_array(qvz_adv, mask = press2D<400)
        
        qc_tend = self.setVarVert(self.qctend[i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window],press2D,mask = 0.005)
        qc_tend = np.ma.masked_array(qc_tend,mask =  press2D>400)
        
        qvh_adv = np.ma.masked_array(qvh_adv, mask = press2D<800)
        
        
        ax.pcolormesh(x2D,press2D, qvz_adv, cmap=plt.get_cmap('BrBG'))
        
        ax.pcolormesh(x2D,press2D,qc_tend, cmap= plt.get_cmap('RdBu_r'))
        
        ax.pcolormesh(x2D,press2D,qvh_adv,cmap = plt.get_cmap('RdBu_r'))
        
        ax.hlines(800, x2D[0,0],x2D[0,-1])

        ax.hlines(400, x2D[0,0],x2D[0,-1])
                        
        ax.plot(x2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical XZ: Adv H Qv (somb inferior),\n  Adv Z Qv (somb medio), Tend. Qi (somb superior)  ")
    
    def eje_formacionLon(self,ax,i):

        y2D, press2D  = self.yp2D()
        
        qvh_adv = self.setVarVert( self.qvh_adv[i,:,self.max_lat-self.window:self.max_lat+self.window, self.max_lon],press2D, mask = 0.05)
        
        ax.pcolormesh(y2D,press2D,qvh_adv, cmap = plt.get_cmap('RdBu_r'), vmin =  -0.2, vmax = 0.2)
                        
        ax.plot(y2D[0,:],self.topo[self.max_lat-self.window:self.max_lat+self.window, self.max_lon],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical YZ: Adv H Qv (somb)  ")
          
    def eje_formacionLat(self,ax,i):
    
        x2D, press2D = self.xp2D()      
        
        qvh_adv = self.setVarVert(self.qvh_adv[i,:,self.max_lat,self.max_lon-self.window:self.max_lon+self.window],press2D, 0.05)
        
        
        self.qvh_plot_1 = ax.pcolormesh(x2D,press2D,qvh_adv,cmap = plt.get_cmap('RdBu_r'), vmin = -0.2, vmax = 0.2)
                        
        ax.plot(x2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title("Corte vertical XZ: Adv H Qv  ")
    
    def eje_wb_1(self,ax,i):
        
        pad = 10
        
        
        x = self.x[pad:-pad,pad:-pad]
        y = self.y[pad:-pad,pad:-pad]
        vmax = np.max(conv_wv[15:-15,15:-15])
        
        #ax.pcolormesh(x,y, conv_wv, cmap= plt.get_cmap('BrBG'), vmax=vmax,vmin = -vmax)
        
        rain = self.rain[i,pad:-pad,pad:-pad]
        
        ax.pcolormesh(x,y, rain, cmap= plt.get_cmap('jet'), vmin =0, vmax= 50)
        
        
        u = self.u[i,0,pad:-pad,1+pad:-pad]
        v = self.v[i,0,1+pad:-pad,pad:-pad]
        
        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(self.x[0,pad:-pad],self.y[pad:-pad,0],u,v,
                       arrowstyle ='->',linewidth=lw, color ='k', density=0.8)

        
        self.bm.drawparallels(np.linspace(-35.,-30,6),labels=[1,0,0,0],color='gray')
        self.bm.drawmeridians(np.linspace(-66,-61,6),labels=[0,0,0,1], color = 'gray')
        
        
        return
        
    def eje_wb_2(self,ax,i):
        
        
        qvtend = self.qvtend[i,0,:,:]
        
        
        vmax = np.max(qvtend[10:-10,10:-10])
        ax.pcolormesh(self.x,self.y,qvtend,cmap= plt.get_cmap('seismic'), vmax = vmax, vmin = -vmax)
                
        u = self.u[i,0,:,1:]
        v = self.v[i,0,1:,:]
        
        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(self.x[0,:],self.y[:,0],u,v,
                       arrowstyle ='->',linewidth=lw, color ='k', density=0.8)

        
        col = 'orange'

        xs,ys = self.get_square()
        
        self.bm.plot(xs, ys, latlon = True)   
        
        plt.scatter(self.max_lat,self.max_lon, marker='o',zorder= 4, c=col)


        
        self.bm.drawparallels(np.linspace(-35.,-30,6),labels=[1,0,0,0],color='gray')
        self.bm.drawmeridians(np.linspace(-66,-61,6),labels=[0,0,0,1], color = 'gray')
        
        return
    
    def eje_wb_3(self,ax,i):
        
        x=self.lon[0,:]
        
        end = 59
        
        press = self.press[i,:,self.max_lat,self.max_lon]
        
        
        x2D, press2D = np.meshgrid(x, press)
        x2D = x2D[:end,:]
        press2D = press2D[:end,:]
        
        qvtend = self.qvtend[i,:,self.max_lat,:]
        qvtend = np.ma.masked_array(qvtend, mask = abs(qvtend)<0.000005)
        qvtend = np.ma.masked_array(qvtend, mask = press2D>70000)
        
        
        qvh_adv = self.qvh_adv[i,:,self.max_lat,:]
        qvh_adv =  np.ma.masked_array(qvh_adv, mask = abs(qvh_adv)<0.5)
        qvh_adv =np.ma.masked_array(qvh_adv,mask = press2D<70000)
        
        
        self.qten_plot = ax.pcolormesh (x2D,press2D,qvtend, cmap = plt.get_cmap('coolwarm'))

        self.qvh_plot = ax.pcolormesh (x2D, press2D, qvh_adv, cmap = plt.get_cmap('Spectral'), vmin=-2, vmax = 2)
        
       

        #self.ascensos = ax.contour(x2D,press2D, self.w[i,:-1,self.max_lat,:])
        
        
        
        """
        clev, clev_bool = self.get_clevs(self.w[i,:-1,self.max_lat,:],ncont=9)
        
        if clev_bool:
            pc = ax.contour(x2D,press2D,self.w[i,:-1,self.max_lat,:],clev, colors = 'k', zorder = 3)
            clab = ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')

        """
        
        plt.ylim(self.pmax,self.pmin)
        
        
        return
    
    def eje_mapa_xy(self,ax,i, viento = True):
        
        print(i)
                
        
        plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left=False,
            labelleft=False,
        )

        self.p_rain = ax.pcolormesh(self.x,self.y,self.rain[i,...],cmap=plt.get_cmap('jet'),vmin = 2,vmax = 40, zorder=3)        
                
            
        ax.contour(self.x,self.y,self.topo, colors = 'black', linewidth = 0.5, zorder = -2)

        col = 'orange'

        
        
 
        
        try: 
            xs,ys = self.get_square()
            self.bm.plot(xs, ys, latlon = True)  
            plt.scatter(self.max_lat,self.max_lon, marker='o',zorder= 4, c=col)
        except:
            print('no se puede mostrar el punto, faltan lat y lon max')
        if viento: 
            try:
                self.plot_viento_xy(i)
            except:
                print('No se puede graficar el viento')
         
        self.bm.drawparallels(np.linspace(-35.,-30,6),labels=[1,0,0,0],color='gray')
        self.bm.drawmeridians(np.linspace(-66,-61,6),labels=[0,0,0,1], color = 'gray')
        
        try:
            plt.title('Trackeo de la celda - %s'%self.get_time(i))
        except:
            plt.title('Trackeo de la celda - tiempo desconocido')
        
        return

    def plot_viento_xy(self,i,):
        u = self.u[i,0,:,1:]
        v = self.v[i,0,1:,:]
        
        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(self.x[self.max_lat,:],self.y[:,self.max_lon],u,v,
                       arrowstyle ='->',linewidth=lw, color ='k', density=0.8)

        return 

    def eje_centrado_xy(self,ax,i,lev_qv = 0,w_level=17):

        lat_ini = self.max_lat-self.window
        lat_end = self.max_lat+self.window

        lon_ini = self.max_lon-self.window
        lon_end = self.max_lon+self.window

        lat = list(range(lat_ini,lat_end,1))
        lon = list(range(lon_ini,lon_end,1))
        
        x = self.x[lat_ini:lat_end,lon_ini:lon_end]
        y = self.y[lat_ini:lat_end,lon_ini:lon_end]
        qv = mapa.carga_sing_var(self,self.data,'QVAPOR',i,lev_qv,lat,lon)

        rain = self.rain [i,self.max_lat-self.window:self.max_lat+self.window, self.max_lon-self.window:self.max_lon+self.window]
        topo = self.topo [lat_ini:lat_end,lon_ini:lon_end]
        w = self.w[i,w_level,lat_ini:lat_end,lon_ini:lon_end]
        u = self.u[i,17,2+lat_ini:lat_end-2,lon_ini+2:lon_end-2]
        v = self.v[i,17,2+lat_ini:lat_end-2,lon_ini+2:lon_end-2]
        
        print(x.shape,y.shape,qv.shape, rain.shape)

        ax.pcolormesh(x ,y ,qv,cmap=plt.get_cmap('BrBG')) 
          
          
        
        self.ps = ax.pcolormesh(x,y,rain, cmap=plt.get_cmap('jet'), vmin =  0, vmax=50, zorder = 2)  

        ax.contour(x,y,topo,colors = 'black', linewidth = 0.5, zorder = 1,alpha = 0.5)

        ax.scatter(x[self.window,self.window],y[self.window,self.window],c= 'orange',zorder = 3)
        

        clev,clev_bool =  self.get_clevs(w)       
        
        if clev_bool:
            pc = ax.contour(x,y,w,clev, cmap = self.wcmap )
            ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')


        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(x[self.window,2:-2],y[2:-2,self.window],u,v,
                                 arrowstyle ='->',linewidth=lw, color ='k', density=0.8)

        

        ax.hlines(y[self.window-1,self.window-1],x[0,0],x[-1,-1], linestyles = 'dashed'
                   ,colors='r',label='Corte vertical')

        plt.scatter(x[self.window-1,self.window-1],y[self.window-1,self.window-1])
        
        
        plt.xticks([])
        plt.yticks([])
        

        plt.title('Mapa centrado en lat: %s lon: %s \n Qv(somb) Precip (somb2) W(cont)'%(self.max_lat,self.max_lon))
        
        return
    
    def eje_yz(self,ax,i):
        
        ### ESTO VA A CAMBIAR CON MAX_LON
        
        y=self.lat[:,self.max_lon]
        
        end = 59
        
        
        press = self.press[i,:,self.max_lat,self.max_lon]
        
        y2D, press2D = np.meshgrid(y, press)
        y2D = y2D[:end,self.max_lat-self.window:self.max_lat+self.window]
        press2D = press2D[:end,self.max_lat-self.window:self.max_lat+self.window]
                
        
        u = self.u[i,:end,self.max_lat-self.window:self.max_lat+self.window,self.max_lon]
        u =  self.mask_topo_simple(u, press2D)
        
        v = self.v[i,:end,self.max_lat-self.window:self.max_lat+self.window,self.max_lon]        
        v =  self.mask_topo(v,press2D,self.max_lat)
        
        w = self.w[i,:end,self.max_lat-self.window:self.max_lat+self.window,self.max_lon]
        w =  self.mask_topo_simple(w, press2D)
        w_somb = np.ma.masked_array(w, abs(w)<1)
        
        titae =self.t[i,:end,self.max_lat-self.window:self.max_lat+self.window,self.max_lon]
        
        titae = self.mask_topo(titae,press2D,self.max_lat)
        #titae = np.ma.array(titae , mask = np.abs(titae) < 0)
        titae = np.ma.masked_array(titae, press2D<= 100) + 290

        red = (70/255, 244/255)
        green = (130/255,164/255)
        blue =(180/255,96/255)

        cmap_per = self.get_cmap(red,green,blue)

         
        q_plot = ax.pcolormesh (y2D,press2D,titae, cmap = cmap_per,vmin = 298, vmax = 312)

        self.ascensos = ax.pcolormesh(y2D,press2D,w_somb,cmap=self.wcmap,vmin = -20,vmax=20)


        clev, clev_bool = self.get_clevs(v,ncont=9)
        
        if clev_bool:
            pc = ax.contour(y2D,press2D,u,clev, colors = 'k', zorder = 3)
            clab = ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')


        a=3 #Separacion de flechas

        #speed = np.sqrt((v[:end,::a]*2) ** 2 + (w[:end,::a]*2) ** 2)
        #lw = speed / speed.max()

        end_v = 30

        M = np.hypot(v[:end_v,::a], w[:end_v,::a])

        plot_u = ax.barbs(y2D[:end_v,3:-3:a],press2D[:end_v,3:-3:a], v[:end_v,3:-3:a],w[:end_v,3:-3:a],
                           color ='#393939')#,  pivot = 'middle',length = 7, linewidth =1 , zorder = 4 )       
        
        
        
        
        ax.plot(y2D[0,:],self.topo[self.max_lat-self.window:self.max_lat+self.window,self.max_lon],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        
        plt.title('Corte vertical YZ: W, THM (somb) - Viento uw y U cont \n Latitud: %s - Longitud: %s'%(str(self.max_lat),str(self.max_lon)))


        
        return
    
    def eje_xz(self,ax,i):
        
        ### ESTO VA A CAMBIAR CON MAX_LAT        
        y_corte = 120
        
        x = self.lon[y_corte,:]
                        
        end = 59
        
        w = self.w[i,:end,self.max_lat,self.max_lon-self.window:self.max_lon+self.window]
        u = self.u[i,:end,self.max_lat,self.max_lon-self.window:self.max_lon+self.window]
        v = self.v[i,:end,y_corte,self.max_lon-self.window:self.max_lon+self.window]
       
        thm =self.t[i,:end,y_corte,self.max_lon-self.window:self.max_lon+self.window]

    
        press = self.press[i,:,y_corte,y_corte]
        
        x2D, press2D = np.meshgrid(x, press)
        x2D = x2D[:end,self.max_lon-self.window:self.max_lon+self.window]
        press2D = press2D[:end,self.max_lon-self.window:self.max_lon+self.window]
                
        titae = self.mask_topo(thm,press2D,y_corte)
        #titae = np.ma.array(titae , mask = np.abs(titae) < 0)
        titae = np.ma.masked_array(titae, press2D<= 100) + 290

        v =  self.mask_topo(v,press2D,y_corte)
        u =  self.mask_topo_simple(u, press2D)
        w =  self.mask_topo_simple(w, press2D)
        w_somb = np.ma.masked_array(w, abs(w)<1)
        
        
        red = (70/255, 244/255)
        green = (130/255,164/255)
        blue =(180/255,96/255)

        cmap_per = self.get_cmap(red,green,blue)

        
        

        self.q_plot = ax.pcolormesh (x2D,press2D,titae, cmap = cmap_per,vmin = 298, vmax = 312)

        self.ascensos = ax.pcolormesh(x2D,press2D,w_somb,cmap=self.wcmap,vmin = -20,vmax=20)


        
        clev, clev_bool = self.get_clevs(v,ncont=9)
        
        if clev_bool:
            pc = ax.contour(x2D,press2D,v,clev, colors = 'k', zorder = 3)
            ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')


        a=4 #Separacion de flechas

        end_v = 30
        ax.barbs(x2D[:end_v,3:-3:a],press2D[:end_v,3:-3:a], u[:end_v,3:-3:a],w[:end_v,3:-3:a],
                           color ='#393939',  pivot = 'middle',length = 7, linewidth =1 , zorder = 4 ) 


        ax.plot(x2D[0,:],self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window],linewidth = 5, color = 'k', zorder= 3 )
        
        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')

        plt.title('Corte vertical XZ: W, THM (somb) - Viento uw y V cont \n - Latitud: %s - Longitud: %s'%(str(self.max_lat),str(self.max_lon)))

        
        
        return     

    def hodografa(self,ax,i):

        ini = 8

        end = 28

        v_c = self.v[i,ini:end,self.max_lat,self.max_lon]
        u_c = self.u[i,ini:end,self.max_lat,self.max_lon]

        cm = plt.cm.get_cmap('viridis')
        colors=[cm(1.*k/(end-ini)) for k in range(end-ini)]

        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('center')
        # Eliminate upper and right axes
        ax.spines['right'].set_color('k')
        ax.spines['top'].set_color('k')
        ax.spines['left'].set_color('k')
        ax.spines['bottom'].set_color('k')
        # Show ticks in the left and lower axes only
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        vmax = np.max(abs(v_c))+2
        umax = np.max(abs(u_c))+2

        plt.xlim(-umax, umax)
        plt.ylim(-vmax,vmax)

        ax.plot(u_c,v_c,color = 'k', linewidth = 2)

        xy = range(ini,end)

        self.hod = ax.scatter(u_c,v_c, c =xy, cmap = cm, zorder = 6 )

        ax.set_facecolor('white')
        ax.grid()

        clevs_hod = np.round(self.press[i,ini:end, self.max_lat,self.max_lon])

        plt.title('Hodografa')

        return



    def buscar_max(self,i,max_lat=-1,max_lon=-1,w_level= 18):
        
        
        wsample = self.w[i,...]
        
        if max_lat == -1 or max_lon == -1: 
            
            maximo = np.argwhere((wsample[w_level,:,:] == wsample[w_level,:,:].max()))[0]

            max_lat_pr = 0
            max_lon_pr = 0
        else:
            
            wminisample = wsample[:,max_lat-self.stride:max_lat+self.stride,max_lon-self.stride_lon:max_lon+self.stride_lon]
            maximo = np.argwhere(wminisample[w_level,:,:] == wminisample[w_level,:,:].max())[0]
            
            max_lat_pr = max_lat - self.stride
            max_lon_pr = max_lon - self.stride_lon 


            

        max_lat =  maximo[0] + max_lat_pr
        max_lon =  maximo[1] + max_lon_pr

        if max_lon + self.window >= self.x.shape[0]: 
            max_lat = 0
            max_lon = 0
            
        self.max_lat = max_lat
        self.max_lon = max_lon
            
            
        
        return    
    
    def mask_topo(self,var,press2D,y_corte):


        for j in range(var.shape[1]):

            top = self.topo[y_corte,j]
            
            var[:,j] = np.ma.array(var[:,j], mask = press2D[:,j] > top)
        
        return var
    
    def mask_topo_simple(self,var, press2D):
    
        var= np.ma.masked_array(var, press2D>self.topo[self.max_lat,self.max_lon-self.window:self.max_lon+self.window])

        return var
    
    def get_cmap(self,red,green,blue):
        
        
        cdict = {
        'blue': [(0,0,1),(0.4,blue[0],blue[0]), (1, 0.95, 1)],
        'green': [(0,0,1),(0.4,green[0],green[0]), (1, 0.99, 0)],
        'red': [(0,0,1),(0.4,red[0],red[0]) , (1,0.99, 1)]}
        cmap_per = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,25)
        
        return cmap_per
    
    def get_clevs(self,variable,ncont = 5):
        
        c_min = np.round(np.min(variable))
        c_max = np.round(np.max(variable))
        c_m = np.max((abs(c_min),abs(c_max)))
        
        clev_bool = c_min != c_max
        
        if clev_bool:
            clev=np.linspace(-c_m,c_m,ncont)
            clev= np.delete(clev,np.argwhere(clev == 0))
            
            return clev,clev_bool
        else: return 0,0
        
    def get_time(self,i):
        
        tiempo = ''.join(str(self.times[i]).replace('b','').replace("'",'').replace('[','').replace(']','').replace('\n','').replace(' ',''))
        return tiempo
               
    def get_square(self):
        
        lower_left = (self.lon[self.max_lat-self.window,self.max_lon-self.window],self.lat[self.max_lat-self.window,self.max_lon-self.window])
        lower_right =(self.lon[self.max_lat-self.window,self.max_lon+self.window],self.lat[self.max_lat-self.window,self.max_lon+self.window])
        upper_left = (self.lon[self.max_lat+self.window,self.max_lon-self.window],self.lat[self.max_lat+self.window,self.max_lon-self.window])
        upper_right =(self.lon[self.max_lat+self.window,self.max_lon+self.window],self.lat[self.max_lat+self.window,self.max_lon+self.window]) 


        xs = [lower_left[0], upper_left[0],
              upper_left[0], upper_right[0],
              upper_right[0], lower_right[0],
              lower_right[0], lower_left[0]]
        ys = [lower_left[1], upper_left[1],
              upper_left[1], upper_right[1],
              upper_right[1], lower_right[1],
              lower_right[1], lower_left[1]]
        
        return xs,ys

    def xp2D(self):
        
        x=self.lon[self.max_lat,self.max_lon-self.window:self.max_lon+self.window]    
        
        press = self.press[self.i,:,self.max_lat,self.max_lon]
        
        x2D, press2D = np.meshgrid(x, press)       
        
        return x2D, press2D
    
    def yp2D(self):

        y=self.lat[self.max_lat-self.window:self.max_lat+self.window, self.max_lon]    
        
        press = self.press[self.i,:,self.max_lat,self.max_lon]
        
        y2D, press2D = np.meshgrid(y, press)   
        
        return y2D, press2D
          
    def setVarVert(self,var,press2D, mask = ''):
        
        var = self.mask_topo_simple(var,press2D)
        
        if mask != '':
            var = np.ma.masked_array(var,mask = abs(var)<mask)
        
        return var

    def set_ik(self,i,k):
        if i == 'a': i = self.i
        if k == 'a': k = self.k

        return i,k
  
