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

def corte_vert(archivo,tiempos = [0], n_var_s = 'no', n_var_c='no',y_corte = 200, n_var_x = 'XLONG', set_topo = True,
              guardar = False, show = True):
    

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
        
        print(time.shape)

        
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
        minimo = -10#np.abs(np.min(ow))
        maximo = +10
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
        wcmap = plt.get_cmap('coolwarm')
        #wcmap.set_bad(color=colmiss)
        #wnorm = Normalize(-xtrm,xtrm)
        for j in range(305):
            
            if set_topo ==True:
                top = topografia[t,i,j]
            else: top = 100000
                

            ##### HAY UN PROBLEMA CON  LA MASCARA
            var_s[:,j] = np.ma.array(var_s[:,j], mask = press2D[:,j] > top)
            var_c[:,j] = np.ma.array(var_c[:,j], mask = press2D[:,j] > top)

        var_s = np.ma.array(var_s , mask = np.abs(var_s) < 1)
        #v = np.ma.array(w , mask = np.abs(w) < 2)
        s = ax.pcolormesh(x2D,press2D,var_s, cmap = wcmap,vmax =maximo,vmin =minimo)
        
        ncont = 5
        c_min, c_max =np.round(np.min(var_c)),np.round(np.max(var_c))
        if c_min != c_max:
            clev=np.linspace(c_min,c_max,ncont)
            clev= np.delete(clev,np.where(clev==0))
            pc = plt.contour(x2D,press2D,var_c,clev, colors = 'k', )

            clab = plt.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')
        
        if set_topo == True: plt.plot(x2D[0,:],topografia[0,i,:],linewidth = 5, color = 'k' )

        
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
            plt.savefig('./img%sZsec/CorteVert%s_%s_%i' %(carpeta,carpeta[1:-1],time,t//4),format='png')
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
        
        print(var_viento)
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
"""def set_var (var):
"""
"""
    if len(var.shape) == 4: var_set = to_np(var[i,...])
    elif len(var.shape) == 3: var_set = to_np(var[i,...])
    else: var_set = to_np(var)
    return var
"""
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
    print(speed.shape)
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
    print(vminr, vmaxr)
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
               lev_rich = 15
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
    
    print('Este gráfico es de una función externa')
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
        print(tiempos[-1])
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
    
    print('Lozada: ',lats[lats == lozada[0]], lons[lons == lozada[1]])

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
        print(tiempos[j] + ' OK')
        
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
        plt.savefig("%s%s-%s-utc" %(savepath,titulo_,tiempos[j]),format='png')
        print('Imagen guardada en %s'%savepath)

    if show == True: plt.show()
    data.close()
    
        
       

