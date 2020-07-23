from netCDF4 import Dataset

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LightSource
from matplotlib.cm import get_cmap
from matplotlib.patches import Rectangle, Polygon


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
def carga_var(data,name_var,time,lev_var,ini = 0,end = 305):
    if name_var != '0':
        if (type (name_var)== str) and (type(time) == int):

            var = getvar(data,name_var, timeidx = time)

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

def carga_v_rain(data, i,lev_viento = 0, lev_rich = 15, richardson = False,
          viento = False,sep_barb = 2, lluvia = False, topografia = False,
          ini =0,end = 305 , var_lluvia = 'RAINNC'): 
    """""""""""""""""""""""""""""""""VIENTO"""""""""""""""""""""""""""""""""
    if viento == True:
        v_x = getvar(data, 'U', timeidx= i ) [lev_viento,ini:end,ini:end]
        v_y = getvar(data, 'V', timeidx= i) [lev_viento,ini:end,ini:end]
        v = [v_x,v_y]
        
    else: v = 0
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""  
    """""""""""""""""""""""""""""LLUVIA"""""""""""""""""""""""""""""
    if lluvia == True:

        raini = getvar(data,var_lluvia,timeidx = i)[ini:end,ini:end]
        rain_ = getvar(data,var_lluvia,timeidx = i-1)[ini:end,ini:end]

        rain = raini-rain_
        
    else: rain = 0

    if richardson == True:
        if viento == False:
            v_x = getvar(data, 'U', timeidx= i ) [lev_viento,ini:end,ini:end]
            v_y = getvar(data, 'V', timeidx= i) [lev_viento,ini:end,ini:end]
        v_x,v_y = to_np(v_x), to_np(v_y)
        cape = np.array(getvar(data,'AFWA_CAPE', timeidx = i))[ini:end,ini:end]
        #cape = np.expand_dims(cape,axis = 0)
        #cape = np.repeat(cape,59,axis = 0) sirve para repetir la matriz en un tercer eje
        v_media = v_x**2+v_y**2
        ri = cape/v_media
        print('Richardson activado')

    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    ri = 0
    return to_np(v),rain, ri

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
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
def set_var (var):
    if len(var.shape) == 4: var_set = to_np(var[i,...])
    elif len(var.shape) == 3: var_set = to_np(var[i,...])
    else: var_set = to_np(var)
    return var

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def sombreado (x,y, var_somb,ax, ini = 0, end = 305, cmap_s = 'coolwarm', zorder = 0, k=1):
    
    var_s = set_var(var_somb)
    ps = ax.contourf(x, y, var_s[ini:end,ini:end]*k, cmap = plt.get_cmap(cmap_s),zorder = zorder)

    return ps

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def contorno (x,y, var_cont,ax, ini = 0, end = 305, cmap_c = 'hot_r', zorder = 1,ncont = 12):
    
    var_c = set_var(var_cont)
    
    c_min, c_max =np.min(var_c),np.max(var_c)

    clev=np.linspace(c_min,c_max,ncont)

    pc = ax.contour(x,y, var_c[ini:end,ini:end],clev, cmap = plt.get_cmap(cmap_c),zorder = zorder)
    clab = plt.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')

    return pc
 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_viento(x,y,v,bm,ini=0,end=305,cmap_viento = 'no',zorder = 2):
                
    
    """""""""""""""CMAP VIENTO"""""""""""""""
    if cmap_viento == 'no':
        cdict = {
        'blue': [(0,1,0.5),(0.4,0.5,0.4) , (1, 1, 1)],
        'green': [(0,0,0.2),(0.4,0.2,0), (1, 0, 0)],
        'red': [(0,0,0) ,(0.4,0,0.4), (1, 1, 1)]}
        cmap_viento = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,10)
    """"""""""""""""""""""""""""""""""""""""""
    x_v = x[0,:]
    y_v = y[:,0]

    v_x = to_np(v[0])[ini:end,ini:end]
    v_y = to_np(v[1])[ini:end,ini:end]
    
    speed = np.sqrt((v_x*2) ** 2 + (v_y*2) ** 2)
    lw = 3*speed / speed.max()

    vi= bm.streamplot(y_v, x_v, v_x, v_y, color = speed , density=1.5
                     ,cmap = cmap_viento , linewidth= lw, arrowsize = 1.5,
                         arrowstyle = '->', zorder = zorder)
        
    
    return vi

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_precip(x,y, rain, bm, ini=0,end = 305, cmap_lluvia='raibow', zorder = 3):
    
    r = np.ma.array(rain, mask= rain < 0.5)
    ra = bm.contourf(x,y,r[ini:end,ini:end], cmap = plt.get_cmap(cmap_lluvia),zorder = zorder)
    return ra    

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##############################################################################
##############################################################################
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_recuadro(x,ax = False, ancho= 0.1,alto=0.1):
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
                   archivo = 'wrfout_d02_2018-11-10_12:00:00',
                   ini = 0, end = 100,rows = 2,cols = 2,richardson = False,
                   var_s = '0',  ts = [0],   lev_s = 0  , us = 0 , cmap_s = 'coolwarm', 
                   var_c1 = '0', lev_c1 = 0 , uc1 = 0, cmap_c1 = 0,
                   var_c2 = '0', lev_c2 = 0 , uc2 = 0, cmap_2 = 0,
                   viento = True, lev_viento = 0,sep_barb = 1,cmap_viento = 'no',
                   titulo = 'Sin titulo definido', guardar = False,
                   topografia = True, topo_mode = 'soft',
                   lluvia = True, cmap_lluvia = 'gist_ncar',
                   lev_rich = 15
                   ):   
    """
    EJEMPLO:
    
    mapa_somb_cont2( archivo='wrfout_d02_2018-11-10_12:00:00' , var_s = 'QVAPOR',lev_s=3,var_c1='TH2',
                    guardar=False,ini = 0, end = 305,rows = 2, cols = 2,
                    ts = [0,11,23,35],tc1 = 10, lluvia= True, topografia=True,
                    titulo='$\Theta$ (cont - K - 2m) - Qv (somb - kg/kg - 2m) V sup - Lluvia (mm/hr)',
                    cmap_s = 'BrBG',cmap_c1='hot_r', viento = True,
                    cmap_lluvia = 'rainbow')
    """               

    """CARGA DE DATOS"""
    data = Dataset(path+archivo,'r')

    lats, lons = getvar(data,'XLAT'), getvar(data,'XLONG')
    topo = np.asarray(getvar(data,'HGT'))
    times =  data.variables['Times']
    
    tiempos = []
    for i in ts:
        tiempos.append( ''.join(str(times[i]).replace('b','').replace("'",'').replace('[','').replace(']','').replace('\n','').replace(' ','')))
    print ('Mapa con límites lat: %f - %f, lon: %f - %f. \n Variables: %s somb, %s cont 1, %s cont 2,  %s. \n Tiempos: '
          %(lats[ini,ini], lats[end,end], lons[ini,ini], lons[end,end], 
            info(var_s,lev_s), info(var_c1,lev_c1),info(var_c2,lev_c2),
            info(viento = viento, lev = lev_viento, lluvia = lluvia, topografia = topografia)), tiempos)
    
    # Cargo el basemap de las latitudes
    bm = get_basemap(lats[ini:end,ini:end])
        
    x, y = bm(to_np(lons)[ini:end,ini:end], to_np(lats)[ini:end,ini:end])
    

    if len(ts)!= rows*cols:
        if np.sqrt(len(ts)) == round(np.sqrt(len(ts))):
            rows = int(np.sqrt(len(ts)))
            cols = rows
        else: 
            rows = int(np.sqrt(len(ts)))+1
            cols = rows
                
    #Guardo los nombres de las variables
    n_var_s, n_var_c1, n_var_c2 = var_s, var_c1, var_c2
    
    lozada = [-31.651943,-64.07947] #Lat Lon de Lozada, Córdoba
    
    lats, lons = np.array(lats),np.array(lons)
    
    print('Lozada: ',lats[lats == lozada[0]], lons[lons == lozada[1]])

    ###############################
    ######  GENERO LA IMAGEN ######
    ###############################
    
    fig = plt.figure(figsize= (6,6))#,constrained_layout=True)
    axx = []
    for i, j in zip(ts,range(len(ts))):
        
        ## Un subplot por cada mapa uso las rows cols y la enumeración j
        ax = fig.add_subplot(int(rows*100+cols*10 + j+1))
        
        
        #Aca guardo los subplots
        axx.append(ax)
        
        """""""""""""""FUNCION DE CARGA DE DATOS"""""""""""""""
        ## Esta carga las cosas raras
        v,rain, ri = carga_v_rain(data = data, i = i,lev_viento = lev_viento, lev_rich=lev_rich, richardson= richardson,
                              viento = viento, sep_barb = sep_barb, lluvia = lluvia, ini =ini,end = end)
        
        ## Esta carga las cosas comunes
        smooth_var_s  = carga_var(data,n_var_s,i ,lev_s,ini,end)
        smooth_var_c1 = carga_var(data,n_var_c1,i,lev_c1,ini,end)
        smooth_var_c2 = carga_var(data,n_var_c2,i, lev_c2,ini,end)
        
        ## Me genero los bools para despues ver que grafico y que no
        if type(smooth_var_s) == int: somb = False
        else: somb = True
        if type(smooth_var_c1) == int: cont1 = False
        else: cont1 = True
        if type(smooth_var_c2) == int: cont2 = False
        else: cont2 = True
        
        

        """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        """""""""ACA EMPIEZA EL MAPA"""""""""
        # DEINO DIA Y HORA
        l=i%24
        d=int(i/24)+12 

        """""""""""""""GRAFICO SOMBREADO"""""""""""""""
        if somb == True:
           
            ps = sombreado(x,y, smooth_var_s,ax,ini =ini, end = end,cmap_s= cmap_s)
                        
            ###TOPOGRAFÍA EN CONTORNO###
        if topografia == True: 
            ax.contour(x,y, topo,cmap = 'gray_r', linewidths = 2, zorder = 0)
                
        """""""""""""""""""""GRAFICO CONTORNO 1"""""""""""""""""""""
        if cont1 == True:
            
            pc1 = contorno(x,y, smooth_var_c1,ax, ini =ini, end = end, cmap_c = cmap_c1, zorder = 1,ncont = 12)
            
        """""""""""""""""""""GRAFICO CONTORNO 2"""""""""""""""""""""
        if cont2 == True:
            
            pc2 = contorno(x,y, smooth_var_c2,ax, ini =ini, end = end, cmap_c = cmap_c2, zorder = 1,ncont = 12)

        """""""""""""""""""""""""""GRAFICO VIENTO"""""""""""""""""""""""""""
        if viento == True:
            
            vi = plot_viento(x,y,v,bm,ini,end,cmap_viento=cmap_viento,zorder= 4)
            
            
        """""""""""""""""""""""""""""LLUVIA"""""""""""""""""""""""""""""
        if lluvia == True:
            
            ra = plot_precip(x,y, rain, bm, ini,end, cmap_lluvia, zorder = 3)

        if richardson == True:

            p_rich = contorno(x,y,ri,ax,cmap_c='viridis')

        
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
        
    if guardar==True: plt.savefig("img/"+titulo_+"%s-utc" %(tiempos[j]),format='png')

    ##################################
    #Esto es para que quede acomodado#
    ##################################
    fig.tight_layout()
    set_size(12,10)
    ##################################
    
    ##################################
    ###"""ACA VAN LAS COLORBARS """###
    ##################################
    if somb == True: 
        caxv = plt.axes([0.98, 0.05, 0.02, 0.85])
        cbars = fig.colorbar(ps, ax=axx,cax=caxv, orientation = 'vertical')
        cbars.ax.set_ylabel(n_var_s + ' (%s)' %data.variables[n_var_s].units,fontsize = 12)
    if viento == True:
        caxh = plt.axes([0.05, -0.1, 0.9, 0.02])
        cbarv = fig.colorbar(vi.lines,ax=axx,cax=caxh, orientation = 'horizontal')    
        cbarv.ax.set_xlabel('Velocidad del viento (%s)' %data.variables['U'].units, fontsize = 12)
    if lluvia == True:
        caxh2 = plt.axes([0.05, -0.02, 0.9, 0.02])
        cbar = fig.colorbar(ra, ax = axx, cax = caxh2 ,orientation = 'horizontal')
        cbar.ax.set_xlabel('Lluvia (%s)' %data.variables['RAINNC'].units, fontsize = 12)
    ###################################
    
    ## El titulo de toda la imagen
    fig.suptitle(titulo+ '\n',fontsize = 16, y = 1)
    plt.show()
    data.close()
    
        
