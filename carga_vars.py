from librerias import *

import carga_datos as carga


def lluvia(data=False,tiempo = 0,lluvia_prev = False, lluviaName = 'RAINNC'):
    """
    Carga la lluvia en self.rain

    Si existe un archivo previo en self.archivoPrevio, 
    intenta cargar la lluvia previa para tener continuidad.
    """
    try:
        rain =  np.array(carga.varIndNC(data,lluviaName, tiempo = 'todo'))
        rain[1:] = rain[1:,...]-rain[:-1,...]
        
        
        if type(lluvia_prev) == bool:    
            print('No hay archivo de lluvia previa')
        else:
            rain[0] = rain[0,...]-lluvia_prev
            
        
        rain = np.ma.masked_array(rain, rain <0.5)
        
    except:
        print('Error al cargar la lluvia, se setea todo en 0')
        rain = np.zeros((305,305))
    try:
        print('carga de lluvia', rain.shape)
        return rain[tiempo]
    except:
        return np.zeros((48,305,305))

def topografia(data,topoName = 'HGT'):
    """
    Carga la topografía en self.topo y hace los cálculos
    para que este en unidades de presión
    """
    try:
        topo = carga.varIndNC(data,topoName)
        topo = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,np.array(topo))))/100
        print('Topografia ok')
        return topo
    except: 
        print('No se pudo cargar la topografía')    
        return
def tiempos(data, timeName = 'Times'):
    """
    Carga los tiempos y los concatena bien
    Los guarda en self.tiempos
    """
    try:    
        tiemposRaw = carga.varIndNC(data,timeName,tiempo='todo')
        tiempos = []
        
        for i in range(len(tiemposRaw)):
            tiempos.append(''.join(str(tiemposRaw[i]).replace('b','').replace("'",'').replace('[','').replace(']','').replace('\n','').replace(' ','')))


        return tiempos
    except:
        print('No se pueden cargar los tiempos')
def press(data, tiempo=0):
    """
    Carga las presiones de base y perturbadas y las suma en hPa
    Tambien genera la pmax y pmin para los gráficos
    Guarda en self.press, self.pmax y self.pmin
    """
    print(np.array(data.variables['P']).shape,type(data))

    pres = carga.varIndNC(data,'P',tiempo=tiempo) 
    presb = carga.varIndNC(data,'PB',tiempo=tiempo) 
    
    press = (pres[:58,:] + presb[:58,:])


    return press    
def z(data,tiempo=0):

    ph = carga.varIndNC(data,'PH',tiempo=tiempo)[:58,:]
    phb = carga.varIndNC(data,'PHB',tiempo = tiempo)[:58,:]

    z = (ph + phb)/9.8

    return z



def coords(data,latName = 'XLAT', lonName ='XLONG'):
    """
    Carga las coordenadas y el basemap.
    """
    try:    
        lat,lon = carga.varIndNC(data,latName), carga.varIndNC(data,lonName)

        bm = get_basemap(lat)
        lat,lon = to_np(lat),to_np(lon)
        
        x,y  = bm(lon,lat)
        return x,y,lat,lon,bm
    except:
        print('No se puede cargar la lat y lon')
def getXY(x,y,mLon,mLat,window=False,corte='xp'):
    
    """
    Setea las coordenadas según el gráfico
    """
    cortes = {
        'xz': setLatVert(x,mLon,mLat,window),
        'yz': setLonVert(x,mLon,mLat,window)
    }

    x = cortes[corte]
    y = y[:,mLat,mLon]


    x, y = np.meshgrid(x, y)   
    
    return x, y

def setLonVert(x, mLon,mLat, window=False):
    """
    Setea la longitud acotada por self.window
    """
    if window is not False:
        try:
            return x[mLat-window:mLat+window,mLon]
        except:
            return x[mLat-window:,mLon]
    else: 
        return x[:,mLon]
def setLatVert(y,mLon,mLat,window=False):
    """
    Setea la Lat acotada por self.window
    """
    if window is not False:
        try:
            return y[mLat,mLon-window:mLon+window]
        except:
            return y[mLat,mLon-window:]
    else: 
        return y[mLat,:]

def setXY(y,mLat,mLon,window=False):

    if not window:
        return y[:]
    else:
        try:
            return y[mLat-window:mLat+window,mLon-window:mLon+window]
        except:
            return y[mLat-window:,mLon-window:]



def setData3D(var,mLat=100,mLon=100,corte ='xp',window=False):

    """
    Setea los datos según las coordenadas
    """

    if corte == 'xz':
        if not window:
            return var[:,mLat,:]
        else:
            try:
                return var[:,mLat,mLon-window:mLon+window]
            except: 
                return var[:,mLat,mLon-window:]

    elif corte =='yz':
        if not window:
            return var[:,:,mLon]
        else:
            try:
                return var[:,mLat-window:mLat+window,mLon]
                
            except:
                return var[mLat-window:,mLon]
    elif corte == 'columna':

        return var [:,mLat,mLon]

    else:
        print('El corte está mal')    



def setData2D(var,mLon=100,mLat=100,corte ='xp',window=False):

    if corte == 'xz':
        return setLatVert(var,mLon,mLat,window)
    elif corte == 'yz':
        return setLonVert(var,mLon,mLat,window)
    else:
        print('El corte está mal')

def setTpotToT(tpot,press):

    tpot +=300

    t =  tpot*(press/1013)**0.2854


    return t

def getVort(u,v):

    [dqu_dx, dqu_dy] = np.gradient(u)
    [dqv_dx, dqv_dy] = np.gradient(v)

    div =dqv_dx- dqu_dy

    return  div
