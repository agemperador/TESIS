from librerias import *

from plots import *


class plots():
    def __init__(self):
        super().__init__()
  

    def especies(self,especies = [],data = False, archivo = '' ):

        """
        Estructura los datos necesarios para el gráfico de especies
        individuales

        Necesita los nombres de las especies
        La Data (ya sea por parámetro o por archivo)
        El nombre del archivo donde esta la data

        Genera la lista de variables  con cargaVarInd
        """

        self.cargaData(archivo,data)
        print('Datos cargados')

        self.cargaVar()

        listaEspecies = []
        for esp in especies:

            listaEspecies.append(self.cargaVarInd(esp))
        

              

        return listaEspecies

    def cargaData(self,archivo, data= False):      

        """
        Carga los datos, ya sea por parametro o por archivo
        Lo mete en self.data 

        """
        print(archivo)
        if data == False:

            try:
                self.data = Dataset(archivo,'r')
                self.archivo = archivo
            except:
                try:
                    self.data = Dataset(self.archivo,'r')
                except:
                    try: 
                        print(self.data)
                    except:
                        print('No se encontro ningun archivo %s'%archivo)
                        return

        else: 
            try:
                self.data = Dataset(archivo,'r')
            except: 
                try:
                    self.data = Dataset(self.archivo)
                except:
                    print('no se encontro ningún archivo')

    def cargaVar(self):
        """
        Carga las variables generales del WRF
        Coordenadas, Tiempos en texto
        Topografía, Lluvia, Presion
        Y las variables de los mapas window y stride

        """
        

        self.cargarCoord()

        self.cargaTiempos()
        self.cargaTopografia()
        self.cargaLluvia()
        self.cargaPress()

        self.window = 20
        self.stride = 11
        self.stride_lon = 16

    def cargaLluvia(self):
        """
        Carga la lluvia en self.rain

        Si existe un archivo previo en self.archivoPrevio, 
        intenta cargar la lluvia previa para tener continuidad.
        """

        try:
            rain =  np.array(self.data.variables['RAINNC'])
            rain[1:] = rain[1:,...]-rain[:-1,...]
            try:    
                data_prev = Dataset(self.archivoPrevio,'r')
                rain[0] = rain[0,...]-np.array(data_prev.variables['RAINNC'][-1,...])
                data_prev.close()
            except:
                print('No hay archivo de lluvia previa')

            self.rain = np.ma.masked_array(rain, rain <2)
            print('lluvia ok')
        except:
            print('Error al cargar la lluvia')
            self.rain = np.zeros((len(self.times),305,305))
    def cargaTopografia(self):
        """
        Carga la topografía en self.topo y hace los cálculos
        para que este en unidades de presión
        """

        topo = np.array(getvar(self.data,'HGT'))
        self.topo = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,topo)))/100
        print('Topografia ok')
    def cargaTiempos(self):
        """
        Carga los tiempos y los concatena bien
        Los guarda en self.tiempos
        """
        tiempos = self.data.variables['Times'][:]
        self.tiempos =[]
        for i in range(len(tiempos)):
            self.tiempos.append(''.join(str(tiempos[i]).replace('b','').replace("'",'').replace('[','').replace(']','').replace('\n','').replace(' ','')))

        print('Times ok')
    def cargaPress(self):
        """
        Carga las presiones de base y perturbadas y las suma en hPa
        Tambien genera la pmax y pmin para los gráficos
        Guarda en self.press, self.pmax y self.pmin
        """
        pres = np.array(self.data.variables['P']) 
        presb = np.array(self.data.variables['PB'])
        self.press = (pres + presb)/100
        self.pmax, self.pmin = 966,100
        print('Press ok')
    def cargarCoord(self):
        """
        Carga las coordenadas y el basemap.
        """
        try:    
            lat,lon = getvar(self.data,'XLAT'), getvar(self.data,'XLONG')
            self.bm = get_basemap(lat)
            self.lat,self.lon = to_np(lat),to_np(lon)
            self.x, self.y  = self.bm(self.lon,self.lat)
        except:
            print('No se puede cargar la lat y lon')

    def cargaVarInd(self, var=''):
        """
        Carga variables segun su nombre

        """

        dataVar = {
            'v': 'QVAPOR',
            'r': 'QRAIN',
            'c': 'QCLOUD',
            'i': 'QICE',
            's': 'QSNOW',
            'w': 'W',
            'u': 'U',
            'v': 'V'
        }
        
        try:
            
            var_name = dataVar[var]
            try:
                if var_name == 'W':
                    variable = Dataset(self.archivo,'r')[var_name][self.i,:-1,...]
                else:

                    variable = Dataset(self.archivo,'r')[var_name][self.i,:,...]
                return variable

            except:
                print('No se puede cargar %s de %s'%(var,self.archivo))

                
        except:  
            print('%s debe ser uno de estos %a'%(var,list(dataVar.keys())))
     

    def mapaEspecies(self,archivo,i, coords = [100,100], save = False, show = True,carpeta = ''):
        """
        Genera el mapa de especies

        """
        
        self.i = i

        listaEspecies = ['r','v','i','c','s']   

        especies = self.especies(especies= listaEspecies, archivo = archivo)

        w = self.cargaVarInd('w')
        u = self.cargaVarInd('u')
        v = self.cargaVarInd('v')
        
        """
        Genera un gráfico múltiple con gs

        """

        gs = [5,2]


        fig = plt.figure(figsize=(10,10))

        plt.title('Especies corte lat(izq) y lon(der) en %i,%i t= %i'%(coords[0],coords[1],i))


        gspec = fig.add_gridspec(gs[0],gs[1])  

        
        for i in range(len(especies)):

            espVar = especies[i]
            espName = listaEspecies[i]

            self.ejeEspecies(espVar,[u,w],v,fig,gs=gspec[i,0],corte='xp',mLat=coords[0],mLon=coords[1], especie = espName)


            self.ejeEspecies(espVar,[v,w],u,fig,gs=gspec[i,1],corte='yp',mLat=coords[0],mLon=coords[1], especie = espName)

        if save== True:
            plt.savefig('%s/especies_collage_%s'%(carpeta,str(self.i)))
            

        if show == True:
            plt.show()
        else:
            plt.close(fig)
        

    def mapaCollage(self,archivo,i, coords = [100,100], save = False, show = True,carpeta = ''):
        """
        Genera el mapa de especies

        """
        """
        Estructura los datos necesarios para el gráfico de especies
        individuales

        Necesita los nombres de las especies
        La Data (ya sea por parámetro o por archivo)
        El nombre del archivo donde esta la data

        Genera la lista de variables  con cargaVarInd
        """

        self.cargaData(archivo)
        print('Datos cargados')

        self.cargaVar()
        
        self.i = i

        w = np.array(self.cargaVarInd('w'))
        w = np.ma.masked_array(w, abs(w) < 1)
        u = self.cargaVarInd('u')
        v = self.cargaVarInd('v')
        
        """
        Genera un gráfico múltiple con gs

        """

        gs = [1,1]


        fig = plt.figure(figsize=(10,10))

        plt.title('Mapas W')


        gspec = fig.add_gridspec(gs[0],gs[1])  

        

        self.ejeCorteVert(somb = w,fig = fig, gs = gspec[0,0],grafico='w')


        if save== True:
            plt.savefig('%s/especies_collage_%s'%(carpeta,str(self.i)))
            

        if show == True:
            plt.show()
        else:
            plt.close(fig)
        

    def ejeCorteVert(self, somb='',viento ='',cont ='',fig='',gs='',mLat =150, mLon =150,corte='xp',grafico='r'):

        ax = fig.add_subplot(gs)

        x2D,press2D = self.setXY(mLat = mLat,mLon = mLon,corte = corte)
        

        if type(somb)!=str:
            
            somb = self.setData(somb,mLat = mLat,mLon = mLon,corte= corte)
            self.graficos(ax,x2D,press2D,somb,plot = 'w')
        if type(viento)!=str:

            u = viento[0]
            w = viento[1]
            w = self.setData(w,mLat = mLat,mLon = mLon,corte= corte)
            u = self.setData(u,mLat = mLat,mLon = mLon,corte= corte)
            self.graficos(ax,x2D,press2D,[u,w],plot = 'uw')

        if type(cont) !=str:

            cont = self.setData(cont,mLat = mLat,mLon = mLon,corte= corte)
            self.graficos(ax,x2D,press2D,cont,plot='cont')



        topo = self.setData(self.topo,mLat = mLat,mLon = mLon,corte= corte)

        
        


        self.plotTopoVert(ax,x2D[0,:],topo)

        plt.ylim((self.pmax,self.pmin))
                


        
        
 


            
    def ejeEspecies(self, somb='',viento ='',cont ='',fig='',gs='',mLat =150, mLon =150,corte='xp',especie='r'):

        ax = fig.add_subplot(gs)

        x2D,press2D = self.setXY(mLat = mLat,mLon = mLon,corte = corte)
        

        if somb!='':
            
            somb = self.setData(somb,mLat = mLat,mLon = mLon,corte= corte)
            self.graficos(ax,x2D,press2D,somb,plot = especie)
        if viento!='':

            u = viento[0]
            w = viento[1]
            w = self.setData(w,mLat = mLat,mLon = mLon,corte= corte)
            u = self.setData(u,mLat = mLat,mLon = mLon,corte= corte)
            self.graficos(ax,x2D,press2D,[u,w],plot = 'uw')

        if cont !='':

            cont = self.setData(cont,mLat = mLat,mLon = mLon,corte= corte)
            self.graficos(ax,x2D,press2D,cont,plot='cont')



        topo = self.setData(self.topo,mLat = mLat,mLon = mLon,corte= corte)

        
        


        self.plotTopoVert(ax,x2D[0,:],topo)

        plt.ylim((self.pmax,self.pmin))
                
    
    def graficos(self,ax,x,y,var,plot='r'):

        """
        Genera los gráficos pertinentes a cada variable
        """

        varName = plot
        
        plots = {
                'r':{'cmap':'Greens','plot':pmesh},
                'v':{'cmap':'Purples','plot':pmesh},
                'i':{'cmap':'Reds','plot':pmesh},
                's':{'cmap':'Greys','plot':pmesh},
                'c':{'cmap':'Blues','plot':pmesh},
                'w':{'cmap': 'seismic','plot':pmesh},
                'uw':{'cmap': 'k','plot':pbarbs},
                'cont':{'cmap':'k','plot':pcont},
                
            

        }
        cmap = plots[varName]['cmap']

        plots[varName]['plot'](ax,x,y,var,cmap)
    


    def setXY(self,mLon=100,mLat=100,corte='xp'):

        """
        Setea las coordenadas según el gráfico
        """
        cortes = {
            'xp': [
                    self.setLonVert(mLon,mLat),
                    self.press[self.i,:,mLon,mLat]
            ],
            'yp':[
                    self.setLatVert(mLon,mLat),
                    self.press[self.i,:,mLon,mLat]
            ]
        }
        x = cortes[corte][0]
        y = cortes[corte][1]


        x, y = np.meshgrid(x, y)   
        
        return x, y
    def setData(self,var,mLon=100,mLat=100,corte ='xp'):

        """
        Setea los datos según las coordenadas
        """

        if corte == 'xp':
            try:
                return var[:,mLon-self.window:mLon+self.window,mLat]
            except:
                return var[mLon-self.window:mLon+self.window,mLat]
            
        elif corte =='yp':
            try:
                return var[:,mLon,mLat-self.window:mLat+self.window]
            except:
                return var[mLon,mLat-self.window:mLat+self.window]
        else:
            print('El corte está mal')    
        

    

    def setLonVert(self, mLon=100,mLat=100):
        """
        Setea la longitud acotada por self.window
        """
        try:
            return self.lon[mLon,mLat-self.window:mLat+self.window]
        except:
            return self.lon[mLon,mLat-self.window:]
    def setLatVert(self,mLon=100,mLat=100):
        """
        Setea la Lat acotada por self.window
        """
        try:
            return self.lat[mLon-self.window:mLon+self.window,mLat]
        except:
            return self.lat[mLon-self.window:,mLat]
            
    def plotTopoVert(self,ax, x,topo):
        #Plot de topografía
        ax.plot(x,topo, linewidth = 5, color='k', zorder = 3)


def buscar_max(data,max_lat=-1,max_lon=-1,time=0,level=18,window=20,stride=11,stride_lon=16):

    w = Dataset(data,'r').variables['W']

    wsample = np.array(w[time,level,...])

    if max_lat == -1 or max_lon == -1: 
        
        maximo = np.argwhere((wsample[:,:] == wsample[:,:].max()))[0]

        max_lat_pr = 0
        max_lon_pr = 0
    else:
        if (max_lat+window>=304):
      
            stride = 0
        if (max_lon+window>=304):
      
            stride_lon=0
        wminisample = wsample[max_lat-stride:max_lat+stride,max_lon-stride_lon:max_lon+stride_lon]
        
        maximo = np.argwhere(wminisample == wminisample.max())[0]

        max_lat_pr = max_lat - stride
        max_lon_pr = max_lon - stride_lon 

        

    max_lat =  maximo[0] + max_lat_pr
    max_lon =  maximo[1] + max_lon_pr

    if max_lon + window >= w.shape[2]: 
        max_lat = 0
        max_lon = 0
        
    
        
    
    return    max_lat,max_lon



        






        

