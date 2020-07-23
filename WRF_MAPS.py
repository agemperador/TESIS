from librerias import *



import mapas

import carga_datos as cargar
import carga_vars as variables

import funciones_wrf as funciones

def sh(script):
    os.system("bash -c '%s'" % script)


#sh("echo $0")
#sh("ls -l")
#sh("echo done")

carpetas = ['control-ndg2', 'compl-ndg' ,'control','control-ndg']


#archivo = 'wrfout_d02_2018-11-10_12:00:00'


track = [   (24,179,121),(25,175,127),(26,171,137),
            (27,170,139),(28,170,146),(29,169,155),(30,169,165),
            (31,170,175),(32,170,185),(33,165,196),
            (34,162,199),(35,172,210),(36,169,225),
            (37,169,231),(38,168,237),(39,161,247),
            (40,175,257),(41,173,264),(42,177,268),
            (43,177,276),(44,174,283),(45,169,283),
            ]

dataTrack = pd.read_csv('/media/agustin/1D36257C5C6244CB/WRF_sim/Datos_Faltantes/System1.trj',sep=' ')
trackedTime = dataTrack['Time-step']
trackedLon = dataTrack['i']
trackedLat = dataTrack['j']


window = 20


#print(data.variables.keys())
#for i,archivo in enumerate(archivos[1:]):
for carpeta in [carpetas[1]]:

    print(carpeta)

    directorio = '/media/agustin/Linux/salidas_wrf/%s/'%carpeta

    #archivos = !ls /media/agustin/Linux/salidas_wrf/control-ndg2/

    archivos = glob.glob(directorio+'wrf*')

    print(archivos)

    #archivosUtiles = [archivos[0],archivos[2],archivos[3]]



    """
    for i, mLat, mLon in track:

        print(i,tiempos[i])

        mapas.convectivo(data,i,x,y,bm,lat,lon,topo,tiempos,mLat,mLon,window,show=False,save=True)
    """
    wMax = []

  

    #archivos = [archivos[3],archivos[2]]
    wMaxPrev = False

    for archivo in [archivos[3]]:

        print(archivo)

        data = cargar.dataNC(archivo)

        w = data.variables['W'][:,:,:,:]  

        #lluviaPrev = np.array(cargar.dataNC(archivos[0]).variables['RAINNC'])[-1,:-1,:-1]

        x,y,lat,lon,bm = variables.coords(data)
        tiempos        = variables.tiempos(data)

        if len(tiempos)< 2:
            continue

        topo           = variables.topografia(data)




        print(archivo)

        print(data.variables.keys())

        stride = 30

        #    i,mLat, mLon = wMaxTLL[0][0],wMaxTLL[1][0],wMaxTLL[2][0]

        #for i,mLat,mLon in track:

        #for i,mLon,mLat  in zip(trackedTime,trackedLon,trackedLat):

        mLat = 240
        mLon = 135

        for i in range(36,47,1):

            lluviaPrev = False

        
            #minCoords = funciones.buscarMin(wData[i,...],mLat,mLon,window,stride)

            #mLatMin = minCoords[1][0][0]
            #mLonMin = minCoords[1][1][0]

            #print(mLat,mLon,mLatMin,mLonMin)

            maxCoords = funciones.buscarMax(w[i,18,...],mLat,mLon,window,stride)

            mLat = maxCoords[1][0][0]
            mLon = maxCoords[1][1][0]

            print(i,mLat,mLon)

            #Centrado en los máximos ascensos
            mapas.convectivo(data,i,x,y,bm,lat,lon,topo,tiempos,mLat,mLon,window,
                                show=False,save=True,lluviaPrevia=lluviaPrev,
                                llueve=True,carpeta=carpeta,wMaxPrev = wMaxPrev,nombre='ascensos',z=False)

            wMaxPrev = (i,mLat,mLon)

            wMax.append(w[i,mLon,mLat])
            #Centrado en los máximos descensos
            """
            mapas.convectivo(data,i,x,y,bm,lat,lon,topo,tiempos,mLatMin,mLonMin,window,
                        show=True,save=False,lluviaPrevia=lluviaPrev,
                        llueve=True,carpeta=carpeta,wMaxPrev = wMaxPrev,nombre = 'descensos',z=False)
            """
        
        plt.plot(wMax)
        plt.savefig ('./img/wMax.png')