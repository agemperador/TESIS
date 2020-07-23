from librerias import *

import mapas

import ejes
import carga_datos as cargar
import carga_vars as variables

def getTrayectorias(otFile='/media/agustin/1D36257C5C6244CB/WRF_sim/Track_Maite/retrayectoriasrelampago/OT1_15min.txt'):
    ot = pd.read_csv(otFile,sep='\t',header=1)
    
    time = ot['Date']
    lon = ot['Centroid Longitude']
    lat = ot['Centroid Latitude']


    return time, lon,lat


if __name__ ==  '__main__':

    track = [   (24,179,121),(25,175,127),(26,171,137),
                (27,170,139),(28,170,146),(29,169,155),(30,169,165),
                (31,170,175),(32,170,185),(33,165,196),
                (34,162,199),(35,172,210),(36,169,225),
                (37,169,231),(38,168,237),(39,161,247),
                (40,175,257),(41,173,264),(42,177,268),
                (43,177,276),(44,174,283),(45,169,283),
                ]

    print('archivo')
    tyTime,tyLon, tyLat = getTrayectorias()

    tyTime = np.asarray([time[:-1].replace(' ','_') for time in tyTime])


    carpeta ='control-ndg'
    directorio = '/media/agustin/Linux/salidas_wrf/%s/'%carpeta


    archivo = glob.glob(directorio+'wrf*')[3]

    data = cargar.dataNC(archivo)

    x,y,lat,lon,bm = variables.coords(data)
    tiempos        = np.asarray(variables.tiempos(data))


    indexTimes = [i for i,j,k in track if tiempos[i] in tyTime] 
    comparartiveTrack = [(i,j,k) for i,j,k in track if tiempos[i] in tyTime] 

    indexLat = [np.argwhere(np.round(lat,4)==j) for j in np.round(tyLat,4) if j in np.round(lat,4)]
    indexLon = [np.argwhere(np.round(lon,4)==j) for j in np.round(tyLon,4) if j in np.round(lon,4)]

    # usar distancia!

    print(indexLon,indexLat)

    mapa = np.zeros_like(lat)
 
    #for j,k in zip(tyLon, tyLat):


    #plt.pcolormesh(lon,lat,mapa)
    #plt.show()

    
    



