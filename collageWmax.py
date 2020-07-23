
from WRF2 import *

from os import listdir
from os.path import isfile, join

directorio = '/media/agustin/Linux/salidas_wrf/control-ndg2/'


archivo = 'wrfout_d02_2018-11-10_12:00:00'

max_lat = 157
max_lon = 294


i = 25

dataPath =  directorio + archivo

data = Dataset(dataPath,'r')

wLevel = 18

wData = np.array(data.variables['W'][:,wLevel,:,:])

wMaxCord = np.where(wData == np.max(wData))
print(wData[wMaxCord],wMaxCord[0],wMaxCord[1])

maxLat, maxLon = wMaxCord[0],wMaxCord[1]

#max_lat,max_lon = buscar_max(data, max_lat=max_lat,max_lon=max_lon,time = i)

mapa = plots()



#mapa.mapaCollage(data,i,coords=[maxLat,maxLon],show=True,save=False, carpeta='./img/control-ndg2/especies/SC1')

