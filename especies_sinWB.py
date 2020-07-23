
from WRF2 import *

from os import listdir
from os.path import isfile, join

directorio = '/media/agustin/Linux/salidas_wrf/control-ndg2/'


archivo = 'wrfout_d02_2018-11-10_12:00:00'

max_lat = 168
max_lon = 62

for i in range(17,37):
    data =  directorio + archivo

    max_lat,max_lon = buscar_max(data, max_lat=max_lat,max_lon=max_lon,time = i)

    mapa = plots()

    mapa.mapaEspecies(data,i,coords=[max_lat,max_lon],show=False,save=True, carpeta='./img/control-ndg2/especies/SC1')




"""
folders = [f for f in listdir(directorio) if not isfile(join(directorio, f))]



for folder in folders:
    
    directorioActual = join(directorio,str(folder+'/'))

    print(directorioActual)
    if folder == 'control' or folder == 'control-ndg2' or folder == 'compl-ndg':
        files = [f for f in listdir(directorioActual)]
                
        for archivo in files:
            print(str(archivo))
   
"""