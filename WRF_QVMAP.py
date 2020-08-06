from librerias import *

import mapas

import carga_datos as cargar
import carga_vars as variables

import funciones_wrf as funciones

import os
def sh(script):
    os.system("bash -c '%s'" % script)

import glob



carpetas = ['control-ndg2', 'compl-ndg' ,'test01','control-ndg']


for carpeta in carpetas:

    print(carpeta)

    directorio = '/media/agustin/Linux/salidas_wrf/%s/'%carpeta

    archivos = glob.glob(directorio+'wrf*')

    print(archivos)

    if carpeta == 'test01': carpeta = 'control'
    for j,archivo in enumerate(archivos):

        
        print(archivo)

        data = cargar.dataNC(archivo)

        #lluviaPrev = np.array(cargar.dataNC(archivos[0]).variables['RAINNC'])[-1,:-1,:-1]
        lluviaPrev = False
        
        x,y,lat,lon,bm = variables.coords(data)
        tiempos        = variables.tiempos(data)
        topo           = variables.topografia(data)

        if len(tiempos)< 2:
            continue

        for i in range(len(tiempos)):

            print(tiempos[i])

            if i ==0:
                dataPrev = cargar.dataNC(archivos[j-1])
                lluviaPrev = data.variables['RAINNC'][-1,:,:] 
                print(lluviaPrev.shape)

            cdict = {
            'blue': [(0,1,0.5),(0.4,0.5,0.4) , (1, 1, 1)],
            'green': [(0,0,0.2),(0.4,0.2,0), (1, 0, 0)],
            'red': [(0,0,0) ,(0.4,0,0.4), (1, 1, 1)]}
            cmapViento = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,10)

            mapas.qvapor(data,i,x,y,bm,lat,lon,topo,tiempos,
                                show=True,save=False,lluviaPrevia=lluviaPrev,
                                llueve=True,carpeta=carpeta,nombre='qvapor',z=False,cmapViento=cmapViento)