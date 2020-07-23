from librerias import *


def buscarMin(var,mLat,mLon,window, stride ):

    if mLat + stride >= var.shape[0]:

        strideLat = var.shape[0]-mLat -1
        if strideLat <window: 
            window = strideLat


    elif mLat -stride <= window:

        strideLon = mLat
        if strideLat < window: 
            window = mLat

    else: strideLat = stride

    if mLon + stride >= var.shape[1]:

        strideLon = var.shape[1]-mLon -1
        if strideLon <window: 
            window = strideLon

    elif mLon -stride <= window:

        strideLon = mLon
        if mLon < window:
            window = mLat
    else: strideLon = stride
        
    minVar = np.min(var[mLat-strideLat:mLat+strideLat,mLon-strideLon:mLon+strideLon])

    varMin = (minVar, np.where(var==minVar))
    

    return varMin

def buscarMax(var,mLat,mLon,window, stride, strideLat = 10 ):

    if mLat + stride >= var.shape[0]:

        print('se pasa 1')
        strideLat = var.shape[0]-mLat -1
        if strideLat <window: 
            window = strideLat


    elif mLat -stride <= window:
        print('se pasa 2')
        strideLon = mLat
        if strideLat < window: 
            window = mLat

    else: strideLat = stride

    if mLon + stride >= var.shape[1]:
        print('se pasa 3')
        strideLon = var.shape[1]-mLon -1
        if strideLon <window: 
            window = strideLon

    elif mLon -stride <= window:
        print('se pasa 4')
        strideLon = mLon
        if mLon < window:
            window = mLat
    else: strideLon = stride
        
    maxVar = np.max(var[mLat-strideLat:mLat+strideLat,mLon-strideLon:mLon+strideLon])

    varMax = (maxVar, np.where(var==maxVar))
    

    return varMax

    
        

