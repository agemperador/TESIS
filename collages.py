from librerias import *

import mapas
import carga_datos as carga
import carga_vars as variables

def convectivo(data,i,x,y,bm,lat,lon,topo,tiempos,mLat,mLon,
                window = False,show=True,save= False, lluviaPrevia = False, 
                llueve=False, carpeta='',wMaxPrev=False, nombre= '', z =False):


    window = 15

    lluvia         = variables.lluvia(data,i,lluvia_prev=lluviaPrevia)

    if z:

        z = variables.z(data,tiempo=i)

        pmin = 0
        pmax = 14000

        press = z

    else:    
        press = variables.press(data,i)
        press = press/100
        pmin = 1000
        pmax = 100 


    w,u,v  = carga.varIndNC(data,'W',tiempo=i,lev=[1,-1]),carga.varIndNC(data,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)
    
    qc,qv,qi= carga.varIndNC(data,'QCLOUD',tiempo=i),carga.varIndNC(data,'QVAPOR',tiempo=i),carga.varIndNC(data,'QICE',tiempo=i)
    qs =carga.varIndNC(data,'QSNOW',tiempo=i)
    qg =carga.varIndNC(data,'QGRAUP',tiempo=i)
    qr =carga.varIndNC(data,'QRAIN',tiempo=i)

    tPot = carga.varIndNC(data,'T',tiempo=i)
    tPot = variables.setTpotToT(tPot,press)

    titaE = carga.varIndNC(data,'THM',tiempo=i)
    
    titaE  = np.ma.masked_array(titaE,press<750)
    



    pmin = 1000
    pmax = 100

    fig = plt.figure(figsize=(20,10))

    gspec = fig.add_gridspec(6,10)



    corte = 'xz'
    gs = gspec[0:2,0:4]

    mapas.convectivo(fig,lat,lon,press,mLat,mLon,window,u,w,tPot,titaE,topo,pmin,pmax,corte,gspec =gs )


    corte = 'yz'
    gs = gspec[2:4,0:4]

    mapas.convectivo(fig,lat,lon,press,mLat,mLon,window,v,w,tPot,titaE,topo,pmin,pmax,corte,gspec =gs )
    
    
    position = [[0,4],[0,6],[0,8],[1,4],[1,6],[1,8]]

    corte = 'xz'
    window = 15

    vpos = 0
    mapas.especies(fig,gspec,lon,lat,position,mLat,mLon,press,topo,w,qv,qc,qi,qr,qs,qg,tPot,corte,vpos,window,pmin,pmax)

    corte = 'yz'
    window = 15

    vpos = 2

    mapas.especies(fig,gspec,lon,lat,position,mLat,mLon,press,topo,w,qv,qc,qi,qr,qs,qg,tPot,corte,vpos,window,pmin,pmax)


    window = 30
    gs = gspec[4:,0:2]




    mapas.xyCentrado (fig,x,y,mLat,mLon,bm,window,topo,qv,lluvia,u,v,w,gspec=gs)


    gs = gspec[4:,2:4]
    mapas.xyTotal(fig,x,y,bm,mLon,mLat,u,v,topo,lluvia,gspec=gs)


    gs = gspec[4:,4:6]

    mapas.mapaHodografa(fig,u,v,mLat,mLon,press,gspec=gs)
    
        
    window = 15

    corte = 'yz'
    gs = gspec[4:,6:8]
    mapas.mapaMultiespecies(fig,lon,lat,press,mLat,mLon,corte,window,w,u,v,qv,qr,qi,qc,qs,qg,topo,pmin,pmax,gspec=gs)



    gs = gspec[4:,8:]

    mapas.mapaMesociclon(fig,x,y,mLat,mLon,window,press,u,v,topo,wMaxPrev,z,gspec=gs)


    try:
        plt.title('%s %s - Lat: %s - Lon: %s'%(carpeta,tiempos[i][4:-2], str(np.round(lat[mLat,mLon],2)), str(np.round(lon[mLat,mLon],2))),fontdict={'fontsize':20  })
    except:
        print(len(tiempos))

    plt.tight_layout()

    if show:
        plt.show()

    if save:
        print(i,'Guardando en ./img/%s/especies/TRACK/final/%s_%s_wmax_%s'%(carpeta,nombre,tiempos[i],carpeta))
        plt.savefig('./img/%s/especies/TRACK/final/%s_%s_%s.png'%(carpeta,carpeta,nombre,tiempos[i]))
        plt.close()


def convectivoLatLon(data,fig,i,press,lat,lon,mLat,mLon,window,topo,pmin,pmax,gs):
    
    w,u,v  = carga.varIndNC(data,'W',tiempo=i,lev=[1,-1]),carga.varIndNC(data   ,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)
       
    tPot = carga.varIndNC(data,'T',tiempo=i)
    tPot = variables.setTpotToT(tPot,press)

    titaE = carga.varIndNC(data,'THM',tiempo=i)
    
    titaE  = np.ma.masked_array(titaE,press<750)


    corte = 'xz'
    gsLon = gs[0]

    mapas.convectivo(fig,lat,lon,press,mLat,mLon,window,u,w,tPot,titaE,topo,pmin,pmax,corte,gspec =gsLon )

    corte = 'yz'
    gsLat = gs[1]

    mapas.convectivo(fig,lat,lon,press,mLat,mLon,window,v,w,tPot,titaE,topo,pmin,pmax,corte,gspec =gsLat )
    

    

def soloxyCentrado(data,fig,i,x,y,mLat,mLon,bm,window,topo,lluvia,press,gs):

    w,u,v  = carga.varIndNC(data,'W',tiempo=i,lev=[1,-1]),carga.varIndNC(data,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)
            
    qv=carga.varIndNC(data,'QVAPOR',tiempo=i)

    tPot = carga.varIndNC(data,'T',tiempo=i)
    tPot = variables.setTpotToT(tPot,press)



    window = 30

    

    mapas.xyCentrado (fig,x,y,mLat,mLon,bm,window,topo,qv,lluvia,u,v,w,gspec=gs,showLat=False)


def mesociclonInd(data,fig,i,x,y,mLat,mLon,window,press,topo,wMaxPrev=False,z=False,gspec=[1,1]):

    u,v  = carga.varIndNC(data,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)


    mapas.mapaMesociclon(fig,x,y,mLat,mLon,window,press,u,v,topo,wMaxPrev,z,gspec)
    

def hodografaInd(data,fig,i,mLat,mLon,press,gspec=[1,1]):

    u,v  = carga.varIndNC(data,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)

    mapas.mapaHodografa(fig,u,v,mLat,mLon,press,gspec=gspec)
    
        