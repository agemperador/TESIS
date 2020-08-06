from librerias import *

import mapas
import collages
import carga_datos as carga
import carga_vars as variables

import funciones_wrf as funciones

import ejes

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


for carpeta in [carpetas[3]]:



    mLat = 240
    mLon = 132



    #Desarrollo compl-ndg
    mLat = 236
    mLon = 142

    #splitting compl-ndg
    SC1COMPLNDG_SPLITTING ={
        'carpeta':'compl-ndg',
        'mLat': 233,
        'mLon': 183,
        'ti': 41
    }

    #sc1 control
    mLat = 160
    mLon = 130 
    ti = 12
    
    #CI compl-ndg
    sc1ComplNdg = {
        'carpeta':'compl-ndg',
        'mLat': 253,
        'mLon': 123,
        'ti':35
    }

    #sc1 control
    sc1Control = {
        'carpeta':'control',
        'mLat': 175,
        'mLon' : 125,
        'ti' : 20
    }

    #sc1 control-ndg
    sc1ControlNdg ={
        'carpeta':'control-ndg',
        'mLat':180,
        'mLon': 120,
        'ti': 23
    }

    sc1ControlNdg2 ={
        'carpeta':'control-ndg2',
        'mLat':179,
        'mLon':121,
        'ti':23
    }


    coords = sc1ControlNdg


    show = False
    save = True

    lluviaPrev = False
    z=False

    qvmap =  True
    conv = False
    meso = False
    hodor = False



    carpeta = coords['carpeta']

    if carpeta == 'control': 
        carpeta = 'test01'
        carpetaName = 'control'
    else:
        carpetaName = carpeta


    print(carpeta)

    directorio = '/media/agustin/Linux/salidas_wrf/%s/'%carpeta

    #archivos = !ls /media/agustin/Linux/salidas_wrf/control-ndg2/

    archivos = glob.glob(directorio+'wrf*')

    print(archivos)

    wMax = []

  
    wMaxPrev = False

    for archivo in [archivos[3]]:

        print(archivo)

        data = carga.dataNC(archivo)

        w = data.variables['W'][:,:,:,:]  

        #lluviaPrev = np.array(cargar.dataNC(archivos[0]).variables['RAINNC'])[-1,:-1,:-1]

        x,y,lat,lon,bm = variables.coords(data)
        tiempos        = variables.tiempos(data)

        if len(tiempos)< 2:
            continue

        topo           = variables.topografia(data)

        print(archivo)

        print(data.variables.keys())

        stride = 20



        ti = coords['ti']
        mLat = coords['mLat']
        mLon = coords['mLon']

        

        for i in range(ti,47,1):

            print(w.shape)
            
            maxCoords = funciones.buscarMax(w[i,18,...],mLat,mLon,window,stride)

            mLat = maxCoords[1][0][0]
            mLon = maxCoords[1][1][0]

            print(i,mLat,mLon)

            lluvia         = variables.lluvia(data,i,lluvia_prev=lluviaPrev)

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


            print(i,mLat,mLon)

            if qvmap:
                
                fig = plt.figure(figsize=(10,10))

                gspec = fig.add_gridspec(1,1)

                gs = gspec[0,0]



                collages.soloxyCentrado(data,fig,i,x,y,mLat,mLon,bm,window,topo,lluvia,press,gs)
                
                plt.title('%s %s - Lat: %s - Lon: %s'%(carpetaName,tiempos[i][4:-2], str(np.round(lat[mLat,mLon],2)), str(np.round(lon[mLat,mLon],2))),fontdict={'fontsize':20  })
            
                plt.tight_layout()

                if show:
                    plt.show()
                
                nombre= 'sc1_xy'


                if save:
                    print(i,'Guardando en ./img/%s/especies/TRACK/final/%s_%s_xy_%s'%(carpetaName,nombre,tiempos[i],carpetaName))
                    plt.savefig('./img/%s/especies/TRACK/final/%s_%s_%s.png'%(carpetaName,carpetaName,nombre,tiempos[i]))
                    plt.close()
            

            if meso: 

                fig = plt.figure(figsize=(9,9))

                #plt.suptitle('%s %s - Lat: %s - Lon: %s'%(carpeta,tiempos[i][4:-2], str(np.round(lat[mLat,mLon],2)), str(np.round(lon[mLat,mLon],2))),fontdict={'fontsize':20  })

                window = 15

                gspec = fig.add_gridspec(1,1)
                gs = gspec[0,0]


                collages.mesociclonInd(data,fig,i,x,y,mLat,mLon,window,press,topo,wMaxPrev,z,gs)
        
            
                plt.tight_layout()

                if show:
                    plt.show()

                nombre= 'sc1_meso'


                if save:
                    print(i,'Guardando en ./img/%s/especies/TRACK/final/%s_%s_wmax_%s'%(carpetaName,nombre,tiempos[i],carpetaName))
                    plt.savefig('./img/%s/especies/TRACK/final/%s_%s_%s.png'%(carpetaName,carpetaName,nombre,tiempos[i]))
                    plt.close()

            if conv: 

                fig = plt.figure(figsize=(12,9))

                #plt.suptitle('%s %s - Lat: %s - Lon: %s'%(carpeta,tiempos[i][4:-2], str(np.round(lat[mLat,mLon],2)), str(np.round(lon[mLat,mLon],2))),fontdict={'fontsize':20  })

                window = 20

                gspec = fig.add_gridspec(2,1)

                gs = (gspec[0,0],gspec[1,0])

                collages.convectivoLatLon(data,fig,i,press,lat,lon,mLat,mLon,window,topo,pmin,pmax,gs)
        
            
                plt.tight_layout()

                if show:
                    plt.show()

                nombre= 'sc1_convectivo'


                if save:
                    print(i,'Guardando en ./img/%s/especies/TRACK/final/%s_%s_wmax_%s'%(carpetaName,nombre,tiempos[i],carpetaName))
                    plt.savefig('./img/%s/especies/TRACK/final/%s_%s_%s.png'%(carpetaName,carpetaName,nombre,tiempos[i]))
                    plt.close()



            if hodor: 

                fig = plt.figure(figsize=(9,9))

                gspec = fig.add_gridspec(1,1)     

                gs = gspec[0,0]           

                collages.hodografaInd(data,fig,i,mLat,mLon,press,gspec=gs)
        
            
                plt.tight_layout()

                if show:
                    plt.show()

                nombre= 'sc1_hodografa'


                if save:
                    print(i,'Guardando en ./img/%s/especies/TRACK/final/%s_%s_%s'%(carpetaName,nombre,tiempos[i],carpetaName))
                    plt.savefig('./img/%s/especies/TRACK/final/%s_%s_%s.png'%(carpetaName,carpetaName,nombre,tiempos[i]))
                    plt.close()


            wMaxPrev = (i,mLat,mLon)

                
            
        break

        """ 
        lluviaPrev = False

        maxCoords = funciones.buscarMax(w[i,18,...],mLat,mLon,window,stride)

        mLat = maxCoords[1][0][0]
        mLon = maxCoords[1][1][0]

        window = 15

        lluvia         = variables.lluvia(data,i,lluvia_prev=lluviaPrev)

        press = variables.press(data,i)
        press = press/100
        pmin = 1000
        pmax = 100 


        w,u,v  = cargar.varIndNC(data,'W',tiempo=i,lev=[1,-1]),cargar.varIndNC(data,'U',tiempo=i),cargar.varIndNC(data,'V',tiempo=i)

        tPot = cargar.varIndNC(data,'T',tiempo=i)
        tPot = variables.setTpotToT(tPot,press)

        titaE = cargar.varIndNC(data,'THM',tiempo=i)
        
        titaE  = np.ma.masked_array(titaE,press<750)

        pmin = 1000
        pmax = 100

        fig = plt.figure(figsize=(10,10))

        gspec = fig.add_gridspec(2,2)

        def convectivo(ax,corte='xz'):
            if corte == 'yz': l = lat
            else: l = lon

            x2D,y2D   = variables.getXY(l,press,mLat = mLat,mLon = mLon,corte = corte,window=window)

            wi    = variables.setData3D(w, mLat, mLon, corte, window)
            ui    = variables.setData3D(u, mLat, mLon, corte, window)
            tPoti  = variables.setData3D(tPot,mLat, mLon, corte, window)

            titaEi = variables.setData3D(titaE,mLat,mLon,corte,window)
            
            topoi = variables.setData2D(topo,mLat = mLat,mLon = mLon,corte= corte,window=window)

            ejes.verticalConvectivo(ax,x2D,y2D,wi,ui,corte,ymin=pmin,ymax=pmax,topo =topoi,viento=True, tpot = tPoti,titaE=titaEi)
        

        ax = fig.add_subplot(gspec[0,0:2])
        corte = 'xz'

        convectivo(ax,corte)


        ax = fig.add_subplot(gspec[1,0:2])
        corte = 'yz'

        convectivo(ax,corte)


        try:
            plt.title('w max %s %s %s'%(str(np.max(w)),carpeta,tiempos[i]))
        except:
            print(len(tiempos))

        plt.tight_layout()

        if show:
            plt.show()

        if save:
            print(i,'Guardando en ./img/%s/especies/TRACK/final/%s_%s_wmax_%s'%(carpeta,nombre,tiempos[i],carpeta))
            plt.savefig('./img/%s/especies/TRACK/final/%s_%s_%s.png'%(carpeta,carpeta,nombre,tiempos[i]))
            plt.close()


    
        """