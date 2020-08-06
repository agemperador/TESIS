from librerias import *

import carga_vars as variables

import carga_datos as carga

import ejes

# Añadir titae , mesociclon




def convectivo(fig,lat,lon,press,mLat,mLon,window,u,w,tPot,titaE,topo,pmin,pmax,corte='xz',gspec = None):
    ax = fig.add_subplot(gspec)

    if corte == 'yz': l = lat
    else: l = lon

    x2D,y2D   = variables.getXY(l,press,mLat = mLat,mLon = mLon,corte = corte,window=window)

    wi    = variables.setData3D(w, mLat, mLon, corte, window)
    ui    = variables.setData3D(u, mLat, mLon, corte, window)
    tPoti  = variables.setData3D(tPot,mLat, mLon, corte, window)

    titaEi = variables.setData3D(titaE,mLat,mLon,corte,window)
    
    topoi = variables.setData2D(topo,mLat = mLat,mLon = mLon,corte= corte,window=window)

    ejes.verticalConvectivo(ax,x2D,y2D,wi,ui,corte,ymin=pmin,ymax=pmax,topo =topoi,viento=True, tpot = tPoti,titaE=titaEi)

#caxW = plt.axes([2/5, 1/3+0.05, 0.005, 2/3-0.1])
#cbarW = fig.colorbar(wPlot, ax = ax, cax = caxW ,orientation = 'vertical')
#cbar.ax.set_yticklabels([vminr,10,20,30,40,vmaxr])
#Defino los limites del gráfico
#ra.set_clim(0, 40)
#cbar.draw_all()
#cbarW.ax.set_ylabel('W', labelpad = -70 , fontsize = 12)



def especies(fig,gspec,lon,lat,position,mLat,mLon,press,topo,w,qv,qc,qi,qr,qs,qg,tPot,corte,vpos,window,pmin,pmax):
    q = []
    print(corte)
    q.append(variables.setData3D(qv, mLat, mLon, corte, window))
    q.append(variables.setData3D(qc, mLat, mLon, corte, window))
    q.append(variables.setData3D(qi, mLat, mLon, corte, window))
    q.append(variables.setData3D(qr, mLat, mLon, corte, window))
    q.append(variables.setData3D(qs, mLat, mLon, corte, window))
    q.append(variables.setData3D(qg, mLat, mLon, corte, window))

    tPoti  = variables.setData3D(tPot,mLat, mLon, corte, window)


    cmaps = ['Purples','Blues','Greys','Greens','Reds','RdPu']
    q_names = ['QV','QC','QI','QR','QS','QG']

    if corte =='xz': l = lon
    elif corte =='yz': l = lat

    x2D,y2D   = variables.getXY(l,press,mLat = mLat,mLon = mLon,corte = corte,window=window)

    topoi = variables.setData2D(topo,mLat = mLat,mLon = mLon,corte= corte,window=window)



    wi = variables.setData3D(w, mLat, mLon, corte, window)


    axEsp1 = []

    for j in range(len(q)):

    
        axEsp1.append(fig.add_subplot(gspec[position[j][0]+vpos,position[j][1]:position[j][1]+2]))

        ejes.verticalEspecies(axEsp1[j],x2D,y2D,q[j],wi,cmaps[j],ymin=pmin,ymax=pmax,topo =topoi,viento=True,tpot=tPoti)

        plt.title(q_names[j])
    




def xyCentrado (fig,x,y,mLat,mLon,bm,window,topo,qv,lluvia,u,v,w,gspec,showLat=False):
    xi = variables.setXY(x,mLat,mLon,window)
    yi = variables.setXY(y,mLat,mLon,window)
    qvi = variables.setXY(qv[0,...] ,mLat,mLon,window)
    r  = variables.setXY(lluvia,mLat,mLon,window)



    ui  = variables.setXY(u[1,...],mLat,mLon,window)
    vi  = variables.setXY(v[1,...],mLat,mLon,window)
    wi  = variables.setXY(w[18,...],mLat,mLon,window)

    topoi = variables.setXY(topo,mLat,mLon,window)

    ax = fig.add_subplot(gspec)

    cdict = {
    'blue': [(0,1,0.5),(0.4,0.5,0.4) , (1, 1, 1)],
    'green': [(0,0,0.2),(0.4,0.2,0), (1, 0, 0)],
    'red': [(0,0,0) ,(0.4,0,0.4), (1, 1, 1)]}
    cmapViento = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,10)


    ejes.xy(ax,xi,yi,qvi,r,bm, u=ui,v=vi,topo = topoi,
    dot = False,w=wi,showLat=showLat,clevs=True,vmin = 0.008,vmax = 0.022,cmapViento=cmapViento,vminViento=0,vmaxViento=60)



def xyTotal(fig,x,y,bm,mLon,mLat,u,v,topo,lluvia,gspec):
    ax = fig.add_subplot(gspec)


    ui  = variables.setXY(u[1,...],mLat,mLon)[:303,:303]
    vi  = variables.setXY(v[1,...],mLat,mLon)[:303,:303]

    x = x[:303,:303]
    y = y[:303,:303]
    topo = topo[:303,:303]
    lluvia = lluvia[:303,:303]


    ejes.xy(ax,x,y,False,lluvia,bm, u=ui,v=vi,topo = topo,viento = True, dot = [mLon,mLat],)


def mapaHodografa(fig,u,v,mLat,mLon,press,gspec):
    ax = fig.add_subplot(gspec)

    ini = 8
    end = 28


    ui  = variables.setData3D(u[ini:end,...],mLat,mLon,corte='columna')
    vi  = variables.setData3D(v[ini:end,...],mLat,mLon, corte='columna')
    pressi = variables.setData3D(press[ini:end,...],mLat,mLon, corte='columna')


    ejes.hodografa(ax,ui,vi,pressi)







def mapaMultiespecies(fig,lon,lat,press,mLat,mLon,corte,window,w,u,v,qv,qr,qi,qc,qs,qg,topo,pmin,pmax,gspec):

    ax = fig.add_subplot(gspec)

    if corte == 'xz': l = lon
    elif corte == 'yz': l = lat


    x2D,y2D   = variables.getXY(l,press,mLat = mLat,mLon = mLon,corte = corte,window=window)

    wi    = variables.setData3D(w, mLat, mLon, corte, window)
    ui    = variables.setData3D(u, mLat, mLon, corte, window)
    vi    = variables.setData3D(v, mLat, mLon, corte, window)

    qvi   = variables.setData3D(qv,mLat,mLon,corte,window)
    qri   = variables.setData3D(qr,mLat,mLon,corte,window)
    qii   = variables.setData3D(qi,mLat,mLon,corte,window)
    qci   = variables.setData3D(qc,mLat,mLon,corte,window)
    qsi   = variables.setData3D(qs,mLat,mLon,corte,window)
    qgi   = variables.setData3D(qg,mLat,mLon,corte,window)



    div = variables.getVort(ui,vi)

    #tPoti  = variables.setData3D(tPot,mLat, mLon, corte, window)

    topoi = variables.setData2D(topo,mLat = mLat,mLon = mLon,corte= corte,window=window)


    ejes.multiEspecies(ax,x2D,y2D,wi,ui,qvi,qri,qci,qii,qgi,qsi,corte,ymin=pmin,ymax=pmax,topo =topoi,viento=True, tpot = False,div=div)


"""
div = variables.getVort(ui,vi)

#div = np.ma.masked_array(div,abs(div)<2)

div = np.ma.masked_array(div,abs(div)<4)

ejes.mesociclon(ax,xi,yi,ui,vi,div,topoi)"""

def mapaMesociclon(fig,x,y,mLat,mLon,window,press,u,v,topo,wMaxPrev=False,z=False,gspec=[1,1]):

    ax = fig.add_subplot(gspec)



    xi = variables.setXY(x,mLat,mLon,window)
    yi = variables.setXY(y,mLat,mLon,window)

    topoi = variables.setXY(topo,mLat,mLon,window)

    if not z : 
        level500 = np.where(abs(press[:,mLat,mLon]-500 )==np.min(abs(press[:,mLat,mLon]-500)))[0][0]

    else: level500=16
    ui  = variables.setXY(u[level500,...],mLat,mLon,window)
    vi  = variables.setXY(v[level500,...],mLat,mLon,window)


    if type(wMaxPrev)!=bool:

        print(wMaxPrev)
        # 0 es t 1 es lat 2 es lon
        uRel = (mLon - wMaxPrev[2]) / 0.25
        vRel = (mLat - wMaxPrev[1]) / 0.25

        print('WMAXPREV')

        ui = ui -  np.mean(ui) # - uRel
        vi = vi -  np.mean(vi) # - vRel

    div = variables.getVort(ui,vi)

    #div = np.ma.masked_array(div,abs(div)<2)

    div = np.ma.masked_array(div,abs(div)<4)

    ejes.mesociclon(ax,xi,yi,ui,vi,div,topoi)

    






def qvapor(data,i,x,y,bm,lat,lon,topo,tiempos,show=True,save= False, lluviaPrevia = False, 
                llueve=False, carpeta='', nombre= '', z =False,cmapViento=False):

    l = 303

    lluviaPrevia = lluviaPrevia[304,304]

    lluvia         = variables.lluvia(data,i)#,lluvia_prev=lluviaPrevia)

    if z:

        z = variables.z(data,tiempo=i)

        pmin = 0
        pmax = 14000

        press = z

        

    else:    
        press          = variables.press(data,i)
        press = press/100
        pmin = 1000
        pmax = 100 


    qv = carga.varIndNC(data,'QVAPOR',tiempo=i)
    titae = carga.varIndNC(data,'THM',tiempo=i)

    w,u,v  = carga.varIndNC(data,'W',tiempo=i,lev=[1,-1]),carga.varIndNC(data,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)


    fig = plt.figure(figsize=(10,10))
  
    gspec = fig.add_gridspec(1,1)

    ax = fig.add_subplot(gspec[0,0])   

    print('Tiempos: ',tiempos[i])



    rain = lluvia[:l,:l]
    
    qvi =  qv[0,:l,:l]
    xi = x[:l,:l]
    yi = y[:l,:l]
    ui = u[0,:l,:l]
    vi = v[0,:l,:l]
    titaei=titae[1,:l,:l]

    titaei = np.ma.masked_array(titaei,abs(titaei)<9)

    topoi=topo[:l,:l]

    ejes.xy(ax,xi,yi,qvi,r=rain,bm=bm,cmapSomb='BrBG',cont = titaei,u=ui,v=vi,topo=topoi,linewidth=3,
            removeContour=[4,5], clevs=True,vmin = 0.008,vmax = 0.022,cmapViento=cmapViento,vminViento=0,vmaxViento=60)


    bm.drawcoastlines(linewidth=0.25)
    bm.drawstates(linewidth=1)
    bm.drawcountries(linewidth=0.25)

    bm.drawparallels(np.linspace(np.min(lat),np.max(lat),6),labels=[1,0,0,0],color='gray')
    bm.drawmeridians(np.linspace(np.min(lon),np.max(lon),6),labels=[0,0,0,1], color = 'gray')


    plt.title(tiempos[i]+'0')



    if show:
        plt.show()

    if save:
        plt.savefig('./img/%s/QVAPOR/new/%s_%s_%s.png'%(carpeta,carpeta,nombre,tiempos[i]))   
        plt.close()




def test(data,i,x,y,bm,lat,lon,topo,tiempos,mLat,mLon,window = False,show=True,save= False, lluviaPrevia = False, llueve=False, carpeta='',wMaxPrev=False):


    lluvia         = variables.lluvia(data,i,lluvia_prev=lluviaPrevia)
    press          = variables.press(data,i)
    pmin = 1000
    pmax = 100 

    

    w,u,v  = carga.varIndNC(data,'W',tiempo=i,lev=[1,-1]),carga.varIndNC(data,'U',tiempo=i),carga.varIndNC(data,'V',tiempo=i)
    
    window = 30

    if type(wMaxPrev)!=bool and False:
        # 0 es t 1 es lat 2 es lon
        uRel = (mLon - wMaxPrev[1][0]) / 0.25
        vRel = (mLat - wMaxPrev[0][0]) / 0.25

        print(wMaxPrev)

        print(uRel,vRel)

        print(np.max(u),np.max(v))

        u = u - uRel 
        v = v - vRel
        
        


    qc,qv,qi= carga.varIndNC(data,'QCLOUD',tiempo=i),carga.varIndNC(data,'QVAPOR',tiempo=i),carga.varIndNC(data,'QICE',tiempo=i)
    qs =carga.varIndNC(data,'QSNOW',tiempo=i)
    qg =carga.varIndNC(data,'QGRAUP',tiempo=i)
    qr =carga.varIndNC(data,'QRAIN',tiempo=i)

    tPot = carga.varIndNC(data,'T',tiempo=1)
    tPot = variables.setTpotToT(tPot,press/100)

    presion = carga.varIndNC(data,'P',tiempo=i)/100

    presion = np.ma.masked_array(presion,abs(presion)<1)

    press = press / 9.8

    pmin = 0
    pmax = 14000


    fig = plt.figure(figsize=(20,10))

    gspec = fig.add_gridspec(6,10)


    xi = variables.setXY(x,mLat,mLon,window)
    yi = variables.setXY(y,mLat,mLon,window)
    topoi = variables.setXY(topo,mLat,mLon,window)
    


    k = 0
    for i,j in enumerate([11,13]):

        pressi = variables.setXY(presion[j,...],mLat,mLon,window)

        
        ui  = variables.setXY(u[j,...],mLat,mLon,window)
        vi  = variables.setXY(v[j,...],mLat,mLon,window)
        
        k= (i*2-i*2%6)//6

        ax = fig.add_subplot(gspec[k*2:k*2  +2,int(i*2)%6:int(i*2)%6+2])

        div = variables.getVort(ui,vi)

        div = np.ma.masked_array(div,abs(div)<4)

        ejes.mesociclon(ax,xi,yi,ui,vi,div,pressi,topoi)

    plt.show()
