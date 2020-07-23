from librerias import *

import carga_vars as variables

import plots


def verticalConvectivo(ax,x,y,w,u, corte = 'xz', window =False, ymin=0,ymax = 100, topo = False, viento = False,tpot=False,
                       titaE=False):

    plt.title('Corte vertical %s'%corte)

    wSomb = np.ma.masked_array(w,abs(w)<0.5)
    
    pm = plots.pmesh(ax,x,y,wSomb,'seismic',vmin = -35,vmax = 35)
    
    if type(tpot)!=bool:

        print('tpot')

        tpot = np.where(tpot < 273.3 ,0,1)

        cs =   ax.contour(x[:,:],y[:,:],tpot[:,:], levels = 0)
        plt.clabel(cs,inline=1, fmt='%1.0f')
        #plots.pcont(ax,x[:-25,:],y[:-25,:],tpot[:-25,:],color='k', especial='t')
        #plots.pcont(ax,x,y,cont,color = 'k', completo = True)

    if type(titaE) != bool:

        print('titae')

        titaE = np.ma.masked_array(titaE, y > topo )

        ax.contour(x[:-40,:],y[:-40,:],titaE[:-40,:],colors= 'indigo', levels= 5)

    if type(topo) is not bool :
        
        try:
            plots.pline(ax,x[0,:],topo)
            
        except:
            print('No se puede graficar la topografìa')

    if viento:

        plots.pbarbs(ax,x,y,[u,w], cmap = 'k')

    plt.ylim((ymin,ymax))

    plt.colorbar(pm, use_gridspec = True)

    return 

def verticalEspecies(ax,x,y,q,w,cmap, corte = 'xz',ymin=0,ymax = 100, topo = False, viento = False,tpot = False):

    #for i in range(len(q)):
    #    plots.pcont(ax,x,y,q[i],color=cmaps[i])

    pm = plots.pmesh(ax,x,y,q,cmap = cmap)

    plt.colorbar(pm,use_gridspec =True)

    plots.pcont(ax,x,y,w,'k')


    if type(tpot)!=bool:

        tpot = np.where(tpot < 273.3 ,0,1)

        cs =   ax.contour(x[:,:],y[:,:],tpot[:,:], levels = 0)
        plt.clabel(cs,inline=1, fmt='%1.0f')

    plt.yticks([])
    plt.xticks()

    if type(topo) is not bool:
        
        try:
            plots.pline(ax,x[0,:],topo)

        except:
            print('No se puede graficar la topografìa')



    plt.ylim((ymin,ymax))



def multiEspecies(ax,x,y,w,u,qv,qr,qc,qi,qg,qs, corte = 'xz', window =False, ymin=0,ymax = 100, topo = False, 
                viento = False,tpot=False, div=False):

    plt.ylim((ymin,ymax))

    plt.title('Corte vertical %s'%corte)

    v = np.ma.masked_array(qv,y<400)
    r = np.ma.masked_array(qr,y<600)
    c = np.ma.masked_array(qc,y>700) 
    c = np.ma.masked_array(c, y<200)
    i = np.ma.masked_array(qi,y>400)
    s = np.ma.masked_array(qs,y>600)
    g = np.ma.masked_array(qg,y>600)

    plots.pmesh(ax,x,y,v,'Purples',alpha=0.5)
    plots.pmesh(ax,x,y,c,'Blues',alpha=0.5)
    plots.pmesh(ax,x,y,i,'Greys',alpha=0.5)

    plots.pcont(ax,x,y,qr,'g')
    plots.pcont(ax,x,y,qg,'indigo')
    plots.pcont(ax,x,y,qs,'r')

    if type(div)!= bool:
        
        plots.pcont(ax,x,y,div,'k')

    if type(topo) is not bool:
        
        try:
            plots.pline(ax,x[0,:],topo)
            
        except:
            print('No se puede graficar la topografìa')






def xy(ax,x,y,somb = False,r = False,bm = False,cmapSomb = 'BrBG',cont=False,u= False,v= False,topo = False, viento = True, 
llueve = True,dot = False,linewidth=2,removeContour=[4,5],clevs=False,vmin = False,vmax = False,cmapViento=False,
vminViento=False,vmaxViento = False):

    if type (somb)!=bool:
        pm = plots.pmesh(ax,x,y,somb,cmap=cmapSomb,vmin = vmin,vmax = vmax)

        #plt.colorbar(pm,use_gridspec = True, shrink=0.5)


    if viento:
        plots.streamplot(ax,x,y,u,v,linewidth=linewidth, cmapViento=cmapViento,vminViento=vminViento,vmaxViento=vmaxViento)

    if llueve:
        print('rain',r.shape)
        pm = plots.pmesh(ax,x,y,r,cmap = 'gist_ncar',zorder = 12 ,vmin=0,vmax=50)
        #if r.shape[0]<300:
        #plt.colorbar(pm,use_gridspec =True,shrink=0.5, pad=0.04)

    if type(cont)!=bool:
        
        plots.pcont(ax,x,y,cont,color='brown',remove=removeContour,clevs=clevs,linewidth=0.7)

    if type( topo)!=bool:
        plots.pcont(ax,x,y,topo,color = 'grey',completo=True,zorder = 0)
        
    if type(dot)!=bool:
        plt.scatter(x[dot[1],dot[0]],y[dot[1],dot[0]],c='orange',zorder = 6)
       
   

    plt.yticks([])
    plt.xticks([])



def hodografa(ax, u,v,press):

    ini = 8
    end = 28

    cm = plt.cm.get_cmap('viridis')
    colors=[cm(1.*k/(end-ini)) for k in range(end-ini)]

    #ax.spines['left'].set_position('center')
    #ax.spines['bottom'].set_position('center')
    
    # Eliminate upper and right axes
    ax.spines['right'].set_color('k')
    ax.spines['top'].set_color('k')
    ax.spines['left'].set_color('k')
    ax.spines['bottom'].set_color('k')
    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')



    vmax = np.max(abs(v))+2
    umax = np.max(abs(u))+2

    plt.xlim(-umax, umax)
    plt.ylim(-vmax,vmax)

    ax.plot(u,v,color = 'k', linewidth = 2)

    xy = range(ini,end)

    hod = ax.scatter(u,v, c =xy, cmap = cm, zorder = 6 )

    plt.colorbar(hod,use_gridspec=True)

    ax.set_facecolor('white')
    ax.grid()

    clevs_hod = np.round(press)

    plt.title('Hodografa')

    return


def mesociclon (ax,x,y,u,v,div,topo=False):

    if type( topo)!=bool:
        print('topo')
        plots.pcont(ax,x,y,topo,color = 'grey',completo=True)
        
    plots.pmesh(ax,x,y,div, cmap='seismic',alpha =0.6)

    #ax.contour(x,y,press,colors='r')
    #plots.pcont(ax,x,y,press,color='r', especial='t')

    plots.streamplot(ax,x,y,u,v)

    
def colorbar(vmin,vmax,cmap, orientation='horizontal', label=False):

    a = np.array([[vmin,vmax]])
    img = plt.imshow(a, cmap=cmap)
    plt.gca().set_visible(False)
    cb = plt.colorbar(orientation=orientation, pad=0.09)
    cb.ax.set_xlabel(label)



        