from librerias import *

def pbarbs(ax,x,y,vecViento, cmap):

    #speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
    #lw = 2*speed / speed.max()
    a=4 #Separacion de flechas    
    b = 2
    end_v = 30
    pad = 3

    x = x[:end_v:b,pad:-pad:a]
    y = y[:end_v:b,pad:-pad:a]
    u = vecViento[0][:end_v:b,pad:-pad:a]
    v = vecViento[1][:end_v:b,pad:-pad:a]

    
    M = np.hypot(u, v)

    ax.barbs(x,y, u,v, color ='#444444',length=6,linewidth=1,alpha=0.7)#,  pivot = 'middle',length = 7, linewidth =1 , zorder = 4 )       

def streamplot(ax,x,y,u,v,linewidth=2,color='k',cmapViento=False,vminViento=False,vmaxViento=False):

    end_v = 30
    pad = 1

    x = x[0,pad:]
    y = y[pad:,0]
    u = u[pad:,pad:]
    v = v[pad:,pad:]



    speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
    lw = 3*linewidth*speed / speed.max()



    if type(cmapViento)== bool:
        s = ax.streamplot(x, y, u,v, color =color,linewidth=lw,zorder = 5, arrowsize = 1.2)
    else: 
        s = ax.streamplot(x, y, u,v, color=speed,linewidth=lw,zorder = 5,cmap=cmapViento, arrowsize=1.2)

    





def pcont(ax,x,y,var,color, completo = False,especial='',remove=[4,5],clevs=False,linewidth=1,zorder = 5, cmap=False):

    m = np.max([abs(np.min(var)),abs(np.max(var))])
    try:
        if especial=='' and type (cmap)==bool:

            if completo == False:

                levels =  np.linspace(-m,m,10)
                cs = ax.contour(x,y,var,levels = levels,colors=color,linewidth=linewidth,zorder=zorder)

                for i in remove:
                    cs.collections[i].remove()
                    

            else: 
                
                cs = ax.contour(x,y,var,colors=color,linewidth=linewidth,zorder=zorder)
        
        elif type (cmap)!=bool:
            print('cmap')
            if completo == False:

                levels =  np.linspace(-m,m,5)
                cs = ax.contour(x,y,var,levels = levels,cmap=cmap,linewidth=linewidth,zorder=zorder)

                for i in remove:
                    cs.collections[i].remove()

            else: 
                
                cs = ax.contour(x,y,var,cmap=cmap,linewidth=linewidth,zorder=zorder)


        elif especial=='t':

            levels =  np.linspace(-m,m,10)

            print(levels.shape)

            cs = ax.contour(x,y,var,levels = levels,colors=color,linewidth=1,zorder=zorder)
            """
            #cs.collections[5].remove()
            #cs.collections[4].remove()
            cs.collections[3].remove()
            cs.collections[6].remove()
            cs.collections[2].remove()
            cs.collections[7].remove()
            cs.collections[1].remove()
            cs.collections[8].remove()
            cs.collections[0].remove()
            cs.collections[9].remove()
            """
        

            #plt.setp( zc, linewidth=5)
        if clevs: plt.clabel(cs,inline=1, fmt='%1.0f')
    except: 
        print('no se pudo graficar el contorno')


def pmesh(ax,x,y,var,cmap='Reds',vmin=False,vmax=False,alpha =1,zorder= 43):
    # pcolormesh
    if vmin == False and vmax == False:
        try:
            vmin = np.min(var)
            vmax = np.max(var)
        except:
            vmin = -2   
            vmax = 2
    try:
        return ax.pcolormesh(x,y,var,cmap = cmap,vmin=vmin,vmax=vmax,alpha=alpha,zorder=0)
    except:
        print('Una variable con cmap %s no se pudo graficar'%cmap)


def pline(ax, x,var,linewidth=5,color='k',zorder=3):
    
    ax.plot(x,var, linewidth = linewidth, color=color, zorder =zorder)