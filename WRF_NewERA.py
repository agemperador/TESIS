

#wb = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqv_d02_2018-11-11_00:00:00'
#normal = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-11_00:00:00'
#normal_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-10_18:00:00'
#wb_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfbqv_d02_2018-11-10_18:00:00'

#data = Dataset(normal,'r')
#data_wb = Dataset(wb,'r')
from WRF import *


class mapa ():
    def __init__(_, normal,normal_prev,wb,wb_prev):
        super().__init__()

        _.data = Dataset(normal,'r')
        _.wb = wb
        _.archivo_prev = normal_prev
        _.archivo_prev2 = wb_prev

    #def cargar_var(_):
        lat,lon = getvar(_.data,'XLAT'), getvar(_.data,'XLONG')
        _.bm = get_basemap(lat)
        _.lat,_.lon = to_np(lat),to_np(lon)
        topo = np.array(getvar(_.data,'HGT'))
        _.topo = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,topo)))


        _.x, _.y  = _.bm(_.lon,_.lat)

        rain =  np.array(_.data.variables['RAINNC'])

        data_prev = Dataset(_.archivo_prev,'r')
        rain[1:] = rain[1:,...]-rain[:-1,...]
        rain[0] = rain[0,...]-np.array(data_prev.variables['RAINNC'][-1,...])

        data_prev.close()

        _.qv =  np.array(_.data.variables['QVAPOR'])
        _.t = np.array(_.data.variables['T'])

        _.rain = np.ma.masked_array(rain, rain <2)

        _.u,_.v,_.w = np.array(_.data.variables['U']),np.array(_.data.variables['V']), np.array(_.data.variables['W'])

        pres = np.array(_.data.variables['P'])
        presb = np.array(_.data.variables['PB'])
        _.press = pres + presb

        _.window = 20
        _.stride = 11
        _.stride_lon = 16


        _.wcmap = plt.get_cmap('seismic')

        _.data.close()

        data_wb = Dataset(_.wb,'r')
        _.qvtend  = to_np(data_wb.variables['QVTTEND'])
        _.qvh_adv = to_np(data_wb.variables['QVH_ADV'])

        data_wb.close()

        return

    def crear_mapa_wb(_,i,gs = [2,2], save = False, carpeta = './'):

        #_.max_lat,_.max_lon = 120,120
        #_.max_lat,_.max_lon = 120,120



        fig = plt.figure(figsize=(10,12))
        _.gs = fig.add_gridspec(gs[0],gs[1])


        ax1 = fig.add_subplot(_.gs[0,0:1])

        mapa.eje_wb_1(_,ax1,i)

        ax2 = fig.add_subplot(_.gs[0,1:2])

        mapa.eje_wb_2(_,ax2,i)

        ax3 = fig.add_subplot(_.gs[1,:])

        mapa.eje_wb_3(_,ax3,i)


        cax_qten = fig.add_axes([0.93, 0.28, 0.02, 0.15])
        cbar_qten = plt.colorbar(_.qten_plot, ax= ax3,cax = cax_qten, orientation = 'vertical')


        cax_qvh = fig.add_axes([0.93, 0.12, 0.02, 0.15])
        cbar_qvh = plt.colorbar(_.qvh_plot, ax = ax3, cax = cax_qvh, orientation = 'vertical')

        plt.tight_layout()

        plt.savefig(carpeta+'collage_wb_%i.png'%i)
        plt.show()

        return

    def eje_wb_1(_,ax,i):

        pad = 1


        print('\033[1;50;40m %s \033[1;0,0m'%str(i))


        #conv_wv = conv_wv[i,pad:-pad,pad:-pad]

        x = _.x[pad:-pad,pad:-pad]
        y = _.y[pad:-pad,pad:-pad]

        rain = _.rain[i,pad:-pad,pad:-pad]

        ax.pcolormesh(x,y, rain, cmap= plt.get_cmap('jet'), vmin =0, vmax= 50)


        u = _.u[i,0,pad:-pad,1+pad:-pad]
        v = _.v[i,0,1+pad:-pad,pad:-pad]

        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(_.x[0,pad:-pad],_.y[pad:-pad,0],u,v,
                       arrowstyle ='->',linewidth=lw, color ='k', density=0.8)


        _.bm.drawparallels(np.linspace(-35.,-30,6),labels=[1,0,0,0],color='gray')
        _.bm.drawmeridians(np.linspace(-66,-61,6),labels=[0,0,0,1], color = 'gray')


        return

    def eje_wb_2(_,ax,i):


        qvtend = _.qvtend[i,0,:,:]


        vmax = np.max(qvtend[10:-10,10:-10])
        ax.pcolormesh(_.x,_.y,qvtend,cmap= plt.get_cmap('seismic'), vmax = vmax, vmin = -vmax)

        u = _.u[i,0,:,1:]
        v = _.v[i,0,1:,:]

        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(_.x[0,:],_.y[:,0],u,v,
                       arrowstyle ='->',linewidth=lw, color ='k', density=0.8)


        col = 'orange'

        xs,ys = mapa.get_square()

        _.bm.plot(xs, ys, latlon = True)

        plt.scatter(_.max_lat,_.max_lon, marker='o',zorder= 4, c=col)



        _.bm.drawparallels(np.linspace(-35.,-30,6),labels=[1,0,0,0],color='gray')
        _.bm.drawmeridians(np.linspace(-66,-61,6),labels=[0,0,0,1], color = 'gray')

        return

    def eje_wb_3(_,ax,i):

        x=_.lon[0,:]

        end = 59

        press = _.press[i,:,_.max_lat,_.max_lon]


        x2D, press2D = np.meshgrid(x, press)
        x2D = x2D[:end,:]
        press2D = press2D[:end,:]

        qvtend = _.qvtend[i,:,_.max_lat,:]
        qvtend = np.ma.masked_array(qvtend, mask = abs(qvtend)<0.000005)
        qvtend = np.ma.masked_array(qvtend, mask = press2D>70000)


        qvh_adv = _.qvh_adv[i,:,_.max_lat,:]
        qvh_adv =  np.ma.masked_array(qvh_adv, mask = abs(qvh_adv)<0.5)
        qvh_adv =np.ma.masked_array(qvh_adv,mask = press2D<70000)


        _.qten_plot = ax.pcolormesh (x2D,press2D,qvtend, cmap = plt.get_cmap('coolwarm'))

        _.qvh_plot = ax.pcolormesh (x2D, press2D, qvh_adv, cmap = plt.get_cmap('Spectral'), vmin=-2, vmax = 2)



        #_.ascensos = ax.contour(x2D,press2D, _.w[i,:-1,_.max_lat,:])



        """
        clev, clev_bool = mapa.get_clevs(_.w[i,:-1,_.max_lat,:],ncont=9)

        if clev_bool:
            pc = ax.contour(x2D,press2D,_.w[i,:-1,_.max_lat,:],clev, colors = 'k', zorder = 3)
            clab = ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')

        """

        plt.ylim(100000,10000)


        return







    def crear_mapa_cv (_,i, gs = [6,2], save = False, carpeta = './'):


        fig = plt.figure(figsize=(10,12))
        _.gs = fig.add_gridspec(gs[0],gs[1])


        ax1 = fig.add_subplot(_.gs[0:2,1])

        #ax_bar_pp = fig.add_axes([0.83, 0.66, 0.02, 0.23])

        mapa.eje_mapa_xy(_,ax1,i)

        ax0 = fig.add_subplot(_.gs[0:2,0])

        mapa.eje_centrado_xy(_,ax0,i,0)

        ax2 = fig.add_subplot(_.gs[2:4,:])

        mapa.eje_xz(_,ax2,i)


        ax3 = fig.add_axes([0.135,0.46,0.155,0.12])

        mapa.hodografa(_,ax3,i)

        ax4 = fig.add_subplot(_.gs[4:6,:])

        mapa.eje_yz(_,ax4,i)


        caxr = fig.add_axes([1., 0.66, 0.02, 0.28])

        cbar_pp = plt.colorbar(_.p_rain, ax = [ax0,ax1],cax = caxr, orientation = 'vertical', pad = 0.03)
        cbar_pp.ax.set_ylabel('Precipitación hr')

        caxq = fig.add_axes([1.0, 0.11, 0.02, 0.45])

        cbar_qv = plt.colorbar(_.q_plot, ax = [ax2,ax4], cax = caxq, orientation = 'vertical' )
        cbar_qv.ax.set_ylabel('qvapor kg/kg')

        caxw = fig.add_axes([0.1,0.06,0.65,0.02])

        cbar_w = plt.colorbar(_.ascensos, ax = ax4, cax = caxw, orientation='horizontal')
        cbar_w.ax.set_xlabel('Ascensos y descensos (m/s)')



        plt.tight_layout()

        if save == True: plt.savefig(carpeta+'collage_track_%s.png'%str(i))

        plt.show()


    def eje_mapa_xy(_,ax,i):



        plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left=False,
            labelleft=False,
        )

        _.p_rain = ax.pcolormesh(_.x,_.y,_.rain[i,...],cmap=plt.get_cmap('jet'), zorder=3)

        ax.contour(_.x,_.y,_.topo, colors = 'black', linewidth = 0.5, zorder = -2)

        col = 'orange'

        xs,ys = mapa.get_square(_)

        _.bm.plot(xs, ys, latlon = True)

        plt.scatter(_.max_lat,_.max_lon, marker='o',zorder= 4, c=col)




        u = _.u[i,0,:,1:]
        v = _.v[i,0,1:,:]

        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(_.x[_.max_lat,:],_.y[:,_.max_lon],u,v,
                       arrowstyle ='->',linewidth=lw, color ='k', density=0.8)


        _.bm.drawparallels(np.linspace(-35.,-30,6),labels=[1,0,0,0],color='gray')
        _.bm.drawmeridians(np.linspace(-66,-61,6),labels=[0,0,0,1], color = 'gray')

        plt.title('Trackeo de la celda')

        return


    def eje_centrado_xy(_,ax,i,lev_qv = 0,w_level=17):

        x = _.x [_.max_lat-_.window:_.max_lat+_.window,_.max_lon-_.window:_.max_lon+_.window]
        y = _.y [_.max_lat-_.window:_.max_lat+_.window,_.max_lon-_.window:_.max_lon+_.window]
        qv = _.qv [i,lev_qv,_.max_lat-_.window:_.max_lat+_.window,_.max_lon-_.window:_.max_lon+_.window]
        rain = _.rain [i,_.max_lat-_.window:_.max_lat+_.window,_.max_lon-_.window:_.max_lon+_.window]
        topo = _.topo [_.max_lat-_.window:_.max_lat+_.window,_.max_lon-_.window:_.max_lon+_.window]
        w = _.w[i,w_level,_.max_lat-_.window:_.max_lat+_.window,_.max_lon-_.window:_.max_lon+_.window]
        u = _.u[i,17,2+_.max_lat-_.window:_.max_lat+_.window-2,2+_.max_lon-_.window:_.max_lon+_.window-2]
        v = _.v[i,17,2+_.max_lat-_.window:_.max_lat+_.window-2,2+_.max_lon-_.window:_.max_lon+_.window-2]


        pqv = ax.pcolormesh(x ,y ,qv,cmap=plt.get_cmap('BrBG'))

        ps = ax.pcolormesh(x,y,rain, cmap=plt.get_cmap('jet'), vmin =  0, vmax=50, zorder = 2)

        ax.contour(x,y,topo,colors = 'black', linewidth = 0.5, zorder = 1,alpha = 0.5)

        ax.scatter(x[_.window,_.window],y[_.window,_.window],c= 'orange',zorder = 3)


        clev,clev_bool =  mapa.get_clevs(_,w)

        if clev_bool:
            pc = ax.contour(x,y,w,clev, cmap = _.wcmap )
            clab = ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')


        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        v_plot0 = plt.streamplot(x[_.window,2:-2],y[2:-2,_.window],u,v,
                                 arrowstyle ='->',linewidth=lw, color ='k', density=0.8)



        ax.hlines(y[_.window-1,_.window-1],x[0,0],x[-1,-1], linestyles = 'dashed'
                   ,colors='r',label='Corte vertical')

        plt.scatter(x[_.window-1,_.window-1],y[_.window-1,_.window-1])


        plt.xticks([])
        plt.yticks([])



        return

    def eje_yz(_,ax,i):

        ### ESTO VA A CAMBIAR CON MAX_LON

        y=_.lat[:,_.max_lon]

        end = 59


        press = _.press[i,:,_.max_lat,_.max_lon]

        y2D, press2D = np.meshgrid(y, press)
        y2D = y2D[:end,_.max_lat-_.window:_.max_lat+_.window]
        press2D = press2D[:end,_.max_lat-_.window:_.max_lat+_.window]


        u = _.u[i,:end,_.max_lat-_.window:_.max_lat+_.window,_.max_lon]
        u =  mapa.mask_topo_simple(_,u, press2D)

        v = _.v[i,:end,_.max_lat-_.window:_.max_lat+_.window,_.max_lon]
        v =  mapa.mask_topo(_,v,press2D,_.max_lat)

        w = _.w[i,:end,_.max_lat-_.window:_.max_lat+_.window,_.max_lon]
        w =  mapa.mask_topo_simple(_,w, press2D)
        w_somb = np.ma.masked_array(w, abs(w)<1)

        titae =_.t[i,:end,_.max_lat-_.window:_.max_lat+_.window,_.max_lon]

        titae = mapa.mask_topo(_,titae,press2D,_.max_lat)
        #titae = np.ma.array(titae , mask = np.abs(titae) < 0)
        titae = np.ma.masked_array(titae, press2D<= 10000) + 290

        red = (70/255, 244/255)
        green = (130/255,164/255)
        blue =(180/255,96/255)

        cmap_per = mapa.get_cmap(_,red,green,blue)


        q_plot = ax.pcolormesh (y2D,press2D,titae, cmap = cmap_per,vmin = 298, vmax = 312)

        ascensos = ax.pcolormesh(y2D,press2D,w_somb,cmap=_.wcmap,vmin = -20,vmax=20)


        clev, clev_bool = mapa.get_clevs(_,v,ncont=9)

        if clev_bool:
            pc = ax.contour(y2D,press2D,u,clev, colors = 'k', zorder = 3)
            clab = ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')


        a=3 #Separacion de flechas

        #speed = np.sqrt((v[:end,::a]*2) ** 2 + (w[:end,::a]*2) ** 2)
        #lw = speed / speed.max()

        end_v = 30

        M = np.hypot(v[:end_v,::a], w[:end_v,::a])

        plot_u = ax.barbs(y2D[:end_v,3:-3:a],press2D[:end_v,3:-3:a], v[:end_v,3:-3:a],w[:end_v,3:-3:a],
                           color ='#393939')#,  pivot = 'middle',length = 7, linewidth =1 , zorder = 4 )




        ax.plot(y2D[0,:],_.topo[_.max_lat-_.window:_.max_lat+_.window,_.max_lon],linewidth = 5, color = 'k', zorder= 3 )

        plt.ylim (96600.,10000)
        plt.ylabel('Presión (hPa)')


        return

    def eje_xz(_,ax,i):

        ### ESTO VA A CAMBIAR CON MAX_LAT
        y_corte = 120

        x = _.lon[y_corte,:]

        end = 59

        w = _.w[i,:end,_.max_lat,_.max_lon-_.window:_.max_lon+_.window]
        u = _.u[i,:end,_.max_lat,_.max_lon-_.window:_.max_lon+_.window]
        v = _.v[i,:end,y_corte,_.max_lon-_.window:_.max_lon+_.window]
        thm =_.t[i,:end,y_corte,_.max_lon-_.window:_.max_lon+_.window]


        press = _.press[i,:,y_corte,y_corte]

        x2D, press2D = np.meshgrid(x, press)
        x2D = x2D[:end,_.max_lon-_.window:_.max_lon+_.window]
        press2D = press2D[:end,_.max_lon-_.window:_.max_lon+_.window]

        titae = mapa.mask_topo(_,thm,press2D,y_corte)
        #titae = np.ma.array(titae , mask = np.abs(titae) < 0)
        titae = np.ma.masked_array(titae, press2D<= 10000) + 290

        v =  mapa.mask_topo(_,v,press2D,y_corte)
        u =  mapa.mask_topo_simple(_,u, press2D)
        w =  mapa.mask_topo_simple(_,w, press2D)
        w_somb = np.ma.masked_array(w, abs(w)<1)


        red = (70/255, 244/255)
        green = (130/255,164/255)
        blue =(180/255,96/255)

        cmap_per = mapa.get_cmap(_,red,green,blue)




        _.q_plot = ax.pcolormesh (x2D,press2D,titae, cmap = cmap_per,vmin = 298, vmax = 312)

        _.ascensos = ax.pcolormesh(x2D,press2D,w_somb,cmap=_.wcmap,vmin = -20,vmax=20)



        clev, clev_bool = mapa.get_clevs(_,v,ncont=9)

        if clev_bool:
            pc = ax.contour(x2D,press2D,v,clev, colors = 'k', zorder = 3)
            clab = ax.clabel(pc,clev,fontsize=12,fmt='%.0f',colors='k')


        a=4 #Separacion de flechas

        print(u.shape,v.shape,w.shape, end,a)

        #speed = np.sqrt((u[:end,::a]*2) ** 2 + (w[:end,::a]*2) ** 2)
        #lw = speed / speed.max()

        end_v = 30

        M = np.hypot(u[:end_v,::a], w[:end_v,::a])

        plot_v = ax.barbs(x2D[:end_v,3:-3:a],press2D[:end_v,3:-3:a], u[:end_v,3:-3:a],w[:end_v,3:-3:a],
                           color ='#393939',  pivot = 'middle',length = 7, linewidth =1 , zorder = 4 )


        ax.plot(x2D[0,:],_.topo[_.max_lat,_.max_lon-_.window:_.max_lon+_.window],linewidth = 5, color = 'k', zorder= 3 )

        plt.ylim (96600.,10000)
        plt.ylabel('Presión (hPa)')

        plt.title('Corte vertical: W, Qv (somb) - Viento uw y V cont - Latitud: %s - Longitud: %s'%(str(_.max_lat),str(_.max_lon)))



        return


    def hodografa(_,ax,i):

        ini = 8

        end = 28

        v_c = _.v[i,ini:end,_.max_lat,_.max_lon]
        u_c = _.u[i,ini:end,_.max_lat,_.max_lon]

        cm = plt.cm.get_cmap('BuPu')
        colors=[cm(1.*k/(end-ini)) for k in range(end-ini)]

        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('center')
        # Eliminate upper and right axes
        ax.spines['right'].set_color('k')
        ax.spines['top'].set_color('k')
        ax.spines['left'].set_color('k')
        ax.spines['bottom'].set_color('k')
        # Show ticks in the left and lower axes only
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        vmax = np.max(abs(v_c))+2
        umax = np.max(abs(u_c))+2

        plt.xlim(-umax, umax)
        plt.ylim(-vmax,vmax)

        ax.plot(u_c,v_c,color = 'k', linewidth = 2)

        xy = range(ini,end)

        sc = ax.scatter(u_c,v_c, c =xy, cmap = cm, zorder = 5 )

        ax.set_facecolor('gray')
        ax.grid()

        clevs_hod = np.round(_.press[i,ini:end, _.max_lat,_.max_lon])


        return

    def buscar_max(_,i,max_lat=-1,max_lon=-1,w_level= 18):


        rsample = _.rain[i,...]
        wsample = _.w[i,...]

        if max_lat == -1 or max_lon == -1:

            maximo = np.argwhere((wsample[w_level,:,:] == wsample[w_level,:,:].max()))[0]

            max_lat_pr = 0
            max_lon_pr = 0
        else:

            wminisample = wsample[:,max_lat-_.stride:max_lat+_.stride,max_lon-_.stride_lon:max_lon+_.stride_lon]
            maximo = np.argwhere(wminisample[w_level,:,:] == wminisample[w_level,:,:].max())[0]

            max_lat_pr = max_lat - _.stride
            max_lon_pr = max_lon - _.stride_lon




        max_lat =  maximo[0] + max_lat_pr
        max_lon =  maximo[1] + max_lon_pr

        if max_lon + _.window >= _.x.shape[0]:
            max_lat = 0
            max_lon = 0

        _.max_lat = max_lat
        _.max_lon = max_lon



        return




    def mask_topo(_,var,press2D,y_corte):


        for j in range(var.shape[1]):

            top = _.topo[y_corte,j]

            var[:,j] = np.ma.array(var[:,j], mask = press2D[:,j] > top)

        return var

    def mask_topo_simple(_,var, press2D):

        var= np.ma.masked_array(var, press2D>_.topo[_.max_lat,_.max_lon-_.window:_.max_lon+_.window])

        return var

    def get_cmap(_,red,green,blue):


        cdict = {
        'blue': [(0,0,1),(0.4,blue[0],blue[0]), (1, 0.95, 1)],
        'green': [(0,0,1),(0.4,green[0],green[0]), (1, 0.99, 0)],
        'red': [(0,0,1),(0.4,red[0],red[0]) , (1,0.99, 1)]}
        cmap_per = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,25)

        return cmap_per

    def get_clevs(_,variable,ncont = 5):

        c_min = np.round(np.min(variable))
        c_max = np.round(np.max(variable))
        c_m = np.max((abs(c_min),abs(c_max)))

        clev_bool = c_min != c_max

        if clev_bool:
            clev=np.linspace(-c_m,c_m,ncont)
            clev= np.delete(clev,np.argwhere(clev == 0))

            return clev,clev_bool
        else: return 0,0



    def get_square(_):

        lower_left = (_.lon[_.max_lat-_.window,_.max_lon-_.window],_.lat[_.max_lat-_.window,_.max_lon-_.window])
        lower_right =(_.lon[_.max_lat-_.window,_.max_lon+_.window],_.lat[_.max_lat-_.window,_.max_lon+_.window])
        upper_left = (_.lon[_.max_lat+_.window,_.max_lon-_.window],_.lat[_.max_lat+_.window,_.max_lon-_.window])
        upper_right =(_.lon[_.max_lat+_.window,_.max_lon+_.window],_.lat[_.max_lat+_.window,_.max_lon+_.window])


        xs = [lower_left[0], upper_left[0],
              upper_left[0], upper_right[0],
              upper_right[0], lower_right[0],
              lower_right[0], lower_left[0]]
        ys = [lower_left[1], upper_left[1],
              upper_left[1], upper_right[1],
              upper_right[1], lower_right[1],
              lower_right[1], lower_left[1]]

        return xs,ys
