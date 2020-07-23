from WRF import *

class especies_wb():
    def __init__(self,archivo,wb_r,wb_i,wb_c,wb_v,wb_s):
        super().__init__()

        self.data = Dataset(archivo,'r')
        self.data_wb_v = Dataset(wb_v,'r')
        self.data_wb_c = Dataset(wb_c,'r')
        self.data_wb_r = Dataset(wb_r,'r')
        self.data_wb_i = Dataset(wb_i,'r')
        self.data_wb_s = Dataset(wb_s,'r')

    def cargar(self,sistema):
        
        self.sistema=sistema

        self.lat = np.array(self.data.variables['lat'])
        self.lon = np.array(self.data.variables['lon'])

        presb = np.array(self.data['p_basebox'])
        presp = np.array(self.data['p_perbox'])
        
        self.press = presb + presp

        self.pmin = 100
        self.pmax = 100000

        self.u = np.array(self.data.variables['uabox'])
        self.v = np.array(self.data.variables['vabox'])
        self.w = np.array(Dataset('/media/agustin/1D36257C5C6244CB/WRF_sim/3D-an_wrfout_d02.nc','r').variables['W'])
        
        folder = '/media/agustin/1D36257C5C6244CB/WRF_sim/'

        rain = np.array(Dataset(folder+'SFC-an_wrfout_d02.nc','r').variables['RAINNC'])

        rain[1:,...] = rain[1:,...]-rain[:-1,...]
        
        self.rain = np.ma.masked_array(rain, rain <2)

        self.medio = 30

        scaleStandar = 0.0000003
        self.data_q = {
            'r': { 
                'data':self.data_wb_r,
                'qhAdvMask':0.0000001,
                'qzAdvMask':0.0000001,
                'cmap': 'Greens',
                'scale': 0.000007,
                'title': 'QRAIN',
                'color':'g'
                },
            'i':{ 
                'data':self.data_wb_i,
                'qhAdvMask':0.00000001,
                'qzAdvMask':0.00000001,
                'cmap':'Greys',
                'scale': scaleStandar,
                'title': 'QICE',
                'color':'grey'
                },

            's': { 
                'data':self.data_wb_s,
                'qhAdvMask':0.00000001,
                'qzAdvMask':0.00000001,
                'cmap':'Reds',
                'scale': 0.000007,
                'title': 'QSNOW',
                'color':'r'
                },
            'c': { 
                'data':self.data_wb_c,
                'qhAdvMask':0.00000001,
                'qzAdvMask':0.00000001,
                'cmap':'Blues',
                'scale': scaleStandar,
                'title': 'QCLOUD',
                'color':'b'
                },
            'v':{ 
                'data':self.data_wb_v,
                'qhAdvMask':0.00000001,
                'qzAdvMask':0.00000001,
                'cmap':'Purples',
                'scale': 0.000003,
                'title': 'QVAPOR',
                'color': 'indigo'
                },
        }



    def mapa_especies (self,i, carpeta='./', save = False, show = True):


        self.i = i
        self.corteFunc = {
            'lat': self.xp2D(),
            'lon': self.yp2D()
        }

        try:
            self.getTopo()
        except:
            print('problemas con la topo')
        gs = [5,4]

        fig = plt.figure(figsize=(15,15))


        self.gs = fig.add_gridspec(gs[0],gs[1])

        qx = ['r','i','s','v','c']        


        for especie in range(4):

     
            ax_lat  = fig.add_subplot(self.gs[especie,0:2])

            self.eje_especies_ind_lat(ax_lat,qx[especie],corte='lat')

            ax_lon  = fig.add_subplot(self.gs[especie,2:4])

            self.eje_especies_ind_lat(ax_lon,qx[especie], corte= 'lon')



        plt.tight_layout()


        if show==True:
            plt.show()
        
        if save == True:
            plt.savefig(carpeta+'/WB_Sys%s_Especies%s.png'%(str(self.sistema),str(self.i)))
            plt.close(fig)


    def mapa_especies_unificado(self,i, carpeta='./', save = False, show = True):


        self.i = i
        self.corteFunc = {
            'lat': self.xp2D(),
            'lon': self.yp2D()
        }

        try:
            self.getTopo()
        except:
            print('problemas con la topo')
        gs = [2,1]

        fig = plt.figure(figsize=(15,15))


        self.gs = fig.add_gridspec(gs[0],gs[1])   

    
        ax_lat  = fig.add_subplot(self.gs[0,0])

        self.eje_especies_unif(ax_lat,corte='lat')

        ax_lon  = fig.add_subplot(self.gs[1,0])

        self.eje_especies_unif(ax_lon, corte= 'lon')



        plt.tight_layout()


        if show==True:
            plt.show()
        
        if save == True:
            plt.savefig(carpeta+'/WB_Sys%s_Especies%s.png'%(str(self.sistema),str(self.i)))
            plt.close(fig)

    def mapa_collage(self,i, carpeta='./', save = False, show = True):


        self.i = i
        self.corteFunc = {
            'lat': self.xp2D(),
            'lon': self.yp2D()
        }

        try:
            self.getTopo()
        except:
            print('problemas con la topo')

        gs = [3,4]

        fig = plt.figure(figsize=(15,15))


        self.gs = fig.add_gridspec(gs[0],gs[1])   

    
        ax_lat  = fig.add_subplot(self.gs[0,0:2])

        self.eje_especies_unif(ax_lat,corte='lat',nombre='Especies total Lat')

        ax_lon  = fig.add_subplot(self.gs[0,2:4])

        self.eje_especies_unif(ax_lon, corte= 'lon',nombre='Especies total Lon')

        
        ax_xy  = fig.add_subplot(self.gs[1:3,:2])

        self.eje_xy(ax_xy)



        plt.tight_layout()


        if show==True:
            plt.show()
        
        if save == True:
            plt.savefig(carpeta+'/WB_Sys%s_Especies%s.png'%(str(self.sistema),str(self.i)))
            plt.close(fig)

    def eje_xy(self,ax):


        
        x = self.lon[0,:,:]
        y = self.lat[0,:,:]
        
        qv = self.data.variables['qvbox'][self.i,0,:,:]

        ifile = self.data.variables['ifile'][self.i]
        jfile = self.data.variables['jfile'][self.i]

        rain = self.rain [self.i,jfile-self.medio-1:jfile+self.medio,ifile-self.medio-1:ifile+self.medio]

        print(self.w.shape,self.v.shape,self.u.shape)
        w = self.w[self.i,jfile-self.medio-1:jfile+self.medio,ifile-self.medio-1:ifile+self.medio]
        u = self.u[self.i,17,:,:]
        v = self.v[self.i,17,:,:]
        

        ax.pcolormesh(x ,y ,qv,cmap=plt.get_cmap('BrBG')) 
          
          
        
        self.ps = ax.pcolormesh(x,y,rain, cmap=plt.get_cmap('jet'), vmin =  0, vmax=50, zorder = 4)  

        try:
            ax.contour(x,y,self.topo,colors = 'black', linewidth = 0.5, zorder = 1,alpha = 0.5)
        except:
            print('problemas con la topo')
        ax.scatter(x[self.medio,self.medio],y[self.medio,self.medio],c= 'orange',zorder = 3)
        
        ax.contour(x,y,w, cmap = plt.get_cmap('seismic') )



        speed = np.sqrt((u*2) ** 2 + (v*2) ** 2)
        lw = 2*speed / speed.max()

        plt.streamplot(x[self.medio,:],y[:,self.medio],u,v,
                                 arrowstyle ='->',linewidth=lw, color ='k', density=0.8)

        

        ax.hlines(y[self.medio-1,self.medio-1],x[0,0],x[-1,-1], linestyles = 'dashed'
                   ,colors='r',label='Corte vertical')

        #plt.scatter(x[self.medio-1,self.medio-1],y[self.medio-1,self.medio-1])
        
        
        plt.xticks([])
        plt.yticks([])
        

        plt.title('Mapa centrado en lat: %s lon: %s \n Qv(somb) Precip (somb2) W(cont)'%(jfile,ifile))
        
        return

    def eje_especies_ind_lat(self,ax,q_name ,corte, nombre =''):


        qvar,qtend,qhadv,qzadv = self.get_q(q_name,corte)

        

        x2D, press2D = self.corteFunc[corte]

        try:     
            if corte =='lat': 
                qzadv = self.mask_topo_lat(qzadv,press2D)
                qhadv = self.mask_topo_lat(qhadv,press2D)
            elif corte == 'lon':
                qzadv = self.mask_topo_lon(qzadv,press2D)
                qhadv = self.mask_topo_lon(qhadv,press2D)
        except:
            print('problemas con  la topo')

        
        self.plot_contour_especies(ax,x2D,press2D,qtend,color='k')

        ax.pcolormesh(x2D,press2D, qvar, cmap=plt.get_cmap(self.data_q[q_name]['cmap']), alpha = 0.5)

        

        gap = 3
        ax.quiver(x2D[gap:,:],press2D[gap:,:],qhadv[gap:,:],qzadv[gap:,:],scale = self.data_q[q_name]['scale'] ,color ='k',units='dots')

        try:
            if corte == 'lat': 
                ax.plot(x2D[0,:],self.topo[self.medio,:],linewidth = 5, color = 'k', zorder= 3 )
            elif corte =='lon':
                ax.plot(x2D[0,:],self.topo[:,self.medio],linewidth = 5, color = 'k', zorder= 3 )
        except: print('sin topo')

        viento =  self.get_viento(corte)

        ax.contour(x2D,press2D,viento, colors = 'grey', zorder = 2,alpha=1)


        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        ax.set_facecolor('white')

        plt.title('%s - %s'%(self.data_q[q_name]['title'],corte))



    def eje_especies_unif(self,ax,corte,  nombre =''):
        
        x2D, press2D = self.corteFunc[corte]

        qx = ['r','i','s','v','c']        

        alpha = 0.5
        
        legend = []
        legendName=[]

        for q_name in qx:
            
            qvar,qtend,qhadv,qzadv = self.get_q(q_name,corte) 

            self.q_plot_somb = ax.pcolormesh(x2D,press2D, qvar, cmap=plt.get_cmap(self.data_q[q_name]['cmap']), alpha = alpha,zorder=2)

            color = self.data_q[q_name]['color']

            l = self.plot_contour_especies(ax=ax, x2D=x2D, press2D=press2D,var=qtend,color = color,zorder = 3)

            legend.append(l[-1])
            legendName.append(self.data_q[q_name]['title'])
            gap = 3
            ax.quiver(x2D[gap:,:],press2D[gap:,:],qhadv[gap:,:],qzadv[gap:,:],scale = self.data_q[q_name]['scale'] ,color =self.data_q[q_name]['color'],units='dots',zorder =3,alpha=0.5)


        ax.legend(legend,legendName)

        try:
            self.plot_topo(ax,x2D,corte)
        except:
            print('problemas con la topo')

        viento =  self.get_viento(corte)
        

        ax.contour(x2D,press2D,viento, colors = 'darkgrey', zorder = 2,alpha=1)

        plt.ylim (self.pmax,self.pmin)
        plt.ylabel('Presión (hPa)')
        ax.set_facecolor('white')

        plt.title(nombre)

    def plot_topo(self,ax,x2D,corte):

        if corte == 'lat': 
            ax.plot(x2D[0,:],self.topo[30,:],linewidth = 5, color = 'k', zorder= 3 )
        elif corte =='lon':
            ax.plot(x2D[0,:],self.topo[:,30],linewidth = 5, color = 'k', zorder= 3 )
        else:
            print('el corte esta mal')



    def get_viento(self,corte):
    
        if corte =='lat':
            viento = self.u[self.i,:,:,self.medio]

        elif corte =='lon':
            viento = self.v[self.i,:,self.medio,:]
        
        return viento

    def get_q(self,q_name, corte):
        qtend_str = q_name + 'ttendschbox'

        qvar_str  = q_name + 'box'

        if q_name =='v': 
            qvar_str  =  'qvbox'
            qtend_str = 'husttendbox'
        

        qh_str = 'hadv_q' + q_name + 'box'

        qz_str = 'zadv_q' + q_name + 'box'

        if corte== 'lat':

            lon = 30
            
            qvar = np.array(self.data.variables[qvar_str][self.i,:,:,lon])
            qtend = np.array(self.data_q[q_name]['data'].variables[qtend_str][self.i,:,:,lon])
            qhadv = np.array(self.data_q[q_name]['data'].variables[qh_str][self.i,:,:,lon])
            qzadv = np.array(self.data_q[q_name]['data'].variables[qz_str][self.i,:,:,lon])
            


        elif corte== 'lon':

            lat = 30

            qvar = np.array(self.data.variables[qvar_str][self.i,:,lat,:])
            qtend = np.array(self.data_q[q_name]['data'].variables[qtend_str][self.i,:,lat,:])
            qhadv = np.array(self.data_q[q_name]['data'].variables[qh_str][self.i,:,lat,:])
            qzadv = np.array(self.data_q[q_name]['data'].variables[qz_str][self.i,:,lat,:])
            
        
        

        qhAdvMask = self.data_q[q_name]['qhAdvMask']
        qzAdvMask = self.data_q[q_name]['qzAdvMask']
        qhadv = np.ma.masked_array(qhadv,abs(qhadv)<qhAdvMask)
        qhadv = np.ma.masked_array(qhadv,abs(qzadv)<qzAdvMask)
        qzadv = np.ma.masked_array(qzadv,abs(qhadv)<qhAdvMask)
        qzadv = np.ma.masked_array(qzadv,abs(qzadv)<qzAdvMask)

        return qvar,qtend,qhadv,qzadv

    def plot_contour_especies(self,ax,x2D,press2D, var, color,zorder = 3):

        def mult(x):
            return x*100000
        var = mult(var)

        
        m = np.max([abs(np.min(var)),abs(np.max(var))])

        levels =  np.linspace(-m,m,10)

        cs = ax.contour(x2D,press2D,var,levels = levels,colors = color, zorder=zorder) 
        cs.collections[5].remove()
        cs.collections[4].remove()

        #plt.setp( zc, linewidth=5)
        #plt.clabel(cs,inline=1, fmt='%1.2f')

        l,_ = cs.legend_elements()
        
        return l

    def xp2D(self):
        
        x=self.lon[0,self.medio,:]
        
        press = self.press[self.i,:,self.medio,self.medio]
        
        x2D, press2D = np.meshgrid(x, press)       
        
        return x2D, press2D

    def yp2D(self):

        y=self.lat[0,:, self.medio]    
        
        press = self.press[self.i,:,self.medio,self.medio]
        
        y2D, press2D = np.meshgrid(y, press)   
        
        return y2D, press2D

    def getTopo(self):

        arch = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-11_00:00:00'        
        i = self.data.variables['ifile'][self.i]
        j = self.data.variables['jfile'][self.i]

        dataTopo = Dataset(arch,'r')

        topo = np.array(getvar(dataTopo,'HGT'))
        self.topo = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,topo)))[j-self.medio-1:j+self.medio,i-self.medio-1:i+self.medio]

    def mask_topo_lat(self,var, press2D):
        
        return  np.ma.masked_array(var, press2D>self.topo[self.medio,:])

    def mask_topo_lon(self,var,press2D):

        return np.ma.masked_array(var, press2D>self.topo[:,0])