from WRF import *

folder = '/media/agustin/1D36257C5C6244CB/WRF_sim/'

tiempos = ["11-10_06:00:00","11-10_12:00:00","11-10_18:00:00","11-11_00:00:00",
           "11-11_06:00:00","11-11_12:00:00"]

data = Dataset(folder+'SFC-an_wrfout_d02.nc','r')

lluvia_raw = np.array(data.variables['RAINNC'][:])
#lluvia = np.array(np.zeros_like(lluvia_raw))
lluvia_h = lluvia_raw[1:]-lluvia_raw[:-1]
lluvia = np.ma.masked_array(lluvia_h,lluvia_h<1)
#lluvia[0,...] = np.ma.masked_array(0,True)
#lluvia[1:,...] = lluvia_ma
ti = data.variables['intTime'][:]

ini = 0
end = 24



for j,t in enumerate(tiempos[1:2]):
    
    dt = 0
    
    wb_v = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqv_d02_2018-%s'%t
    wb_c = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqc_d02_2018-%s'%t
    wb_r = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqr_d02_2018-%s'%t
    wb_i = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqi_d02_2018-%s'%t
    wb_s = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqs_d02_2018-%s'%t
    if t[-6:-4] == '06': 
        t[-6:-4]='00'
        dt = 24
    elif t[-6:-4] == '18':
        t[-6:-4] ='12'
        st = 24
    normal = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-%s'%t
    normal_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-10_12:00:00'
    

    start_time = time.time()
    mapas = mapa(normal,normal_prev,wb_v,wb_c,wb_r,wb_i,wb_s)
    print('Cargando Variables - %s'%(time.time() - start_time))
    mapas.cargar_var()
    
    print('Variables Cargadas - %s'%(time.time() - start_time))

    
    
    for i in range(144):
        
        fig = plt.figure(figsize = (10,10))

        mapas.rain = lluvia
        ax = fig.add_subplot(212)
        mapas.eje_mapa_xy(ax,i,viento= False)

        plt.colorbar(mapas.p_rain)
        plt.title('Tiempo %s'%str(ti[i]))

        plt.savefig('./img/FF/WRF_FF_Corregido/Lluvia/lluvia_track%s'%str(ti[i]))
        


    
        