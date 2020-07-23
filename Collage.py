from WRF import *

#Los tiempos
tiempos = ["11-10_06:00:00","11-10_12:00:00","11-10_18:00:00","11-11_00:00:00",
           "11-11_06:00:00","11-11_12:00:00"]

#coordenadas cercanas al máximo para este tiempo
max_lat = 125
max_lon = 97

ini = 0
end = 24

for j,t in enumerate(tiempos[1:]):
    
    dt = 0
    
    wb_v = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqv_d02_2018-%s'%t
    wb_c = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqc_d02_2018-%s'%t
    wb_r = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqr_d02_2018-%s'%t
    wb_i = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqi_d02_2018-%s'%t
    wb_s = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqs_d02_2018-%s'%t

    #Hago esto porque los datos de wb estan cada 6 hs y los otros cada 12
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

    #Cargo las variables iniciales
    mapas.cargar_var()
    
    print('Variables Cargadas - %s'%(time.time() - start_time))

    #Esto es para el trackeo
    mapas.window = 12
    mapas.stride_lon = 17
    mapas.stride = 18

    
    #Arranco en 14 porque en este tiempo hay máximos cerca de esas coordenadas
    for i in range(14,end,1):
        
        #Track
        mapas.buscar_max(i+dt,max_lat=max_lat,max_lon = max_lon)
      
        #Este mapa hacia el gráfico de las especies
        #mapas.crear_mapa_especies(i+dt,i,save = False,carpeta = './img/FF/WB')

        mapas.crear_mapa_WB(i+dt,i,save = False,carpeta = './img/FF/WB')
        #Este mapa hace el collage 
        #mapas.crear_mapa_SC_CI(i+dt,i,save = False,carpeta = './img/FF/WB')

        
        print("--- %s seconds ---" % (time.time() - start_time))

        break
        
    
        
    ini = 0
    max_lat = mapas.max_lat
    max_lon = mapas.max_lon

    break