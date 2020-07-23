from WRF import *

from WRF_NewERA import mapa


wb = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqv_d02_2018-11-11_00:00:00'
normal = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-11_00:00:00'
normal_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-10_18:00:00'
wb_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfbqv_d02_2018-11-10_18:00:00'
    
data = Dataset(normal,'r')
data_wb = Dataset(wb,'r')


mapa = mapa(normal,normal_prev,wb,wb_prev)
mapa.cargar_var()



i = 2
max_lat = 200
max_lon = 200


mapa.buscar_max(i,max_lat=max_lat,max_lon=max_lon)
print(mapa.max_lat, mapa.max_lon)
max_lat = mapa.max_lat
max_lon = mapa.max_lon
mapa.crear_mapa_wb(i, save = True,carpeta='/home/agustin/TESIS/Scripts/img/FF/WB/')