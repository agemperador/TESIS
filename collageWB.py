from WRF import *

from WRF_NewERA import mapa


wb = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfwbqv_d02_2018-11-11_00:00:00'
normal = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-11_00:00:00'
normal_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfout_d02_2018-11-10_18:00:00'
wb_prev = '/media/agustin/1D36257C5C6244CB/WRF_sim/wrfbqv_d02_2018-11-10_18:00:00'

data = Dataset(normal,'r')
data_wb = Dataset(wb,'r')


mapa = mapa(normal,normal_prev,wb,wb_prev)
#self.cargar_var()

print(' \033[1;32;40m Variables cargadas \033[1;0,0m')


i = 2
max_lat = 200
max_lon = 200


self.buscar_max(i,max_lat=max_lat,max_lon=max_lon)
print(self.max_lat, self.max_lon)
max_lat = self.max_lat
max_lon = self.max_lon
self.crear_mapa_wb(i, save = True,carpeta='/home/agustin/TESIS/Scripts/img/FF/WB/')
