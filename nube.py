from WRF import corte_vert

path = '/media/agustin/Linux/salidas_wrf/compl-ndg/'
archivo ='wrfout_d02_2018-11-10_12:00:00'
#a =!ls {path}/compl-ndg/*wrfout_d02_2018*
    
archivo= path+archivo
#for archivo in a:
t = 16
corte_vert(archivo,tiempos= t, n_var_s ='QICE', n_var_c = 'V',minimo=0, maximo=0.0016,mask_s=0.0001, cmap = 'Blues', ini=30,end=60, guardar = True)