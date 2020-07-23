from WRF import *

folder = '/media/agustin/1D36257C5C6244CB/WRF_sim/tracked_files/'

for sist in range(1,5):
    data_3Dpath = folder + 'System'+str(sist) + '_3D-wrfout_traj_B61-R30.nc'
    data_rpath = folder+'System'+str(sist) + '_qr_traj_B61-R30.nc'
    data_ipath = folder+'System'+str(sist) + '_qi_traj_B61-R30.nc'
    data_spath = folder+'System'+str(sist) + '_qs_traj_B61-R30.nc'
    data_vpath = folder+'System'+str(sist) + '_qv_traj_B61-R30.nc'
    data_cpath = folder+'System'+str(sist) + '_qc_traj_B61-R30.nc'

    datas = [data_rpath,data_ipath, data_spath,data_vpath, data_cpath]

    extension = 'ttendschbox'
    var_wb = 'r'


    from w_budget_functions import especies_wb as esp_wb

    mapa_esp = esp_wb(data_3Dpath, wb_s=data_spath,
                wb_c=data_cpath,wb_i=data_ipath,
                wb_r=data_rpath,wb_v=data_vpath)
    try:
        mapa_esp.cargar(sist)
    except:
        continue
    i = 5
    for i in range(29):
        mapa_esp.mapa_especies(i, save = True,show=False, carpeta='/home/agustin/TESIS/Scripts/img/FF/WRF_FF_Corregido/Especies')
        print('mapa %s sistema %s'%(str(i),str(sist)))





