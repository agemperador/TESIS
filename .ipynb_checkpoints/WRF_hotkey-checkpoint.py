from WRF import *

def mapa_sup(a, rango = [0,48,4] ,todoT = False, topografia = False, lluvia= True, viento =True
            , guardar = False, show = True):

    lluvia_previa = 0
    carpeta_add = ''
    if todoT == True:
        rango = [0,48,1]
        carpeta_add = 'TodoT/'
        
        
    
    
    for arch in a:
        for t in range(rango[0],rango[1],rango[2]):
            if arch[32:-30] == '/test01/': carpeta = '/control/'
            else: carpeta = arch[32:-30]
            if carpeta== 'WRF_sim/': carpeta = '/FF/'
            archivo =  str(arch[32:])
            if t == 23:
                lluvia_previa = to_np(Dataset(arch).variables['RAINNC'])[t,:-1,:-1]
            elif t ==0:
                lluvia_inicial = to_np(Dataset(arch).variables['RAINNC'][t,:-1,:-1])-lluvia_previa

            elif t == 1:
                lluvia_inicial = -1

            mapa_somb_cont2(path= str(arch[:32]),archivo=str(arch[32:]), 
                            var_s='Q2',vmin = 0.008,vmax = 0.022,cmap_s = 'BrBG',
                            var_c1='TH2',cmap_c1='hot_r',
                            guardar=guardar,savepath = './img%sQVAPOR/%s'%(carpeta,carpeta_add), show = show,ini = 0, end = 304,rows = 1, cols = 1,
                            ts = [t],  modo = 'normal', #lev_rich=k,
                            titulo="Q 2m (s) - $\Theta$ 2m (c) -  V sup (slines) - Lluvia (mm/hr)",

                            lluvia= lluvia, cmap_lluvia = 'gist_ncar', vminr = 0.5, vmaxr=50,
                            topografia=topografia,
                            viento = viento, lluvia_inicial = lluvia_inicial
                            )

