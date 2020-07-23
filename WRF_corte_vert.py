from WRF import *

def set_corte_vert(data,y_corte, t, topo,var_s,var_c,press,mask_s=1,cmap = 'coolwarm'):    

    
    var_s = np.array(var_s)
    var_c = np.array(var_c)
    
    x = data.variables['XLONG'][t,y_corte,:]
    
    press = np.array(press[t,:,y_corte,y_corte]) 

    var_s = var_s[t,:,y_corte,:]
    var_c = var_c[t,:,y_corte,:]
    
    x2D, press2D = np.meshgrid(x, press)

    colmiss='#AAAAAA'
    wcmap = plt.get_cmap(cmap)

    for j in range(305):

        top = topo[t,y_corte,j]

        var_s[:,j] = np.ma.array(var_s[:,j], mask = press2D[:,j] > top)
        var_c[:,j] = np.ma.array(var_c[:,j], mask = press2D[:,j] > top)

    var_s = np.ma.array(var_s , mask = np.abs(var_s) < mask_s)
    
    return x2D, press2D,var_s,var_c

def get_var_corte_vert(archivos,index, variables):
    
    if type(archivos) == list:
        archivo = archivos[index]
    elif type(archivos) == str:
            archivo = archivos
    else: 
        print('El archivo debe ser una lista de strings o una sola string')
        return 0
    print(archivo)
    
    data = Dataset(archivo,'r')
    
    lat = getvar(data,'XLAT')
    lon = getvar(data,'XLONG')
    bm = get_basemap(lon)

    x,y = bm(to_np(lon),to_np(lat))
    
    tiempos = data.variables['Times']

    topo = data.variables['HGT']
    topo = np.asarray(list(map(lambda x: -1.225*9.8*x +101300 ,topo)))

    
    pres = np.array(data.variables['P'])
    presb = np.array(data.variables['PB'])

    press = pres + presb
    
    varlist = []
    
    for var in variables:
        
        
        if var == 'RAINNC':
            
                r = np.array(data.variables[var])
                
                r_horaria = r[1:,...]-r[:-1,...]
    
                r_0 = Dataset(archivos[index-1],'r').variables['RAINNC'][-1,...]

                r_00 = r[0,...]-r_0

                r_h = r_00

                r_horaria = np.concatenate((np.reshape(r_h,(1,305,305)),r_horaria), axis = 0 )
                
                varlist.append(r_horaria)
                
        else: varlist.append(np.array(data.variables[var]))

    return data,x,y,bm, topo, press,tiempos, varlist

    
    
    
    
    