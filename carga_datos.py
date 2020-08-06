from librerias import *



def dataNC(archivo, data= False):      

    """
    Carga los datos de netcdf, ya sea por parametro o por archivo
    Lo mete en self.data 

    """
    if not data:

        try:
            data = Dataset(archivo,'r')
            return data

        except:
            try: 
                print('No se encontro ningun archivo %s'%archivo)
            except:
                print('Falta ingresar el archivo por parámetro')
    else: 
        try:
            print(data.variabes.keys())
            return data
        except:
            print('Los datos no son de un archivo nc')


def varIndNC(data, var, tiempo = 0,lev = False,lat = False,lon=False):
    """
    Carga variables segun su nombre

    """
    

    try:
        if str(type(data))== "<class 'netCDF4._netCDF4.Dataset'>":

            if not lat: lat = [0,-1]            
            if not lon: lon = [0,-1]
            if not lev: lev = [0,-1]
            
            if var =='XLAT' or var=='XLONG':
            
                variable = getvar(data,var)
                return variable

            elif var == 'HGT':
                variable = getvar(data,var) 
                return variable

            if len(data.variables[var].shape)==4:     
                
                if tiempo == 'todo':
                    variable = data.variables[var][:,lev[0]:lev[1],lat[0]:lat[1],lon[0]:lon[1]]
                else: 
                    variable = data.variables[var][tiempo,lev[0]:lev[1],lat[0]:lat[1],lon[0]:lon[1]]

            elif len(data.variables[var].shape) == 3:
                if tiempo == 'todo':
                    variable = data.variables[var][:,lat[0]:lat[1],lon[0]:lon[1]]
                else: 
                    variable = data.variables[var][tiempo,lat[0]:lat[1],lon[0]:lon[1]]       

            elif len(data.variables[var].shape)==2:
                variable = data.variables[var][lat[0]:lat[1],lon[0]:lon[1]]

            elif len(data.variables[var].shape) == 1:
                if tiempo == 'todo':
                    variable = data.variables[var][:]
                else: 
                    variable = data.variables[var][tiempo]
                    
            else:
                print('No se puede cargar %s'%(var))
                return 
            
            return np.array(variable)
    except:
        assert('El archivo no es netcdf')

