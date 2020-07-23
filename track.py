from WRF import *
import pandas as pd


data = Dataset('/media/agustin/1D36257C5C6244CB/WRF_sim/3D-an_wrfout_d02.nc')

w = np.array( data.variables['W'])

window = 20
stride = 15
stride_lon = 16

def buscar_max(w,max_lat=-1,max_lon=-1,time=0,window=20,stride=11,stride_lon=16):

    wsample = w[time,...]
    
    if max_lat == -1 or max_lon == -1: 
        
        maximo = np.argwhere((wsample[:,:] == wsample[:,:].max()))[0]

        max_lat_pr = 0
        max_lon_pr = 0
    else:
        if (max_lat+window>=304):
            stride = 0
        if (max_lon+window>=304):
            stride_lon=0
        wminisample = wsample[max_lat-stride:max_lat+stride,max_lon-stride_lon:max_lon+stride_lon]
        maximo = np.argwhere(wminisample[:,:] == wminisample[:,:].max())[0]
        
        max_lat_pr = max_lat - stride
        max_lon_pr = max_lon - stride_lon 

        

    max_lat =  maximo[0] + max_lat_pr
    max_lon =  maximo[1] + max_lon_pr

    if max_lon + window >= w.shape[1]: 
        max_lat = 0
        max_lon = 0
        
    max_lat = max_lat
    max_lon = max_lon
        
        
    
    return    max_lat,max_lon

coords = []
max_lat,max_lon=230,210
for t in range(95,110,1):
    max_lat,max_lon = buscar_max(w, max_lat=max_lat,max_lon=max_lon,time = t)
    coords.append((t,max_lat,max_lon))
    w_s = np.ma.masked_array(w[t,...],abs(w[t,...])<1)
    plt.figure()
    plt.imshow(w_s)
    plt.scatter(max_lon,max_lat,marker='o',c='orange')
    plt.title('Tiempo: %i Lat: %i  Lon: %i'%(t,max_lat,max_lon))
    plt.show()
    print(t, max_lat,max_lon)

cord = pd.DataFrame(coords,columns=['Time','Lat','Lon'])

cord.to_csv('Coordenadas_5.csv',index=False)
