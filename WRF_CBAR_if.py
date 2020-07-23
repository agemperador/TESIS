from librerias import *
import ejes
import mapas


cdict = {
'blue': [(0,1,0.5),(0.4,0.5,0.4) , (1, 1, 1)],
'green': [(0,0,0.2),(0.4,0.2,0), (1, 0, 0)],
'red': [(0,0,0) ,(0.4,0,0.4), (1, 1, 1)]}
cmapViento = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,100)


fig = plt.figure(figsize=(20,5))


ejes.colorbar(0,30,cmapViento,label= 'Velocidad del viento m/s')

ejes.colorbar(0.008,0.022,'BrBG',label = 'QVAPOR gr/kg')

ejes.colorbar(0,50,'gist_ncar',label ='lluvia mm')

plt.savefig('./img/colorbars/colorbarQvapor.png')
plt.show()