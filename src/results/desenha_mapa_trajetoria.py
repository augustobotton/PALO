import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt



def desenha_mapa_trajetoria(ori, traj):
    # Funcao para desenhar um mapa em projecao, junto de uma trajetoria fornecida
    # Entradas
    # or = [delta0,lamb0,h0]: origem do mapa. delta0 (graus) latitude. lamb0 (graus) longitude. h0 (m) altitude.
    # traj=[lat long], matriz onde a primeira coluna sao as latitudes e a segunda as longitudes de uma trajetoria, em graus

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)

    m = Basemap(projection='cyl', resolution='c', ax=ax, llcrnrlon=-180, llcrnrlat=-90, urcrnrlon=180, urcrnrlat=90)

    m.drawcoastlines(color='green')
    m.drawcountries(color='green')
    m.drawrivers(color='blue')
    m.fillcontinents(color='lightgreen')

    # Plotting trajectory
    lats = traj[:, 0]
    longs = traj[:, 1]
    x, y = m(longs, lats)
    m.plot(x, y, marker='o', color='red')
    m.plot(x[0], y[0], marker='*', color='red')

    plt.show()