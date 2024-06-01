import numpy as np
from mayavi import mlab

def desenhar_orbita(longitude, latitude, altitude, texture_path):
    # Converter latitude, longitude e altitude para coordenadas cartesianas
    earth_radius = 6371  # Raio médio da Terra em quilômetros
    x = (earth_radius + altitude) * np.cos(latitude) * np.cos(longitude)
    y = (earth_radius + altitude) * np.cos(latitude) * np.sin(longitude)
    z = (earth_radius + altitude) * np.sin(latitude)

    # Criar uma figura
    fig = mlab.figure(size=(800, 600), bgcolor=(0.1, 0.1, 0.1))

    # Adicionar a textura da Terra ao globo
    surface = mlab.mesh(x, y, z, colormap='gist_earth', scalars=z, figure=fig)
    surface.actor.actor.texture = texture_path

    # Adicionar a trajetória do veículo espacial
    mlab.plot3d(x, y, z, color=(1, 0, 0), tube_radius=0.1, figure=fig)

    # Configurar a visualização do globo
    mlab.view(azimuth=0, elevation=90, distance='auto', focalpoint=(0, 0, 0))
    
    # Vetores de longitude, latitude e altitude em radianos
    longitude = long  # Vetor de longitudes em radianos
    latitude = delta  # Vetor de latitudes em radianos
    altitude = h  # Vetor de altitudes
    texture_path = "G:\.shortcut-targets-by-id\1-35xUduX2YD_n1NJU0ku_vtiEkY6R714\trab2\aula_29_certo\ultimofinalmesmo29\earth.jpg"  # Caminho para a imagem de textura da Terra

    desenhar_orbita(long, delta, h, texture_path)


    # Mostrar a figura
    mlab.show()

