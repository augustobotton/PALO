from mayavi import mlab
import numpy as np
from tvtk.api import tvtk

# Definir o tamanho da malha
phi, theta = np.mgrid[0:np.pi:101j, 0:2*np.pi:101j]

# Coordenadas esféricas
r = 1
x = r * np.sin(phi) * np.cos(theta)
y = r * np.sin(phi) * np.sin(theta)
z = r * np.cos(phi)

# Carregar a imagem
img = tvtk.JPEGReader(file_name=r"G:\.shortcut-targets-by-id\1-35xUduX2YD_n1NJU0ku_vtiEkY6R714\trab2\aula_29_certo\ultimofinalmesmo29\earth.jpg")
texture = tvtk.Texture(input_connection=img.output_port, interpolate=1)

# Criar a esfera
mesh = mlab.mesh(x, y, z, representation='surface', color=(1, 1, 1))

# Aplicar a textura
mesh.actor.enable_texture = True
mesh.actor.tcoord_generator_mode = 'plane'
mesh.actor.actor.texture = texture

mlab.show()
