import astropy.units as u
import plotly.graph_objects as go
from astropy.coordinates import EarthLocation, CartesianRepresentation


def cartesian_to_geographic(x, y, z):
    # Converter de coordenadas cartesianas para locação na Terra
    cart_rep = CartesianRepresentation(x * u.km, y * u.km, z * u.km)
    earth_loc = EarthLocation.from_geocentric(cart_rep.x, cart_rep.y, cart_rep.z)

    return earth_loc.lat.deg, earth_loc.lon.deg, earth_loc.height.to(u.km)


def ground_track(t, trajetorias):
    # Constante de rotação da Terra (graus por segundo)
    taxa_rotacao_terra = 360 / 86400  # 360 graus em 86400 segundos (24 horas)

    # Listas para armazenar latitudes, longitudes e altitudes
    lats = []
    lons = []
    alts = []

    # Iterar sobre cada conjunto de posições na matriz de trajetórias
    for i in range(trajetorias.shape[0]):
        x = trajetorias[i, 0]
        y = trajetorias[i, 1]
        z = trajetorias[i, 2]
        lat, lon, alt = cartesian_to_geographic(x, y, z)

        # Ajustar a longitude de acordo com a rotação da Terra
        tempo_passado = t[i]  # Tempo em segundos desde o início da propagação
        ajuste_longitude = taxa_rotacao_terra * tempo_passado
        lon_ajustada = (lon + ajuste_longitude) % 360  # Garante que a longitude fique entre 0 e 360

        lats.append(lat)
        lons.append(lon_ajustada)
        alts.append(alt)

    # Plotar o ground track usando Plotly
    fig = go.Figure(go.Scattergeo(
        lon=lons,
        lat=lats,
        mode='lines+markers',  # Linhas e marcadores para visualizar o caminho
        line=dict(width=2, color='blue'),
        marker=dict(size=2, color='red'),
    ))

    fig.update_layout(
        title='Ground Track da Trajetória Orbital',
        geo=dict(
            projection_type='orthographic',
            showland=True,
            showcountries=True,
            showsubunits = True,
            landcolor="rgb(212, 212, 212)",
            countrycolor="rgb(204, 204, 204)",
        ),
    )

    fig.show()