import numpy as np
import plotly.graph_objects as go

from src.domain.modelos.planeta.ConstrutorPlaneta import terra


def inercial_para_geografica(t, trajetorias, planeta=terra):
    """
    Converte posições inerciais (ECI) em latitude e longitude fixas no planeta.

    A longitude fixa no planeta é a longitude inercial menos o ângulo de rotação
    acumulado do meridiano de referência: ωE·(t − t_g), onde t_g é o instante em
    que a longitude celeste do meridiano de referência é nula.

    :param t: Vetor de tempos [s].
    :param trajetorias: Matriz de estados (N x >=3), posições em km nas 3 primeiras colunas.
    :param planeta: Planeta com velocidade_inercial_de_rotacao [rad/s] e tempo_longitude_celeste_nula [s].
    :return: Latitude geocêntrica [deg], longitude em [-180, 180) [deg], raio [km].
    """
    t = np.asarray(t, dtype=float)
    x, y, z = (np.asarray(trajetorias)[:, i] for i in range(3))

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    lat = np.degrees(np.arcsin(z / r))

    angulo_rotacao = planeta.velocidade_inercial_de_rotacao * (t - planeta.tempo_longitude_celeste_nula)
    lon = np.degrees(np.arctan2(y, x) - angulo_rotacao)
    lon = (lon + 180.0) % 360.0 - 180.0  # Garante que a longitude fique em [-180, 180)

    return lat, lon, r


def ground_track(t, trajetorias, planeta=terra):
    lats, lons, _ = inercial_para_geografica(t, trajetorias, planeta)

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
