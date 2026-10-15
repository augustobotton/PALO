"""
Gera os dados exibidos na página inicial (docs/index.html).

Todos os números e trajetórias da página saem daqui, calculados pela própria biblioteca PALO.
Execute a partir da raiz do repositório:

    python docs/gerar_dados_pagina.py

Saída: docs/assets/dados.js (define window.PALO_DADOS). A simulação de lançamento leva alguns minutos.
"""
import json
import os
import sys
import tempfile
import time

import matplotlib

matplotlib.use("Agg")
import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)

from src.domain.modelos.Simulacao import Simulacao
from src.domain.modelos.foguete.aerodinamica.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.estrutura.ConstrutorModeloEstrutural import ConstrutorModeloEstrutural
from src.domain.modelos.foguete.propulsao.ConstrutorModeloPropulsivo import ConstrutorModeloPropulsivo
from src.domain.modelos.foguete.veiculo_lancador.ConstrutorFoguete import ConstrutorDeFoguete
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_cowell import propaga_cowell
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_encke import propagacao_encke
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.modelos.orbitas.utilidades.vetor_de_estados_para_coe import vetor_de_estados_para_coe
from src.domain.modelos.planeta.ConstrutorBaseDeLancamento import ConstrutorBaseDeLancamento
from src.domain.modelos.planeta.ConstrutorPlaneta import ConstrutorPlaneta, terra as terra_km
from src.domain.modelos.planeta.ModeloAtmosferico import ModeloAtmosferico
from src.domain.modelos.planeta.ModeloGravidadePlanetaAxisimetrico import calcular_gravidade_axisimetrico

MU_KM = 398600.4418
RE_KM = 6378.1370
OMEGA_TERRA = 7.2921150e-5
JSON_ATM = os.path.join(RAIZ, "src", "domain", "modelos", "planeta", "dados_JSON_planetas",
                        "dados_atmosfericos_terra.json")


def r(valores, casas=4):
    """Arredonda para reduzir o tamanho do arquivo de dados."""
    return [round(float(v), casas) for v in valores]


def terra_em_metros():
    return (ConstrutorPlaneta()
            .com_delta_temperatura_atm(10)
            .com_raio_equatorial(RE_KM * 1e3)
            .com_velocidade_inercial_de_rotacao(OMEGA_TERRA)
            .com_gravidade(9.80665)
            .com_mut(MU_KM * 1e9)
            .com_J2(0.00108263)
            .com_J3(-0.00000254)
            .com_J4(-0.00000161)
            .com_tempo_longitude_celeste_nula(0)
            .com_modelo_atmosferico(ModeloAtmosferico(JSON_ATM))
            .construir())


def elementos_orbitais():
    """Casos de referência para a calculadora de elementos orbitais da página."""
    casos = {
        "leo": ("Órbita baixa inclinada (exemplo de Encke)", [-2384.46, 5729.01, 3050.46], [-7.36138, -2.98997, 1.64354]),
        "decaimento": ("Satélite de 100 kg (exemplo de Cowell)", [5872.94259799039, -662.622151725639, 3007.48704280510],
                       [-2.89355195466701, 4.09603530867803, 6.14446573555145]),
        "equatorial": ("Órbita equatorial elíptica", [8200, 8200, 0], [-4.5, 4.5, 0]),
        "hiperbolica": ("Trajetória hiperbólica (Tewari, ex. 5.2)", [-5000, 0, 12500], [5, -8, 0]),
    }
    saida = {}
    for chave, (nome, rr, vv) in casos.items():
        o = Orbita.criar_pelo_vetor_de_estado(np.array(rr, float), np.array(vv, float), MU_KM)
        saida[chave] = {
            "nome": nome, "r": rr, "v": vv,
            "palo": {"a": o.semi_eixo_maior, "e": o.excentricidade, "i": o.inclinacao, "raan": o.raan,
                     "w": o.arg_periastro, "theta": o.anomalia_verdadeira, "h": o.quantidade_momento_angular},
        }
    return saida


def lancamento():
    """Exemplo adquire_gso_exemplo.py: veículo de 3 estágios saindo de Alcântara rumo à órbita geoestacionária."""
    terra = terra_em_metros()
    isp, mp, queima = [251, 271, 315], [77366, 11058, 225.33], [62, 64.62, 279]
    me = [7750, 1367, 64.7544]
    propulsivo = (ConstrutorModeloPropulsivo().com_impulso_especifico(isp).com_massa_propelente_estagios(mp)
                  .com_duracao_queima_estagios(queima).com_tempo_primeira_queima_terceiro_estagio(213)
                  .com_tempo_espera_separacao([2, 2, 2]).com_tempo_espera_ignicao([5, 580]).com_h0(0.0)
                  .com_planeta(terra).construir())
    estrutural = (ConstrutorModeloEstrutural().com_massa_estrutural_por_estagio(me).com_massa_de_carga_util(13)
                  .com_area_secao_transversal_1_estagio(4.6 * 5 / 3).com_area_secao_transversal_2_estagio(1.5)
                  .com_area_secao_transversal_3_estagio(1.5).com_area_secao_transversal_carga_util(1.5)
                  .com_comprimento_carga_util(1).com_comprimento_total_do_foguete(7.33 + 7.1 + 6.28)
                  .com_comprimento_sem_1_estagio(7.1 + 6.28).com_comprimento_sem_2_estagio(6.28).construir())
    foguete = (ConstrutorDeFoguete().com_modelo_propulsivo(propulsivo).com_modelo_estrutural(estrutural)
               .com_modelo_aerodinamico(ModeloAerodinamico()).construir())
    alcantara = (ConstrutorBaseDeLancamento().com_altitude_base(0).com_latitude_inicial().com_longitude_inicial()
                 .com_comprimento_trilho(estrutural.comprimento_total_do_foguete).construir())

    orbita_alvo = Orbita.circular(42.164140e6, np.deg2rad(8))
    simulacao = Simulacao(terra, alcantara, foguete, [120000, 1, np.deg2rad(80), orbita_alvo])

    # Simulacao.simular() grava um .pkl no diretório atual; roda em um diretório temporário para não sujar o repositório
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        try:
            t, y = simulacao.simular()
        finally:
            os.chdir(cwd)

    # Amostragem: fina na subida, esparsa na costa até o apogeu e na órbita final
    idx = np.unique(np.concatenate([
        np.searchsorted(t, np.arange(0, 1200, 2)),
        np.searchsorted(t, np.arange(1200, t[-1], 60)),
        [len(t) - 1],
    ]))
    idx = idx[idx < len(t)]
    tt, yy = t[idx], y[:, idx]

    massa = []
    for k in range(len(tt)):
        _, m, _, _ = propulsivo.propulsao_n_estagios(tt[k], yy[:, k], estrutural)
        massa.append(float(np.squeeze(m)))

    altitude_km = (yy[3] - terra.raio_equatorial) / 1e3
    lon_inercial = yy[5] + OMEGA_TERRA * tt  # longitude fixa à Terra + rotação acumulada

    eventos = {
        "ignicao": r(propulsivo.tempos_de_ignicao, 1),
        "fim_queima": r(propulsivo.tempos_de_fim_de_queima, 1),
        "separacao": r(propulsivo.tempos_de_separacao, 1),
    }
    return {
        "t": r(tt, 1),
        "altitude_km": r(altitude_km, 2),
        "velocidade_ms": r(yy[0], 1),
        "elevacao_graus": r(np.rad2deg(yy[2]), 2),
        "latitude_graus": r(np.rad2deg(yy[4]), 3),
        "longitude_graus": r(np.rad2deg(yy[5]), 3),
        "longitude_inercial_rad": r(lon_inercial, 5),
        "massa_kg": r(massa, 1),
        "eventos": eventos,
        "veiculo": {"isp_s": isp, "propelente_kg": mp, "estrutura_kg": me, "queima_s": queima, "carga_util_kg": 13,
                    "massa_inicial_kg": round(float(sum(mp) + sum(me) + 13), 1)},
        "alvo": {"raio_km": 42164.140, "inclinacao_graus": 8},
        "azimute_final_graus": round(float(np.rad2deg(simulacao.azimute_final)), 3),
        "base": {"nome": "Centro de Lançamento de Alcântara",
                 "latitude_graus": round(float(np.rad2deg(alcantara.latitude_base_de_lancamento)), 4),
                 "longitude_graus": round(float(np.rad2deg(alcantara.longitude_base_de_lancamento)), 4)},
    }


def perturbacao_j2():
    """Exemplo de propagação de Encke: 48 h com perturbação J2."""
    r0 = np.array([-2384.46, 5729.01, 3050.46])
    v0 = np.array([-7.36138, -2.98997, 1.64354])
    orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, MU_KM)
    t, y = propagacao_encke(0, 48 * 3600, orbita)
    t = np.array(t)

    raan, w = [], []
    for k in range(len(t)):
        coe = vetor_de_estados_para_coe(y[k, :3], y[k, 3:], MU_KM)
        raan.append(coe[3])
        w.append(coe[4])
    raan = np.unwrap(np.array(raan))
    w = np.unwrap(np.array(w))

    passo = max(1, len(t) // 1500)
    sel = slice(None, None, passo)

    # Traço no solo: longitude fixa à Terra = longitude inercial - rotação da Terra (GMST nulo em t = 0)
    raio = np.linalg.norm(y[:, :3], axis=1)
    lat = np.rad2deg(np.arcsin(y[:, 2] / raio))
    lon = np.rad2deg(np.arctan2(y[:, 1], y[:, 0]) - OMEGA_TERRA * t)
    lon = (lon + 180) % 360 - 180

    periodo = calcular_periodo_orbital(orbita.semi_eixo_maior, MU_KM)
    return {
        "t_h": r(t[sel] / 3600, 4),
        "x": r(y[sel, 0], 1), "y": r(y[sel, 1], 1), "z": r(y[sel, 2], 1),
        "lat": r(lat[sel], 3), "lon": r(lon[sel], 3),
        "delta_raan_graus": r(np.rad2deg(raan - raan[0])[sel], 5),
        "delta_w_graus": r(np.rad2deg(w - w[0])[sel], 5),
        "periodo_min": round(periodo / 60, 2),
        "inclinacao_graus": round(float(np.rad2deg(orbita.inclinacao)), 3),
        "a_km": round(float(orbita.semi_eixo_maior), 2),
        "e": round(float(orbita.excentricidade), 5),
    }


def decaimento_por_arrasto():
    """Exemplo de Cowell: satélite de 100 kg, CD 2,2, área de 0,785 m², até a altitude crítica de 200 km."""
    r0 = np.array([5872.94259799039, -662.622151725639, 3007.48704280510])
    v0 = np.array([-2.89355195466701, 4.09603530867803, 6.14446573555145])
    orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, MU_KM)
    t, y = propaga_cowell(0, 120 * 86400, orbita, [terra_km, 100, np.pi / 4, 2.2, 200])

    altitude = np.linalg.norm(y[:, :3], axis=1) - RE_KM
    # Apogeu e perigeu de cada revolução (janelas de um período)
    periodo = calcular_periodo_orbital(orbita.semi_eixo_maior, MU_KM)
    dias, apogeus, perigeus = [], [], []
    inicio = 0.0
    while inicio < t[-1]:
        janela = (t >= inicio) & (t < inicio + periodo)
        if janela.sum() > 3:
            dias.append(inicio / 86400)
            apogeus.append(altitude[janela].max())
            perigeus.append(altitude[janela].min())
        inicio += periodo
    passo = max(1, len(dias) // 600)
    return {
        "dias": r(dias[::passo], 3), "apogeu_km": r(apogeus[::passo], 2), "perigeu_km": r(perigeus[::passo], 2),
        "dias_ate_200km": round(float(t[-1] / 86400), 2),
        "apogeu_inicial_km": round(float(apogeus[0]), 1), "perigeu_inicial_km": round(float(perigeus[0]), 1),
        "revolucoes": len(dias),
    }


def modelos_planeta():
    atm = ModeloAtmosferico(JSON_ATM)
    alts = np.concatenate([np.arange(0, 100, 2), np.arange(100, 1001, 10)])
    temperatura, densidade = [], []
    for h in alts:
        res = atm.calcula(h * 1e3, 0, 1, 0)
        temperatura.append(float(np.squeeze(res[0])))
        densidade.append(float(np.squeeze(res[3])))

    terra = terra_em_metros()
    lats = np.linspace(-90, 90, 181)
    g_radial, g_lat = [], []
    for lat in lats:
        gc, gd = calcular_gravidade_axisimetrico(terra.raio_equatorial, np.deg2rad(lat), terra)
        g_radial.append(float(gc))
        g_lat.append(float(gd))
    g_esferica = terra.mut / terra.raio_equatorial ** 2
    return {
        "atmosfera": {"altitude_km": r(alts, 0), "temperatura_K": r(temperatura, 2),
                      "densidade": [float(f"{d:.4e}") for d in densidade]},
        "gravidade": {"latitude_graus": r(lats, 0), "radial": r(g_radial, 6), "latitudinal": r(g_lat, 7),
                      "esferica": round(float(g_esferica), 6)},
    }


def main():
    dados = {"gerado_em": time.strftime("%Y-%m-%d")}
    etapas = [("elementos", elementos_orbitais), ("modelos", modelos_planeta), ("j2", perturbacao_j2),
              ("arrasto", decaimento_por_arrasto), ("lancamento", lancamento)]
    for chave, funcao in etapas:
        inicio = time.time()
        print(f"[{chave}] calculando...")
        dados[chave] = funcao()
        print(f"[{chave}] ok em {time.time() - inicio:.0f} s")

    destino = os.path.join(RAIZ, "docs", "assets", "dados.js")
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    with open(destino, "w", encoding="utf-8") as f:
        f.write("// Arquivo gerado por docs/gerar_dados_pagina.py. Não edite à mão.\n")
        f.write("window.PALO_DADOS = ")
        json.dump(dados, f, ensure_ascii=False, separators=(",", ":"))
        f.write(";\n")
    print(f"Dados gravados em {destino} ({os.path.getsize(destino) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
