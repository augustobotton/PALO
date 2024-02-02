import numpy as np


class ModeloPropulsivo:
    def __init__(self):
        self.tempo_da_segunda_queima_3_estagio = None
        self.velocidade_de_exaustao = None
        self.impulso_especico_por_estagio = np.array([251, 271, 315])  # s
        self.massa_propelente_estagios_1_2 = np.array([55290, 11058, 0, 0])
        self.massa_propelente_terceiro_estagio = 224.53
        self.duracao_de_queima_primeiro_estagio = 62.0233
        self.duracao_de_queima_segundo_estagio = 64.6105
        self.duracao_total_de_queima_do_terceiro_estagio = 277.5325

        self.tempo_de_espera_separacao_1_2 = 2
        self.tempo_de_espera_separacao_2_3 = 2
        self.tempo_de_espera_separacao_3_4 = 2
        self.tempo_espera_ignicao_2_estagio = 5  # s
        self.tempo_espera_ignicao_3_estagio = 300  # s
        self.tempo_da_primeira_queima_3_estagio = 222.8  # 222.8

    def calculacoisas(self):
        self.tempo_da_segunda_queima_3_estagio = self.duracao_total_de_queima_do_terceiro_estagio - self.tempo_da_primeira_queima_3_estagio
        tempos_de_ignicao = np.array([0, 0, 0, 0])
        instantes_de_fim_de_queima = np.array([0, 0, 0, 0])
        instante_de_separacao = np.array([0, 0, 0, 0])
        instantes_de_fim_de_queima[0] = tempos_de_ignicao[0] + self.duracao_de_queima_primeiro_estagio
        instante_de_separacao[0] = instantes_de_fim_de_queima[0] + self.tempo_de_espera_separacao_1_2
        tempos_de_ignicao[1] = instante_de_separacao[0] + self.tempo_espera_ignicao_2_estagio
        instantes_de_fim_de_queima[1] = tempos_de_ignicao[1] + self.duracao_de_queima_segundo_estagio
        instante_de_separacao[1] = instantes_de_fim_de_queima[1] + self.tempo_de_espera_separacao_2_3
        tempos_de_ignicao[2] = instante_de_separacao[1] + self.tempo_espera_ignicao_3_estagio
        instantes_de_fim_de_queima[2] = tempos_de_ignicao[2] + self.tempo_da_primeira_queima_3_estagio
        tempos_de_ignicao[3] = 1e10
        instantes_de_fim_de_queima[3] = tempos_de_ignicao[3] + self.tempo_da_segunda_queima_3_estagio
        instante_de_separacao[2] = instantes_de_fim_de_queima[3] + self.tempo_de_espera_separacao_3_4

        self.massa_propelente_estagios_1_2[
            2] = self.massa_propelente_terceiro_estagio * self.tempo_da_primeira_queima_3_estagio / self.duracao_total_de_queima_do_terceiro_estagio
        self.massa_propelente_estagios_1_2[
            3] = self.massa_propelente_terceiro_estagio * self.tempo_da_segunda_queima_3_estagio / self.duracao_total_de_queima_do_terceiro_estagio
        self.massa_inicial_do_foguete = np.sum(self.massa_propelente_estagios_1_2) + np.sum(
            self.massa_estrutural_por_estagio) + self.massa_de_carga_util
        self.distancia_radial_inicial = planeta.raio_equatorial_terrestre + h0

    def calculate_thrust_mass(self, t, ti, tq, ts, Isp, mp, ms, m0, g, estagio, mL=None):
        """
        Calculate thrust and mass for a given stage of the rocket.
        :param t: Current time
        :param ti: Ignition times
        :param tq: Burnout times
        :param ts: Stage separation times
        :param Isp: Specific impulse
        :param mp: Propellant mass
        :param ms: Structural mass
        :param m0: Initial mass
        :param g: Gravity
        :param estagio: Stage number (1-indexed)
        :param mL: Payload mass (optional)
        :return: Tuple of thrust and mass
        """
        estagio -= 1  # Adjust for 0-indexed arrays
        if t <= ti[estagio]:
            m = m0 if estagio == 0 else calculate_mass_after_separation(self, ti, tq, ts, mp, ms, m0, estagio)
            ft = 0
        elif t <= tq[estagio]:
            md = -mp[estagio] / (tq[estagio] - ti[estagio])
            m = m0 + md * (t - ti[estagio]) if estagio == 0 else calculate_mass_after_separation(ti, tq, ts, mp, ms, m0,
                                                                                                 estagio) + md * (
                                                                         t - ti[estagio])
            ft = -g * Isp[estagio] * md
        elif t <= ts[estagio]:
            m = calculate_mass_after_separation(ti, tq, ts, mp, ms, m0, estagio)
            ft = 0
        else:
            m = self.calculate_mass_after_separation(ti, tq, ts, mp, ms, m0, estagio + 1)
            m = mL if mL and estagio == len(ti) - 1 else m
            ft = 0
        return ft, m

    def calculate_mass_after_separation(self, ti, tq, ts, mp, ms, m0, estagio):
        """
        Calculate rocket mass after separation of a given stage.
        :param ti: Ignition times
        :param tq: Burnout times
        :param ts: Stage separation times
        :param mp: Propellant mass
        :param ms: Structural mass
        :param m0: Initial mass
        :param estagio: Stage number (1-indexed)
        :return: Mass after stage separation
        """
        mass = m0
        for i in range(estagio):
            mass -= mp[i] + ms[i]
        return mass

    def propulsao_N_estagios(self, tempo, vetor_de_estados):

        N = len(ti)

        # Determine the current stage based on time
        current_stage = sum([tempo > ts[i] for i in range(N)])

        # Calculate thrust and mass for the current stage
        ft, m = calculate_thrust_mass(tempo, ti, tq, ts, Isp, mp, ms, m0, g, current_stage, mL)

        # Additional calculations
        V, A, phi, r, delta = vetor_de_estados
        h = r - Re
        if h < 200e3:
            epsl = 0
            mu = 0
        else:
            _, phii, Ai = Vrel2Vine(V, phi, A, we, r, delta)
            mu = np.arcsin(np.cos(A) * np.cos(phii) * np.sin(Ai) - np.sin(A) * np.cos(phii) * np.cos(Ai))
            epsl = -np.arctan2(-np.cos(phi) * np.sin(phii) + np.sin(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) +
                               np.sin(phi) * np.cos(A) * np.cos(phii) * np.cos(Ai), np.sin(phi) * np.sin(phii) +
                               np.cos(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) + np.cos(phi) * np.cos(A) * np.cos(
                phii) * np.cos(Ai))

        return ft, m, mu, epsl
