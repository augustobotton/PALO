from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita


class ConstrutorOrbita:
    def __init__(self):
        self._semi_eixo_maior = None
        self._excentricidade = None
        self._inclinacao = None
        self._raan = None
        self._arg_periastro = None
        self._anomalia_verdadeira = None
        self._tempo_periastro = None
        self._parametro = None

    def com_semi_eixo_maior(self, semi_eixo_maior: float):
        self._semi_eixo_maior = semi_eixo_maior
        return self

    def com_excentricidade(self, excentricidade: float):
        self._excentricidade = excentricidade
        return self

    def com_inclinacao(self, inclinacao: float):
        self._inclinacao = inclinacao
        return self

    def com_raan(self, raan: float):
        self._raan = raan
        return self

    def com_arg_periastro(self, arg_periastro: float):
        self._arg_periastro = arg_periastro
        return self

    def com_anomalia_verdadeira(self, anomalia_verdadeira: float):
        self._anomalia_verdadeira = anomalia_verdadeira
        return self

    def com_tempo_de_periastro(self, tempo_periastro: float):
        self._tempo_periastro = tempo_periastro
        return self

    def com_parametro(self, parametro: float):
        self._parametro = parametro
        return self

    def com_elementos_keplerianos(self, elementos):
        self._semi_eixo_maior = elementos[0]
        self._excentricidade = elementos[1]
        self._inclinacao = elementos[2]
        self._raan = elementos[3]
        self._arg_periastro = elementos[4]
        self._anomalia_verdadeira = elementos[5]
        return self

    def construir(self):
        return Orbita(
            self._semi_eixo_maior,
            self._excentricidade,
            self._inclinacao,
            self._raan,
            self._arg_periastro,
            self._anomalia_verdadeira,
            self._tempo_periastro,
            self._parametro
        )
