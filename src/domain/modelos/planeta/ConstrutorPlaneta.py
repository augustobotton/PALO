from src.domain.modelos.planeta.ModeloPlaneta import ModeloPlaneta


class ConstrutorPlaneta:
    def __init__(self):
        # Inicializa os atributos do construtor
        self.delta_temperatura_atm = None
        self.raio_equatorial = None
        self.velocidade_inercial_de_rotacao = None
        self.gravidade = None
        self.mut = None
        self.J2 = None
        self.J3 = None
        self.J4 = None
        self.tempo_longitude_celeste_nula = None
        self.modelo_atmosferico = None

    def com_delta_temperatura_atm(self, valor: float):
        self.delta_temperatura_atm = valor
        return self

    def com_raio_equatorial(self, valor: float):
        self.raio_equatorial = valor
        return self

    def com_velocidade_inercial_de_rotacao(self, valor: float):
        self.velocidade_inercial_de_rotacao = valor
        return self

    def com_gravidade(self, valor: float):
        self.gravidade = valor
        return self

    def com_mut(self, valor: float):
        self.mut = valor
        return self

    def com_J2(self, valor: float):
        self.J2 = valor
        return self

    def com_J3(self, valor: float):
        self.J3 = valor
        return self

    def com_J4(self, valor: float):
        self.J4 = valor
        return self

    def com_tempo_longitude_celeste_nula(self, valor: float):
        self.tempo_longitude_celeste_nula = valor
        return self

    def com_modelo_atmosferico(self, valor: float):
        self.modelo_atmosferico = valor
        return self

    def construir(self) -> ModeloPlaneta:
        # Valida os campos obrigatórios
        if None in (self.delta_temperatura_atm, self.raio_equatorial, self.velocidade_inercial_de_rotacao,
                    self.gravidade, self.mut, self.J2, self.J3, self.J4,
                    self.tempo_longitude_celeste_nula):
            raise ValueError(
                "Todos os parâmetros devem ser fornecidos antes de construir o modelo do planeta.")

        # Retorna a instância de ModeloPlaneta configurada
        return ModeloPlaneta(
            self.delta_temperatura_atm,
            self.raio_equatorial,
            self.velocidade_inercial_de_rotacao,
            self.gravidade,
            self.mut,
            self.J2,
            self.J3,
            self.J4,
            self.tempo_longitude_celeste_nula,
			self.modelo_atmosferico
        )