from src.domain.modelos.planeta.Planeta import ModeloPlaneta


class ConstrutorPlaneta:
    """
    Uma classe para construir modelos de planetas com características atmosféricas e físicas específicas.

    Esta classe permite a construção passo a passo de um modelo de planeta, definindo vários parâmetros,
    como a variação de temperatura da atmosfera, raio equatorial, velocidade de rotação inercial, gravidade
    e modelo atmosférico, entre outros. Ela garante que todos os parâmetros necessários sejam definidos antes
    da construção do modelo do planeta.
    """

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
        self.massa = None

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

    def com_massa(self, valor: float):
        self.massa = valor
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
            self.modelo_atmosferico,
            self.massa
        )


terra = (ConstrutorPlaneta()
         .com_delta_temperatura_atm(10)
         .com_raio_equatorial(6378.1370)
         .com_velocidade_inercial_de_rotacao(7.2921150e-5)
         .com_gravidade(9.80665)
         .com_mut(398600.4418)
         .com_J2(0.00108263)
         .com_J3(-0.00000254)
         .com_J4(-0.00000161)
         .com_tempo_longitude_celeste_nula(0)
         .com_massa(5.972e24)
         .construir())
