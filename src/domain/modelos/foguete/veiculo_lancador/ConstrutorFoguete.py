from src.domain.modelos.foguete.aerodinamica.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.estrutura.ModeloEstrutural import ModeloEstrutural
from src.domain.modelos.foguete.propulsao.ModeloPropulsivo import ModeloPropulsivo
from src.domain.modelos.foguete.veiculo_lancador.ModeloFoguete import ModeloFoguete


class ConstrutorDeFoguete:
    def __init__(self):
        self.modelo_propulsivo = None
        self.modelo_estrutural = None
        self.modelo_aerodinamico = None

    def com_modelo_propulsivo(self, modelo_propulsivo: ModeloPropulsivo):
        self.modelo_propulsivo = modelo_propulsivo
        return self

    def com_modelo_estrutural(self, modelo_estrutural: ModeloEstrutural):
        self.modelo_estrutural = modelo_estrutural
        return self

    def com_modelo_aerodinamico(self, modelo_aerodinamico: ModeloAerodinamico):
        self.modelo_aerodinamico = modelo_aerodinamico
        return self

    def construir(self):
        if not self.modelo_propulsivo or not self.modelo_estrutural or not self.modelo_aerodinamico:
            raise ValueError("Todos os modelos (propulsivo, estrutural, aerodinâmico) devem ser fornecidos.")
        return ModeloFoguete(self.modelo_propulsivo, self.modelo_estrutural, self.modelo_aerodinamico)