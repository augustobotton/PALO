from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita


class ModeloOrbitaEliptica(Orbita):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)