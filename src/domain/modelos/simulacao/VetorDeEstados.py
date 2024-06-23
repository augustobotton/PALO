import numpy as np
class VetorEstado:
    def __init__(self, V=0, A=0, phi=0, r=0, delta=0):
        self.V = V
        self.A = A
        self.phi = phi
        self.r = r
        self.delta = delta

    def __str__(self):
        return f"Velocidade (m/s): {self.V}, Azimute (rad): {self.A}, Elevação (rad): {self.phi}, Distância radial (m): {self.r}, Latitude (rad): {self.delta}"

    def atualizar(self, V=None, A=None, phi=None, r=None, delta=None):
        """
        Atualiza os atributos do vetor de estado com novos valores fornecidos.

        Parâmetros:
        V (float, opcional): Novo módulo do vetor velocidade.
        A (float, opcional): Novo ângulo de azimute.
        phi (float, opcional): Novo ângulo de elevação.
        r (float, opcional): Nova distância radial.
        delta (float, opcional): Nova latitude.
        """
        if V is not None:
            self.V = V
        if A is not None:
            self.A = A
        if phi is not None:
            self.phi = phi
        if r is not None:
            self.r = r
        if delta is not None:
            self.delta = delta

    def como_array(self):
        """
        Retorna os atributos do vetor de estado como uma array NumPy.
        """
        return np.array([self.V, self.A, self.phi, self.r, self.delta])