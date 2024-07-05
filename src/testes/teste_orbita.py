import unittest
import numpy as np
from src.domain.modelos.orbitas.Orbita import Orbita


class OrbitaTestes(unittest.TestCase):

    def test_orbita_inicializacao_com_parametros_validos(self):
        orbita = Orbita(semi_eixo_maior=7000e3, excentricidade=0.01, inclinacao=0, raan=0,
                        argumento_periastro=0, anomalia_verdadeira=0)
        self.assertEqual(orbita.semi_eixo_maior, 7000e3)
        self.assertEqual(orbita.excentricidade, 0.01)
        self.assertEqual(orbita.inclinacao, 0)
        self.assertEqual(orbita.raan, 0)
        self.assertEqual(orbita.arg_periastro, 0)
        self.assertEqual(orbita.anomalia_verdadeira, 0)

    def test_orbita_define_mu_altera_constante_gravitacional(self):
        orbita = Orbita()
        novo_mu = 42828.375214  # Valor de mu para Marte
        orbita.define_mu(novo_mu)
        self.assertEqual(orbita.mu, novo_mu)

    def test_orbita_calcula_parametro_orbital_corretamente(self):
        orbita = Orbita(semi_eixo_maior=7000e3, excentricidade=0.01)
        parametro_esperado = 7000e3 * (1 - 0.01 ** 2)
        self.assertAlmostEqual(orbita.calcular_parametro_orbital(), parametro_esperado)

    def test_orbita_calcula_apoastro_e_periastro_corretamente(self):
        orbita = Orbita(semi_eixo_maior=7000e3, excentricidade=0.01)
        apoastro_esperado = 7000e3 * (1 + 0.01)
        periastro_esperado = 7000e3 * (1 - 0.01)
        self.assertAlmostEqual(orbita.calcula_apoastro(), apoastro_esperado)
        self.assertAlmostEqual(orbita.calcula_periastro(), periastro_esperado)

    def test_orbita_circular_cria_orbita_com_excentricidade_zero(self):
        orbita = Orbita.circular(semi_eixo_maior=7000e3, inclinacao=0)
        self.assertEqual(orbita.excentricidade, 0)

    def test_orbita_criar_pelo_vetor_de_estado_calcula_parametros_corretamente(self):
        posicao = np.array([7000, 0, 0])
        velocidade = np.array([0, 7.5, 0])
        mu = 398600.4418  # Constante gravitacional da Terra
        orbita = Orbita.criar_pelo_vetor_de_estado(posicao, velocidade, mu)
        self.assertAlmostEqual(orbita.semi_eixo_maior, 7000, delta=1e3)
        self.assertAlmostEqual(orbita.excentricidade, 0, delta=0.1)
        self.assertAlmostEqual(orbita.inclinacao, 0, delta=0.1)

    def test_orbita_calcula_vetor_estado_retorna_posicao_e_velocidade_corretas(self):
        orbita = Orbita(semi_eixo_maior=7000, excentricidade=0.01, inclinacao=0, raan=0,
                        argumento_periastro=0, anomalia_verdadeira=np.pi / 2)
        posicao, velocidade = orbita.calcular_vetor_estado()
        self.assertTrue(np.allclose(posicao, np.array([0, 7000 * (1 - 0.01 ** 2) / (1 + 0.01), 0]), atol=1e3))
        self.assertTrue(np.allclose(velocidade, np.array([-7.5, 0, 0]), atol=1e2))


if __name__ == '__main__':
    unittest.main()
