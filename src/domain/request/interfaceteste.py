import tkinter as tk
from tkinter import messagebox

from src.domain.modelos.foguete.ModeloEstrutural import ModeloEstruturalBuilder
from src.domain.modelos.foguete.ModeloFoguete import *
from src.domain.request.simula import *


def iniciar_simulacao():
	try:
		planeta = (ConstrutorDePlanetas()
		           .com_delta_temperatura_atm(10)
		           .com_raio_equatorial(6378.1370e3)
		           .com_velocidade_inercial_de_rotacao(7.2921150e-5)
		           .com_gravidade_padrao_nivel_do_mar(9.80665)
		           .com_mut(3.986004418e14)
		           .com_J2(0.00108263)
		           .com_J3(-0.00000254)
		           .com_J4(-0.00000161)
		           .com_tempo_longitude_celeste_nula(0)
		           .construir())

		base_de_lancamento = (ConstrutorBaseDeLancamento()
		                      .com_altitude_base(0)
		                      .com_comprimento_trilho(50)
		                      .construir())
		modelo_estrutural = ModeloEstruturalBuilder().com_massa_estrutural_por_estagio([0.5, 0.5, 0.5,
		                                                                                0.5]).com_massa_de_carga_util(
			0.5).com_comprimento_carga_util(1).com_comprimento_total_do_foguete(50

		).com_comprimento_sem_1_estagio(30).com_comprimento_sem_2_estagio(15).construir(

		)

		modelo_propulsivo = ConstrutorModeloPropulsivo()
		foguete = ConstrutorDeFoguete().com_modelo_estrutural(modelo_estrutural).com_modelo_propulsivo(

		).com_modelo_aerodinamico().construir()

		tempo_simulacao = int(entry_tempo_simulacao.get())
		velocidade_inicial = float(entry_velocidade_inicial.get())
		angulo_elevacao_inicial = float(entry_angulo_elevacao_inicial.get())
		inclinacao_orbita = float(entry_inclinacao_orbita.get())
		altitude_geo_sincrona = float(entry_altitude_geo_sincrona.get())

		simulacao = Simulacao(planeta, base_de_lancamento, foguete, tempo_simulacao,
		                      velocidade_inicial, angulo_elevacao_inicial, inclinacao_orbita,
		                      altitude_geo_sincrona)
		tempos, estados = simulacao.simular()
		print(estados)
		messagebox.showinfo("Simulação Completa", "A simulação foi concluída com sucesso!")
	except Exception as e:
		messagebox.showerror("Erro na Simulação", str(e))


# Configuração da interface gráfica
root = tk.Tk()
root.title("Simulação de Lançamento de Foguete")

tk.Label(root, text="Tempo de Simulação (s)").grid(row=0)
tk.Label(root, text="Velocidade Inicial (m/s)").grid(row=1)
tk.Label(root, text="Ângulo de Elevação Inicial (graus)").grid(row=2)
tk.Label(root, text="Inclinação da Órbita (graus)").grid(row=3)
tk.Label(root, text="Altitude Geo-Síncrona (m)").grid(row=4)

entry_tempo_simulacao = tk.Entry(root)
entry_velocidade_inicial = tk.Entry(root)
entry_angulo_elevacao_inicial = tk.Entry(root)
entry_inclinacao_orbita = tk.Entry(root)
entry_altitude_geo_sincrona = tk.Entry(root)

entry_tempo_simulacao.grid(row=0, column=1)
entry_velocidade_inicial.grid(row=1, column=1)
entry_angulo_elevacao_inicial.grid(row=2, column=1)
entry_inclinacao_orbita.grid(row=3, column=1)
entry_altitude_geo_sincrona.grid(row=4, column=1)

tk.Button(root, text='Iniciar Simulação', command=iniciar_simulacao).grid(row=5, column=0, columnspan=2,
                                                                          pady=10)

root.mainloop()
