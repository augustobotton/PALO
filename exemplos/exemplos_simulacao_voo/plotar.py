import pickle

from src.domain.construtorderesultados.resultados_voo_ascendente import plotaresultados

with open('resposta_simulacao.pkl', 'rb') as f:
    loaded_resposta = pickle.load(f)

with open('simulacao.pkl', 'rb') as f:
    simulacaoa = pickle.load(f)

print(loaded_resposta)
t = loaded_resposta.t
X = loaded_resposta.y
resp = t, X
plotaresultados(resp, simulacaoa)