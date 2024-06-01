import numpy as np

def exemplo_8p4():
    # Dados
    m0 = np.array([1.094266769684211e+05, 5.284772960000001e+04, 19806.1060, 6046.9099])
    mp = np.array([55000, 2.8978709948e+04, 1.279605244274179e+04, 4568.81457])
    mL = 1.124811598501497e+03
    ve = 1.0e+03 * np.array([2.363318181818182, 2.844900000000000, 2.844900000000000, 4.463550000000001])

    # Determina as massas estruturais de cada estagio
    N = len(m0)
    ms = np.zeros(N)
    ms[-1] = m0[-1] - mL - mp[-1]
    ms[:-1] = m0[:-1] - m0[1:] - mp[:-1]

    # Massas ao final da queima de cada estagio
    mf = m0 - mp

    # Mostra as massas
    print('Massa da carga util (kg):', mL)
    print('Massas no inicio da queima de cada estagio:', m0)
    print('Massas no final da queima de cada estagio:', mf)
    print('Massas estruturais de cada estagio:', ms)
    print('Massas de propelente de cada estagio:', mp)

    # Calculo das derivadas de massa de carga util com respeito as massas estruturais
    dmLds = der_mL_dms(mp, ms, mL, ve)
    print('Vetor de derivadas da carga util com respeito as massas estruturais:', dmLds)

    # Calculo das derivadas de massa de carga util com respeito as massas de propelente
    dmLdp = der_mL_dmp(mp, ms, mL, ve)
    print('Vetor de derivadas da carga util com respeito as massas de propelente:', dmLdp)

def der_mL_dms(mp, ms, mL, ve):
    # Numero de estagios
    N = len(mp)

    # Vetores de massa final e inicial por estagio
    mf = np.zeros(N)
    m0 = np.zeros(N)
    for i in range(N):
        m0[i] = mL + np.sum(mp[i:]) + np.sum(ms[i:])
        mf[i] = m0[i] - mp[i]

    # Derivadas da massa de carga util com respeito as massas estruturais
    dmLds = np.zeros(N)
    for k in range(N):
        dmLds[k] = -np.sum(ve[:k+1]*(1/m0[:k+1]-1/mf[:k+1]))/np.sum(ve*(1/m0-1/mf))
    
    return dmLds

def der_mL_dmp(mp, ms, mL, ve):
    # Numero de estagios
    N = len(mp)

    # Vetores de massa final e inicial por estagio
    mf = np.zeros(N)
    m0 = np.zeros(N)
    for i in range(N):
        m0[i] = mL + np.sum(mp[i:]) + np.sum(ms[i:])
        mf[i] = m0[i] - mp[i]

    # Derivadas da massa de carga util com respeito as massas propelente
    dmLdp = np.zeros(N)
    for k in range(N):
        dmLdp[k] = -(ve[k]/m0[k] + np.sum(ve[:k]*(1/m0[:k]-1/mf[:k])))/np.sum(ve*(1/m0-1/mf))

    return dmLdp

# Exemplo de uso da função
exemplo_8p4()