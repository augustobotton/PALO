from setuptools import setup, find_packages

setup(
    name="palo",
    version="0.1.0",
    author="Augusto Botton Pozzebon",
    author_email="bottonaugusto@gmail.com",
    description="Pacote de Análise de Lançamento e Órbitas (PALO): uma biblioteca para cálculos orbitais e simulação de voo ascendente de foguetes.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/augustobotton/PALO",
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        'tqdm',
        'plotly',
        'astropy'
    ],
)