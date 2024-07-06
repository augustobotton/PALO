from setuptools import setup, find_packages

setup(
    name="Pacote de Análise de Lançamento e Órbita (PALO)",
    version="0.1.0",
    author="Augusto Botton Pozzebon",
    author_email="bottonaugusto@gmail.com",
    description="Uma biblioteca para cálculos orbitais e simulação de voo ascendente de foguetes.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/your_username/your_repository",  # Replace with your repository URL
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        'tqdm',

    ],
)