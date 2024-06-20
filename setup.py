from setuptools import setup, find_packages

setup(
    name="your_package_name",  # Replace with your package name
    version="0.1.0",
    author="Augusto Botton Pozzebon",
    author_email="bottonaugusto@gmail.com",
    description="A Python library for orbital mechanics and rocket launch simulations",
    long_description=open('README.md').read(),  # Ensure you have a README.md file
    long_description_content_type='text/markdown',
    url="https://github.com/your_username/your_repository",  # Replace with your repository URL
    packages=find_packages(),  # Automatically finds all the packages in your source directory
    python_requires='>=3.6',  # Replace with your minimum Python version requirement
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",

    ],
)