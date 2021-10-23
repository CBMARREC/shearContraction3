# Instalación y setup

⚠️ Esta guía está escrita suponiendo que usas Mac. Los mismos comandos pueden no funcionar en otros sistemas operativos.

Para poder instalar Fenics sin mucho dolor de cabeza lo mejor es usar [Anaconda](https://www.anaconda.com/products/individual).
1. Descarga e instala Anaconda usando el link anterior.
2. Instala y crea un virtual environment para fenics.

    `conda create -n fenics -c conda-forge fenics`

3. Activa el virtual environment.

    `conda activate fenics`
4. Instala el resto de dependencias.

    ```
    conda install -c conda-forge mshr
    conda install hdf5
    pip install scipy
    pip install matplotlib
   ```
  
## IDE

Para mayor facilidad a la hora de leer el código recomiendo usar [Pycharm](https://www.jetbrains.com/pycharm/). Con el link ese puedes descargar e instalar la versión de la comunidad que es gratuita.
