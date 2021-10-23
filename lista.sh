#!/bin/bash
# Lista de programas para correr
#
# Standard commands block:
#
# mkdir ./results/vtkFiles/estudioResolucion/res07
# python 4real7.py
# mv ./results/vtkFiles/detF* ./results/vtkFiles/estudioResolucion/res07
# mv ./results/vtkFiles/dens* ./results/vtkFiles/estudioResolucion/res07
# mv ./results/vtkFiles/u* ./results/vtkFiles/estudioResolucion/res07
# cp -r ./results/vtkFiles/saved_functions ./results/vtkFiles/estudioResolucion/res07
#
mkdir ./results/vtkFiles/noContract/pasos001
python shearCellContract001.py
mv ./results/vtkFiles/detF* ./results/vtkFiles/noContract/pasos001
mv ./results/vtkFiles/dens* ./results/vtkFiles/noContract/pasos001
mv ./results/vtkFiles/u* ./results/vtkFiles/noContract/pasos001
cp -r ./results/vtkFiles/saved_functions ./results/vtkFiles/noContract/pasos001
#
mkdir ./results/vtkFiles/noContract/pasos0005
python shearCellContract0005.py
mv ./results/vtkFiles/detF* ./results/vtkFiles/noContract/pasos0005
mv ./results/vtkFiles/dens* ./results/vtkFiles/noContract/pasos0005
mv ./results/vtkFiles/u* ./results/vtkFiles/noContract/pasos0005
cp -r ./results/vtkFiles/saved_functions ./results/vtkFiles/noContract/pasos0005
#
mkdir ./results/vtkFiles/noContract/pasos0001
python shearCellContract0001.py
mv ./results/vtkFiles/detF* ./results/vtkFiles/noContract/pasos0001
mv ./results/vtkFiles/dens* ./results/vtkFiles/noContract/pasos0001
mv ./results/vtkFiles/u* ./results/vtkFiles/noContract/pasos0001
cp -r ./results/vtkFiles/saved_functions ./results/vtkFiles/noContract/pasos0001
