from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
import numpy
import copy
import os

from fastaReader import fastaReader

import pandas as pd 

if __name__ == "__main__":
    numeroDeBacterias = 10
    numRandomBacteria = 1
    iteraciones = 10
    tumbo = 10  # numero de gaps a insertar
    nado = 3
    secuencias = list()

    secuencias = fastaReader().seqs
    names = fastaReader().names

    # Hacemos que todas las secuencias sean listas de caracteres
    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])

    globalNFE = 0  # Numero de evaluaciones de la funcion objetivo

    dAttr = 0.1
    wAttr = 0.002
    hRep = dAttr
    wRep = 0.001

    manager = Manager()
    numSec = len(secuencias)
    print("numSec: ", numSec)

    poblacion = manager.list(range(numeroDeBacterias))
    names = manager.list(names)
    NFE = manager.list(range(numeroDeBacterias))

    def poblacionInicial():  # lineal
        # Crece la poblacion al numero de bacterias
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)

    def printPoblacion():
        for i in range(numeroDeBacterias):
            print(poblacion[i])

    operadorBacterial = bacteria(numeroDeBacterias)
    veryBest = [None, None, None]  # indice, fitness, secuencias

    # Registra el tiempo de inicio
    start_time = time.time()

    print("poblacion inicial ...")
    poblacionInicial()

    tasa_mutacion = 0.01  # Probabilidad de mutación para cada carácter

    for it in range(iteraciones):
        print("poblacion inicial creada - Tumbo ...")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)
        
        print("Tumbo Realizado - Aplicando Mutación ...")
        operadorBacterial.mutacion(poblacion, tasa_mutacion)
        
        print("Mutación Aplicada - Cuadrando ...")
        operadorBacterial.cuadra(numSec, poblacion)
        
        print("poblacion inicial cuadrada - Creando granLista de Pares...")
        operadorBacterial.creaGranListaPares(poblacion)
        
        print("granLista creada - Evaluando Blosum Parallel")
        operadorBacterial.evaluaBlosum()  # paralelo
        
        print("blosum evaluado - creando Tablas Atract Parallel...")
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)

        operadorBacterial.creaTablaInteraction()
        print("tabla Interaction creada - creando tabla Fitness")
        operadorBacterial.creaTablaFitness()
        print("tabla Fitness creada ")
        
        globalNFE += operadorBacterial.getNFE()
        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
        
        if (veryBest[0] == None) or (bestFitness > veryBest[1]):  # Remplaza el mejor
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = copy.deepcopy(poblacion[bestIdx])
        
        operadorBacterial.replaceWorst(poblacion, veryBest[0])
        operadorBacterial.resetListas(numeroDeBacterias)

    print("Very Best: ", veryBest)
    # Calcula el tiempo de ejecución
    execution_time = time.time() - start_time
    print("--- %s seconds ---" % execution_time)

    # ---- Código añadido para guardar resultados en hojjas de Excel ----

    # Define los nombres de archivo de Excel
    results_file_name = "resultados.xlsx"
    very_best_file_name = "very_best_results.xlsx"

    # Guarda los resultados de las iteraciones
    new_data = pd.DataFrame({
        "Bacteria": [f"Bacteria {i+1}" for i in range(len(veryBest[2]))],
        f"Iteración {iteraciones}": ["".join(seq) for seq in veryBest[2]]
    })

    if os.path.exists(results_file_name):
        df_existing = pd.read_excel(results_file_name, engine='openpyxl')
        df_updated = pd.concat([df_existing, new_data.set_index('Bacteria')], axis=1)
        df_updated.to_excel(results_file_name, index=True, engine='openpyxl')
    else:
        new_data.to_excel(results_file_name, index=False, engine='openpyxl')

    # Guarda los resultados de Very Best
    best_data = pd.DataFrame({
        "Iteración": [iteraciones],
        "Best Index": [veryBest[0]],
        "Best Fitness": [veryBest[1]],
        "Execution Time (seconds)": [execution_time],
        "BLOSUM Score": [operadorBacterial.blosumScore[veryBest[0]]],
        "Interaction": [operadorBacterial.tablaInteraction[veryBest[0]]],
        "NFE": [operadorBacterial.getNFE()]
    })

    if os.path.exists(very_best_file_name):
        df_existing_best = pd.read_excel(very_best_file_name, engine='openpyxl')
        df_updated_best = pd.concat([df_existing_best, best_data], axis=0)
        df_updated_best.to_excel(very_best_file_name, index=False, engine='openpyxl')
    else:
        best_data.to_excel(very_best_file_name, index=False, engine='openpyxl')

