import osmnx as ox
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from shapely import wkt
from shapely.geometry import Point, MultiPolygon
import networkx as nx
import random
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from pybib.api import *
import time
from datetime import datetime
import copy
from collections import OrderedDict

def merge_edges(removed_edges):
    
    c = 1
    while(c != 0):
        c = 0
        colunas = np.array(list(removed_edges.keys()))
        for i in colunas[::-1]:
            if(removed_edges[i] == 24):
                break
            for j in removed_edges.keys():
                if(len(removed_edges[j])!=24):
                    if(len(removed_edges[j]) + len(removed_edges[i]) <= 24):
                        if(removed_edges[j][-1][1] == removed_edges[i][0][0]):
                            removed_edges[j] = removed_edges[j] + removed_edges[i]
                            del removed_edges[i]
                            c += 1
                            colunas = np.array(list(removed_edges.keys()))
                            break
                        if(removed_edges[i][-1][1] == removed_edges[j][0][0]):
                            removed_edges[j] = removed_edges[i] + removed_edges[j]
                            del removed_edges[i]
                            c += 1
                            colunas = np.array(list(removed_edges.keys()))
                            break
        print('volta!',c)
    return removed_edges

def caminhar_e_remover(Grafo, max_passos=24,seed = 42):
    """
    Percorre o grafo G removendo arestas de forma aleatória até completar
    max_passos ou até não ser possível continuar. Se completar max_passos,
    continua a partir do último nó. Caso contrário, remove nós com grau zero
    e reinicia o processo a partir de um novo nó aleatório.
    
    Args:
        G (networkx.DiGraph): Grafo direcionado.
        max_passos (int): Número máximo de passos por percurso.
        
    Returns:
        list: Lista de arestas removidas na ordem.
    """
    G = Grafo.copy()
    arestas_removidas = {}
    
    # Verifica se o grafo tem arestas
    if G.number_of_edges() == 0:
        return arestas_removidas
    
    # Inicializa escolhendo um nó aleatório
    random.seed(seed)
    current_node = random.choice(list(G.nodes))
    count =0
    while True:
        removidas_atual = []
        passos = 0
        
        while passos < max_passos:
            arestas_saindo = list(G.out_edges(current_node))
            
            # Se não houver arestas saindo, interrompe o percurso
            if not arestas_saindo:
                break
            
            # Escolhe uma aresta de saída aleatória
            aresta_escolhida = random.choice(arestas_saindo)
            
            # Remove a aresta do grafo
            G.remove_edge(*aresta_escolhida)
            removidas_atual.append(aresta_escolhida)
            
            # Move para o próximo nó
            current_node = aresta_escolhida[1]
            passos += 1
        
        # Adiciona as arestas removidas ao total
        arestas_removidas[count] = removidas_atual
        count += 1
        if passos == max_passos:
            # Continua a partir do current_node
            if G.number_of_edges() == 0:
                break
            continue
        else:
            # Remove nós que ficaram com grau zero
            nos_para_remover = [no for no in G.nodes if G.in_degree(no) == 0 and G.out_degree(no) == 0]
            G.remove_nodes_from(nos_para_remover)
            
            # Verifica se ainda há arestas no grafo
            if G.number_of_edges() == 0:
                break
            
            # Escolhe um novo nó aleatório para continuar
            current_node = random.choice(list(G.nodes))
    
    return arestas_removidas

def verificar_e_carregar_grafo(nome_cidade, caminho_arquivo):
    
    filepath = f'./input/shapefile/{caminho_arquivo}_graph_shapefile'
    if os.path.exists(filepath):
        print(f"Carregando o grafo de '{nome_cidade}' a partir do arquivo '{caminho_arquivo}'.")
        G = ox.load_graphml(f'input/{caminho_arquivo}.graphml')
    else:
        print(f"Arquivo '{caminho_arquivo}' não encontrado. Baixando dados e criando o grafo para '{nome_cidade}'.")
        G = ox.graph_from_place(nome_cidade, network_type='drive')
        ox.save_graph_geopackage(G, filepath)
        ox.save_graphml(G, f'input/{caminho_arquivo}.graphml')
    return G

def largest_strongly_connected_component(graph):
    # Encontra todos os componentes fortemente conectados
    scc = list(nx.strongly_connected_components(graph))
    # Encontra o maior componente fortemente conectado
    largest_scc = max(scc, key=len)
    # Extrai o subgrafo correspondente ao maior componente fortemente conectado
    return graph.subgraph(largest_scc).copy()

# Definindo o nome da cidade e o caminho do arquivo
nome_cidade = "Fortaleza, Ceará ,Brasil"
caminho_arquivo = "fortaleza"

# Chama a função para verificar e carregar (ou criar) o grafo
G = verificar_e_carregar_grafo(nome_cidade, caminho_arquivo)

print(f'Antes da Remoção: {len(G.edges())}')

Grafo = largest_strongly_connected_component(G)
Grafo.remove_edges_from(nx.selfloop_edges(Grafo))

print(f'Depos da Remoção: {len(Grafo.edges())}')
removed_edges = caminhar_e_remover(Grafo)
removed_edges = {i: removed_edges[i] for i in removed_edges if(len(removed_edges[i]) != 0)}
dicionario_ordenado = OrderedDict(
    sorted(
        removed_edges.items(),
        key=lambda item: len(item[1]),
        reverse=True
    )
)
removed_edges = dict(dicionario_ordenado)
dicionario_copia = removed_edges.copy()

dicionario_copia = copy.deepcopy(removed_edges)
dic = merge_edges(dicionario_copia)
dic = {i: dic[i] for i in dic if(len(dic[i]) != 0)}
dicionario_ordenado = OrderedDict(
    sorted(
        dic.items(),
        key=lambda item: len(item[1]),
        reverse=True
    )
)

dic = dict(dicionario_ordenado)

with open("token.txt", "r") as file:
    token = file.read().strip()
count = 0
while True:
    agora = datetime.now()
    if agora.hour == 12 and agora.minute == 0:
        tempo = 0
        for j in dic:
            removed_edge = dic[j]
            inicio = time.time()
            r = request_api(dic,0,f'data_{int(agora.hour)}_{int(agora.minute)}.txt')
            fim = time.time()
            tempo += inicio - fim
            if(r == -1):
                continue
            if(r == 0):
                break

            if((count+1)%290 == 0):
                print(f'{count+1}/{len(dic.keys())}')
                if(tempo < 60):
                    time.sleep(60)
                tempo -= 60
            count += 1

        break