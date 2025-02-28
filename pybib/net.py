import networkx as nx
import osmnx as ox
import os
import numpy as np
import random

def dict_copy(dicionario):
    copia = {}
    for i in dicionario:
        copia[i] = dicionario[i].copy()
    return copia


def merge_edges(removed_edges):
    
    dicionario = dict_copy(removed_edges)
    removed_edges_after = dict_copy(removed_edges)
    colunas = np.array(list(dicionario.keys()))
    
    for i in colunas[::-1]:
        if(dicionario[i] == 23):
            break
        for j in colunas:
            if(i > j):
                if(len(dicionario[j]) + len(dicionario[i]) <= 23):
                    
                    if(dicionario[j][-1][1] == dicionario[i][0][0]):
                        removed_edges_after[j] = dicionario[j] + dicionario[i] 
                        del removed_edges_after[i]
                        return removed_edges_after,1
                    else:
                        if(dicionario[j][0][0] == dicionario[i][-1][1]):
                            removed_edges_after[j] = dicionario[i] + dicionario[j]
                            del removed_edges_after[i]
                            return removed_edges_after,1
            else:
                break
    return dicionario,0


def remove_isolated_nodes(G):
    """
    Remove todos os nós sem ligações de um grafo direcionado G.

    Parâmetros:
    G (networkx.DiGraph): O grafo direcionado do qual os nós isolados serão removidos.

    Retorna:
    networkx.DiGraph: O grafo modificado sem os nós isolados.
    """
    isolated_nodes = [node for node, degree in G.degree() if degree == 0]
    G.remove_nodes_from(isolated_nodes)
    return G
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
    Group = {}
    
    # Inicializa escolhendo um nó aleatório
    random.seed(seed)
    list_nodes = list(G.nodes)
    current_node = random.choice(list_nodes)
    count = 0
    while True:
        Nodes = [current_node]
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

            if G.out_degree(current_node) == 0:
                list_nodes.remove(current_node)
            
            # Move para o próximo nó
            current_node = aresta_escolhida[1]
            Nodes.append(current_node)

            passos += 1
        print(G.number_of_edges())
        if(len(Nodes) > 1):
            Group[count] = Nodes.copy()
            count += 1
        if G.number_of_edges() == 0:
            break
        if(G.out_degree(current_node) == 0):
            current_node = random.choice(list_nodes)
        
            
    
    return Group
def largest_strongly_connected_component(graph):
    # Encontra todos os componentes fortemente conectados
    scc = list(nx.strongly_connected_components(graph))
    # Encontra o maior componente fortemente conectado
    largest_scc = max(scc, key=len)
    # Extrai o subgrafo correspondente ao maior componente fortemente conectado
    return graph.subgraph(largest_scc).copy()
def remove_residential_edges(G):
    """
    Remove todas as arestas do grafo G que têm o atributo 'highway' igual a 'residential'.
    
    Parâmetros:
    G (networkx.Graph): O grafo do qual as arestas serão removidas.
    
    Retorna:
    networkx.Graph: O grafo modificado sem as arestas 'residential'.
    """
    edges_to_remove = []
    for u, v, k, data in G.edges(keys=True, data=True):
        highway = data.get('highway', [])
        if isinstance(highway, list):
            if 'residential' in highway:
                edges_to_remove.append((u, v, k))
        elif highway == 'residential':
            edges_to_remove.append((u, v, k))
    
    G.remove_edges_from(edges_to_remove)
    return G
