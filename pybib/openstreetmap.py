import osmnx as ox 
import os

def verificar_e_carregar_grafo(nome_cidade, caminho_arquivo):
    
    # Verifica se o arquivo já existe
    filepath = f'./input/shapefile/{caminho_arquivo}_graph_shapefile'
    if os.path.exists(filepath):
        print(f"Carregando o grafo de '{nome_cidade}' a partir do arquivo '{caminho_arquivo}'.")
        G = ox.load_graphml(f'input/{caminho_arquivo}.graphml')
    else:
        print(f"Arquivo '{caminho_arquivo}' não encontrado. Baixando dados e criando o grafo para '{nome_cidade}'.")
        # Baixa os dados da cidade especificada e cria o grafo
        G = ox.graph_from_place(nome_cidade, network_type='drive')
        # Salva o grafo em um arquivo para uso futuro
        ox.save_graph_geopackage(G, filepath)
        ox.save_graphml(G, f'input/{caminho_arquivo}.graphml')
    return G