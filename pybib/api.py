import requests
import numpy as np
import time
import pandas as pd
from datetime import datetime
# Função para carregar a chave da API de um arquivo
def carregar_api_key(file_path):
    with open(file_path, 'r') as f:
        return f.readline().strip()

# Função para calcular o tempo entre múltiplas coordenadas (origem e destino) em tempo real
def calcular_tempos_viagem(origens, destinos, chave_api):
    url = "https://maps.googleapis.com/maps/api/distancematrix/json"
    origem_str = "|".join([f"{lat},{lon}" for lat, lon in origens])
    destino_str = "|".join([f"{lat},{lon}" for lat, lon in destinos])
    
    params = {
        "origins": origem_str,
        "destinations": destino_str,
        "key": chave_api,
        "departure_time": "now",  # Tempo real
        "mode": "driving"  # Modo de transporte
    }

    resposta = requests.get(url, params=params)
    dados = resposta.json()

    if dados["status"] == "OK":
        tempos_resultados = []
        for i in range(len(origens)):
            linha_tempos = []
            for j in range(len(destinos)):
                elemento = dados["rows"][i]["elements"][j]
                if elemento["status"] == "OK":
                    tempo_estimado = elemento["duration"]["value"]  # Tempo em segundos
                    linha_tempos.append(tempo_estimado)
                else:
                    linha_tempos.append(np.nan)  # Valor indefinido para casos de erro
            tempos_resultados.append(linha_tempos)
        
        return np.diagonal(np.array(tempos_resultados))
    else:
        raise ValueError("Erro na resposta da API")

def get_time(hora,minuto,Grafo,removed_edges):

    with open("token.txt", "r") as file:
        token = file.read().strip()

    while(True):
        agora = datetime.now()
        if agora.hour == hora and agora.minute == minuto:
            removed_edge = removed_edges[i]
            nodes = [removed_edge[0][0]] + np.array(removed_edge)[:,1].tolist()
            nodes = [(Grafo.nodes[node]['x'],Grafo.nodes[node]['y']) for node in nodes]
            edges= []
            for i in range(len(nodes)-1):
                edges.append(nodes[i]+nodes[i+1])
            edges = np.array(edges)
            if request_api(nodes,edges, token) == 0:
                break
            print(i)
            if((i+1)%290 == 0):
                print('Tempo')
                time.sleep(62)
            break
def request_api(removed_edge,coord,filename):
    nodes = [removed_edge[0][0]] + np.array(removed_edge)[:,1].tolist()
    if(coord ==0 ):
        nodes = [(Grafo.nodes[node]['x'],Grafo.nodes[node]['y']) for node in nodes]
    else:
        nodes = [(Grafo.nodes[node]['y'],Grafo.nodes[node]['x']) for node in nodes]
    edges= []
    for i in range(len(nodes)-1):
        edges.append(nodes[i]+nodes[i+1])
    edges = np.array(edges)
    coordinates_str = ";".join([f"{lon},{lat}" for lon, lat in nodes])

    # URL da API Directions do Mapbox
    url = f"https://api.mapbox.com/directions/v5/mapbox/driving/{coordinates_str}?geometries=geojson&access_token={token}"

    # Fazendo a requisição
    response = requests.get(url)

    if response.status_code == 200:
        matrix_data = response.json()
        if((matrix_data['code'] == 'NoSegment') or (matrix_data['code'] == 'NoRoute')):
            print('Deu erro na coleta!')
            return -1
        else:
            durations = pd.DataFrame(matrix_data['routes'][0]['legs'])['duration'].values
            tempo = np.hstack((edges,durations.reshape(-1,1)))
            with open(filename, 'a') as f:
                np.savetxt(f, tempo, fmt='%f', delimiter=',')
            return 1
    else:
        print("Erro na requisição:", response.status_code, response.text)
        return 0