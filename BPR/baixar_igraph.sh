#!/bin/bash

# Script para baixar, compilar e instalar a versão mais recente da biblioteca igraph em C.
# Este script é destrutivo: ele removerá qualquer pasta 'igraph' existente no diretório atual.

# -- Funções de ajuda --
# Função para verificar se um comando existe
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# -- Passo 0: Verificar Pré-requisitos --
echo "Verificando se as ferramentas necessárias estão instaladas..."
PREREQS=("curl" "jq" "wget" "tar" "cmake" "gcc" "g++")
MISSING_PREREQS=()

for tool in "${PREREQS[@]}"; do
    if ! command_exists "$tool"; then
        MISSING_PREREQS+=("$tool")
    fi
done

if [ ${#MISSING_PREREQS[@]} -ne 0 ]; then
    echo "Erro: As seguintes ferramentas não foram encontradas: ${MISSING_PREREQS[*]}"
    echo "Por favor, instale-as para continuar."
    echo "Em Debian/Ubuntu: sudo apt install build-essential curl jq wget tar cmake"
    echo "Em Fedora/CentOS: sudo dnf install gcc-c++ curl jq wget tar cmake"
    exit 1
fi
echo "Todas as ferramentas necessárias foram encontradas."

# --- Configuração ---
REPO="igraph/igraph"
NEW_DIR_NAME="igraph"
echo "--------------------------------------------------------"
echo "Iniciando o processo para obter a última versão de $REPO"
echo "--------------------------------------------------------"

# --- Passo 1: Encontrar e Baixar a Versão Mais Recente ---
echo "Buscando o URL de download da versão mais recente..."
DOWNLOAD_URL=$(curl -s "https://api.github.com/repos/$REPO/releases/latest" | \
               jq -r '.assets[] | select(.name | endswith(".tar.gz")) | .browser_download_url')

if [ -z "$DOWNLOAD_URL" ]; then
    echo "Erro: Não foi possível encontrar o URL de download. Verifique a API do GitHub ou sua conexão."
    exit 1
fi
echo "URL encontrado: $DOWNLOAD_URL"

FILENAME=$(basename "$DOWNLOAD_URL")
echo "Baixando $FILENAME..."
wget -q --show-progress "$DOWNLOAD_URL"
if [ $? -ne 0 ]; then echo "Erro no download."; exit 1; fi

# --- Passo 2: Extrair e Renomear ---
echo "Extraindo arquivo..."
tar -xzf "$FILENAME"
if [ $? -ne 0 ]; then echo "Erro na extração."; exit 1; fi

ORIGINAL_DIR_NAME=${FILENAME%.tar.gz}

echo "Limpando diretórios antigos (se a pasta '$NEW_DIR_NAME' já existir, será removida)..."
if [ -d "$NEW_DIR_NAME" ]; then
    rm -rf "$NEW_DIR_NAME"
fi

echo "Renomeando '$ORIGINAL_DIR_NAME' para '$NEW_DIR_NAME'..."
mv "$ORIGINAL_DIR_NAME" "$NEW_DIR_NAME"
if [ $? -ne 0 ]; then echo "Erro ao renomear o diretório."; exit 1; fi

# --- Passo 3: Compilar e Instalar ---
echo "--------------------------------------------------------"
echo "Iniciando compilação e instalação automatizada..."
echo "--------------------------------------------------------"

# Entra no diretório principal do igraph
cd "$NEW_DIR_NAME" || { echo "Erro: não foi possível entrar no diretório '$NEW_DIR_NAME'"; exit 1; }

echo "[1/5] Criando diretório de compilação 'build'..."
mkdir build && cd build || { echo "Erro ao criar ou entrar no diretório 'build'"; exit 1; }

echo "[2/5] Configurando o projeto com CMake..."
cmake ..
if [ $? -ne 0 ]; then echo "Erro na configuração com CMake."; exit 1; fi

echo "[3/5] Compilando a biblioteca (isso pode levar alguns minutos)..."
cmake --build .
if [ $? -ne 0 ]; then echo "Erro durante a compilação."; exit 1; fi

echo "[4/5] Executando os testes de verificação..."
cmake --build . --target check
if [ $? -ne 0 ]; then echo "Atenção: Os testes falharam. A instalação pode estar instável."; fi # Não sai, mas avisa

echo "[5/5] Instalando a biblioteca no sistema (será solicitada a senha de sudo)..."
sudo cmake --install .
if [ $? -ne 0 ]; then echo "Erro durante a instalação com sudo."; exit 1; fi

# --- Passo 4: Limpeza Final ---
echo "Atualizando o cache de bibliotecas compartilhadas do sistema..."
sudo ldconfig

echo "Limpando arquivos de download..."
# Volta para o diretório original para remover o .tar.gz
cd ../..
rm "$FILENAME"

echo ""
echo "--------------------------------------------------------"
echo "SUCESSO! A versão mais recente do igraph foi instalada."
echo "--------------------------------------------------------"