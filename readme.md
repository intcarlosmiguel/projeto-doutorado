Professor, segue aqui a descrição do problema:

Na minha pesquisa do doutorado eu preciso coletar dados da API do Gooogle do tempo de viagem dentro de Fortaleza, contudo isso tem custo, cerca de 0.005 centavos por requisição no instante de tempo atual. Essa requisição é feita, da forma mais barata, na forma de caminho, eu seleciono dois pontos na rede com no máximo 24 pontos de parada, totalizando 25 arestas no máximo, e mando para a API e ela me devolve o tempo para cada aresta. Então para essa API se eu mandar uma aresta só e um caminho no grafo com no máximo 25 arestas eu pago o mesmo preço. Portanto para otimizar o gasto precisaria não repetir arestas e eu enviar o máximo de caminhos com 25 arestas possível. Consequentemente um problema muito parecido com esse na literatura é encontrar um caminho Euleriano no grafo, pois ao encontrar um caminho que não se repete no grafo, eu poderia contar quantas arestas têm nesse caminho e dividir em pacotes de 25 e mandar para API.

Antes de entrar com mais detalhes do problema vou falar um pouco do grafo de Fortaleza, apesar de que eu possa estudar outras cidades. O grafo de fortaleza é um grafo direcionado e para que ele seja euleriano ele tem que satisfazer as seguintes restrições:

1. Grafo deve ser conexo.
2. Seja Δ o a diferença do número de ligações saindo com número de ligações entrando, ou todos os nós o saldo é zero ou existe apenas um nó com saldo positivo e um nó com saldo negativo.

O grafo da rede tem um total de 96.029 arestas e 36.630 nós para tornar o grafo conexo houve uma redução para 95.370 arestas e 36.456 nós. Nesse caso seriam necessários 3815 requisições, no caso ideal, para cobrir toda a rede. Por fim segue um grafico da distribuição de Δ, quase 90% do grafo tem grau zero e tem 5% de saldo positivo e 5% de saldo negativo e alguns poucos tem saldo 2 e -2. 

Dado essas características eu pesquisei e um problema na qual eu acho mais próximo do que o meu na literatura é o problema de encontar o maior subgrafo euleriano, na qual ja foi mostrado que é NP completo (Ref: W. R. Pulleyblank, A note on graphs spanned by Eulerian graphs, J. Graph Theory 3, 1979, pp. 309–310), mas também isso não é tão problemático para mim, pois com uma situação aproximativa, em vez de passar 3.800 requisições passar 4.000 ou 5.000 não é tão problemático para mim, contudo atualmente estou passando mais de 12.000 requisições e isso já limita muito o que eu posso fazer. Atualmente tenho feito a seguinte heurística: eu inicio em um ponto aleatório do grafo e faço um caminho toda vez que passo em uma aresta eu a retiro do grafo e e quando eu paro em algum nó eu reinicio em outro aleatoriamente, com isso eu conseguir os 12.000 requisições.

Quando eu fui pesquisar sobre esse problema os estudos eram muito mais focados na parte de matematica do que exatamente resolve-lo, mas também não tenho muita experiência nessa área, então estou um pouco perdido e também não sei se eu devo tentar resolver como algoritmo de otimização, mas também não sei exatamente como modelar, enfim não tenho muita experiência na área e preciso de um norte para onde devo começar a tratar.

Att,
Carlos Miguel