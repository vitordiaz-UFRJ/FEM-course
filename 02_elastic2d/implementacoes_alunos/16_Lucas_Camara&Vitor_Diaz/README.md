1. Comparação de Métodos de Armazenamento Matricial em MEF: Skyline vs CSR


2. DESCRIÇÃO

Este trabalho implementa e compara dois métodos de armazenamento esparso da matriz de rigidez global K no contexto do Método dos Elementos Finitos: o Skyline (armazenamento em perfil de coluna) e o CSR simétrico (Compressed Sparse Row, triângulo superior).

Ambos os programas são derivados de um código base em Fortran 90, que resolve problemas de elasticidade plana 2D. O solver utilizado nos dois métodos é o Gradiente Conjugado Precondicionado (PCG) com precondicionador diagonal de Jacobi, mantido idêntico para que qualquer diferença de desempenho seja atribuída exclusivamente ao método de armazenamento.

Os resultados mostram que o CSR apresenta vantagem crescente em memória e tempo à medida que o problema cresce, armazenando apenas os coeficientes não-nulos da matriz, enquanto o Skyline guarda toda a faixa do perfil de coluna, incluindo zeros estruturais internos.


3. FORMULAÇÃO TEÓRICA

O problema de elasticidade linear em MEF resulta no sistema:

    K · u = f

onde K é a matriz de rigidez global (simétrica, positiva definida e esparsa), u é o vetor de deslocamentos nodais e f é o vetor de forças nodais. A dimensão de K é neq × neq, onde neq é o número de graus de liberdade livres.

Armazenamento Skyline

Para cada coluna j, armazena todos os coeficientes entre o primeiro não-zero acima da diagonal e a própria diagonal (o "perfil"). Usa dois vetores:

    am(ncs)     — valores reais (real*8, 8 bytes cada)
    jdiag(neq)  — ponteiros da diagonal (integer*4, 4 bytes cada)

    Memória = ncs x 8 + neq x 4 bytes

O tamanho do perfil ncs depende da numeração dos nós — numerações ruins inflam o perfil com zeros desnecessários.

Armazenamento CSR Simétrico

Armazena apenas os coeficientes realmente não-nulos do triângulo superior. Usa três vetores:

    val(nnz)       — valores não-nulos (real*8, 8 bytes cada)
    col_ind(nnz)   — índice de coluna de cada valor (integer*4, 4 bytes cada)
    row_ptr(neq+1) — ponteiro de início de cada linha (integer*4, 4 bytes cada)

    Memória = nnz x 8 + (nnz + neq + 1) x 4 bytes

A fase simbólica (geração de row_ptr e col_ind) é feita em duas passadas sobre a conectividade da malha. Cada linha é ordenada por coluna crescente para permitir busca binária durante a montagem.

Solver PCG

O sistema é resolvido iterativamente. A operação central de cada iteração é o produto matriz-vetor K·d, que no Skyline percorre ncs valores e no CSR percorre apenas nnz valores — essa diferença explica diretamente a diferença de tempo no solver.

4. COMO COMPILAR E RODAR

Requisitos: gfortran

Compilação:

    gfortran -std=legacy Program_Skyline.f90 -o Skyline
    gfortran -std=legacy Program_CSR.f90 -o CSR

A flag -std=legacy suprime avisos sobre construções de Fortran 77 presentes no código base que são válidas mas marcadas como obsoletas no padrão Fortran 2018.

Execução (Windows):

    .\Skyline
    .\CSR

O programa solicita três entradas interativamente:

    Arquivo de dados:
    malha_q4_70x70.dat

    Arquivo de saida:
    resultado.txt

    Arquivo VTK de saida:
    resultado.vtk

5. CASO TESTE

Problema físico:

    Domínio:              chapa quadrada 1,0 m x 1,0 m
    Material:             EPD, E = 1,0, v = 0,3
    Condição de contorno: engaste na borda esquerda (x = 0), ux = uy = 0
    Carregamento:         q = 0,20 kN/m na borda superior (y = 1,0 m)

Malhas testadas:

    Malha    Elementos    Nos      Equacoes
    10x10    100          121      220
    20x20    400          441      840
    30x30    900          961      1.860
    40x40    1.600        1.681    3.280
    50x50    2.500        2.601    5.100
    60x60    3.600        3.721    7.320
    70x70    4.900        5.041    9.940

6. VERIFICAÇÃO

Os dois métodos foram verificados confirmando que produzem exatamente a mesma solução. A norma de energia e os deslocamentos nodais são idênticos em todos os casos testados:

    Malha    Norma de energia (Skyline)    Norma de energia (CSR)
    70x70    19,30283679                   19,30283679

Os arquivos de saída foram comparados linha a linha com o comando diff — zero diferenças em todos os casos. O caso Ex1_quad_25e.dat (malha 5x5) foi rodado nos dois programas e os resultados conferem com a saída original, confirmando que a implementação do CSR preservou integralmente a lógica do código base.


7. REFERÊNCIAS

BATHE, K. J. Finite Element Procedures. 2. ed. Prentice Hall, 2014.
HUGHES, T. J. R. The Finite Element Method. Dover Publications, 2000.
SAAD, Y. Iterative Methods for Sparse Linear Systems. 2. ed. SIAM, 2003.
DAVIS, T. A. Direct Methods for Sparse Linear Systems. SIAM, 2006.


