! ======================================================================
!     Programa MEF - Metodo dos Elementos Finitos
!     Elasticidade Plana 2D (Estado Plano de Tensoes e Deformacoes)
!
!     Elementos disponiveis:
!       1 - T3  Estado Plano de Deformacoes (EPD)
!       2 - T3  Estado Plano de Tensoes    (EPT)
!       3 - Q4  Estado Plano de Deformacoes (EPD)
!       4 - Q4  Estado Plano de Tensoes    (EPT)
!       5 - T3  Problema de potencial
!       6 - T3  Axissimetrico
!
!     Solver: Gradiente Conjugado Precondicionado (PCG)
!             com precondicionador diagonal (Jacobi)
!
!     Armazenamento da matriz de rigidez: CSR (Compressed Sparse Row)
!     Apenas o TRIANGULO SUPERIOR e armazenado (matriz simetrica),
!     para comparacao justa de memoria com a versao Skyline (que
!     tambem guarda so "metade" da matriz, dentro do perfil).
!
!     >>> VERSAO CSR (derivada de programMEF_f90.f90) <<<
!     Todas as rotinas de leitura de dados, elementos finitos
!     (elmt01..elmt01_axi), solver PCG e escrita de resultados
!     SAO IDENTICAS a versao Skyline. As UNICAS rotinas alteradas
!     foram as relacionadas ao armazenamento da matriz global:
!       profil   -> profil_csr   (fase simbolica: gera row_ptr/col_ind)
!       addstf   -> addstf_csr   (montagem: busca binaria na linha)
!       matvec   -> matvec_csr   (produto matriz-vetor via CSR)
!     e a subrotina contr, que precisou ajustar a alocacao de
!     memoria (agora usa row_ptr, col_ind, val em vez de jdiag, a).
! ======================================================================
      program mef
      implicit none
      integer, parameter :: npos = 1000000
      integer :: m(npos)
      integer :: max
      common m
      common /size/ max
      real*8 :: a(1)
      equivalence (m(1),a(1))
      max = npos/2
      call contr(m,a)
      stop
      end program mef

! ======================================================================
      subroutine contr(m,a)
!     Subrotina controladora principal.
!     Orquestra todas as etapas do calculo MEF.
! ======================================================================
      implicit none
      integer :: max
      common /size/ max
      integer :: m(*)
      real*8  :: a(1)
      real*8  :: dot
      character*80 :: fname
      external dot
      integer :: nin, nout, nnode, numel, numat, nen, ndf, ndm
      integer :: nen1, i1, i2, i3, i4, i5, i6
      integer :: j7, j8, j9, j10, j11, j12, j13, j14
      integer :: j15, j16, j17, j18, j19, j20, jtmp
      integer :: nst, neq, nnz, nvtk
      real*8  :: t0, t1, t_simbolico, t_montagem, t_solver
      real*8  :: mem_reais, mem_inteiros
      real*8  :: mem_bytes, mem_kb, mem_mb
!
!     Abertura de arquivos:
!
      nin  = 1
      nout = 2
      print*, 'Arquivo de dados:'
      read(*,'(a)') fname
      open(nin,file=fname)
      print*, 'Arquivo de saida:'
      read(*,'(a)') fname
      open(nout,file=fname)
!
!     Leitura das variaveis principais:
!       nnode = numero de nos
!       numel = numero de elementos
!       numat = numero de materiais
!       nen   = nos por elemento (3=triangulo, 4=quadrilatero)
!       ndf   = graus de liberdade por no (2 para 2D)
!       ndm   = dimensao do problema (2 para 2D)
!
      read(nin,*) nnode,numel,numat,nen,ndf,ndm
!
!     Verificacao basica de consistencia dos parametros de entrada:
!
      if (nnode .le. 0 .or. numel .le. 0) then
        print*, '*** ERRO: nnode e numel devem ser positivos.'
        stop
      endif
      if (nen .ne. 3 .and. nen .ne. 4) then
        print*, '*** ERRO: nen deve ser 3 (T3) ou 4 (Q4). Lido:', nen
        stop
      endif
      if (ndm .ne. 2) then
        print*, '*** ERRO: programa valido apenas para 2D (ndm=2)'
        stop
      endif
!
!     Alocacao de memoria (vetor unico m/a):
!
!      1       i1       i2       i3       i4     i5       i6
!      ---------------------------------------------------------
!      |   e   |   ie   |   ix   |   id   |   x   |   f   |
!      ---------------------------------------------------------
!
      nen1 = nen+1
      i1 = ( 1 + numat*10 - 1)*2 + 1
      i2 =  i1 + numat
      i3 =  i2 + numel*nen1
      i4 = (i3 + nnode*ndf)/2 + 1
      i5 =  i4 + nnode*ndm
      i6 =  i5 + nnode*ndf
      call mem(i6)
!
!     Leitura dos dados:
!
      call rdata(a,m(i1),m(i2),m(i3),a(i4),a(i5),nnode,numel,numat, &
                 nen,ndm,ndf,nin)
!
!     Numeracao das equacoes livres (nos nao restritos):
!
      call numeq(m(i3),nnode,ndf,neq)
      print*, 'Numero de equacoes (graus de liberdade livres):', neq
!
!     Alocacao de memoria para u e row_ptr (CSR):
!       row_ptr tem neq+1 posicoes (convencao CSR padrao)
!
!      i4    i5    i6    j7
!      ------------------------------
!      |  x  |  f  |  u  | row_ptr  |
!      ------------------------------
!
      j7 = (i6 + neq - 1)*2 + 1
!
!     Fase simbolica CSR (1a passada): gera row_ptr e descobre nnz
!     (numero de coeficientes nao-nulos no triangulo superior).
!     m(jtmp) e um vetor de marcacao TEMPORARIO de tamanho neq,
!     alocado em uma regiao PROPRIA (apos row_ptr), que so sera
!     reaproveitada para outra coisa depois que soubermos nnz.
!
      jtmp = j7 + neq + 1
      call mem(jtmp + neq - 1)
      call cpu_time(t0)
      call profil_csr(m(i2),m(i3),m(j7),m(jtmp),numel,nen,ndf,neq,nnz)
      print*, 'Numero de nao-nulos no triangulo superior (nnz):', nnz
!
!     Alocacao de memoria para col_ind e val (CSR) e arrays locais:
!     col_ind e alocado numa posicao NOVA, apos o vetor de marcacao
!     temporario (que nao e mais reaproveitado, por seguranca).
!
!      j7        jtmp           j8        j9      j10   j11  j12   j13   j14
!      -------------------------------------------------------------------
!      | row_ptr | marca(temp) | col_ind |  val  |  xl  |ul  | fl  | sl  | ld |
!      -------------------------------------------------------------------
!
      j8  = jtmp + neq
      nst = nen*ndf
      j9  = (j8  + nnz)/2 + 1
      j10 =  j9  + nnz
      j11 =  j10 + nen*ndm
      j12 =  j11 + nst
      j13 =  j12 + nst
      j14 = (j13 + nst*nst - 1)*2 + 1
      j15 = (j14 + nst)/2 + 1
      call mem(j15)
!
!     Preenchimento de col_ind (2a passada da fase simbolica):
!
      call fillcol_csr(m(i2),m(i3),m(j7),m(j8),m(jtmp),numel,nen,ndf, &
                       neq,nnz)
      call cpu_time(t1)
      t_simbolico = t1 - t0
!
!     Alocacao de memoria para vetores do solver PCG:
!       x=solucao, r=residuo, z=precond*r, prec=precond, d=direcao
!
!      j15   j16   j17   j18   j19   j20
!      -----------------------------------
!      |  u  |  x  |  r  |  z  | prec| d  |
!      -----------------------------------
!
      j16 = j15 + neq
      j17 = j16 + neq
      j18 = j17 + neq
      j19 = j18 + neq
      j20 = j19 + neq
      call mem(j20 + neq)
!
!     Forcas nodais equivalentes:
!       Transfere forcas prescritas para o vetor global u
!
      call pload(m(i3),a(i5),a(i6),nnode,ndf)
!
!     Montagem da matriz de rigidez global K (formato CSR) e
!     correcao do vetor de forcas:
!       isw=2  => calcula rigidez e forcas internas do elemento
!       afl=.true. => corrige o vetor de forcas (deslocamentos prescritos)
!       bfl=.true. => monta a matriz de rigidez global
!
      call cpu_time(t0)
      call pform_csr(a,m(i1),m(i2),m(i3),a(i4),a(i5),a(i6), &
                 m(j7),m(j8),a(j9),a(j10),a(j11),a(j12),a(j13),m(j14), &
                 numel,ndm,ndf,nen,nst,2,.true.,.true.,neq)
      call cpu_time(t1)
      t_montagem = t1 - t0
!
!     Resolucao do sistema de equacoes K*u = f
!     via Gradiente Conjugado Precondicionado (PCG):
!
      print*, 'Iniciando solver PCG...'
      call cpu_time(t0)
      call pcg_csr(a(j9),m(j7),m(j8),a(i6),a(j16),a(j17),a(j18), &
               a(j19),a(j20),neq)
      call cpu_time(t1)
      t_solver = t1 - t0
!
!     Relatorio de desempenho (tempo e memoria) - ARMAZENAMENTO CSR:
!       Memoria da matriz:
!         - nnz valores real*8 (8 bytes cada) = vetor val (coeficientes)
!         - (neq+1 + nnz) valores integer*4 (4 bytes cada)
!             neq+1 = vetor row_ptr (ponteiros de linha)
!             nnz   = vetor col_ind (indices de coluna)
!
      mem_reais    = dble(nnz)
      mem_inteiros = dble(neq+1) + dble(nnz)
      mem_bytes = mem_reais * 8.d0 + mem_inteiros * 4.d0
      mem_kb    = mem_bytes / 1024.d0
      mem_mb    = mem_kb   / 1024.d0
      print*, '============================================='
      print*, ' RELATORIO DE DESEMPENHO - CSR (triang. sup.)'
      print*, '============================================='
      print*, ' Numero de equacoes (neq)        :', neq
      print*, ' Nao-nulos armazenados (nnz)     :', nnz
      print*, ' --- Memoria da matriz K ---'
      print*, ' Reais real*8  : val     =', nnz, 'x 8 bytes =', &
              nint(mem_reais*8.d0), 'bytes'
      print*, ' Inteiros int*4: row_ptr =', neq+1, 'x 4 bytes'
      print*, ' Inteiros int*4: col_ind =', nnz,   'x 4 bytes'
      print*, '   Total inteiros        =', nint(mem_inteiros), &
              'x 4 bytes =', nint(mem_inteiros*4.d0), 'bytes'
      print*, ' Total memoria matriz K        :', nint(mem_bytes), 'bytes'
      write(*,'(a,f12.3,a)') '                                  ', &
              mem_kb, ' KB'
      write(*,'(a,f12.6,a)') '                                  ', &
              mem_mb, ' MB'
      print*, ' --- Tempo ---'
      print*, ' Tempo fase simbolica (profil+fillcol):', t_simbolico, 's'
      print*, ' Tempo montagem      (pform_csr)  :', t_montagem,  's'
      print*, ' Tempo solver PCG    (pcg_csr)    :', t_solver,    's'
      print*, ' Tempo total (simb+mont+solver)   :', &
              t_simbolico+t_montagem+t_solver, 's'
      print*, '============================================='
!
!     Impressao dos resultados (formato texto):
!
      call wdata(m(i2),m(i3),a(i4),a(i5),a(i6),nnode,numel,ndm,nen, &
                 ndf,nout)
!
!     Impressao dos resultados (formato VTK para ParaView):
!
      nvtk = 3
      print*, 'Arquivo VTK de saida (ex: resultado.vtk):'
      read(*,'(a)') fname
      open(nvtk,file=fname,status='replace',action='write')
      call wvtk(m(i2),m(i3),a(i4),a(i5),a(i6),nnode,numel,ndm,nen, &
                ndf,nvtk)
      close(nvtk)
      print*, 'Arquivo VTK gerado com sucesso.'
      print*, 'Calculo concluido. Resultados gravados com sucesso.'
      return
      end subroutine contr

! ======================================================================
      subroutine rdata(e,ie,ix,id,x,f,nnode,numel,numat,nen,ndm,ndf, &
                        nin)
!     Leitura dos dados de entrada:
!       e  = propriedades dos materiais  (real*8, 10 x numat)
!       ie = tipo de elemento por material (integer, numat)
!       ix = conectividade + material    (integer, nen+1 x numel)
!       id = condicoes de contorno       (integer, ndf x nnode)
!       x  = coordenadas nodais          (real*8, ndm x nnode)
!       f  = forcas nodais               (real*8, ndf x nnode)
! ======================================================================
      implicit none
      integer :: nnode,numel,numat,ndm,nen,ndf,nin
      integer :: ie(*),ix(nen+1,*),id(ndf,*)
      integer :: iaux(6)
      real*8 :: e(10,*),x(ndm,*),f(ndf,*),aux(6),dum(1)
      integer :: i,j,k,n,ma,iel
!
!     Inicializa id (todos livres) e f (sem forcas):
!
      do 5 n = 1, nnode
        do 5 j = 1, ndf
          id(j,n) = 0
          f(j,n)  = 0.d0
    5 continue
!
!     Leitura das propriedades dos materiais:
!       Linha: ma iel
!       Em seguida, a rotina do elemento le suas constantes fisicas.
!
      do 100 i = 1, numat
        read(nin,*) ma,iel
        if (ma .lt. 1 .or. ma .gt. numat) then
          print*, '*** ERRO rdata: indice de material invalido:', ma
          stop
        endif
        call elmlib(e(1,ma),dum,dum,dum,dum,1,iel,1,1,1,nin)
        ie(ma) = iel
  100 continue
!
!     Leitura das coordenadas nodais:
!       Linha: no  x  y  (z e lido mas ignorado para 2D)
!
      do 200 i = 1, nnode
        read(nin,*) k, (x(j,k),j=1,ndm)
  200 continue
!
!     Leitura da conectividade dos elementos:
!       Linha: nel  n1 n2 n3 [n4]  ma
!       Importante: nos em ordem anti-horaria!
!
      do 300 i = 1, numel
        read(nin,*) k, (ix(j,k),j=1,nen+1)
  300 continue
!
!     Leitura das condicoes de contorno (deslocamentos prescritos):
!       Linha: no  id1  id2   (id=1 => grau de liberdade restrito)
!       Termina com: 0  0  0
!
  400 continue
        read(nin,*) k, (iaux(j), j = 1, ndf)
        if (k .le. 0) goto 500
        do 410 j = 1, ndf
          id(j,k) = iaux(j)
  410   continue
      goto 400
  500 continue
!
!     Leitura das forcas nodais aplicadas:
!       Linha: no  fx  fy
!       Termina com: 0  0.0  0.0
!
  510 continue
        read(nin,*) k, (aux(j), j = 1, ndf)
        if (k .le. 0) goto 600
        do 520 j = 1, ndf
          f(j,k) = aux(j)
  520   continue
      goto 510
  600 continue
      return
      end subroutine rdata

! ======================================================================
      subroutine numeq(id,nnode,ndf,neq)
!     Numeracao sequencial das equacoes livres.
!     Nos livres  (id=0) recebem numero de equacao: 1, 2, 3, ...
!     Nos restritos (id=1) recebem id=0 (sem equacao associada).
! ======================================================================
      implicit none
      integer :: nnode,ndf,neq,id(ndf,*)
      integer :: n,i,j
      neq = 0
      do 100 n = 1, nnode
      do 100 i = 1, ndf
        j = id(i,n)
        if (j .eq. 0) then
          neq = neq + 1
          id(i,n) = neq
        else
          id(i,n) = 0
        endif
  100 continue
      return
      end subroutine numeq

! ======================================================================
      subroutine profil_csr(ix,id,row_ptr,marca,numel,nen,ndf,neq,nnz)
!     Fase simbolica CSR (1a passada): determina, para cada linha
!     (equacao) i, quantas colunas j >= i ela acopla (incluindo a
!     propria diagonal), com base na conectividade dos elementos.
!     So o TRIANGULO SUPERIOR e considerado (matriz simetrica).
!
!     Resultado: row_ptr(1:neq+1) no formato CSR padrao, onde
!     row_ptr(i) = posicao inicial da linha i em val/col_ind, e
!     row_ptr(neq+1)-1 = nnz (numero total de nao-nulos guardados).
!
!     marca(1:neq) e um vetor de marcacao TEMPORARIO (zerado e
!     reusado a cada linha) para nao contar a mesma coluna duas
!     vezes quando dois elementos compartilham a mesma aresta/no.
! ======================================================================
      implicit none
      integer :: numel,nen,ndf,neq,nnz,ix(nen+1,*),id(ndf,*)
      integer :: row_ptr(*),marca(*)
      integer :: nel,i,j,k,l,noi,noj,kk,ll,lin,col,cont
      call mzero(row_ptr,1,neq+1)
      call mzero(marca,1,neq)
!
!     Conta, para cada linha "lin", quantas colunas distintas
!     col >= lin aparecem nos elementos que tocam "lin". Usa marca()
!     com um "carimbo" (numero da propria linha) para saber se
!     aquela coluna ja foi contada nesta linha, sem precisar zerar
!     marca() inteiro a cada iteracao (custaria O(neq) por linha).
!
      do 50 lin = 1, neq
        cont = 0
        do 400 nel = 1, numel
        do 400 i = 1, nen
          noi = ix(i,nel)
          if (noi .eq. 0) goto 400
          do 300 k = 1, ndf
            kk = id(k,noi)
            if (kk .ne. lin) goto 300
            do 200 j = 1, nen
              noj = ix(j,nel)
              if (noj .eq. 0) goto 200
              do 100 l = 1, ndf
                ll = id(l,noj)
                if (ll .eq. 0) goto 100
                col = max0(kk,ll)
                if (col .lt. lin) goto 100
                if (marca(col) .eq. lin) goto 100
                marca(col) = lin
                cont = cont + 1
  100         continue
  200       continue
  300     continue
  400   continue
        row_ptr(lin) = cont
   50 continue
!
!     Soma cumulativa: converte contagem por linha em ponteiros
!     absolutos no vetor val/col_ind (convencao CSR: row_ptr(1)=1).
!
      row_ptr(neq+1) = 1
      do 500 i = 1, neq
        cont = row_ptr(i)
        row_ptr(i) = row_ptr(neq+1)
        row_ptr(neq+1) = row_ptr(neq+1) + cont
  500 continue
      nnz = row_ptr(neq+1) - 1
      return
      end subroutine profil_csr

! ======================================================================
      subroutine fillcol_csr(ix,id,row_ptr,col_ind,marca,numel,nen, &
                              ndf,neq,nnz)
!     Fase simbolica CSR (2a passada): preenche col_ind com as
!     colunas de cada linha, usando row_ptr ja calculado por
!     profil_csr. Ao final, ordena cada linha por coluna crescente
!     (necessario para a busca binaria usada na montagem).
!
!     marca(1:neq) e o MESMO vetor de marcacao temporario usado em
!     profil_csr (passado explicitamente como argumento, em vez de
!     array automatico local, por seguranca/clareza).
! ======================================================================
      implicit none
      integer :: numel,nen,ndf,neq,nnz,ix(nen+1,*),id(ndf,*)
      integer :: row_ptr(*),col_ind(*),marca(*)
      integer :: nel,i,j,k,l,noi,noj,kk,ll,lin,col,pos,ult
      call mzero(marca,1,neq)
!
!     Preenchimento: para cada linha, percorre os elementos de novo
!     e escreve as colunas distintas a partir de row_ptr(lin).
!
      do 50 lin = 1, neq
        pos = row_ptr(lin)
        do 400 nel = 1, numel
        do 400 i = 1, nen
          noi = ix(i,nel)
          if (noi .eq. 0) goto 400
          do 300 k = 1, ndf
            kk = id(k,noi)
            if (kk .ne. lin) goto 300
            do 200 j = 1, nen
              noj = ix(j,nel)
              if (noj .eq. 0) goto 200
              do 100 l = 1, ndf
                ll = id(l,noj)
                if (ll .eq. 0) goto 100
                col = max0(kk,ll)
                if (col .lt. lin) goto 100
                if (marca(col) .eq. lin) goto 100
                marca(col) = lin
                col_ind(pos) = col
                pos = pos + 1
  100         continue
  200       continue
  300     continue
  400   continue
   50 continue
!
!     Ordenacao de cada linha por coluna crescente (insertion sort,
!     eficiente pois cada linha costuma ter poucos elementos):
!
      do 700 lin = 1, neq
        ult = row_ptr(lin+1) - 1
        do 600 pos = row_ptr(lin)+1, ult
          col = col_ind(pos)
          k = pos - 1
  550     continue
          if (k .lt. row_ptr(lin)) goto 590
          if (col_ind(k) .le. col) goto 590
          col_ind(k+1) = col_ind(k)
          k = k - 1
          goto 550
  590     continue
          col_ind(k+1) = col
  600   continue
  700 continue
      return
      end subroutine fillcol_csr

! ======================================================================
      subroutine pload(id,f,u,nnode,ndf)
!     Montagem do vetor de forcas nodais globais.
!     Transfere as forcas do array local f(ndf,nnode) para o
!     vetor global u(neq), usando a numeracao de equacoes id.
! ======================================================================
      implicit none
      integer :: nnode, ndf, id(ndf,*)
      real*8 :: f(ndf,*), u(*)
      integer :: i,j,k
      do 100 i = 1, nnode
      do 100 j = 1, ndf
        k = id(j,i)
        if (k .gt. 0) u(k) = f(j,i)
  100 continue
      return
      end subroutine pload

! ======================================================================
      subroutine pform_csr(e,ie,ix,id,x,f,u,row_ptr,col_ind,val, &
                        xl,ul,fl,sl,ld,numel,ndm,ndf,nen,nst,isw, &
                        afl,bfl,neq)
!     Loop nos elementos: monta a matriz de rigidez global (CSR,
!     triangulo superior) e corrige o vetor de forcas para
!     deslocamentos prescritos.
!       isw = codigo de instrucao para a rotina de elemento
!       afl = .true. => corrige o vetor de forcas
!       bfl = .true. => monta a matriz de rigidez global
!     Identico ao pform original, exceto a chamada final que usa
!     addstf_csr (montagem via busca binaria) em vez de addstf.
! ======================================================================
      implicit none
      integer :: numel,ndm,ndf,nen,nst,isw,neq
      integer :: ie(*),ix(nen+1,*),id(ndf,*),ld(ndf,*)
      integer :: row_ptr(*),col_ind(*)
      real*8  :: e(10,*),x(ndm,*),f(ndf,*),u(*),val(*),xl(ndm,*),ul(ndf,*)
      real*8  :: fl(*),sl(*)
      logical :: afl, bfl
      integer :: nel,i,j,no,ma,iel,k
!
!     Loop nos elementos:
!
      do 700 nel = 1, numel
!
!       Loop nos nos do elemento: monta vetores locais xl, ul, ld
!
        do 600 i = 1, nen
          no = ix(i,nel)
          if (no .gt. 0) goto 300
!
!         No = 0: zera vetores locais
!
          do 100 j = 1, ndm
            xl(j,i) = 0.d0
  100     continue
          do 200 j = 1, ndf
            ul(j,i) = 0.d0
            ld(j,i) = 0
  200     continue
          goto 600
!
!         No > 0: copia coordenadas e deslocamentos/forcas prescritas
!
  300     continue
          do 400 j = 1, ndm
            xl(j,i) = x(j,no)
  400     continue
          do 500 j = 1, ndf
!           Numeracao global da equacao do grau de liberdade (j,no):
            k = id(j,no)
            ld(j,i) = k
!           Deslocamento prescrito (k=0 indica no restrito):
            ul(j,i) = 0.d0
            if (k .eq. 0) ul(j,i) = f(j,no)
  500     continue
  600   continue
        ma  = ix(nen+1,nel)
        iel = ie(ma)
!
!       Calcula rigidez e forcas locais do elemento:
!
        call elmlib(e(1,ma),xl,ul,fl,sl,nel,iel,ndm,nst,isw,1)
!
!       Acumula contribuicoes local -> global (formato CSR):
!
        call addstf_csr(row_ptr,col_ind,val,u,fl,sl,ld,nst,afl,bfl,neq)
  700 continue
      return
      end subroutine pform_csr

! ======================================================================
      subroutine addstf_csr(row_ptr,col_ind,val,b,p,s,ld,nst,afl,bfl, &
                             neq)
!     Montagem local->global da matriz de rigidez e vetor de forcas.
!     Armazenamento CSR: apenas o TRIANGULO SUPERIOR e armazenado.
!     Para cada par (linha,coluna) com linha<=coluna, busca a
!     posicao correspondente em col_ind(row_ptr(linha):row_ptr(linha+1)-1)
!     usando busca binaria (cada linha esta ordenada, ver fillcol_csr).
!       p = vetor de forcas locais do elemento
!       s = matriz de rigidez local do elemento
!       ld = numeracao global das equacoes do elemento
! ======================================================================
      implicit none
      integer :: row_ptr(*),col_ind(*),ld(*),nst,neq
      real*8 :: val(*),b(*),p(*),s(nst,*)
      logical :: afl,bfl
      integer :: j,k,i,m,pos
      do 200 j = 1, nst
        k = ld(j)
        if (k .eq. 0) goto 200
!
!       Correcao do vetor de forcas para deslocamentos prescritos:
!       b(k) -= s(i,j) * u_prescrito  (contribuicao de cada gdl)
!
        if (afl) b(k) = b(k) - p(j)
        if (.not. bfl) goto 200
!
!       Montagem na matriz de rigidez (so guarda quando m<=k, que e
!       exatamente a mesma condicao de filtro do addstf original;
!       a diferenca e que agora (m,k) com m<=k JA E o par
!       (linha,coluna) do triangulo superior, sem espelhamento):
!
        do 100 i = 1, nst
          m = ld(i)
          if (m .gt. k .or. m .eq. 0) goto 100
          call busca_csr(row_ptr,col_ind,m,k,neq,pos)
          val(pos) = val(pos) + s(i,j)
  100   continue
  200 continue
      return
      end subroutine addstf_csr

! ======================================================================
      subroutine busca_csr(row_ptr,col_ind,lin,col,neq,pos)
!     Busca binaria da coluna "col" dentro da linha "lin" do CSR.
!     Retorna em "pos" o indice em val()/col_ind() correspondente
!     ao par (lin,col). Para se nao encontrar (erro de logica).
! ======================================================================
      implicit none
      integer :: row_ptr(*),col_ind(*),lin,col,neq,pos
      integer :: lo,hi,mid
      lo = row_ptr(lin)
      hi = row_ptr(lin+1) - 1
   10 continue
      if (lo .gt. hi) goto 900
      mid = (lo+hi)/2
      if (col_ind(mid) .eq. col) then
        pos = mid
        return
      else if (col_ind(mid) .lt. col) then
        lo = mid + 1
      else
        hi = mid - 1
      endif
      goto 10
  900 continue
      print*, '*** ERRO busca_csr: posicao (',lin,',',col,')'
      print*, '    nao encontrada na estrutura CSR.'
      print*, '    Isto indica erro na fase simbolica (profil_csr).'
      stop
      end subroutine busca_csr

! ======================================================================
      subroutine elmlib(e,xl,ul,fl,sl,nel,iel,ndm,nst,isw,nin)
!     Biblioteca de elementos: despachante para cada tipo de elemento.
!       iel = 1 => T3 Estado Plano de Deformacoes (EPD)
!       iel = 2 => T3 Estado Plano de Tensoes    (EPT)
!       iel = 3 => Q4 Estado Plano de Deformacoes (EPD)
!       iel = 4 => Q4 Estado Plano de Tensoes    (EPT)
!       iel = 5 => T3 Conducao de calor 2D
!       iel = 6 => T3 Axissimetrico
! ======================================================================
      implicit none
      integer :: nel,iel,ndm,nst,isw,nin
      real*8 :: e(*),xl(*),ul(*),fl(*),sl(*)
      goto (100,200,300,400,500,600) iel
      print*, '*** ERRO elmlib: tipo de elemento ',iel, &
              ' nao existe. Elemento: ',nel
      stop
  100 call elmt01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  200 call elmt02(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  300 call elmt03(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  400 call elmt04(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  500 call elmt01_t(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  600 call elmt01_axi(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
      end subroutine elmlib

! ======================================================================
      subroutine elmt01(e,x,u,p,s,nel,ndm,nst,isw,nin)
!     Elemento triangular linear T3
!     Hipotese: Estado Plano de DEFORMACOES (EPD)
!     Parametros do material: E (modulo de Young), nu (Poisson)
!
!     Matriz constitutiva EPD:
!       c = E(1-nu)/[(1+nu)(1-2nu)]
!       D = | c          c*nu/(1-nu)  0            |
!           | c*nu/(1-nu)  c          0            |
!           | 0            0    E/[2(1+nu)]        |
! ======================================================================
      implicit none
      integer :: nel,ndm,nst,isw,nin
      real*8  :: e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  :: det,xj11,xj12,xj21,xj22,hx(3),hy(3)
      real*8  :: xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  :: my,nu,a,b,c,wt
      integer :: i,j,k,l
      goto (100,200) isw
!
!     isw=1: leitura das constantes fisicas do material
!
  100 continue
      read(nin,*) e(1),e(2)
      return
!
!     isw=2: calculo da matriz de rigidez
!
!     Matriz Jacobiana:
!       J = | x1-x3  y1-y3 |    det(J) = 2 * Area do triangulo
!           | x2-x3  y2-y3 |
!
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
!
!     Inversa da matriz Jacobiana:
!
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
!
!     Derivadas das funcoes de interpolacao em relacao a x e y:
!       hx(i) = dN_i/dx,  hy(i) = dN_i/dy
!
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
!
!     Matriz constitutiva D (Estado Plano de Deformacoes):
!
      my  = e(1)
      nu  = e(2)
      a   = 1.d0+nu
      b   = a*(1.d0-2.d0*nu)
      c   = my*(1.d0-nu)/b
      d11 = c
      d12 = my*nu/b
      d21 = d12
      d22 = c
      d33 = my/(2.d0*a)
!
!     Matriz de rigidez: Ke = B^T D B * det/2
!       wt = Area do triangulo = det/2
!
      wt = 0.5d0*det
      do 220 j = 1, 3
        k = (j-1)*2+1
        do 210 i = 1, 3
          l = (i-1)*2+1
!         s(l,k)     += (dNi/dx*d11*dNj/dx + dNi/dy*d33*dNj/dy)*A
          s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
!         s(l,k+1)   += (dNi/dx*d12*dNj/dy + dNi/dy*d33*dNj/dx)*A
          s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
!         s(l+1,k)   += (dNi/dy*d21*dNj/dx + dNi/dx*d33*dNj/dy)*A
          s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
!         s(l+1,k+1) += (dNi/dy*d22*dNj/dy + dNi/dx*d33*dNj/dx)*A
          s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
  210   continue
  220 continue
!
!     Vetor de forcas internas: p = Ke * u
!
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt01: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end subroutine elmt01

! ======================================================================
      subroutine elmt02(e,x,u,p,s,nel,ndm,nst,isw,nin)
!     Elemento triangular linear T3
!     Hipotese: Estado Plano de TENSOES (EPT)
!     Parametros do material: E, nu, espessura t
!
!     Matriz constitutiva EPT:
!       a = E/(1-nu^2)
!       D = | a      a*nu      0          |
!           | a*nu   a         0          |
!           | 0      0    E/[2(1+nu)]     |
! ======================================================================
      implicit none
      integer :: nel,ndm,nst,isw,nin
      real*8  :: e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  :: det,xj11,xj12,xj21,xj22,hx(3),hy(3)
      real*8  :: xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  :: my,nu,thic,a,wt
      integer :: i,j,k,l
      goto (100,200) isw
!
!     isw=1: leitura das constantes fisicas: E, nu, espessura
!
  100 continue
      read(nin,*) e(1),e(2),e(3)
      return
!
!     isw=2: calculo da matriz de rigidez
!
!     Matriz Jacobiana:
!
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
!
!     Inversa da Jacobiana:
!
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
!
!     Derivadas das funcoes de interpolacao:
!
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
!
!     Matriz constitutiva D (Estado Plano de Tensoes):
!
      my   = e(1)
      nu   = e(2)
      thic = e(3)
      a    = my/(1.d0-nu*nu)
      d11  = a
      d12  = a*nu
      d21  = d12
      d22  = a
      d33  = my/(2.d0*(1.d0+nu))
!
!     Matriz de rigidez: Ke = B^T D B * (det/2) * t
!       wt = Area * espessura
!
      wt = 0.5d0*det*thic
      do 220 j = 1, 3
        k = (j-1)*2+1
        do 210 i = 1, 3
          l = (i-1)*2+1
          s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
          s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
          s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
          s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
  210   continue
  220 continue
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt02: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end subroutine elmt02

! ======================================================================
      subroutine elmt03(e,x,u,p,s,nel,ndm,nst,isw,nin)
!     Elemento quadrilatero bilinear Q4
!     Hipotese: Estado Plano de DEFORMACOES (EPD)
!     Parametros do material: E, nu
!
!     Integracao numerica: Gauss 2x2 (4 pontos)
!     Pontos de Gauss: (r,s) = (+/-1/sqrt(3), +/-1/sqrt(3))
!     Pesos: wi = wj = 1.0 (todos iguais, ja incorporados)
!
!     Funcoes de interpolacao bilineares (coordenadas isoparametricas):
!       N1=(1+r)(1+s)/4, N2=(1-r)(1+s)/4
!       N3=(1-r)(1-s)/4, N4=(1+r)(1-s)/4
! ======================================================================
      implicit none
      integer :: nel,ndm,nst,isw,i,j,k,l,m,n,nin
      real*8  :: e(*),x(ndm,*),u(nst),p(nst),s(nst,nst),xj(ndm,2)
      real*8  :: det,hx(4),hy(4),xji(ndm,2),hr(4),hs(4)
      real*8  :: rr,ss,d11,d12,d21,d22,d33
      real*8  :: my,nu,a,b,c,wt
      goto (100,200) isw
!
!     isw=1: leitura das constantes fisicas: E, nu
!
  100 continue
      read(nin,*) e(1),e(2)
      return
!
!     isw=2: calculo da matriz de rigidez
!
!     Matriz constitutiva D (Estado Plano de Deformacoes):
!
  200 continue
      my  = e(1)
      nu  = e(2)
      a   = 1.d0+nu
      b   = a*(1.d0-2.d0*nu)
      c   = my*(1.d0-nu)/b
      d11 = c
      d12 = my*nu/b
      d21 = d12
      d22 = c
      d33 = my/(2.d0*a)
!
!     Zeragem da matriz de rigidez local:
!
      do 20 i = 1, nst
      do 20 j = 1, nst
        s(i,j) = 0.d0
   20 continue
!
!     Loop nos 4 pontos de Gauss (2x2):
!     Coordenadas: +/-1/sqrt(3) = +/-0.577350269189626
!     Pesos: todos 1.0
!
      do 600 i = 1, 2
      do 600 j = 1, 2
        if (i .eq. 1 .and. j .eq. 1) then
          rr =  0.577350269189626d0
          ss = -0.577350269189626d0
        else if (i .eq. 1 .and. j .eq. 2) then
          rr =  0.577350269189626d0
          ss =  0.577350269189626d0
        else if (i .eq. 2 .and. j .eq. 1) then
          rr = -0.577350269189626d0
          ss =  0.577350269189626d0
        else
          rr = -0.577350269189626d0
          ss = -0.577350269189626d0
        endif
!
!       Derivadas das funcoes de interpolacao em r e s:
!         hr(i) = dNi/dr,  hs(i) = dNi/ds
!
        hr(1) =   (1.d0+ss) / 4.d0
        hr(2) = - (1.d0+ss) / 4.d0
        hr(3) = - (1.d0-ss) / 4.d0
        hr(4) =   (1.d0-ss) / 4.d0
        hs(1) =   (1.d0+rr) / 4.d0
        hs(2) =   (1.d0-rr) / 4.d0
        hs(3) = - (1.d0-rr) / 4.d0
        hs(4) = - (1.d0+rr) / 4.d0
!
!       Matriz Jacobiana: J = dX/d(r,s)
!         J(1,1)=sum(hr*x), J(2,1)=sum(hs*x)
!         J(1,2)=sum(hr*y), J(2,2)=sum(hs*y)
!
        xj(1,1) = hr(1)*x(1,1)+hr(2)*x(1,2)+hr(3)*x(1,3)+hr(4)*x(1,4)
        xj(2,1) = hs(1)*x(1,1)+hs(2)*x(1,2)+hs(3)*x(1,3)+hs(4)*x(1,4)
        xj(1,2) = hr(1)*x(2,1)+hr(2)*x(2,2)+hr(3)*x(2,3)+hr(4)*x(2,4)
        xj(2,2) = hs(1)*x(2,1)+hs(2)*x(2,2)+hs(3)*x(2,3)+hs(4)*x(2,4)
!
!       Determinante da Jacobiana: det(J) = J11*J22 - J21*J12
!       det > 0 garante elemento valido (numeracao anti-horaria)
!
        det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
        if (det .le. 0.d0) goto 1000
!
!       Inversa da Jacobiana: J^-1 = (1/det)*adj(J)
!
        xji(1,1) =  xj(2,2) / det
        xji(1,2) = -xj(1,2) / det
        xji(2,1) = -xj(2,1) / det
        xji(2,2) =  xj(1,1) / det
!
!       Derivadas das funcoes de interpolacao em x e y:
!         hx(k) = dNk/dx = J^-1(1,1)*hr(k) + J^-1(1,2)*hs(k)
!         hy(k) = dNk/dy = J^-1(2,1)*hr(k) + J^-1(2,2)*hs(k)
!
        do 300 k = 1, 4
          hx(k) = xji(1,1)*hr(k) + xji(1,2)*hs(k)
          hy(k) = xji(2,1)*hr(k) + xji(2,2)*hs(k)
  300   continue
!
!       Contribuicao ao Ke: Ke += B^T*D*B * det(J) * wi*wj  (wi=wj=1)
!
        wt = det
        do 500 m = 1, 4
          k = (m-1)*2+1
          do 400 n = 1, 4
            l = (n-1)*2+1
            s(l,k)     = s(l,k)    +(hx(n)*d11*hx(m)+hy(n)*d33*hy(m))*wt
            s(l,k+1)   = s(l,k+1)  +(hx(n)*d12*hy(m)+hy(n)*d33*hx(m))*wt
            s(l+1,k)   = s(l+1,k)  +(hy(n)*d12*hx(m)+hx(n)*d33*hy(m))*wt
            s(l+1,k+1) = s(l+1,k+1)+(hy(n)*d22*hy(m)+hx(n)*d33*hx(m))*wt
  400     continue
  500   continue
  600 continue
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt03: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end subroutine elmt03

! ======================================================================
      subroutine elmt04(e,x,u,p,s,nel,ndm,nst,isw,nin)
!     Elemento quadrilatero bilinear Q4
!     Hipotese: Estado Plano de TENSOES (EPT)
!     Parametros do material: E, nu, espessura t
!
!     Integracao numerica: Gauss 2x2 (4 pontos)
!     Igual ao elmt03, com matriz D de EPT e fator espessura.
! ======================================================================
      implicit none
      integer :: nel,ndm,nst,isw,i,j,k,l,m,n,nin
      real*8  :: e(*),x(ndm,*),u(nst),p(nst),s(nst,nst),xj(ndm,2)
      real*8  :: det,hx(4),hy(4),xji(ndm,2),hr(4),hs(4)
      real*8  :: rr,ss,d11,d12,d21,d22,d33
      real*8  :: my,nu,a,thic,wt
      goto (100,200) isw
!
!     isw=1: leitura das constantes fisicas: E, nu, espessura
!
  100 continue
      read(nin,*) e(1),e(2),e(3)
      return
!
!     isw=2: calculo da matriz de rigidez
!
!     Matriz constitutiva D (Estado Plano de Tensoes):
!
  200 continue
      my   = e(1)
      nu   = e(2)
      thic = e(3)
      a    = my/(1.d0-nu*nu)
      d11  = a
      d12  = a*nu
      d21  = d12
      d22  = a
      d33  = my/(2.d0*(1.d0+nu))
!
!     Zeragem da matriz de rigidez local:
!
      do 20 i = 1, nst
      do 20 j = 1, nst
        s(i,j) = 0.d0
   20 continue
!
!     Loop nos 4 pontos de Gauss (2x2):
!
      do 600 i = 1, 2
      do 600 j = 1, 2
        if (i .eq. 1 .and. j .eq. 1) then
          rr =  0.577350269189626d0
          ss = -0.577350269189626d0
        else if (i .eq. 1 .and. j .eq. 2) then
          rr =  0.577350269189626d0
          ss =  0.577350269189626d0
        else if (i .eq. 2 .and. j .eq. 1) then
          rr = -0.577350269189626d0
          ss =  0.577350269189626d0
        else
          rr = -0.577350269189626d0
          ss = -0.577350269189626d0
        endif
!
!       Derivadas das funcoes de interpolacao em r e s:
!
        hr(1) =   (1.d0+ss) / 4.d0
        hr(2) = - (1.d0+ss) / 4.d0
        hr(3) = - (1.d0-ss) / 4.d0
        hr(4) =   (1.d0-ss) / 4.d0
        hs(1) =   (1.d0+rr) / 4.d0
        hs(2) =   (1.d0-rr) / 4.d0
        hs(3) = - (1.d0-rr) / 4.d0
        hs(4) = - (1.d0+rr) / 4.d0
!
!       Matriz Jacobiana:
!
        xj(1,1) = hr(1)*x(1,1)+hr(2)*x(1,2)+hr(3)*x(1,3)+hr(4)*x(1,4)
        xj(2,1) = hs(1)*x(1,1)+hs(2)*x(1,2)+hs(3)*x(1,3)+hs(4)*x(1,4)
        xj(1,2) = hr(1)*x(2,1)+hr(2)*x(2,2)+hr(3)*x(2,3)+hr(4)*x(2,4)
        xj(2,2) = hs(1)*x(2,1)+hs(2)*x(2,2)+hs(3)*x(2,3)+hs(4)*x(2,4)
!
!       Determinante da Jacobiana:
!
        det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
        if (det .le. 0.d0) goto 1000
!
!       Inversa da Jacobiana:
!
        xji(1,1) =  xj(2,2) / det
        xji(1,2) = -xj(1,2) / det
        xji(2,1) = -xj(2,1) / det
        xji(2,2) =  xj(1,1) / det
!
!       Derivadas das funcoes de interpolacao em x e y:
!
        do 300 k = 1, 4
          hx(k) = xji(1,1)*hr(k) + xji(1,2)*hs(k)
          hy(k) = xji(2,1)*hr(k) + xji(2,2)*hs(k)
  300   continue
!
!       Contribuicao ao Ke: Ke += B^T*D*B * det(J) * t * wi*wj
!
        wt = det * thic
        do 500 m = 1, 4
          k = (m-1)*2+1
          do 400 n = 1, 4
            l = (n-1)*2+1
            s(l,k)     = s(l,k)    +(hx(n)*d11*hx(m)+hy(n)*d33*hy(m))*wt
            s(l,k+1)   = s(l,k+1)  +(hx(n)*d12*hy(m)+hy(n)*d33*hx(m))*wt
            s(l+1,k)   = s(l+1,k)  +(hy(n)*d12*hx(m)+hx(n)*d33*hy(m))*wt
            s(l+1,k+1) = s(l+1,k+1)+(hy(n)*d22*hy(m)+hx(n)*d33*hx(m))*wt
  400     continue
  500   continue
  600 continue
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt04: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end subroutine elmt04

! ======================================================================
      subroutine elmt01_t(e,x,u,p,s,nel,ndm,nst,isw,nin)
!     Elemento triangular linear T3
!     Hipotese: Conducao de calor 2D (equacao do calor)
!     Parametros do material: k ( condutividade termica)
!
! ======================================================================
      implicit none
      integer :: nel,ndm,nst,isw,nin
      real*8  :: e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  :: det,xj11,xj12,xj21,xj22,hx(3),hy(3)
      real*8  :: xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  :: difus,wt
      integer :: i,j
      goto (100,200) isw
!
!     isw=1: leitura das constantes fisicas do material
!
  100 continue
      read(nin,*) e(3)
      return
!
!     isw=2: calculo da matriz de rigidez
!
!     Matriz Jacobiana:
!       J = | x1-x3  y1-y3 |    det(J) = 2 * Area do triangulo
!           | x2-x3  y2-y3 |
!
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
!
!     Inversa da matriz Jacobiana:
!
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
!
!     Derivadas das funcoes de interpolacao em relacao a x e y:
!       hx(i) = dN_i/dx,  hy(i) = dN_i/dy
!
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
!
      difus  = e(3)
!
!       wt = Area do triangulo = det/2
!
      wt = 0.5d0*det
      do 220 j = 1, 3
        do 210 i = 1, 3
!         s(i,j)     += difus*(dNi/dx*dNj/dx + dNi/dy*dNj/dy)*A
          s(i,j)     = ( hx(i)*difus*hx(j) + hy(i)*difus*hy(j) ) * wt
  210   continue
  220 continue
!
!     Vetor de forcas internas: p = Ke * u
!
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt01_t: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end subroutine elmt01_t

! ======================================================================
      subroutine elmt01_axi(e,x,u,p,s,nel,ndm,nst,isw,nin)
!     Elemento triangular linear T3
!     Hipotese: Axissimetrico
!     Parametros do material: E (modulo de Young), nu (Poisson)
!
!     Matriz constitutiva EPD:
!       c = E(1-nu)/[(1+nu)(1-2nu)]
!       D = | c          c*nu/(1-nu)  0           c*nu/(1-nu) |     |
!           | c*nu/(1-nu)  c          0           c*nu/(1-nu) |
!           | 0            0    E/[2(1+nu)]             0     |
!           | c*nu/(1-nu)  c*nu/(1-nu)  0               c     |
! ======================================================================
      implicit none
      integer :: nel,ndm,nst,isw,nin
      real*8  :: e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  :: det,xj11,xj12,xj21,xj22,hx(4),hy(4)
      real*8  :: xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  :: my,nu,a,b,c,wt,raio
      integer :: i,j,k,l
      real*8  :: PI
      parameter (PI = 3.1415926535)
      goto (100,200) isw
!
!     isw=1: leitura das constantes fisicas do material
!
  100 continue
      read(nin,*) e(1),e(2)
      return
!
!     isw=2: calculo da matriz de rigidez
!
!     Matriz Jacobiana:
!       J = | x1-x3  y1-y3 |    det(J) = 2 * Area do triangulo
!           | x2-x3  y2-y3 |
!
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      raio = 1/3.d0*(x(1,1)+x(1,2)+x(1,3))
      if (det .le. 0.d0) goto 1000
!
!     Inversa da matriz Jacobiana:
!
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
!
!     Derivadas das funcoes de interpolacao em relacao a r e z:
!       hx(i) = dN_i/dr,  hy(i) = dN_i/dz
!
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
!
!     Matriz constitutiva D (problema axissimetrico):
!
      my  = e(1)
      nu  = e(2)
      a   = 1.d0+nu
      b   = a*(1.d0-2.d0*nu)
      c   = my*(1.d0-nu)/b
      d11 = c
      d12 = my*nu/b
      d21 = d12
      d22 = c
      d33 = my/(2.d0*a)
!
!     Matriz de rigidez: Ke = 2 * pi * B^T D B * r * det/2
!       wt = Area do triangulo = det/2
!
      wt = 2.d0*PI*0.5d0*raio*det
      do 220 j = 1, 3
        k = (j-1)*2+1
        do 210 i = 1, 3
          l = (i-1)*2+1
!         s(l,k)     += (dNi/dx*d11*dNj/dx + dNi/dy*d33*dNj/dy)*A
          s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
!         s(l,k+1)   += (dNi/dx*d12*dNj/dy + dNi/dy*d33*dNj/dx)*A
          s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
!         s(l+1,k)   += (dNi/dy*d21*dNj/dx + dNi/dx*d33*dNj/dy)*A
          s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
!         s(l+1,k+1) += (dNi/dy*d22*dNj/dy + dNi/dx*d33*dNj/dx)*A
          s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
  210   continue
  220 continue
!
!     Vetor de forcas internas: p = Ke * u
!
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt01: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end subroutine elmt01_axi

! ======================================================================
      subroutine actcol(a,b,jdiag,neq,afac,back)
! ======================================================================
!     Solver direto: decomposicao LtDL com armazenamento skyline.
!     Valido somente para matrizes simetricas.
!
!     Parametros de entrada:
!       a     = coeficientes da matriz (armazenados por altura de coluna)
!       b     = vetor independente
!       jdiag = vetor apontador do armazenamento skyline
!       neq   = numero de equacoes
!       afac  = .true. => fatoriza a matriz
!       back  = .true. => retrosubstitui
!
!     Parametros de saida:
!       a     = coeficientes da matriz triangularizada
!       b     = vetor solucao
! ======================================================================
      implicit real*8 (a-h,o-z)
      common/engys/ aengy
      dimension a(*),b(*),jdiag(*)
      logical afac,back
      aengy = 0.0d0
      jr = 0
      do 600 j = 1,neq
        jd = jdiag(j)
        jh = jd - jr
        is = j - jh + 2
        if(jh-2) 600,300,100
  100   if(.not.afac) goto 500
        ie = j - 1
        k  = jr + 2
        id = jdiag(is - 1)
        do 200 i = is, ie
          ir = id
          id = jdiag(i)
          ih = min0(id-ir-1,i-is+1)
          if(ih.gt.0) a(k) = a(k) - dot(a(k-ih),a(id-ih),ih)
  200   k = k + 1
  300   if(.not.afac) goto 500
        ir = jr+1
        ie = jd - 1
        k  = j - jd
        do 400 i = ir, ie
          id = jdiag(k+i)
          if(a(id).eq.0.0d0) goto 400
          d    = a(i)
          a(i) = a(i)/a(id)
          a(jd) = a(jd) - d*a(i)
  400   continue
  500   if(back) b(j) = b(j) - dot(a(jr+1),b(is-1),jh-1)
  600 jr = jd
      if(.not.back) return
      do 700 i = 1,neq
        id = jdiag(i)
        if(a(id).ne.0.0d0) b(i) = b(i)/a(id)
  700 aengy = aengy + b(i)*b(i)*a(id)
      j  = neq
      jd = jdiag(j)
  800 d  = b(j)
      j  = j - 1
      if(j.le.0) return
      jr = jdiag(j)
      if(jd-jr.le.1) goto 1000
      is = j - jd + jr + 2
      k  = jr - is + 1
      do 900 i = is,j
  900 b(i) = b(i) - a(i+k)*d
 1000 jd = jr
      goto 800
      end subroutine actcol

! ======================================================================
      subroutine lku(s,u,p,nst)
!     Produto matriz-vetor local: p = s * u
!     Calcula o vetor de forcas internas do elemento.
! ======================================================================
      implicit none
      integer :: nst
      real*8 :: s(nst,nst),u(nst),p(nst)
      integer :: i,j
      do 200 i = 1, nst
        p(i) = 0.d0
        do 100 j = 1, nst
          p(i) = p(i) + s(i,j) * u(j)
  100   continue
  200 continue
      return
      end subroutine lku

! ======================================================================
      real*8 function dot(a,b,n)
!     Produto interno (escalar) de dois vetores: dot = a . b
! ======================================================================
      implicit none
      integer :: n
      real*8 :: a(*), b(*)
      integer :: i
      dot = 0.d0
      do 100 i = 1, n
        dot = dot + a(i) * b(i)
  100 continue
      return
      end function dot

! ======================================================================
      subroutine mem(npos)
!     Verifica se ha memoria suficiente no vetor global.
!     Para se a posicao requerida ultrapassar o limite maximo.
! ======================================================================
      implicit none
      integer :: max
      common /size/ max
      integer :: npos
      if ( (npos-1) .gt. max ) then
        print*, '*** ERRO: Memoria insuficiente!'
        print*, '    Posicao requerida:', npos-1, &
                '  Maximo disponivel:', max
        print*, '    Aumente o parametro npos no programa principal.'
        stop
      endif
      return
      end subroutine mem

! ======================================================================
      subroutine mzero(m,i1,i2)
!     Zeragem de vetor inteiro: m(i1:i2) = 0
! ======================================================================
      implicit none
      integer :: i1,i2,m(*)
      integer :: i
      do 100 i = 1, i2-i1+1
        m(i) = 0
  100 continue
      return
      end subroutine mzero

! ======================================================================
      subroutine azero(a,i1,i2)
!     Zeragem de vetor real*8: a(i1:i2) = 0.d0
! ======================================================================
      implicit none
      integer :: i1,i2
      real*8 :: a(*)
      integer :: i
      do 100 i = 1, i2-i1+1
        a(i) = 0.d0
  100 continue
      return
      end subroutine azero

! ======================================================================
      subroutine wdata(ix,id,x,f,u,nnode,numel,ndm,nen,ndf,nout)
!     Escrita dos resultados no arquivo de saida.
!     Formato compativel com visualizadores de pos-processamento.
!
!     Blocos gravados:
!       coor  => coordenadas nodais
!       elem  => conectividade dos elementos
!       nosc  => deslocamentos nodais
!       nvec  => campo vetorial com z=0 (para visualizacao 3D)
!       end   => marcador de fim de arquivo
! ======================================================================
      implicit none
      integer :: nnode,numel,ndm,nen,ndf,ix(nen+1,*),id(ndf,*),nout
      real*8  :: x(ndm,*),f(ndf,*),u(*),aux(6)
      integer :: i,j,k
      write(nout,'(a,i5)') 'coor ', nnode
      do 100 i = 1, nnode
        write(nout,'(i10,3e15.5)') i,(x(j,i),j=1,2),0.0
  100 continue
      write(nout,'(a,i10)') 'elem ', numel
      do 200 i = 1, numel
        write(nout,'(10i10)') i,nen,(ix(j,i),j=1,nen)
  200 continue
!
!     Deslocamentos nodais:
!       Nos livres  (id>0): deslocamento calculado u(id)
!       Nos restritos (id=0): deslocamento prescrito f(j,i)
!
      write(nout,'(a,i2)') 'nosc ',ndf
      do 300 i = 1, nnode
        do 310 j = 1, ndf
          aux(j) = f(j,i)
          k = id(j,i)
          if(k .gt. 0) aux(j) = u(k)
  310   continue
        write(nout,'(i10,6e15.5e3)') i,(aux(j),j=1,ndf)
  300 continue
!
!     Campo vetorial (adiciona z=0 para pos-processamento 3D):
!
      write(nout,'(a)') 'nvec '
      do 400 i = 1, nnode
        do 410 j = 1, ndf
          aux(j) = f(j,i)
          k = id(j,i)
          if(k .gt. 0) aux(j) = u(k)
  410   continue
        write(nout,'(i10,6e15.5e3)') i,(aux(j),j=1,ndf),0.0d0
  400 continue
      write(nout,'(a)') 'end '
      return
      end subroutine wdata

! ======================================================================
      subroutine pcg_csr(val,row_ptr,col_ind,b,x,r,z,prec,d,neq)
!     Solver iterativo: Gradiente Conjugado Precondicionado (PCG)
!     Precondicionador: diagonal de K (precondicionador de Jacobi)
!
!     Resolve: K * x = b
!     Identico ao pcg original (skyline), exceto que a matriz e
!     passada em formato CSR (val,row_ptr,col_ind) e a chamada
!     interna usa matvec_csr em vez de matvec.
!
!     Parametros:
!       val      = valores nao-nulos da matriz (triangulo superior, CSR)
!       row_ptr  = ponteiros de linha (CSR)
!       col_ind  = colunas de cada valor (CSR)
!       b    = vetor de forcas (entrada) / deslocamentos (saida)
!       x    = vetor solucao (trabalho interno)
!       r    = vetor residuo (trabalho)
!       z    = vetor precondicionado (trabalho)
!       prec = precondicionador diagonal M = diag(K)
!       d    = direcao de descida conjugada (trabalho)
!       neq  = numero de equacoes
! ======================================================================
      implicit none
      integer :: neq,maxit,i,iter
      integer :: row_ptr(*),col_ind(*)
      real*8  :: val(*),prec(*),x(*),r(*),z(*),b(*),d(*)
      real*8  :: dnew,dold,alpha,beta,tol,norma_energia,dot
      external dot
      tol   = 1.0d-07
      maxit = 1000
!
!     Precondicionador diagonal: prec(i) = K(i,i) = val(row_ptr(i))
!     (a diagonal e sempre o primeiro elemento de cada linha, pois
!     col_ind esta ordenado e a menor coluna possivel da linha i,
!     no triangulo superior, e a propria diagonal i):
!
      do 10 i = 1, neq
        prec(i) = val(row_ptr(i))
        if (dabs(prec(i)) .lt. 1.d-30) then
          print*, '*** AVISO pcg_csr: diagonal nula na equacao', i
          print*, '    Verifique as condicoes de contorno.'
          prec(i) = 1.d-30
        endif
   10 continue
!
!     Inicializacao: x=0, r=b-K*0=b, d=M^-1*r
!
      do 20 i = 1, neq
        x(i) = 0.d0
   20 continue
      call matvec_csr(neq,val,row_ptr,col_ind,x,z)
      do 30 i = 1, neq
        r(i) = b(i) - z(i)
        d(i) = r(i) / prec(i)
   30 continue
      dnew = dot(r,d,neq)
      tol  = tol*tol*(dabs(dnew))
!
!     Loop principal do PCG:
!
      do 200 iter = 1, maxit
        call matvec_csr(neq,val,row_ptr,col_ind,d,z)
        alpha = dnew / dot(d,z,neq)
        do 100 i = 1, neq
          x(i) = x(i) + alpha * d(i)
          r(i) = r(i) - alpha * z(i)
          z(i) = r(i) / prec(i)
  100   continue
        dold = dnew
        beta = dot(r,z,neq) / dold
        dnew = dot(r,z,neq)
        do 150 i = 1, neq
          d(i) = z(i) + beta * d(i)
  150   continue
!
!       Criterio de convergencia: ||r||_M < tol
!
        if (dabs(dnew) .lt. tol) then
          print*, 'PCG convergiu em', iter, 'iteracoes.'
          goto 300
        endif
  200 continue
!
!     Nao convergiu:
!
      print*, '*** AVISO pcg_csr: nao convergiu em', maxit, 'iteracoes.'
      print*, '    Residuo final (dnew):', dnew
      print*, '    Tolerancia  (tol)   :', tol
!
!     Norma de energia: ||u||_K = sqrt(u^T K u)
!
  300 continue
      call matvec_csr(neq,val,row_ptr,col_ind,x,z)
      norma_energia = dsqrt(dabs(dot(x,z,neq)))
      write(*,'(a,i6)')    'Numero de equacoes    :', neq
      write(*,'(a,f15.8)') 'Norma de energia      :', norma_energia
!
!     Copia solucao de volta para b (vetor de deslocamentos):
!
      do 400 i = 1, neq
        b(i) = x(i)
  400 continue
      return
      end subroutine pcg_csr

! ======================================================================
      subroutine matvec_csr(neq,val,row_ptr,col_ind,x,y)
!     Produto matriz-vetor para matriz simetrica no formato CSR
!     (triangulo superior). Calcula: y = K * x
!
!     Aproveita a simetria: para cada elemento val(k) armazenado na
!     linha "lin", coluna "col" (col>=lin), contribui para y(lin)
!     (sempre) e, quando col<>lin, tambem para y(col), explorando
!     K(lin,col) = K(col,lin). Mesmo princiio do matvec skyline,
!     adaptado para a estrutura de linhas/colunas explicitas do CSR.
! ======================================================================
      implicit none
      integer :: neq,i,lin,col,k
      integer :: row_ptr(*),col_ind(*)
      real*8  :: val(*),x(*),y(*)
      do 50 i = 1, neq
        y(i) = 0.d0
   50 continue
      do 200 lin = 1, neq
        do 100 k = row_ptr(lin), row_ptr(lin+1)-1
          col = col_ind(k)
          if (col .eq. lin) then
!
!           Elemento diagonal: contribui so uma vez para y(lin)
!
            y(lin) = y(lin) + val(k)*x(lin)
          else
!
!           Elemento fora da diagonal: contribui para y(lin) e,
!           por simetria, tambem para y(col)
!
            y(lin) = y(lin) + val(k)*x(col)
            y(col) = y(col) + val(k)*x(lin)
          endif
  100   continue
  200 continue
      return
      end subroutine matvec_csr

! ======================================================================
      subroutine wvtk(ix,id,x,f,u,nnode,numel,ndm,nen,ndf,nvtk)
!     Escrita dos resultados em formato VTK para visualizacao no ParaView.
!
!     Formato: VTK Unstructured Grid ASCII
!       - Coordenadas nodais (com z=0 para 2D)
!       - Conectividade dos elementos
!       - Dados nodais: deslocamentos (u_x, u_y, u_z)
! ======================================================================
      implicit none
      integer :: nnode,numel,ndm,nen,ndf,ix(nen+1,*),id(ndf,*),nvtk
      real*8  :: x(ndm,*),f(ndf,*),u(*),aux(3)
      integer :: i,j,k
!
!     Cabecalho VTK:
!
      write(nvtk,'(a)') '# vtk DataFile Version 2.0'
      write(nvtk,'(a)') 'MEF - Elasticidade Plana 2D'
      write(nvtk,'(a)') 'ASCII'
      write(nvtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
!
!     Bloco de coordenadas nodais:
!
      write(nvtk,'(a,i8,a)') 'POINTS ', nnode, ' float'
      do 100 i = 1, nnode
        write(nvtk,'(3e16.8)') x(1,i), x(2,i), 0.0d0
  100 continue
!
!     Bloco de conectividade dos elementos:
!     Total de valores = numel*(nen+1)  [numero de nos + tipo]
!     Tipo VTK: 5 para triangulos (T3), 9 para quadrilateros (Q4)
!
      write(nvtk,'(a,2i8)') 'CELLS ', numel, numel*(nen+1)
      do 200 i = 1, numel
        write(nvtk,'(10i8)') nen, (ix(j,i)-1, j=1,nen)
  200 continue
!
!     Tipos de celulas VTK:
!     5 = VTK_TRIANGLE (T3)
!     9 = VTK_QUAD (Q4)
!
      write(nvtk,'(a,i8)') 'CELL_TYPES ', numel
      if (nen .eq. 3) then
        do 300 i = 1, numel
          write(nvtk,'(i2)') 5
  300   continue
      else if (nen .eq. 4) then
        do 350 i = 1, numel
          write(nvtk,'(i2)') 9
  350   continue
      endif
      if (ndf .eq. 2) then
!
!     Dados nodais: deslocamentos (u_x, u_y, u_z com z=0)
!
         write(nvtk,'(a)') 'VECTORS Deslocamentos float'
         write(nvtk,'(a,i8)') 'POINT_DATA ', nnode
         do 400 i = 1, nnode
           aux(1) = f(1,i)
           aux(2) = f(2,i)
           aux(3) = 0.0d0
           k = id(1,i)
           if (k .gt. 0) aux(1) = u(k)
           k = id(2,i)
           if (k .gt. 0) aux(2) = u(k)
           write(nvtk,'(3e16.8)') aux(1), aux(2), aux(3)
  400    continue
         return
      elseif (ndf .eq. 1) then
!
!     Dados nodais: temperatura
!
         write(nvtk,'(a,i8)') 'POINT_DATA ', nnode
         write(nvtk,'(a)') 'SCALARS Temperatura float'
         write(nvtk,'(a)') 'LOOKUP_TABLE default'
         do 500 i = 1, nnode
           aux(1) = f(1,i)
           k = id(1,i)
           if (k .gt. 0) aux(1) = u(k)
           write(nvtk,'(1e16.8)') aux(1)
  500    continue
         return
      endif
      end subroutine wvtk
