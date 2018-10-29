      
c      nnoel=30 !no máximo 30 nos por elemento apesar do programa poder ser geral
      not=nnoel
      nglel=3*nnoel
      nzexp=(nglel**2)*nel*2
      nzexp1=Int((2*nzexp+3*n+1)*1.2)
      
      ALLOCATE (prsi0(nel,nnoel,ng3,3)) !variável de estado (mede o tamanho da superfície de plastificaçã)
      ALLOCATE (ep(nel,nnoel,ng3,3,6)) !deformação plástica acumulada para cada ponto de integração
      allocate (itagnos(nnos)) !identifica se o nó é de extremidade 

      allocate (sigmaabaixo(nnos,3),sigmaacima(nnos,3),momentof(nnos)) !tensoes e momento fletor

      allocate (anx(nnos),any(nnos))
      ALLOCATE (eeq(nel,nnoel,ng3,3)) !deformação equivalente
      ALLOCATE (temc(nel))      !temperatura total acima da barra
      ALLOCATE (temb(nel))      !temperatura total abaixo da barra
      ALLOCATE (dtemc(nel))     !acréscimo por passo de tempo ou de carga
      ALLOCATE (dtemb(nel))     !acréscimo por passo de tempo ou de carga
      ALLOCATE (calt(nel))      !coeficiente de dilatação térmica do material por elemento
      
      prsi0=0.
      ep=0.
      eeq=0.
      temc=0.
      temb=0.
      dtemc=0.
      dtemb=0.
      calt=0.
      allocate(alfaglob(not))   !vetor global com os angulos dos elementos de domínio

      ALLOCATE (rmz(nel,not))   !esforço solicitante momento fletor
      ALLOCATE (rvy(nel,not))   !esforço solicitante força cortante
      ALLOCATE (rnx(nel,not))   !esforço solicitante força normal
      
      ALLOCATE (fd(nel,3))      !posição da fibra da seção transversal (3 fibras por enquanto)
      ALLOCATE (brd(nel,3))     !largura da fibra
      ALLOCATE (yb(nel))        !posição do cg (módulo de engenharia)
      ALLOCATE (rib(nel))       !momento de inércia
      ALLOCATE (hb(nel,3))      !altura da fibra
      
      ALLOCATE (qx(nel,nnoel))  !força distribuída horizontal
      ALLOCATE (qy(nel,nnoel))  !força distribuída vertical   
      ALLOCATE (qn(nel,nnoel))  !força não conservativa ortogonal   
      ALLOCATE (qt(nel,nnoel))  !força não conservativa tangencial
      
      ALLOCATE (ro(nel))        !densidade (massa) do material
      ALLOCATE (ram(nel))       !amortecimento viscoso
      ALLOCATE (rmp(not,not,nel)) !matriz de massa compacta
      ALLOCATE (rmg(nglel,nglel,nel)) !matriz de massa ampliada
      ALLOCATE (sx(nel,nnoel))  !tensao x
      ALLOCATE (sy(nel,nnoel))  !tensao y
      ALLOCATE (sxy(nel,nnoel)) !tensao xy
      ALLOCATE (in(noel,noel)) !incidência dos elementos
      ALLOCATE (rhook(nel,6,6)) !lei constitutiva (elástica) por elemento
      ALLOCATE (indx(n)) ! versão casca (ainda não apagar)
      ALLOCATE (itip(100,4)) ! (versão casca (ainda não apagar)
      ALLOCATE (inct(55,3))     !(versão casca)
      ALLOCATE (ko(n))          !restrição graus de liberdade / aplicação de recalque
      ALLOCATE (ko1(n))         !restrição parcial para controle do impacto
      ALLOCATE (IRN(nzexp))     !variável de controle da coluna para solver
      ALLOCATE (ICN(nzexp))     !idem
      ALLOCATE (IRNa(nzexp))    !controle para a linha
      ALLOCATE (ICNa(nzexp))    !idem
      ALLOCATE (IW(nzexp1))     !controle da ordem da matriz vetor
      ALLOCATE (IW1(2*n))       !idem
      ALLOCATE (IKEEP(n,3))     !idem
      
      rrjac=1.
      rfat=1.
      in=0
      indx=0
      itip=0
      inct=0
      ko=0
      ko1=0
      irn=0
      icn=0
      iw=0
      iw1=0
      ikeep=0
      
      ALLOCATE (pl0(nel,6,nnoel)) !posição inicial no elemento
      ALLOCATE (pl(nel,6,nnoel)) !posição atual no elemento
      
      ALLOCATE (rsl(nel,6,nnoel)) !resíduo dinâmico do método de Newmark
      ALLOCATE (qsl(nel,6,nnoel)) !resíduo dinâmico do método de Newmark
      
      ALLOCATE (drl1(nnoel))    !variável para cálculo de derivada de função de forma qualquer em função de variável adimensional
      ALLOCATE (drl2(nnoel))    !idem casca (não apagar)
      ALLOCATE (rh(nel))        !casca
      ALLOCATE (rf(nel))        !casca
      ALLOCATE (rp(nel))        !casca
      ALLOCATE (drex(6,nnoel))  !derivada da deformação longitudinal em relação aos graus de liberdade
      ALLOCATE (drey(6,nnoel))  !derivada da deformação longitudinal em relação aos graus de liberdade
      ALLOCATE (drez(6,nnoel))  !casca
      ALLOCATE (drexy(6,nnoel)) !derivada da distorcao em relação aos graus de liberdade
      ALLOCATE (drexz(6,nnoel)) !casca
      ALLOCATE (dreyz(6,nnoel)) !casca
      
      ALLOCATE (aurxy(6,nnoel)) !auxiliares para derivada de tensor de alongamento de Cauchy Green
      ALLOCATE (aur1xy(6,nnoel)) !idem
      ALLOCATE (aurxz(6,nnoel)) !idem
      ALLOCATE (aur1xz(6,nnoel)) !idem
      ALLOCATE (auryz(6,nnoel)) !idem
      ALLOCATE (aur1yz(6,nnoel)) !idem
      
      ALLOCATE (p1(3,6,nnoel))  !derivada da mudança de configuração em relação ao parâmetro nodal (direcao x)
      ALLOCATE (p2(3,6,nnoel))  !derivada da mudança de configuração em relação ao parâmetro nodal (direção y)
      ALLOCATE (p3(3,6,nnoel))  !casca
      ALLOCATE (daf(3,3,nnoel,6)) !escrita da derivada acima em forma matricial para inicio de montagem da matriz hessiana
      ALLOCATE (dch(3,3,6,nnoel)) !derivada do alongamento de cauchy green em relação aos parâmetros nodais
      

      aurxy=0.
      aur1xy=0.
      aurxz=0.
      aur1xz=0.
      auryz=0.
      aur1yz=0.
      
      
      dch=0.
      daf=0.
      p1=0.
      p2=0.
      p3=0.
      
      pl0=0.
      pl=0.
      pl00=0.
      drl1=0.
      drl2=0.
      rh=0.
      rf=0.
      rp=0.
      drex=0.
      drey=0.
      drez=0.
      derexy=0.
      drexz=0.
      

      ALLOCATE (rl1(nnoel))     !auxilia na montagem da função de forma de ordem qualquer
      ALLOCATE (rl2(nnoel))     !idem casca
      ALLOCATE (xs1(nnoel))     !coordenada nodal em variáveis adimensionais para montagem das funções de forma generalizadas
      ALLOCATE (xs2(nnoel))     !idem casca
      ALLOCATE (qsi(20))        !estádi de gaus no comprimento
      ALLOCATE (qsi3(20))       !estadio de gauss na altura
      ALLOCATE (w(20))          !peso de gauss no comprimento
      ALLOCATE (w3(20))         !peso de gaus na altura
      ALLOCATE (vint(n))        !vetor de forças internas e inerciais (global)
      ALLOCATE (reacao(n))      !vetor de forças reativas (global)
      ALLOCATE (vtl(nglel,nel)) !vetor de fornças internas e inerciais (local)
      
      ALLOCATE (frib(nglel))    !força de inércia base (local)
      ALLOCATE (fri(nglel))     !força de inércia total (local)
      ALLOCATE (rg(nel))        ! módulo de elasticidade transversal
      ALLOCATE (dp(n))          !incremento de posição
      ALLOCATE (df(n))          !incremento de forã externa
      ALLOCATE (v(n))           !vetor independente para solução do sistema
      
	
      ALLOCATE (vs(n))          !velocidade nodal global
      ALLOCATE (as(n))          !aceleração nodal global
      ALLOCATE (as1(n))         !aceleração nodal global (passado)
      ALLOCATE (rs(n))          !resíduo newmark (conhecido passado)
      ALLOCATE (qs(n))          !resíduo newmark (conhecido passado)
      
      
      ALLOCATE (rhs(n))         !vetor de carga para solver
      ALLOCATE (h(26,n))        !conjunto de linhas para matriz hessiana (até 8 nos)
      ALLOCATE (fi(nnoel))      !funções de forma 
      ALLOCATE (dfi(nnoel,2))   !derivada das funções de forma
      ALLOCATE (f(n))           !forças externas (global)
      ALLOCATE (aux(n))         !auxiliar geral
      ALLOCATE (vv(n))          !auxiliar
      ALLOCATE (p(n))           !valores nodais da posição atual dos nós
      ALLOCATE (solint(n))      !posição atual dos nós interpolada
      
      
      ALLOCATE (rkls(nglel,nglel)) !matriz hessiana local (não usada)
      ALLOCATE (re(nel))        !modulo de elasticidade longitudinal
      ALLOCATE (ri(nel))        !coeficiente de poisson
      ALLOCATE (rc(nel))        !adicional nlg
      ALLOCATE (resp(nel))      !casca
      ALLOCATE (ra(nel))        !casca
      ALLOCATE (p0(n))          !posição inicial dos nós
      ALLOCATE (p0i(n))         !posição passada para imposição de impacto geométrico
      ALLOCATE (rkl(nglel,nglel,nel)) !Hessiana local usada
      ALLOCATE (A(nzexp))       !controle solver
      ALLOCATE (Aag(nzexp))     !idem
      ALLOCATE (W27(nzexp1))    !idem
      ALLOCATE (valor(nel,nnoel)) !preparação de pos processamento
      ALLOCATE (fno(n))         !preparação de pós processamento
      valor=0.
      fno=0.
      
      rl1=0.
      rl2=0.
      xs1=0.
      xs2=0.
      qsi=0.
      w=0.
      vit=0.
      vtl=0.
      rg=0.
      dp=0.
      df=0.
      v=0.
      h=0.
      fi=0.
      dfi=0.
      f=0.
      aux=0.
      vint=0.
      reacao=0.

      vv=0.
      p=0.
      rkls=0.
      re=0.
      ri=0.
      rc=0.
      resp=0.
      ra=0.
      p0=0.
      rkl=0.
      a=0.
      w27=0.
      rhs=0.
      vs=0.
      qs=0.
      rs=0.
      as=0.
      as1=0.
      rsl=0.
      qsl=0.
      rmp=0.
      rmg=0.
      frib=0.
      fri=0.
      
      fd=0.
      brd=0.
      yb=0.
      hb=0.
      rib=0.
      
      
      qx=0.
      qy=0.
      qn=0.
      qt=0.

      rmz=0.
      rqy=0.
      rnx=0.
      
      write(5,*) 'grandezas para alocacao dinamica'
      write(5,*) 'n,noel,nnoel,nglel,nzexp,nnos,nel,nzexp1'
      write(5,*) n,noel,nnoel,nglel,nzexp,nnos,nel,nzexp1

      allocate(coordint(2,nnos),confini(2,nnos)) 
    

c     allocate(norma(normax,normay))


      allocate(norma(nnos,2))
      
      allocate(Nx(nnos))
      allocate(Ny(nnos))

      allocate(velint(2,nnos))

      

      

