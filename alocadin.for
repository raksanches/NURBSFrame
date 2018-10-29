      
c      nnoel=30 !no m�ximo 30 nos por elemento apesar do programa poder ser geral
      not=nnoel
      nglel=3*nnoel
      nzexp=(nglel**2)*nel*2
      nzexp1=Int((2*nzexp+3*n+1)*1.2)
      
      ALLOCATE (prsi0(nel,nnoel,ng3,3)) !vari�vel de estado (mede o tamanho da superf�cie de plastifica��)
      ALLOCATE (ep(nel,nnoel,ng3,3,6)) !deforma��o pl�stica acumulada para cada ponto de integra��o
      allocate (itagnos(nnos)) !identifica se o n� � de extremidade 

      allocate (sigmaabaixo(nnos,3),sigmaacima(nnos,3),momentof(nnos)) !tensoes e momento fletor

      allocate (anx(nnos),any(nnos))
      ALLOCATE (eeq(nel,nnoel,ng3,3)) !deforma��o equivalente
      ALLOCATE (temc(nel))      !temperatura total acima da barra
      ALLOCATE (temb(nel))      !temperatura total abaixo da barra
      ALLOCATE (dtemc(nel))     !acr�scimo por passo de tempo ou de carga
      ALLOCATE (dtemb(nel))     !acr�scimo por passo de tempo ou de carga
      ALLOCATE (calt(nel))      !coeficiente de dilata��o t�rmica do material por elemento
      
      prsi0=0.
      ep=0.
      eeq=0.
      temc=0.
      temb=0.
      dtemc=0.
      dtemb=0.
      calt=0.
      allocate(alfaglob(not))   !vetor global com os angulos dos elementos de dom�nio

      ALLOCATE (rmz(nel,not))   !esfor�o solicitante momento fletor
      ALLOCATE (rvy(nel,not))   !esfor�o solicitante for�a cortante
      ALLOCATE (rnx(nel,not))   !esfor�o solicitante for�a normal
      
      ALLOCATE (fd(nel,3))      !posi��o da fibra da se��o transversal (3 fibras por enquanto)
      ALLOCATE (brd(nel,3))     !largura da fibra
      ALLOCATE (yb(nel))        !posi��o do cg (m�dulo de engenharia)
      ALLOCATE (rib(nel))       !momento de in�rcia
      ALLOCATE (hb(nel,3))      !altura da fibra
      
      ALLOCATE (qx(nel,nnoel))  !for�a distribu�da horizontal
      ALLOCATE (qy(nel,nnoel))  !for�a distribu�da vertical   
      ALLOCATE (qn(nel,nnoel))  !for�a n�o conservativa ortogonal   
      ALLOCATE (qt(nel,nnoel))  !for�a n�o conservativa tangencial
      
      ALLOCATE (ro(nel))        !densidade (massa) do material
      ALLOCATE (ram(nel))       !amortecimento viscoso
      ALLOCATE (rmp(not,not,nel)) !matriz de massa compacta
      ALLOCATE (rmg(nglel,nglel,nel)) !matriz de massa ampliada
      ALLOCATE (sx(nel,nnoel))  !tensao x
      ALLOCATE (sy(nel,nnoel))  !tensao y
      ALLOCATE (sxy(nel,nnoel)) !tensao xy
      ALLOCATE (in(noel,noel)) !incid�ncia dos elementos
      ALLOCATE (rhook(nel,6,6)) !lei constitutiva (el�stica) por elemento
      ALLOCATE (indx(n)) ! vers�o casca (ainda n�o apagar)
      ALLOCATE (itip(100,4)) ! (vers�o casca (ainda n�o apagar)
      ALLOCATE (inct(55,3))     !(vers�o casca)
      ALLOCATE (ko(n))          !restri��o graus de liberdade / aplica��o de recalque
      ALLOCATE (ko1(n))         !restri��o parcial para controle do impacto
      ALLOCATE (IRN(nzexp))     !vari�vel de controle da coluna para solver
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
      
      ALLOCATE (pl0(nel,6,nnoel)) !posi��o inicial no elemento
      ALLOCATE (pl(nel,6,nnoel)) !posi��o atual no elemento
      
      ALLOCATE (rsl(nel,6,nnoel)) !res�duo din�mico do m�todo de Newmark
      ALLOCATE (qsl(nel,6,nnoel)) !res�duo din�mico do m�todo de Newmark
      
      ALLOCATE (drl1(nnoel))    !vari�vel para c�lculo de derivada de fun��o de forma qualquer em fun��o de vari�vel adimensional
      ALLOCATE (drl2(nnoel))    !idem casca (n�o apagar)
      ALLOCATE (rh(nel))        !casca
      ALLOCATE (rf(nel))        !casca
      ALLOCATE (rp(nel))        !casca
      ALLOCATE (drex(6,nnoel))  !derivada da deforma��o longitudinal em rela��o aos graus de liberdade
      ALLOCATE (drey(6,nnoel))  !derivada da deforma��o longitudinal em rela��o aos graus de liberdade
      ALLOCATE (drez(6,nnoel))  !casca
      ALLOCATE (drexy(6,nnoel)) !derivada da distorcao em rela��o aos graus de liberdade
      ALLOCATE (drexz(6,nnoel)) !casca
      ALLOCATE (dreyz(6,nnoel)) !casca
      
      ALLOCATE (aurxy(6,nnoel)) !auxiliares para derivada de tensor de alongamento de Cauchy Green
      ALLOCATE (aur1xy(6,nnoel)) !idem
      ALLOCATE (aurxz(6,nnoel)) !idem
      ALLOCATE (aur1xz(6,nnoel)) !idem
      ALLOCATE (auryz(6,nnoel)) !idem
      ALLOCATE (aur1yz(6,nnoel)) !idem
      
      ALLOCATE (p1(3,6,nnoel))  !derivada da mudan�a de configura��o em rela��o ao par�metro nodal (direcao x)
      ALLOCATE (p2(3,6,nnoel))  !derivada da mudan�a de configura��o em rela��o ao par�metro nodal (dire��o y)
      ALLOCATE (p3(3,6,nnoel))  !casca
      ALLOCATE (daf(3,3,nnoel,6)) !escrita da derivada acima em forma matricial para inicio de montagem da matriz hessiana
      ALLOCATE (dch(3,3,6,nnoel)) !derivada do alongamento de cauchy green em rela��o aos par�metros nodais
      

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
      

      ALLOCATE (rl1(nnoel))     !auxilia na montagem da fun��o de forma de ordem qualquer
      ALLOCATE (rl2(nnoel))     !idem casca
      ALLOCATE (xs1(nnoel))     !coordenada nodal em vari�veis adimensionais para montagem das fun��es de forma generalizadas
      ALLOCATE (xs2(nnoel))     !idem casca
      ALLOCATE (qsi(20))        !est�di de gaus no comprimento
      ALLOCATE (qsi3(20))       !estadio de gauss na altura
      ALLOCATE (w(20))          !peso de gauss no comprimento
      ALLOCATE (w3(20))         !peso de gaus na altura
      ALLOCATE (vint(n))        !vetor de for�as internas e inerciais (global)
      ALLOCATE (reacao(n))      !vetor de for�as reativas (global)
      ALLOCATE (vtl(nglel,nel)) !vetor de forn�as internas e inerciais (local)
      
      ALLOCATE (frib(nglel))    !for�a de in�rcia base (local)
      ALLOCATE (fri(nglel))     !for�a de in�rcia total (local)
      ALLOCATE (rg(nel))        ! m�dulo de elasticidade transversal
      ALLOCATE (dp(n))          !incremento de posi��o
      ALLOCATE (df(n))          !incremento de for� externa
      ALLOCATE (v(n))           !vetor independente para solu��o do sistema
      
	
      ALLOCATE (vs(n))          !velocidade nodal global
      ALLOCATE (as(n))          !acelera��o nodal global
      ALLOCATE (as1(n))         !acelera��o nodal global (passado)
      ALLOCATE (rs(n))          !res�duo newmark (conhecido passado)
      ALLOCATE (qs(n))          !res�duo newmark (conhecido passado)
      
      
      ALLOCATE (rhs(n))         !vetor de carga para solver
      ALLOCATE (h(26,n))        !conjunto de linhas para matriz hessiana (at� 8 nos)
      ALLOCATE (fi(nnoel))      !fun��es de forma 
      ALLOCATE (dfi(nnoel,2))   !derivada das fun��es de forma
      ALLOCATE (f(n))           !for�as externas (global)
      ALLOCATE (aux(n))         !auxiliar geral
      ALLOCATE (vv(n))          !auxiliar
      ALLOCATE (p(n))           !valores nodais da posi��o atual dos n�s
      ALLOCATE (solint(n))      !posi��o atual dos n�s interpolada
      
      
      ALLOCATE (rkls(nglel,nglel)) !matriz hessiana local (n�o usada)
      ALLOCATE (re(nel))        !modulo de elasticidade longitudinal
      ALLOCATE (ri(nel))        !coeficiente de poisson
      ALLOCATE (rc(nel))        !adicional nlg
      ALLOCATE (resp(nel))      !casca
      ALLOCATE (ra(nel))        !casca
      ALLOCATE (p0(n))          !posi��o inicial dos n�s
      ALLOCATE (p0i(n))         !posi��o passada para imposi��o de impacto geom�trico
      ALLOCATE (rkl(nglel,nglel,nel)) !Hessiana local usada
      ALLOCATE (A(nzexp))       !controle solver
      ALLOCATE (Aag(nzexp))     !idem
      ALLOCATE (W27(nzexp1))    !idem
      ALLOCATE (valor(nel,nnoel)) !prepara��o de pos processamento
      ALLOCATE (fno(n))         !prepara��o de p�s processamento
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

      

      

