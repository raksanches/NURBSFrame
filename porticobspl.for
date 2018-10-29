c     ***************************************************
c      Pprograma de analise n�o linear geom�trica e f�sica
c      de p�rtico plano incluindo impacto
c      Formula��o posicional - Lagrangiana totl com
c      Cinem�tica exata
c      Autores:   Humberto Breves Coda
c                 Rodrigo Ribeiro Paccola
c     Institui��o EESC - USP
c     Departamento de Engenharia de Estruturas

c     Adapta��es: Ana Kelle Malaguti 
c                 Rodolfo A. K. Sanches
c     Universidade Tecnol�gica Federal do Paran�
c     Campus Pato Branco - 2013
c      ********************************************
c      
c      *******************************************
c      Inicio do programa principal
c      *******************************************   

c     ********************************************
c     O m�dulo BIG_STUFF serve para se poder definir aloca��o din�mica
c     de vari�veis �nicas (globais) para todo o programa
c     sem a necessidade de se utilizar lista de argumentos
c     nas chamadas de subrotinas
c     Nas subrotinas fazemos a chamada do modulo 
c     e adicionamos o include do subprograma alocadin.for  
c     que aloca o tamanho das vari�veis
c     O include do bloclagatbig.for faz o mesmo papel 
c     para as vari�veis est�ticas 
c    *************************************************


c     ******************************************
c     O significado de cada vari�vel din�mica est� 
c     listado no subprograma alocadin.for
c     *******************************************
      MODULE BIG_STUFF

      integer*4,allocatable:: IRN(:),ICN(:),IW(:),IKEEP(:,:),IW1(:)
     $     ,ic(:,:),ko(:),ko1(:),inct(:,:),kn(:),im(:,:),itip(:,:)
     $     ,indx(:),in(:,:),irna(:),icna(:),glgl(:,:) ,itr(:,:),ngl(:)
     $     ,iglv(:,:),notl(:),itipoel(:),itagnos(:)

      REAL*8,ALLOCATABLE::A(:),w27(:),p0(:),p0i(:),rkl(:,:,:),Aag(:)
     $     ,re(:),ri(:),rc(:),resp(:),ra(:),f(:),aux(:),vv(:),p(:),df(:)
     $     ,rkls(:,:),fi(:),dfi(:,:),dp(:),v(:),h(:,:),vs(:),as(:)
     $     ,as1(:), rhook(:,:,:),sx(:,:),sy(:,:),sxy(:,:),ro(:),ram(:)
     $     ,rs(:),qs(:), qsi(:),w(:),vint(:),vtl(:,:),rg(:),qsi3(:)
     $     ,w3(:),reacao(:), rl1(:),rl2(:),xs1(:),xs2(:),pl0(:,:,:),pl(:
     $     ,:,:),qsl(:,:,:), drl1(:),drl2(:),rh(:),rf(:),rp(:),RHS(:),
     $     drex(:,:),drey(:,:) ,drez(:,:),drexy(:,:),drexz(:,:),dreyz(:
     $     ,:), p1(:,:,:),p2(:,:,:) ,p3(:,:,:),daf(:,:,:,:),dch(:,:,:,:)
     $     , aurxy(:,:),aur1xy(:,:) ,aurxz(:,:),aur1xz(:,:),rsl(:,:,:)
     $     ,auryz(:,:),aur1yz(:,:) ,valor(:,:),fno(:),rmp(:,:,:),rmg(:,:
     $     ,:), fri(:),frib(:),rmfi(: ,:),fd(:,:),brd(:,:),yb(:),hb(:,:)
     $     ,rib(:) ,qx(:,:),qy(:,:),qn(: ,:),qt(:,:),rmz(:,:),rvy(:,:)
     $     ,rnx(:,:),alfaglob(:),solint(:),Anx(:),Any(:),sigmaabaixo(:
     $     ,:),sigmaacima(:,:),momentof(:),coordint(:,:),confini(:,:)
     $     ,norma(:,:),Nx(:),Ny(:),velint(:,:),pesoNURBS(:)


c      *************************************************************
c      para plasticidade e temperatura temos as seguintes vari�veis
c      seus significados tamb�m est�o no subprograma alocadina
c      *************************************************************
      REAL*8,ALLOCATABLE::prsi0(:,:,:,:),ep(:,:,:,:,:),eeq(:,:,:,:)
     $     ,defl(:,:),rhdg(:,:),rhvg(:,:),temc(:),temb(:),dtemc(:)
     $     ,dtemb(:) ,calt(:)

     
      CHARACTER*30 Entrada      !nome do arquivo de entrada

      END MODULE BIG_STUFF
      
c     USE MSFLIB
c     use msimsl
      USE BIG_STUFF
      
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      
c     **********************************************************
c     Este subprograma serve para as vari�veis menores, aloca��o est�tica
c     **************************************************************
      include 'bloclagatbig.for'

      INTEGER(2) IHRI,IMINI,ISECI,I100THI,IHRF,IMINF,ISECF,I100THF !so uso para marcar tempo
      
c     **********************************************************
c     Abertura dos arquivos de dados e saida
c     **************************************************************
      write(*,*) 'nome do arquivo de entrada'
      READ(*,*)Entrada

      open(unit=3,file=entrada,status='unknown')
      open(unit=5,file='saida',status='unknown')
      open(unit=11,file='grafico1.dat',status='unknown')
      open(unit=15,file='grafico2.dat',status='unknown')
      open(unit=12,file='contrdesl1.dat',status='unknown')
      open(unit=16,file='contrdesl2.dat',status='unknown')
      open(unit=7,file='desloc',status='unknown')
      open(unit=8,file='desnov',status='unknown')
      open(unit=66,file='ma26coments',status='unknown')
      open(unit=9,file='tensaoacima',status='unknown')
      open(unit=19,file='tensaoabaixo',status='unknown')
      open(unit=90,file='esforcos',status='unknown')
      open(unit=4,file='entrplast.txt',status='unknown')
      open(unit=50,file='auxten.txt',status='unknown')
      open(unit=51,file='noatual.txt',status='unknown')
      open(unit=123,file='graficosdesloc.dat',status='unknown')
      open(unit=38,file='derivadas.txt',status='unknown')
c$$$      open(unit=43,file='deslocno3.txt',status='unknown') 
c$$$      open(unit=44,file='deslocno4.txt',status='unknown') 
c$$$      open(unit=45,file='deslocno5.txt',status='unknown') 
c$$$      open(unit=46,file='deslocno6.txt',status='unknown') 
      open(unit=48,file='velocidade8.txt',status='unknown') 
      open(unit=49,file='deslocno12.txt',status='unknown') 
      open(unit=2,file='analitica.txt',status='unknown')  
      open(unit=1,file='deslocno4.txt',status='unknown')
      open(unit=43,file='Pxu.txt',status='unknown') 
      open(unit=44,file='Pxv.txt',status='unknown') 

      
c      Capturando instante de inicio do processamento
c      CALL GETTIM (IHRI, IMINI, ISECI, I100THI)
      write(*,*) 'inicio',ihri,imini,iseci,i100thi

      pi=2.*datan2(1.d00,0.d00) !definido a grandeza trigonom�trica pi

      call dadosiniciais ! Camada da subrotina que faz a leitura dos dados iniciais
      call lerplast             ! Leitura de informa��es sobre plasticidade
      call matrizmassa          !subrotina que gera a matriz de massa que n�o varia em todo o processo

c     **********************************************
c     inicio do processo din�mico (ou incremental) / iterativo (impl�cito)
c     **********************************************
      ipt=0 !zerando o passo de tempo

       do ifc=1,nfc
c      do ifc=1,nfc !La�o para numero de fases de carga (comportamento temporal de condi��es de contorno)
       read(3,*) !nova lista para cada fase de carga
       call dadosdin !Entrada dos dados que podem variar com o tempo


c     coletando pontos de gauss pra impressao da geometria
       call gauss(ng,qsi,w)
c$$$c     encontra e imprime coordenadas da estrutra para pontos de gauss
c$$$       call imprege
c$$$c     encontra vetor normal unit�rio no n�
c$$$       call vetnorm


c     dividindo em passos de carga:
c$$$
c$$$          if (ifc.lt.240) then
c$$$          deltapx=0.d0
c$$$          deltapy=-50.d0
c$$$          deltapz=0.d0
c$$$          else 
c$$$          deltapx=0.d0
c$$$          deltapy=-3.d0
c$$$          deltapz=0.d0  
c$$$          endif
c$$$         
c$$$c          do i=1,nnc
c$$$             df(glgl(ncarreg,1))=df(glgl(ncarreg,1))+deltapx
c$$$             df(glgl(ncarreg,2))=df(glgl(ncarreg,2))+deltapy
c$$$             df(glgl(ncarreg,3))=df(glgl(ncarreg,3))+deltapz
c$$$c          enddo
c$$$         
        

       do ipc=1,npt !Passos de carga ou tempo
 
          
          ipt=ipt+1             !Atualiza��o do tempo
          qs=p/rbn/(dt*dt)+vs/rbn/dt+(1.d0/2.d0/rbn-1.d0)*as !Res�duo do M�todo de Newmark
          rs=vs+dt*(1.d0-rgn)*as  !Res�duo d� M�todo de Newmark
          ia=0                  !Contador de itera��es
          
c     variacao temporal da temperatura, por enquanto segue o mesmo padr�o da carga externa
          temb=dtemb*           !Abaixo, din�mico
     #(ca+cb*ipc*dt+cc*(ipc*dt)**2+cd*sin(ce*ipc*dt)
     #+cf*cos(cg*ipc*dt)
     #+crh*exp(cri*ipc*dt)+cj*exp(ck*ipc*dt))
          
        temc=dtemc*  !Acima din�mico
     #  (ca+cb*ipc*dt+cc*(ipc*dt)**2+cd*sin(ce*ipc*dt)
     #  +cf*cos(cg*ipc*dt)
     #  +crh*exp(cri*ipc*dt)+cj*exp(ck*ipc*dt))

c     variacao temporal da carga incluindo distribu�da conservativa j� alocada no vetor de for�as
      
        do i=1,n                !N�mero total de graus de liberdade
           if(ko(i).eq.0) f(i)=df(i)* !din�mico
     #(ca+cb*ipc*dt+cc*(ipc*dt)**2
     #+cd*dsin(ce*ipc*dt)+cf*dcos(cg*ipc*dt)
     #+crh*dexp(cri*ipc*dt)+cj*dexp(ck*ipc*dt))
        enddo


        f(68) = -1000.*ipc
c     variacao temporal da posi��o para est�tico
	  do i=1,n   !N�mero total de graus de liberdade
         if(ko(i).eq.1) p(i)=p0(i)+dp(i)*  !din�mico
     #   (ba+bb*ipc*dt+bc*(ipc*dt)**2
     #   +bd*sin(be*ipc*dt)+bf*cos(bg*ipc*dt)
     #   +bh*exp(bri*ipc*dt)+bj*exp(bk*ipc*dt))
        enddo


        call alocanco  !Aloca��o do carregamento ortogonal (n�o conservativo)

10      continue !La�o de itera��es
        if (ia.le.6) call destroca1 !Atualiza��o das vari�veis ap�s solu��o do sistema
        if (ia.gt.6) call destroca1 !idem com bissec��o (n�o utilizado)

        call matriz !Chamada da cria��o das matrizes e vetores do problema
 
        call condcon !Aplica��o da condi��o de contorno

        call ma27 !Solver de dom�nio P�blico HSL (UK)
        call calculanorma       !Calculo da norma de parada das itera��es
        ia=ia+1                 !contador de itera��es
        if (ia.gt.1005) goto 20 !Limite de itera��es
        if ((rnorma.gt.tol).or.(ia.lt.2)) goto 10 !Teste de parada
 20     continue                !Fim das itera��es
        
        call destroca1 !Utilizando a �ltima itera��o ainda n�o aproveitada
        as1=p/rbn/(dt*dt)-qs      !Atualizando acelera��o presente para ser usada
        vs=vs+dt*(1.d0-rgn)*as+rgn*dt*as1 !Acertando velocidade passado (desnecess�rio guardar presente)
        as=as1                  !acelera��o passado sendo atualizada
        write(*,*) 'numero de iteracoes',ia,ipt,ipt



           write(*,*) 'Interpolando solu��o' 
           
c           rnxint=0.d0
c           rnyint=0.d0
           
           solint=0.d0
           confini=0.d0 
           coordint=0.d0
           itagnos=0
           velint=0.d0
           do jel=1,nel
              j=jel
              testefunc=0.d0

c     write(*,*) jel,itipoel(jel)
              if((itipoel(jel).ne.2).and.(itipoel(jel).ne.3))then
                 do inolocal=2,3
                    i=ic(jel,inolocal)
                    
                    if(itagnos(i).eq.0)then
                       xsi1=2.d0*inolocal-5.d0
                       
                       do ifunclocal=1,4
                          
                          ifunc=ic(jel,ifunclocal)
                                                   
                          call formaederi0
                          
                          testefunc=testefunc+fi(ifunclocal)
c                          rnxint(i)=rnxint(i)+fi(ifunclocal)*RNx(i)
 
                         do idof=1,3
                             
                             solint(glgl(i,idof))=solint(glgl(i,idof))
     $                            +p(glgl(ifunc,idof))*fi(ifunclocal)
                             if (idof .lt. 3) then
                                coordint(idof,i)=coordint(idof
     $                               ,i)+fi(ifunclocal)*p(glgl(ifunc
     $                               ,idof))
                              
                               

                                confini(idof,i)=confini(idof,i)
     $                               +fi(ifunclocal)*p0(glgl(ifunc
     $                               ,idof))
                              
c                               obter velocidade do n� interpolada
                                velint(idof,i)=velint(idof,i)
     $                               +fi(ifunclocal)*vs(glgl(ifunc
     $                               ,idof))
                               
                                                                  
                             endif 
                             
                          enddo !idof
                          
                       enddo    !ifunclocal
                       
                       itagnos(i)=1
c$$$  
c$$$  else
c$$$  
                    endif
                    
                 enddo          !inolocal
                 
              endif             !itipoel
              
              
              if(itipoel(jel).eq.2)then
                 
                 do inolocal=1,3,2
                    i=ic(jel,inolocal)
                    if(itagnos(i).eq.0)then
                       do ifunclocal=1,4
                          ifunc=ic(jel,ifunclocal)
                          xsi1=1.d0*inolocal-2.d0
                          
                          call formaederi0
                          
c     write(*,*)jel,xsi1                  
c     read(*,*)
                          
                          do idof=1,3
                             
                             solint(glgl(i,idof))=solint(glgl(i,idof))
     $                            +p(glgl(ifunc,idof))*fi(ifunclocal)
                             if (idof .lt. 3) then
c obter a coordenada do n� interpolada
                                coordint(idof,i)=coordint(idof
     $                               ,i)+fi(ifunclocal)*p(glgl(ifunc
     $                               ,idof))

                                confini(idof,i)=confini(idof,i)
     $                               +fi(ifunclocal)*p0(glgl(ifunc
     $                               ,idof))

c     obter velocidade do n� interpolada
                                velint(idof,i)=velint(idof,i)
     $                               +fi(ifunclocal)*vs(glgl(ifunc
     $                               ,idof))
                                
                             endif 
                          enddo
                          
                       enddo
                       itagnos(i)=1
                    else
                    endif
                    
                 enddo
              endif
              
      
              if(itipoel(jel).eq.3)then
                 do inolocal=2,4,2
                    i=ic(jel,inolocal)
                    if(itagnos(i).eq.0)then
                       do ifunclocal=1,4
                          ifunc=ic(jel,ifunclocal)
                          xsi1=1.d0*inolocal-3.d0
                          
                          call formaederi0
                          
c     write(*,*)jel,xsi1                  
c     read(*,*)
                 
                         
         
                          do idof=1,3
                             
                             solint(glgl(i,idof))=solint(glgl(i,idof))
     $                            +p(glgl(ifunc,idof))*fi(ifunclocal)
                             if (idof .lt. 3) then
                                coordint(idof,i)=coordint(idof
     $                               ,i)+fi(ifunclocal)*p(glgl(ifunc
     $                               ,idof))
                                confini(idof,i)=confini(idof,i)
     $                               +fi(ifunclocal)*p0(glgl(ifunc
     $                               ,idof))
c     obter velocidade do n� interpolada
                                velint(idof,i)=velint(idof,i)
     $                               +fi(ifunclocal)*vs(glgl(ifunc
     $                               ,idof))
                                
                             endif 
                          enddo
                          
                       enddo
                       itagnos(i)=1
                    else
                    endif
                    
                 enddo
              endif
              
              
c              write(*,*)testefunc,jel
              
c$$$              
c$$$
c$$$              if(itipoel(jel).eq.3)then
c$$$
c$$$                 do inolocal=2,3
c$$$
c$$$                    i=ic(jel,inolocal)
c$$$
c$$$
c$$$                    if(itagnos(i).eq.0)then
c$$$
c$$$                    do ifunclocal=1,3
c$$$
c$$$                       ifunc=ic(jel,ifunclocal)
c$$$                       
c$$$                       xsi1=2.d0*inolocal-5.d0
c$$$                       
c$$$                       call formaederi0
c$$$                       
c$$$                       do idof=1,3
c$$$                          
c$$$                          solint(glgl(i,idof))=solint(glgl(i,idof))
c$$$     $                        +p(glgl(ifunc,idof))*fi(ifunclocal)
c$$$                       enddo
c$$$
c$$$                    enddo
c$$$                    itagnos(i)=1
c$$$                    else
c$$$                    endif
c$$$
c$$$                 enddo
c$$$              endif
              
           enddo
           
c     Montagem do gr�fico de deslocamento dos n�s

c$$$      write(43,*) dt*ipt,(coordint(1,3)- confini(1,3))
c$$$      write(44,*) dt*ipt,(coordint(1,4)- confini(1,4))
c$$$      write(45,*) dt*ipt,(coordint(1,5)- confini(1,5))
c$$$      write(46,*) dt*ipt,(coordint(1,6)- confini(1,6))          
c$$$      write(48,*) dt*ipt,(coordint(1,8)- confini(1,8)) 
!          write(49,*) dt*ipt,(coordint(1,12)- confini(1,12)) 
!         write(1,*) dt*ipt,(coordint(1,4)- confini(1,4)) 
c$$$c     write(48,*) dt*ipt,(coordint(1,8)- confini(1,8)) 
c$$$c     write(48,*) dt*ipt,(coordint(1,8)- confini(1,8)) 

!        write(48,*) dt*ipt,(velint(1,8))
           write(43,*) (-100. + p(glgl(23
     $          ,1)))/100., 1000.*ipc*10000./(1000000000.)
           write(44,*) ( p(glgl(23
     $          ,2)))/100.,1000.*ipc*10000./(1000000000.)
       
          
      write(5,*)'deslocamentos verticais'

           do i=1,nnos  
c             write(5,*)p0(glgl(i,1)),p(glgl(i,2)),solint(glgl(i,2))

              write(123,*)solint(glgl(i,1)),solint(glgl(i,2))
     $             ,solint(glgl(i,3))

           enddo

           

        rteste=int(ipt/1.)-ipt/1. !Frequencia de impressao para posproc
        if ((dabs(rteste).lt.0.01).or.(ipt.eq.1)) then
           if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     $          .lt.1.0d-6)) then
           call impreten
           call tennodal
           call saidaposparaview

           endif
c       pos-processamentos gerais
           call saida           !Lista completa de posi��es e tempo
c           call saidaposnov     !Pos processamento do deslocamento
c           call tensaoacima     !Calcula as tens�es nominais tensao nominal
c           call saidaposnovtensaoacima !Pos processa tensoes nominais
c           call tensaoabaixo    !Cacula tensoe principais
c           call saidaposnovtensaoabaixo !Pos-processa tensoe principais
c           call momento         !Calcula os esforcos solicitantes
c     call saidaposmomento !Pos processa os esforcospor n�
c           call momentonovo     !pos processa os esforcos por elemento	
        else
        endif
      enddo

      enddo
      write(*,*)'fim'
c     Calculando o tempo de processamento
c     CALL GETTIM (IHRF, IMINF, ISECF, I100THF) !Hora final de processamento
      WRITE (*,*) 'TEMPO DE EXECUCAO [s]:',
     #3600*(IHRF-IHRI)+60*(IMINF-IMINI)+(ISECF-ISECI)
      
      WRITE (5,*) 'TEMPO DE EXECUCAO [s]:',
     #3600*(IHRF-IHRI)+60*(IMINF-IMINI)+(ISECF-ISECI)
      stop
      end																
      
c     *******************************************
c     Etradas dos dados iniciais
c     deve-se observar que algumas vari�veis s�o calculadas
c     para os elementos finitos
c     *******************************************   
      subroutine dadosiniciais  !Comentado anteriormente (n�o passamos argumentos vari�veis globais)
      USE BIG_STUFF             !Defini��o das vari�veis din�micas
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for' !defini��o das vari�veis est�ticas
      
c     Novidade
c     Para problemas din�micos podemos encontrar a posi��o est�tica de equil�brio     de determinada carga permanente usando a op��o abaixo
      read(3,*)                 !Dizer se deve ler posicao atual diferente da inicial 1 sim, 0 nao
      read(3,*) iposic
      
c     nnos = n�mero de n�s
c     nel = n�mero de elementos
c     nfc = n�mero de fases de carga
c     tol = toler�ncia de parada (posi��o)
c     ngf = n� para gr�fico 1
c     ngf1 = n� para gr�fico 2
      read(3,*)
      write(5,*)'nnos,ne,nfc,tol,ngf,ngf1'
      read(3,*) nnos,nel,nfc,tol,ngf,ngf1
      write(5,20) nnos,nel,nfc,tol,ngf,ngf1
 20   format(3(1x,I7),1x,f20.7,2(1x,I7))
c     *************************************************
      read(3,*)
      read(3,*) ng3,nprint      !N�mero de pontos de Gauss na altura e
      write(5,30)
 30   format(2(1x,I7))

c     frequ�ncia para impress�o de resultados
c     O n�mero de pontos de integra��o na longitudinal �
c     igual ao n�mero de n�s do Elemento finito
      
      ALLOCATE (notl(nel),itipoel(nel),pesoNURBS(nnos)) !Esta � uma das aloca��es din�micas feita antes
      
c     Novidade
      read(3,*)                 ! Elemento, n�mero de n�s do elemento
      write(5,*) ' Elemento, n�mero de n�s do elemento'
      do j=1,nel
         read(3,*) k,notl(k)
         write(5,30) k,notl(k)
      enddo
      
      k=0
      do i=1,nel
         if(notl(i).gt.k) k=notl(i)
      enddo
      nnoel=k+1                 !+1 por caussa do solver (esparsidade)
      
      read(3,*)                 !tipo elemento (1 dom�nio 2 contorno esq 3 contorno dir)
      write(5,*) 'tipo elemento quanto ao contorno'
      do i=1,nel
         read(3,*)itipoel(i)
         write(5,*)itipoel(i)
      enddo
      read(3,*)                 !peso das NURBS
      write(5,*) 'peso das NURBS'
      do i=1,nnos
         read(3,*)pesoNURBS(i)
         write(5,*)pesoNURBS(i)
      enddo
      
c     Lembrando que o c�lculo do n�mero de graus de liberdade 
c     depende das r�tulas precisamos chamar uma subrotina que vai
c     alocar dinamicamente algumas vari�veis separadamente da alocadin,
c     fazer a leitura da incid�ncia dos elementos, antes das coordenadas dos nos
c     bem como da condi��o de rotula��o dos nos das barras
      call rotular

      include 'alocadin.for'    !Aloca��o din�mica das demais vari�veis
      
c     coordenadas dos nos n�o leio giro
c     no,x,y,z
      read (3,*)
      read (3,*)	  
      
      write(5,*) 'coordenadas dos nos'
      write(5,*) 'no,x,y,teta'
      
      do i=1,nnos
         read(3,*) k,p0(glgl(k,1)),p0(glgl(k,2))
         write(5,10) k,p0(glgl(k,1)),p0(glgl(k,2)),p0(glgl(k,3))
      enddo
      
 10	  format(1x,i7,3(1x,f27.10))
      
c     Propriedades f�sicas dos elementos
c     Elemento,e,mi,altura,area,ybarra,momin,ro,cm
c     Ainda n�o estamos entrando com as fibras no elemento diretamente pelas
c     limita��es dos usu�rios atuais
      read (3,*)
      read (3,*)	  
      
      write(5,*) 'elemento,e,mi,altura,area,ybarra,momin,ro,cm'
      
      do j=1,nel
         read(3,*) k,re(k),ri(k),resp(k),ra(k),yb(k),rib(k),ro(k),ram(k)
         write(5,12) k,re(k),ri(k),resp(k),ra(k),yb(k),rib(k),ro(k)
     $        ,ram(k)
         rg(k)=re(k)/(1.d0+ri(k))/2.d0 ! M�dulo de elasticidade transversal
      enddo
      
c     Esta subrotina reduz o n�mero de fibras a tr�s apensa,
c     n�o tirando a generalidade do programa, 
c     mas apenas o reduzindo �s necessidades
c     atuais dos usu�rios do SET-EESC-USP 
      call secoes 
      
c     ************************************************
c     subrotina que define a lei constitutiva
c     ****************************************************
      call definelei

      read(3,*)                 !numero de nos com velocidade inicial
      write(5,*) 'numero de nos com velocidade inicial'
      read(3,*) nnv
      write(5,*) nnv
      read(3,*)                 !vx,vy
      write(5,*) 'no,vx,vy'
      do i=1,nnv
         read(3,*) k,vs(glgl(k,1)),vs(glgl(k,2))
         write(5,12) k,vs(glgl(k,1)),vs(glgl(k,2))
      enddo
      
      
      read(3,*)                 ! Impacto sem atrito, regi�o quadrada
      read(3,*)                 ! Entrar com coordenadas da caixa de impacto e coeficiente de atrito (0-1)
      read(3,*)                 ! xp,yp,xc,yc,rimp
      read(3,*)  xp,yp,xc,yc,rimp
      write(5,*)'coord. da caixa de impacto e coef. de atrito'
      write(5,13) xp,yp,xc,yc,rimp
      
c     Chamada da subrotina que prepara as grandezas iniciais lidas para os n�o
c     em grandezas nos elementos finitos
c     Caractcter�stica do m�todo posicional
      call geometriainicial
      

 12   format(1x,i7,8(1x,f27.10))
      
 13   format(5(1x,f27.10))
 15   format(1x,i7,6(1x,f27.10))
      
      return
      end

c      
c      *******************************************
c      Novidade
c      Nesta subrotina os valores necessarios para 
c      gerar as funcoes de forma sao calculados
c      *******************************************   
      subroutine suportefuncoesdeforma
      USE BIG_STUFF
	 implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	 include 'bloclagatbig.for'
c     C�lculo da posi��o dos pontos nodais no espa�o de Gauss
c     para elemento de qualquer ordem
      do is=1,notl(j)
c$$$         if(itipoel(j).eq.2)then !esquerda           
c$$$            xsi1(1)=0
c$$$            xsi1(2)=1
c$$$
c$$$         endif
c$$$         if(itipoel(j).eq.1)then !esquerda           
c$$$            xsi1(2)=0
c$$$            xsi1(3)=1
c$$$         endif
c$$$         if(itipoel(j).eq.1)then !esquerda           
c$$$            xsi1(3)=0
c$$$            xsi1(4)=1
c$$$         endif

         
         xs1(is)=(2.*is-(1.+notl(j)))/(notl(j)-1.) 
      enddo
      return
      end



c      *******************************************
c      Etrada de dados que podem mudar com o tempo
c      Na realidade para fases de carga
c      que definem intervalos de tempo cujas condi��es de contorno n�o se alteram
c      *******************************************   
      subroutine dadosdin
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      read(3,*)                 !numero de nos restritos, numero de nos carregados, num. de dts
      write(5,*)'num. de nos restr., num. de nos carreg., num. de dts'
      read(3,*) nnr,nnc,npt
      write(5,13) nnr,nnc,npt
      
      read(3,*)                 !nelcar(num. elementos carregados qx,qy),nelnco(num. elementos carregados qn,qt),elm. com temperatura      
      write(5,*) 'Num.El.com: qx,qy; qn,qt; T1,T2'
      read(3,*) nelcar,nelnco,ntem
      write(5,13)nelcar,nelnco,ntem
      
c     Novidade (par�metros de Newmark estabelecidos pelo usu�rio
c     Estruturas isoest�ticas e hiperest�ticas Beta=0.25, Gama=0.5
c     M�quinas Beta=0.33, Gama=0.55
c     Impacto Beta=0.5, Gama=1.5 
c     Problemas est�ticos Beta=1.0, Gama=1.0 
c     Leitura dos parametros do integrador temporal de newmark
c     beta, gama, dt (intervalo de tempo)
      read(3,*)
      read(3,*)
      read(3,*) rbn,rgn,dt
      write(5,*) 'beta,gma,dt'
      write(5,16) rbn,rgn,dt
      
      read (3,*)
      read (3,*)	  
c     restricao de graus de liberdade
c     no,ko1,ko2,ko3
      write(5,*) 'restricao de graus de liberdade'
      write(5,*) 'no,ko1,ko2,ko3'
      
      do i=1,nnr
         read(3,*) k,ko(glgl(k,1)),ko(glgl(k,2)),ko(glgl(k,3))
         write(5,13) k,ko(glgl(k,1)),ko(glgl(k,2)),ko(glgl(k,3))
      enddo
      
      read (3,*)
      read (3,*)	  
c     leitura dos valores dos recalques
c     no,px,py,pz,a,b
      write(5,*) 'leitura dos valores dos recalques'
      write(5,*) 'no,px,py,pz'
      
      do i=1,nnr
         read(3,*) k,dp(glgl(k,1)),dp(glgl(k,2)),dp(glgl(k,3))
         write(5,15) k,dp(glgl(k,1)),dp(glgl(k,2)),dp(glgl(k,3))
      enddo
      
      read(3,*)
c     comportamento temporal do movimento de base
      read(3,*) ba,bb,bc,bd,be,bf,bg,bh,bri,bj,bk
      write(5,*) 'movimento de base: ba,bb,bc,bd,be,bf,bg,bh,bri,bj,bk'
      write(5,16) ba,bb,bc,bd,be,bf,bg,bh,bri,bj,bk
      
      read (3,*)
      read (3,*)	  
c     nos com carregamento prescrito
c     no,f1,f2,m3
      write(5,*) 'nos com carregamento prescrito'
      write(5,*) 'no,f1,f2,f3,ma'
      
      do i=1,nnc
         read(3,*) k,df(glgl(k,1)),df(glgl(k,2)),df(glgl(k,3))
         write(5,15) k,df(glgl(k,1)),df(glgl(k,2)),df(glgl(k,3))
      ncarreg=k
      enddo
      
c     valor das cargas distribuidas conservativas nas extremidades
c     calculo para os outros nos (comportamento linear do caregamento)
      read(3,*)                 !cargas nas extremidades
      read(3,*)                 ! elemeto, qx no inicial, qx no final, qy no inicial, qy no final
      write(5,*) 'cargas nas extremidades'
      write(5,*) 'elemeto, qx no inicial, qx no final, 
     #qy no inicial, qy no final'
      
      do j=1,nelcar
         read(3,*) k,qx(k,1),qx(k,notl(k)),qy(k,1),qy(k,notl(k))
         write(5,12) k,qx(k,1),qx(k,notl(k)),qy(k,1),qy(k,notl(k))
         dqx=(qx(k,notl(k))-qx(k,1))/(notl(k)-1.)
         dqy=(qy(k,notl(k))-qy(k,1))/(notl(k)-1.)
         do jk=2,notl(k)-1
            qx(k,jk)=qx(k,jk-1)+dqx
            qy(k,jk)=qy(k,jk-1)+dqy
         enddo
      enddo
      
c     valor das cargas distribuidas n�o conservativas nas extremidades
c     calculo para os outros nos (comportamento linear do caregamento)
      read(3,*)                 !cargas nas extremidades
      read(3,*)                 ! elemeto, qn no inicial, qn no final, qt no inicial, qt no final
      write(5,*) 'cargas nas extremidades'
      write(5,*) 'elemeto, qn no inicial, qn no final, 
     #qt no inicial, qt no final'
      do j=1,nelnco
         read(3,*) k,qn(k,1),qn(k,notl(k)),qt(k,1),qt(k,notl(k))
         write(5,12) k,qn(k,1),qn(k,notl(k)),qt(k,1),qt(k,notl(k))
         dqn=(qn(k,notl(k))-qn(k,1))/(notl(k)-1.)
         dqt=(qt(k,notl(k))-qt(k,1))/(notl(k)-1.)
         do jk=2,notl(k)-1
            qn(k,jk)=qn(k,jk-1)+dqn
            qt(k,jk)=qt(k,jk-1)+dqt
         enddo
      enddo
      
c     valor das temperaturas constante ao longo do elemento
      read(3,*)                 !temperatura e coeficiente de dilata��o t�rmica
      read(3,*)                 ! elemeto, temp.acima, temp. abaixo, dil. termica
      write(5,*) 'temperatura'
      write(5,*) 'elemeto,  temp.acima, temp. abaixo, dilata��o t�rmica'
      
      do j=1,ntem
         read(3,*) k,dtemc(k),dtemb(k),calt(k)
         write(5,12) k,dtemc(k),dtemb(k),calt(k)
      enddo
      
c     Colocar no carregamento nodal a carga distribuida conservativa
c     A n�o conservativa deve ser feita para cada passo de tempo e n�o para cada fase de carga
      call alocadistr
      
c     comportamento temporal da carga
      read(3,*)
      read(3,*) ca,cb,cc,cd,ce,cf,cg,crh,cri,cj,ck
      write(5,*)'compotamento carga: ca,cb,cc,cd,ce,cf,cg,crh,cri,cj,ck'
      write(5,16) ca,cb,cc,cd,ce,cf,cg,crh,cri,cj,ck
  
12	format(1x,i7,7(1x,f27.10))
13	format(1x,i7,6(1x,i7))
14	format(2(1x,i7),6(1x,f27.10))
15	format(1x,i7,6(1x,f27.10))
16	format(11(1x,f27.10))
      return
      end

c$$$
c$$$      subroutine geometriainicial
c$$$      USE BIG_STUFF	 
c$$$	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
c$$$	include 'bloclagatbig.for'
c$$$
c$$$!     A posi��o angular para cada elemento � calculada
c$$$!     A posi��o de transla��o � a mesma do elemento global
c$$$      do j=1,nel
c$$$	 do k=1,notl(j)
c$$$	  nog=ic(j,k)
c$$$        do igl=1,2 !igl � o grau de liberdade x,y
c$$$         pl0(j,igl,k)=p0(glgl(nog,igl)) !posi��o inicial no elemento
c$$$	  enddo
c$$$       enddo 
c$$$
c$$$c      calcular a derivada em cada ponto do elemento
c$$$       do is1=1,notl(j)
c$$$        kk=is1
c$$$	  xsi1=(2.*is1-(1.+notl(j)))/(notl(j)-1.)
c$$$
c$$$c       ***********************
c$$$c       variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c$$$        call suportefuncoesdeforma !coordenadas de gaus para os pontos nodais
c$$$        call formaederi0 !Novidade (fun��es de forma e derivadas para qualquer ordem de elemento
c$$$        dpl=0.
c$$$        do k=1,notl(j)
c$$$	   dpl(1,1)=dpl(1,1)+dfi(k,1)*pl0(j,1,k)
c$$$	   dpl(2,1)=dpl(2,1)+dfi(k,1)*pl0(j,2,k)
c$$$        enddo
c$$$
c$$$c       Calculando a direcao ortogonal      
c$$$        vn2=dpl(1,1)
c$$$        vn1=-dpl(2,1)
c$$$        rmod=dsqrt(vn1**2+vn2**2)
c$$$ 	  vn1=vn1/rmod
c$$$	  vn2=vn2/rmod
c$$$
c$$$c       �ngulo
c$$$        alfrad=datan2(vn2,vn1) 
c$$$
c$$$c       Acerto do �ngulo para situa��es cr�ticas
c$$$        if (dabs(alfrad+pi).lt.1.0d-8) alfrad=pi
c$$$         pl0(j,3,kk)=alfrad !atribuindo posi��o angular
c$$$	 enddo
c$$$c      Outras situa��es cr�ticas
c$$$       do kk =2,notl(j)
c$$$	  if((pl0(j,3,kk)-pl0(j,3,1)).gt.pi) 
c$$$     #  pl0(j,3,kk)=pl0(j,3,kk)-2*pi
c$$$	  if((pl0(j,3,1)-pl0(j,3,kk)).gt.pi) 
c$$$     #  pl0(j,3,kk)=pl0(j,3,kk)+2*pi
c$$$       enddo
c$$$	enddo
c$$$
c$$$c     Novidade
c$$$c     Caso a posi��o no in�cio de um proceso (uauslamente din�mico)
c$$$c     seja fruto de um equil�brio de carregamento est�tico
c$$$c     permanente, este j� foi processado anteriormente e o arquivo
c$$$c     contendo estas informa��es � lido na posi��o atual
c$$$      if (iposic.eq.1) then 
c$$$       read (51,*)
c$$$       read (51,*)	  
c$$$c      coordenadas dos nos
c$$$c      no,x,y,z
c$$$       write(5,*) 'coordenadas dos nos diferente da inicial'
c$$$       write(5,*) 'no,x,y,z'
c$$$       do i=1,nnos
c$$$        read(51,*) k,p(glgl(k,1)),p(glgl(k,2)),p(glgl(k,3))
c$$$        write(5,10) k,p(glgl(k,1)),p(glgl(k,2)),p(glgl(k,3))
c$$$	 enddo
c$$$	 else
c$$$       p=p0  ! aqui a primeira tentativa de solu��o � a posi��o inicial 
c$$$      endif
c$$$       
c$$$ 10	format(1x,i7,3(1x,f27.10))
c$$$      return
c$$$      end
c$$$
c$$$

c      
c      *******************************************
c      Nesta subrotina a posi��o inicial nodal � atribuida �s vari�veis
c      locais, no elemento
c      calculamos os cossenos diretores iniciais para se tornarem as posi��es angulares
c      iniciais nos elementos (estas s�o diferentes para cada elemento, mesmo conectados)
c      *******************************************   
      subroutine geometriainicial
      USE BIG_STUFF	 
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
!     A posi��o angular para cada elemento � calculada
!     A posi��o de transla��o � a mesma do elemento global
      do j=1,nel
         do k=1,notl(j)
            nog=ic(j,k)
            do igl=1,2          !igl � o grau de liberdade x,y
               pl0(j,igl,k)=p0(glgl(nog,igl)) !posi��o inicial no elemento
            enddo
         enddo 
         
c     calcular a derivada em cada ponto do elemento
         if((itipoel(j).eq.2).or.(itipoel(j).eq.3))then
            numerointerno=3
         else
            numerointerno=2
         endif

         do is1=1,numerointerno             !notl(j)
            nog=ic(j,is1+1)
            xsi1=2.d0*is1-3.d0

            if(itipoel(j).eq.2)then
               nog=ic(j,is1)
               xsi1=1.d0*is1-2.d0
            endif
c     if(itipoel(j).eq.2)nog=ic(j,is1)
            if(itipoel(j).eq.3)then
               nog=ic(j,is1+1)
               xsi1=1.d0*is1-3.d0
            endif

           
            
            kk=is1



c     xsi1=(2.*is1-(1.+notl(j)))/(notl(j)-1.)
            
c     ***********************
c     variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c     call suportefuncoesdeforma !coordenadas de gaus para os pontos nodais
            call formaederi0    !Novidade (fun��es de forma e derivadas
                                !bspline nao uniformes cubicas
            
            dpl=0.d0
            do k=1,notl(j)
               dpl(1,1)=dpl(1,1)+dfi(k,1)*pl0(j,1,k) !dx_dxsi
               dpl(2,1)=dpl(2,1)+dfi(k,1)*pl0(j,2,k) !dy_dxsi
            enddo
            
c     Calculando a direcao ortogonal      
            vn2=dpl(1,1)
            vn1=-dpl(2,1)
            rmod=dsqrt(vn1**2+vn2**2)
            vn1=vn1/rmod
            vn2=vn2/rmod
            
c     �ngulo
            alfrad=datan2(vn2,vn1) 
c     calcula alfa global (preciso para as b-splines)
 
      if (dabs(alfrad+pi).lt.1.0d-8) alfrad=pi
           alfaglob(nog)=alfrad
         enddo

      enddo
         
         do jel=1,nel
            do kk=1,notl(jel)
               nog=ic(jel,kk)
               pl0(jel,3,kk)=alfaglob(nog)
            enddo
            
c     Outras situa��es cr�ticas
            j=jel
            do kk =2,notl(j)
               if((pl0(j,3,kk)-pl0(j,3,1)).gt.pi) 
     #pl0(j,3,kk)=pl0(j,3,kk)-2*pi
               if((pl0(j,3,1)-pl0(j,3,kk)).gt.pi) 
     #pl0(j,3,kk)=pl0(j,3,kk)+2*pi
            enddo

c            do kk=1,notl(j)
c               write(*,*)'j,kk,pl0',j,kk,pl0(j,3,kk)
c               read(*,*)
c            enddo
         enddo






         
c     precisa casos especiais para contorno n�o suave (refazer itipoel esq e dir)
      
c     Novidade
c     Caso a posi��o no in�cio de um proceso (uauslamente din�mico)
c     seja fruto de um equil�brio de carregamento est�tico
c     permanente, este j� foi processado anteriormente e o arquivo
c     contendo estas informa��es � lido na posi��o atual
      if(iposic.eq.1) then 
         read (51,*)
         read (51,*)	  
c     coordenadas dos nos
c     no,x,y,z
         write(5,*) 'coordenadas dos nos diferente da inicial'
         write(5,*) 'no,x,y,z'
         do i=1,nnos
            read(51,*) k,p(glgl(k,1)),p(glgl(k,2)),p(glgl(k,3))
            write(5,10) k,p(glgl(k,1)),p(glgl(k,2)),p(glgl(k,3))
         enddo
      else
         p=p0                   ! aqui a primeira tentativa de solu��o � a posi��o inicial 
      endif
      
 10   format(1x,i7,3(1x,f27.10))
      return
      end

c     *******************************************
c     Novidade
c     Nesta subrotina, ap�s uma solu��o de cada itera��o
c     as posi��es locais s�o atualizadas
c     nos elementos para c�culo de nova
c     for�a interna e matriz Hessiana
c      *******************************************   
      subroutine geometriaatual
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      do k=1,notl(j)
         nog=ic(j,k)
         do igl=1,2             !igl � o grau de liberdade x,y
            pl(j,igl,k)=p(glgl(nog,igl)) !posi��o atual no elemento
            rsl(j,igl,k)=rs(glgl(nog,igl)) !correcao dinamica atual no elemento
            qsl(j,igl,k)=qs(glgl(nog,igl)) !correcao dinamica no elemento
         enddo
      enddo
      
c     La�o para giro, devemos somar os giros nas posi��es angulares iniciais 
c     para calcularmos as deforma��es nos elementos
      do k=1,notl(j)
         nog=ic(j,k)
         pl(j,3,k)=pl0(j,3,k)+p(glgl(nog,iglv(nog,j)))
      enddo 
      
      return
      end
      
c     
c     *******************************************
c     Subrotina para calcular a matriz de massa para todos os elementos
c     a matriz de massa � uma grandeza lagrangeana e n�o muda
c     (conserva��o de masss) � calculada uma s�vez em todo o processo
c
c    n�o inclu�mos in�rcia de rota��o pois as barras s�o esbeltas
c     *******************************************   
      subroutine matrizmassa
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      do j=1,nel
         call massa             !c�lculo efetivo para cada elemento
      enddo

      return
      end
c     
c     *******************************************
c     Subrotina para calcular a matriz hessiana
c     For�as internas inercias de (algoritmo de Newmark)
c     For�as internas elastopl�sticas
c     *******************************************   
      subroutine matriz
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      rkl=0.d0                    !hessiana local
      vtl=0.d0                    !vetor forca int local
      vv=0.d0                     ! vetor independente
      vint=0.d0                   !vetor de forca interna
      
c     nel � o numero de elementos
      do j=1,nel
         call elemento
      enddo
      
c     montagen do vetor Hessiano (novidade)
c     de rigidez a(valores diferentes de zero) ma27
      kkk=0
      
      do ib=1,nnos
         h=0.d0
         do jj=1,kn(ib)         !jj vai at� o n�mero de elementos que o n� participa
            j=im(ib,jj)         ! j � o numero de cada elemento (incid�ncia inversa)


            call montalinha(ib,jj) !montando linhas do n� ib
         enddo


         
         do i=1,ngl(ib)
            ii=(ib-1)*3+i
            ii=glgl(ib,i)
            do jj=ii,n
               if (dabs(h(i,jj)).gt.1.d-16) then !no vetor-matriz s� entra valores n�o nulos
                  kkk=kkk+1
                  irn(kkk)=ii
                  icn(kkk)=jj
                  a(kkk)=h(i,jj)
               else
               endif
            enddo
         enddo
      enddo
      
      nz=kkk                    !controle de tamanho do vetor matriz para ma27
      n27=n                     !n�mero de graus de liberdade informado ao ma27
      
      return
      end

c      
c     *******************************************
c     subrotina que cria a matriz Hessiana (rigidez)
c     e for�as internas para cada elemento finito
c     *******************************************   
      subroutine elemento
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      call leielast             !defini��o da lei constitutiva el�stica para cada elemento
      
      fri=0.d0                    !for�a residual inercial
      frib=0.d0                   !for�a residual de amortecimento
      ue=0.d0                     !energia de deforma��o (curiosidade)

      ng=notl(j) !notl(j)
c     atualizar configura��o nos elementos
      call geometriaatual 

c      residuo inercial
      do k=1,notl(j)
         iloc=(k-1)*3
         do igl=1,3
            i1=iloc+igl
            frib(i1)=ram(j)*rsl(j,igl,k)+(1.+rgn*dt*ram(j))* (pl(j,igl
     $           ,k)/rbn/(dt*dt)-qsl(j,igl,k))
         enddo
      enddo
      
      do i1=1,3*notl(j)
         do j1=1,3*notl(j)
            fri(i1)=fri(i1)+rmg(i1,j1,j)*frib(j1)
         enddo
      enddo
c     Fim do res�duo inercial

c     Aquisi��o das coordenadas e pesos de Gauss
      call gauss(ng,qsi,w)
      call gauss(ng3,qsi3,w3)
c     call dgqrul(ng,1,0.d0,0.d0,0,qsi,qsi,w)
c     call dgqrul(ng3,1,0.d0,0.d0,0,qsi3,qsi3,w3)
      
      do ig1=1,ng               !La�o nos pontos de Gauss ao longo do elemento
         xsi1=qsi(ig1)
c     **********************************************************************
c     Novidade
c     variaveis necess�rias para fun��o de forma e derivadas de ordem qualquer
c         call suportefuncoesdeforma

         call formaederi0      
         
         do igf=1,3             !La�o para fibras da se��o transversal (novidade)
            do ig3=1,ng3        !La�o para Pontos de Gauss em cada fibra
               up=0.d0            !Dire��o local da fibra na configura��o inicial (ao longo do elemento)
               vp=0.d0            !Dire��o local da fibra na configura��o inicial (transversal ao elemento)
               wp=0.d0            !Casca
c     Para a terceira dire��o a fun��o de forma � o pr�prio qsi
               xsi3=qsi3(ig3)   !Ponto de Gauss
               
c     Novidade (conceito posicional)
c     Chamada da subrotina que calcula o gradiente da 
c     Mudan�a de configura��o fict�cia, do espa�o adimensional de Gauss
c     Para a configura��o inicial do corpo
               call calc_A0     !matriz gradiente do adimensional para inicial      
               up(1)=dcos(salfa-pi/2.d0)
               up(2)=dsin(salfa-pi/2.d0)
               vp(1)=dcos(salfa)
               vp(2)=dsin(salfa)
               
c     Novidade (conceito posicional)
c     Chamada da subrotina que calcula o gradiente da 
c     Mudan�a de configura��o fict�cia, do espa�o adimensional de Gauss
c     Para a configura��o atual do corpo
               call calc_Af     !gradiente da transformacao total (dps)
               call cauchy      !Tensor de alongamento de Cauchy Green
               acp=0.d0
               aval=0.d0
               
c     Op��o de trabalho em dire��o principal do alongamento de cauchy-green
c     call devcsf(2,ch,2,aval,acp,2)
c     up(1)=acp(1,1)
c     up(2)=acp(2,1)
c     vp(1)=acp(1,2)
c     vp(2)=acp(2,2)

c        up(1)=1. !op��o de coordenadas globais
c	   vp(2)=1.

c        Novidade
c        A partir do tensor de alongamento de cauchy Green calculam-se
c        as deforma��es de engenharia
c        Os softwreas tradicionais n�o usam medidas de engenharia
         call deformacoes !deformacoes longitudinais, distorcoes e alongamentos
         call energia !c�lculo da energia de deforma��o, opcional
         drex=0.d0
         drey=0.d0
         drez=0.d0
         drexy=0.d0
         drexz=0.d0
         dreyz=0.d0
         p1=0.d0
         p2=0.d0
         p3=0.d0
         daf=0.d0
         dch=0.d0
         
         do li =1,notl(j) 
            call derivadps      !derivadas dos dps 
c     Derivada do gradiente da mudan�a de configura��o
c     fict�cia (Gauss-atual) A1 em rela��o aos graus de liberdade
            call dcauchy        !derivada do tensor decauchy para os gl do n� li
            call deridef        !derivadas das deforma��es em rela��o ao gl local para
            call derivener      !vetor de forcas internas ou derivada da energia de
c     deforma��o em rela��o aos par�metros nodais gls
         enddo
         
         do li=1,notl(j)
            do lj=1,notl(j)         
               pp1=0.d0
               pp2=0.d0
               pp3=0.d0
               d2rex=0.d0
               d2rey=0.d0
               d2rexy=0.d0
               d2ch=0.d0
               call d2cauchy    !Segunda derivada do tensor de Cauchy-Green em rela��o aos gls
               call deri2def    !Segunda derivada da deforma��o em rela��o aos gls
c     Determina��o da matriz hessiana
c     ou segunda derivada da energia de deforma��o em rela��o aos par�metros nodais, gls
               call segderene
            enddo
         enddo
      enddo                     !ng1
      enddo                     !igf=1,3 s�o as tr�s faixas
      enddo                     !ng
      
      if ((rgn.ne.1.).or.(rgn.ne.1.).or.(rgn.ne.1.)) then !din�mico
         do i1=1,3*notl(j)
            vtl(i1,j)=vtl(i1,j)+fri(i1)
            do j1=1,3*notl(j)
               rkl(i1,j1,j)=rkl(i1,j1,j)+rmg(i1,j1,j)*
     #(1.d0+rgn*dt*ram(j))/rbn/(dt*dt)
            enddo
         enddo
      else
c     est�tico
      endif
      
      return
      end

c     *******************************************
c     Subrotina que monta a matriz de rigidez global
c     monta-se apenas o n�mero de linhas referente aos graus
c     de liberdade de cada n� utilizando a incid�ncia inversa
c     pega-se este conjunto e leva-se os valores n�o nulos
c     para a Matriz-Vetor que ser� utilizada no solver MA27
c     Novidade (sistema de montagem muito econ�mico e r�pido, 
c     al�m de simples, tem�po de processamento similar
c     ao dos comerciais internacionais)
c     *******************************************   
      subroutine montalinha(ib,j1)
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'



      nol=ib                    !no global
c     j1 � o jj fora da  subrotina, n�mero local do elemento na incid�ncia inversa
      do ie=1,notl(j)           !n�mero de n� do elemento j
         if (ic(j,ie).eq.nol) then !identifica n� local do elemento a ser considerado
            i=ie
            goto 100
         else 
         endif
      enddo
 100  continue
      
c     o i � o n� local do elemnto j global, relacionado ao no global nol e elemento loca jj
      nol=ic(j,i)
      nl0=(nol-1)*3
      
c     montagen do vetor de forcas internas
      do k=1,3
         ill=(i-1)*3+k
         ilg=glgl(ib,k)         !ib e nol sao iguais
         if(k.eq.3) ilg=glgl(ib,iglv(ib,j)) !r�tulas
         vint(ilg)=vint(ilg)+vtl(ill,j) !vtl forca interna
      enddo
      
c     montagem da hessiana
      do k=1,3
         ill=(i-1)*3+k
         ilg=glgl(ib,k)
         if(k.eq.3) ilg=glgl(ib,iglv(ib,j)) !r�tulas
         if(ko(ilg).eq.0) then
            do jj=1,notl(j)
               noj=ic(j,jj)
               nc0=(noj-1)*3
               do kj=1,3
                  icl=(jj-1)*3+kj
c	    icg=nc0+kj
                  icg=glgl(noj,kj)
                  if(kj.eq.3) icg=glgl(noj,iglv(noj,j)) !r�tulas
                  if (ko(icg).eq.0) then
                     kd=k
                     if(k.eq.3) kd=iglv(ib,j) !r�tula
                     h(kd,icg)=rkl(ill,icl,j)+h(kd,icg) !k � sempre de 1 ate o numero de gl(3)
                  else
                  endif
               enddo
            enddo
c     o conjunto de linhas est� montado
         else 
         endif
      enddo
      
      return 
      end
	

c    ************************************************
c    subrotina que calcula a matriz gradiente da
c    mudan�a de configura��o fict�cia do epa�o adimensional de
c    Gauss para a configura��o inicial (a0)
c    ****************************************************
      Subroutine calc_A0
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	include 'bloclagatbig.for'

c     a1=resp(j)/2. !altura de uma se��o inteira (abolido)
      a1=hb(j,igf)/2.d0 !altura de uma fibra

      alfa=0.d0 !�ngulo em um ponto de gauss
      dalfa=0.d0 !derivada do mesmo em rela��o � coordenada xsi
c     eta � o xsi3       
      pxsi=0.d0                   !derivada da posi��o em rela��o � xsi (longitudinal)
      peta=0.d0                   !derivada da posi��o em rela��o � eta (transversal)
      
      do il=1,notl(j)
         alfa=alfa+fi(il)*pl0(j,3,il)
         dalfa=dalfa+dfi(il,1)*pl0(j,3,il)
      enddo
      
      salfa=alfa                !�ngulo que � usado para calcular a dire��o local desejada
      do il=1,notl(j)
         pxsi(1)=pxsi(1)+dfi(il,1)*pl0(j,1,il) !(1) posi��o x (na linha refer�ncia)
         pxsi(2)=pxsi(2)+dfi(il,1)*pl0(j,2,il) !(2) posi��o y (na linha refer�ncia)
      enddo
      
      pxsi(1)=pxsi(1)-dsin(alfa)*dalfa*(a1*xsi3+fd(j,igf)) !na posi��o efetiva
      pxsi(2)=pxsi(2)+dcos(alfa)*dalfa*(a1*xsi3+fd(j,igf)) !na posi��o efetiva
      
      peta(1)=dcos(alfa)*a1      !derivada em rela��o a xsi3
      peta(2)=dsin(alfa)*a1      !derivada em rela��o a xsi3
      
c     matriz gradiente para configura��o inicial
      a0(1,1)=pxsi(1)
      a0(1,2)=peta(1)
      a0(2,1)=pxsi(2)
      a0(2,2)=peta(2)
      
      rjac0=a0(1,1)*a0(2,2)-a0(1,2)*a0(2,1) !jacobiano para integra��o
      rjac1=dsqrt(pxsi(1)*pxsi(1)+pxsi(2)*pxsi(2)) !para usar no c�culo da forca distribuida
      
      rjac0=rjac0*brd(j,igf)    !no volume total
      rjac2=a1*brd(j,igf)       !para usar no c�lculo dos esfor�os solicitante (na �rea)
      return
      end
      

c     ************************************************
c     subrotina que calcula a matriz gradiente da
c     mudan�a de configura��o fict�cia do epa�o adimensional de
c     Gauss para a configura��o final (af)
c     ****************************************************
      Subroutine calc_Af
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     a1=resp(j)/2.
      a1=hb(j,igf)/2.d0           !considera-se a espessura constante
      
c     repete-se o c�lculo de a0 trocando pl0 por pl
      alfa=0.d0
      dalfa=0.d0
      
      pxsi=0.d0
      peta=0.d0
      
c     eta � o xsi3       
      
      do il=1,notl(j)
         alfa=alfa+fi(il)*pl(j,3,il)
         dalfa=dalfa+dfi(il,1)*pl(j,3,il)
      enddo
      
      do il=1,notl(j)
         pxsi(1)=pxsi(1)+dfi(il,1)*pl(j,1,il)
         pxsi(2)=pxsi(2)+dfi(il,1)*pl(j,2,il)
      enddo
      
      pxsi(1)=pxsi(1)-dsin(alfa)*dalfa*(a1*xsi3+fd(j,igf))
      pxsi(2)=pxsi(2)+dcos(alfa)*dalfa*(a1*xsi3+fd(j,igf))
      
      peta(1)=dcos(alfa)*a1
      peta(2)=dsin(alfa)*a1
      
c     Matriz gradiente da mudan�a de configura��o (gauss - final)
      af(1,1)=pxsi(1)
      af(1,2)=peta(1)
      af(2,1)=pxsi(2)
      af(2,2)=peta(2)
      return
      end
      

c     ************************************************
c     Subrotina que calcula o tensor de cauchy e
c     o gradientes da mudan�a de confihura��o 
c     da configura��o inicial para final
c     ****************************************************
      Subroutine cauchy
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	include 'bloclagatbig.for'

      call inverse(2,a0,a0inv)  ! inversao da matriz a0

c     Produto de af pela inversa de a0
      aaux=0.d0
      do ii=1,2
         do jj=1,2
            do kk=1,2
               aaux(ii,jj)=aaux(ii,jj)+af(ii,kk)*a0inv(kk,jj)
            enddo
         enddo
      enddo
      
c     produto da transposta do resultado anterior pelo resultado anterior
c     alongamento de Cauchy green
      ch=0.d0
      do ii=1,2
         do jj=1,2
            do kk=1,2
               ch(ii,jj)=ch(ii,jj)+aaux(kk,ii)*aaux(kk,jj)
            enddo
         enddo
      enddo
      
      return
      end

c     ************************************************
c     Subrotina que calcula deforma��ess, distor��es e alongamentos
c     Faz uso do tensor alongamento de cauchy-green (ch)
c     ****************************************************
      Subroutine deformacoes
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      aup=0.d0
      rlex=0.d0                   !alongamento na dire��o up 
      do ii=1,2
         do jj=1,2
            aup(ii)=aup(ii)+ch(ii,jj)*up(jj)
         enddo
         rlex=rlex+aup(ii)*up(ii)
      enddo
      rlex=dsqrt(rlex)
      
      aup=0.d0
      rley=0.d0                   !alongamento na dire��o vp
      do ii=1,2
         do jj=1,2
            aup(ii)=aup(ii)+ch(ii,jj)*vp(jj)
         enddo
         rley=rley+aup(ii)*vp(ii)
      enddo
      rley=dsqrt(rley)
      
      goto 10                   !2d
c     casca retirado
      aup=0.d0
      rlez=0.d0
      do ii=1,2
         do jj=1,2
            aup(ii)=aup(ii)+ch(ii,jj)*wp(jj)
         enddo
         rlez=rlez+aup(ii)*wp(ii)
      enddo
      rlez=dsqrt(rlez)
 10   continue                  !fim casca
      
c     Deforma��es de engenharia
      rex=rlex-1.d0
      rey=rley-1.d0
c     rez=rlez-1. !casca
      rez=0.d0
     
      aup=0.d0
      rlexy=0.d0 !'alongamento distorcional'
      do ii=1,2
         do jj=1,2
            aup(ii)=aup(ii)+ch(ii,jj)*vp(jj)
         enddo
         rlexy=rlexy+aup(ii)*up(ii)
      enddo
      rexy=(pi/2.d0-dacos(rlexy/rlex/rley))/2.d0 !deforma��o de distor��o
      rexy=rexy*2.d0              !distor��o de engenharia

c     Casca
      goto 20                   !2d
      aup=0.d0
      rlexz=0.d0
      do ii=1,2
         do jj=1,2
            aup(ii)=aup(ii)+ch(ii,jj)*wp(jj)
         enddo
         rlexz=rlexz+aup(ii)*up(ii)
      enddo
      rexz=(pi/2.d0-dacos(rlexz/rlex/rlez))/2.d0
      rexz=rexz*2.d0 !gama

      aup=0.d0
      rleyz=0.d0
      do ii=1,2
         do jj=1,2
            aup(ii)=aup(ii)+ch(ii,jj)*wp(jj)
         enddo
         rleyz=rleyz+aup(ii)*vp(ii)
      enddo
      reyz=(pi/2.d0-dacos(rleyz/rley/rlez))/2.d0
      reyz=reyz*2.d0              !gama

 20   continue                  !fim casca
      rexz=0.d0                   !2d
      reyz=0.d0

      return
      end

c     ************************************************
c     Calculou-se af que foi utilizada para calcular ch
c     Ainda teremos que calcular a derivada de ch em rela��o 
c     Aos par�metros inc�gnitos: posi��es nodais atuais
c     Para tanto devemos derivar as derivadas das posi��es
c     em rela��o �s coordenadas adimensionais
c     agora em rela��o �s posi��es nodais
c     ****************************************************
      Subroutine derivadps
      USE BIG_STUFF	 
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     Na express�o pi(jj,kk,li) entende-se
c     i=1 se falamos da p1 (posi��o x) (mapa real)
c     i=2 de falamos da p2 (posi��o y) (mapa real)

c     jj=1 para xsi longitudinal
c     jj=2 para eta (xsi3 transversal)

c     kk=1 derivando no gl1 (posi��o nodal x)
c     kk=2 derivando no gl2 (posi��o noda y)
c     kk=3 derivando no gl3 (posi��onodal teta)

c     li para o n� a que se refere a derivada que se calcula
 

c      a1=resp(j)/2.
      a1=hb(j,igf)/2.d0

c     larg=ra(j)
      larg=brd(j,igf)
      alfa=0.d0
      dalfa=0.d0

      pxsi=0.d0
      peta=0.d0
c     eta � o xsi3       

      do il=1,notl(j)
         alfa=alfa+fi(il)*pl(j,3,il)
         dalfa=dalfa+dfi(il,1)*pl(j,3,il)
      enddo

      p1(1,1,li)=dfi(li,1) 
      p1(2,1,li)=0.d0

c     zero para p1 eta e delta derivados em p2 e p3

      p1(1,3,li)=(a1*xsi3+fd(j,igf))*
     #(-dcos(alfa)*fi(li)*dalfa-dsin(alfa)*dfi(li,1))
      p1(2,3,li)=-dsin(alfa)*a1*fi(li)


      p2(1,2,li)=dfi(li,1) 
      p2(2,2,li)=0

c     zero para p1 eta e delta derivados em p2 e p3
      
      p2(1,3,li)=(a1*xsi3+fd(j,igf))*
     #(-dsin(alfa)*fi(li)*dalfa+dcos(alfa)*dfi(li,1))
      p2(2,3,li)=dcos(alfa)*a1*fi(li)
      
      return
      end

c     ************************************************
c     Subrotina que calcula derivadas do tensor de 
c     alongamento de Cauchy-green para cada grau de liberdade
c     de cada n� do elemento finito
c     ****************************************************
      Subroutine dcauchy
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	include 'bloclagatbig.for'

c     mudando a nota��o da subrotina anterior para tensor derivada de af
      do k2=1,3
         daf(1,1,li,k2)=p1(1,k2,li) !k2 � o gl refernte ao n� li
         daf(1,2,li,k2)=p1(2,k2,li)
         daf(2,1,li,k2)=p2(1,k2,li)
         daf(2,2,li,k2)=p2(2,k2,li)
         aaux=0.d0
         do ii=1,2
            do jj=1,2
	   do kk=1,2
        aaux(ii,jj)=aaux(ii,jj)+af(ii,kk)*a0inv(kk,jj)
      enddo
      enddo
      enddo
      
      aaux1=0.
       do ii=1,2
	  do jj=1,2
	   do kk=1,2
	    aaux1(ii,jj)=aaux1(ii,jj)+daf(ii,kk,li,k2)*a0inv(kk,jj)
         enddo
	  enddo
	 enddo

       do ii=1,2
	  do jj=1,2
	   dch(ii,jj,k2,li)=0.d0  !zerando a derivada do alongamento de Cauxhy-Grenn
        enddo
	 enddo

       do ii=1,2
	  do jj=1,2
	   do kk=1,2
	    dch(ii,jj,k2,li)=dch(ii,jj,k2,li)+aaux1(kk,ii)*aaux(kk,jj)+
     #    aaux(kk,ii)*aaux1(kk,jj)
         enddo
	  enddo
	 enddo
      enddo

      return
	end

c     ************************************************
c     Subrotina que calcula derivadas das deforma��es
c     em fun��o das posi��es nodais (gls) locais
c     para cada no li=> d�o origem �s for�as nodais equivalentes
c     ****************************************************
      Subroutine deridef
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

c     precisamos registrar qual � o n� do elemento que estamos considerando
c     para usar na segunda derivada (Hessiana)
      
      do k1=1,3
c     derivada da deformacao ex
         aup=0.d0
         do i1=1,2
            do j1=1,2
               aup(i1)=aup(i1)+dch(i1,j1,k1,li)*up(j1)
            enddo
         enddo
         do i1=1,2
            drex(k1,li)=drex(k1,li)+aup(i1)*up(i1)/rlex/2.d0
         enddo
         
c     derivada da deformacao ey
         aup=0.
         do i1=1,2
            do j1=1,2
               aup(i1)=aup(i1)+dch(i1,j1,k1,li)*vp(j1)
            enddo
         enddo
         do i1=1,2
            drey(k1,li)=drey(k1,li)+aup(i1)*vp(i1)/rley/2.d0
         enddo
         
c     casca
         goto 10                ! 2d
c     derivada da deformacao ez
         aup=0.
         do i1=1,2
            do j1=1,2
               aup(i1)=aup(i1)+dch(i1,j1,k1,li)*wp(j1)
            enddo
         enddo
         do i1=1,2
            drez(k1,li)=drez(k1,li)+aup(i1)*wp(i1)/rlez/2.d0
         enddo
 10      continue
         
c      derivada da deformacao exy
         aup=0.d0
         aup1=0.d0
         do i1=1,2
            do j1=1,2
               aup(i1)=aup(i1)+dch(i1,j1,k1,li)*vp(j1)
               aup1(i1)=aup1(i1)+ch(i1,j1)*vp(j1)
            enddo
         enddo
         aur=0.d0
         aur1=0.d0
         do i1=1,2
            aur=aur+aup(i1)*up(i1)
            aur1=aur1+aup1(i1)*up(i1)
         enddo
         aurxy(k1,li)=aur
         rlexy=aur1

         drexy(k1,li)=(aur*rlex*rley-aur1*(rley*drex(k1,li) +rlex
     $        *drey(k1,li))) /(2.d0*rlex*rley*dsqrt((rlex*rley)*(rlex
     $        *rley)-aur1*aur1))

         drexy(k1,li)=drexy(k1,li)*2.d0 !derivada da distor��o

c     casca
         goto 20                !2d
c     derivada da deformacao exz
         aup=0.d0
         aup1=0.d0

         do i1=1,2
            do j1=1,2
               aup(i1)=aup(i1)+dch(i1,j1,k1,li)*wp(j1)
               aup1(i1)=aup1(i1)+ch(i1,j1)*wp(j1)
            enddo
         enddo

         aur=0.d0
         aur1=0.d0
         do i1=1,3
            aur=aur+aup(i1)*up(i1)
            aur1=aur1+aup1(i1)*up(i1)
         enddo

         aurxz(k1,li)=aur
         rlexz=aur1

         drexz(k1,li)=(aur*rlex*rlez-aur1*(rlez*drex(k1,li) +rlex
     $        *drez(k1,li))) /(2.d0*rlex*rlez*dsqrt((rlex*rlez)*(rlex
     $        *rlez)-aur1 *aur1))

         drexz(k1,li)=drexz(k1,li)*2.d0 !gama

c      derivada da deformacao eyz
         aup=0.d0
         aup1=0.d0

         do i1=1,2
            do j1=1,2
               aup(i1)=aup(i1)+dch(i1,j1,k1,li)*wp(j1)
               aup1(i1)=aup1(i1)+ch(i1,j1)*wp(j1)
            enddo
         enddo

         aur=0.d0
         aur1=0.d0
         
         do i1=1,2
            aur=aur+aup(i1)*vp(i1)
            aur1=aur1+aup1(i1)*vp(i1)
         enddo

         auryz(k1,li)=aur
         rleyz=aur1
         
         dreyz(k1,li)=(aur*rley*rlez-aur1*(rlez*drey(k1,li)+rley*drez(k1
     $        ,li)))/(2.d0*rley*rlez*dsqrt((rley*rlez)*(rley*rlez)-aur1
     $        *aur1))

         dreyz(k1,li)=dreyz(k1,li)*2.d0 !gama
 20      continue               !fim casca
         
      enddo
      return
      end
c     ************************************************
c     Subrotina opcional calcula a energia de deform�o
c     ****************************************************
      Subroutine energia
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
            rrk=re(j)/(1.d0+ri(j))/(1.d0-2.d0*ri(j))
      rri=ri(j)
      rrg=rg(j)
      rue=rrk/2.d0*((1.d0-rri)*(rex*rex+rey*rey)+ 2.d0*rri*(rex*rey))+
     $     rrg*(rexy*rexy) 
      return
      end

c     ************************************************
c     Subrotina que calcula as tens�es de engenharia
c     Chama subrotina de retorno (plasticidade)
c     e calcula as for�as internas
c     ****************************************************
      Subroutine derivener
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
c     c�lculo do valor da temperatura para o ponto de integra��o
      tpos=(a1*xsi3+fd(j,igf))*(temc(j)-temb(j))/2.d0+ (temc(j)+temb(j))
     $     /2.d0
      
c     CHAMADA DA SUBROTINA DE RETORNO.
      et(1)=rex-tpos*calt(j)    !deforma��es totais
      et(2)=rey-tpos*0.d0
      et(3)=rez                 !casca=0 aqui
      et(4)=rexy
      et(5)=rexz                !casca=0 aqui
      et(6)=reyz                !casca=0 aqui
      
c     chamada do modelo constitutivo elastopl�stico
      if (ia.ge.1) call retorno !precau��o para problemas controle delsocamento
c     call retorno
c     vamos corrigir o valor das tens�es calculadas
      
      do iep=1,6
         ee(iep)=et(iep)-ep(j,ig1,ig3,igf,iep) !deforma��o pl�stica acumulada
      enddo
      s=0.d0
      
      do i=1,6
         do jr=1,6
            s(i)=s(i)+ceret(i,jr)*ee(jr) !lei constitutiva elastopl�stica
         enddo
      enddo
      
      duex=s(1)                 !tens�es armazenadas em vari�veis que ser�o utilizadas (melhorar escrita para pr�xima vers�o)
      duey=s(2)
      duez=s(3)
      duexy=s(4)
      duexz=s(5)
      dueyz=s(6)
      
c     due � a derivada da energia em rela��o a cada gl
c     de um determindo n� li de um elemento j
c     para um ponto de Gauss
      
      iloc=(li-1)*3             !pequena incid�ncia local
      
      do k1=1,3
         due(k1)=duex*drex(k1,li)+duey*drey(k1,li)+ (duexy*drexy(k1,li)) !gama
             
c     aqui vamos fazer uma pequena montagen do vetor de for�as
c     temos o due(3) que devera se transformar no vtl(3*notl(j))
         i1=iloc+k1
c     vtl(i1,j)=vtl(i1,j)+due(k1)*w(ig1)*w(ig2)*w3(ig3)*rjac0 !casca
         vtl(i1,j)=vtl(i1,j)+due(k1)*w(ig1)*rjac0*w3(ig3)
      enddo
      
      return
      end


c     ************************************************
c     Subrotina que calcula a matriz hessiana local
c     � uma matriz semi tangente fisicamente pois se usa 
c     a tes�o real nos termos nlg e el�sticas nos termos
c     geom�tricamente lineares (� mais segura doque a tangente total)
c     pois evita singularidades, apesar de taxa de convergencia ligeiramente inferior
c     ****************************************************
      Subroutine segderene
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	include 'bloclagatbig.for'

      d2ue=0.d0
      s=0.d0 !recuperando as tens�es
 	do i=1,6
	 do jr=1,6
	  s(i)=s(i)+ceret(i,jr)*ee(jr)
       enddo
	enddo

      duex=s(1) !escrita a ser atualizada
      duey=s(2)
      duez=s(3)
      duexy=s(4)
      duexz=s(5)
      dueyz=s(6)

      iloc=(li-1)*3 !incid�ncia de montagem local
      jloc=(lj-1)*3

      do k1=1,3
       i1=iloc+k1
       do k2=1,3
        j1=jloc+k2
        d2ue(k1,k2)=
     #  rhook(j,1,1)*drex(k1,li)*drex(k2,lj)+
     #  rhook(j,2,2)*drey(k1,li)*drey(k2,lj)

        d2ue(k1,k2)=d2ue(k1,k2)+
     #  (rhook(j,3,3)*drexy(k1,li)*drexy(k2,lj))*1. !?

        d2ue(k1,k2)=d2ue(k1,k2)+
     #  rhook(j,1,2)*drex(k1,li)*drey(k2,lj)+
     #  rhook(j,2,1)*drex(k2,lj)*drey(k1,li)

        d2ue(k1,k2)=d2ue(k1,k2)+duex*d2rex(k1,k2)+duey*d2rey(k1,k2)+
     #  duexy*d2rexy(k1,k2)

        rkl(i1,j1,j)=rkl(i1,j1,j)+d2ue(k1,k2)*w(ig1)*w3(ig3)*rjac0
       enddo

      enddo
      return
      end

c     *****************************************
c     Subrotina que aplica for�as internas e externas
c     respeitando as condi��es de contorno
c     ***************************************
      subroutine condcon
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      v=0.
      do i=1,n
         if (ko(i).eq.1) then   !ser� usado no impacto entre duas estruturas (futura vers�o)
         else
            v(i)=f(i)-vint(i)  
         endif
      enddo
      
      return 
      end
      
c     *****************************************
c     Subrotina que destroca vetores com bisseccao !desativada
c     ***************************************
      subroutine destroca
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      rrrf=1.d0
      rrr1=0.d0
      rrr2=0.d0

      do j=1,n
	 if(ko(j).ne.1) rrr1=rrr1+v(j)**2
	 if(ko(j).ne.1) rrr2=rrr2+p0(j)**2
      enddo

      rrnorma=dsqrt(rrr1/rrr2)
      rrrfat=rrnorma/tol
      
	if (rrrfat.gt.1.) rrrf=rrrfat

	do j=1,n
       if (ko(j).eq.1) then
        f(j)=v(j)/rrrf !V � O DELTA P
        v(j)=0.d0
        else
	  p(j)=p(j)+v(j)/rrrf
       endif
      enddo

      if (ia.gt.20) ria=2.d0
      if (ia.gt.30) ria=3.d0
      if (ia.gt.40) ria=2.d0
      if (ia.gt.50) ria=3.d0
      if (ia.gt.60) ria=2.d0

      if(ia.gt.20) then
       if(int(ia/ria+.5d0).eq.(ia/ria+.5d0)) p=(p+p0i)/2.d0
	 p=(p+p0i)/2.d0
	 else
	endif
      p0i=p
      v=0.d0

      call box
      return 
      end

c     *****************************************
c     Subrotina que ap�s a solu��o do sistem
c     Atualiza a posi��o dos pontos livres
c     Garda trajet�ria para impacto
c     Chama rotina de impacto
c     ***************************************
      subroutine destroca1
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      p0i=p !guardando valor anterior para c�lculo de impacto

	do j=1,n
       if ((ko(j).eq.1)) then
        f(j)=vint(j)
        reacao(j)=vint(j)
        else
        p(j)=p(j)+v(j)
       endif
c      if ((ko1(j).eq.1)) vint(j)=0. !lembrete para aplica��o da for�a interna no corpo alvo
      enddo

      ko1=0
      call box !impacto com atrito
      v=0.d0
      return 
      end


c     *****************************************
c     Impacto com atrito
c     por enquanto trata regi�o retangular
c     para limite de an�lise
c     ***************************************
      subroutine box
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      xb0=xp
	xb1=xc
      yb0=yp
	yb1=yc

	do j=1,nnos

c     marcando o n� impactante      
      if (p(glgl(j,1)).ge.xb1) then
	ko1(glgl(j,1))=1
	ko1(glgl(j,2))=1
      endif

      if (p(glgl(j,1)).le.xb0) then
	ko1(glgl(j,1))=1
	ko1(glgl(j,2))=1
      endif

      if (p(glgl(j,2)).ge.yb1) then
	ko1(glgl(j,1))=1
	ko1(glgl(j,2))=1
	endif

      if (p(glgl(j,2)).le.yb0) then
	ko1(glgl(j,1))=1
	ko1(glgl(j,2))=1
      endif

c     Novidade (procedimento geom�trico de impacto)
c     ao se corrigir a posi��o do n� impactante
c     for�a-se o c�lculo de um for�a interna
c     que ao final das intera��es
c     resulta na for�a externa aplicada no n� impactante
c     ou seja, se iguala � for�a de contato nodal equivalente

      if (p(glgl(j,1)).ge.xb1) then
	if (dabs(p(glgl(j,1))-p0i(glgl(j,1))).gt.1.d-13) then
	a2=(p(glgl(j,2))-p0i(glgl(j,2)))/(p0i(glgl(j,1))-p(glgl(j,1)))
	b2=p0i(glgl(j,2))-a2*p0i(glgl(j,1))
	ym=a2*xb1+b2
	p(glgl(j,2))=p(glgl(j,2))*(1.-rimp)+ym*rimp
	p(glgl(j,1))=xb1
	else
	p(glgl(j,1))=xb1
	endif
	endif

      if (p(glgl(j,1)).le.xb0) then
	if (dabs(p(glgl(j,1))-p0i(glgl(j,1))).gt.1.d-13) then
	a2=(p(glgl(j,2))-p0i(glgl(j,2)))/(p0i(glgl(j,1))-p(glgl(j,1)))
	b2=(p0i(glgl(j,2))-a2*p0i(glgl(j,1)))
	ym=a2*xb0+b2
	p(glgl(j,2))=p(glgl(j,2))*(1.-rimp)+ym*rimp
	p(glgl(j,1))=xb0
	else
	p(glgl(j,1))=xb0
	endif
      endif

      if (p(glgl(j,2)).ge.yb1) then
      if (dabs(p(glgl(j,2))-p0i(glgl(j,2))).gt.1.d-13) then
	a2=(p(glgl(j,2))-p0i(glgl(j,2)))/(p0i(glgl(j,1))-p(glgl(j,1)))
	b2=p0i(glgl(j,2))-a2*p0i(glgl(j,1))
      xm=(yb1-b2)/a2
      p(glgl(j,1))=p(glgl(j,1))*(1.-rimp)+xm*rimp
	p(glgl(j,2))=yb1
	else
	p(glgl(j,2))=yb1
	endif
	endif

      if (p(glgl(j,2)).le.yb0) then
      if (dabs(p(glgl(j,2))-p0i(glgl(j,2))).gt.1.d-13) then
	a2=(p(glgl(j,2))-p0i(glgl(j,2)))/(p0i(glgl(j,1))-p(glgl(j,1)))
	b2=p0i(glgl(j,2))-a2*p0i(glgl(j,1))
      xm=(yb0-b2)/a2

      p(glgl(j,1))=p(glgl(j,1))*(1.-rimp)+xm*rimp
	p(glgl(j,2))=yb0
	else
	p(glgl(j,2))=yb0
	endif
	endif

      enddo

      return 
      end



c     *****************************************
c     Subrotina que calcula norma
c     crit�rio de parada
c     ***************************************
      subroutine calculanorma
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      r1=0.
      r2=0.
      do j=1,n
	 if((ko(j).ne.1).and.(ko1(j).ne.1)) r1=r1+v(j)**2
	 if(ko(j).ne.1) r2=r2+p0(j)**2
      enddo
      write(*,*) 'erro',dsqrt(r1/r2)
      rnorma=dsqrt(dabs(r1/r2))
      return 
      end


c     *******************************************
c     Saida de resultados gerais
c     ******************************************* 
      subroutine saida  
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     # .lt.1.0d-6)) then

      write(5,*)
      write(5,*)'posicao dos nos, passo, numero de iteracoes',ipt,ia
      write(5,*) 'no,x,y'
c     posi��o dos nos
c     no,x,y
      do i=1,nnos
       write(5,10) i,(p(glgl(i,kk)),kk=1,ngl(i))
      enddo


      write(5,*)
      write(5,*)'For�as reativas',ipt,ia
      write(5,*) 'no,fx,fy,mz'
c     posi��o dos nos
c     no,x,y
      do i=1,nnos
      km=0
	do kk=1,ngl(i)
	if(ko(glgl(i,kk)).eq.0) goto 25
       if(km.eq.0) write(5,10) i,(reacao(glgl(i,kk1)),kk1=1,ngl(i))
      km=1
 25   continue
      enddo
      enddo

c     escrita do primeiro no escolhido para gr�fico
      nprint1=glgl(ngf,1)
      nprint2=glgl(ngf,2)
      nprint3=glgl(ngf,3)

c	write(11,11) ipt*dt
c     #,(p(nprint1)-p0(nprint1)),(p(nprint2)-p0(nprint2))
c     #,(p(nprint3)-p0(nprint3)),f(nprint1),f(nprint2),f(nprint3)

c     colocando a velocidade no grafico1
	write(11,11) ipt*dt
     #,(p(nprint1)-p0(nprint1)),(p(nprint2)-p0(nprint2))
     #,(p(nprint3)-p0(nprint3)),vs(nprint1),vs(nprint2),vs(nprint3)

	write(12,11) ipt*dt,
     #(p(nprint1)-p0(nprint1)),(p(nprint2)-p0(nprint2))
     #,(p(nprint3)-p0(nprint3)),vint(nprint1)
     #,vint(nprint2),vint(nprint3)

c     escrita do segundo no escolhido para gr�fico
      nprint1=glgl(ngf1,1)
      nprint2=glgl(ngf1,2)
      nprint3=glgl(ngf1,3)

	write(15,11) ipt*dt,
     #(p(nprint1)-p0(nprint1)),(p(nprint2)-p0(nprint2))
     #,(p(nprint3)-p0(nprint3)),f(nprint1),f(nprint2),f(nprint3)

	write(16,11) ipt*dt,
     #(p(nprint1)-p0(nprint1)),(p(nprint2)-p0(nprint2))
     #,(p(nprint3)-p0(nprint3)),vint(nprint1)
     #,vint(nprint2),vint(nprint3)
      else
	endif
 10	format(1x,i5,6(1x,f17.8))
 11	format(1x,7e27.4)
   
      return
      end



c      *******************************************
c      Saida para pos processamento posi��es e deslocamentos
c      ******************************************* 
      subroutine saidaposnov
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	include 'bloclagatbig.for'

      if(ipt.eq.1) then
       write(8,200)
	 write(8,*)
	 write(8,*)'teste-geral-padrao - txt'
	 write(8,*)
	 write(8,*)'ites de linhas em branco ou com textos,'
	 write(8,*)'a diferente de #'
	 write(8,*)
       write(8,*)'nnos nelem nlistas'
	 write(8,201)
	 write(8,*) nnos,nel,3*(npt/nprint)
       write(8,*)
	 write(8,*)'coordx coordy coordz deslx delsy deslz'
	 write(8,201)
	 do i=1,nnos
	  write(8,111) p0(glgl(i,1)),p0(glgl(i,2)),0.,
     #  p(glgl(i,1))-p0(glgl(i,1)),p(glgl(i,2))-p0(glgl(i,2))
     #,p(glgl(i,3))-p0(glgl(i,3)) !giro do primeiro elemento
	 enddo
	 write(8,*)
	 write(8,*)'tipoetc     *grauaprox n�1 n�2...non'
	 write(8,201)
       do j=1,nel
        write(8,20) 1,(notl(j)-1),(ic(j,k),k=1,notl(j))      
	 enddo
       else
	endif

 120  format(18i6)
 111  format (6(1x,e17.6))
 112  format (4(1x,e17.6))


 200  format(20hpegar comentarios em)
 201  format(1h#)

c     inicio das listas    

      if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     # .lt.1.0d-6)) then
	 write(8,*)
	 write(8,*)'t�tulo da lista e dados dos n�s'
	 write(8,201)
       write(8,*) 'DESL X',ipt
	 do ibet=1,nnos
	  write(8,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),0.,p(glgl(ibet,1))
     #-p0(glgl(ibet,1))
	 enddo
 	 write(8,*)
	 write(8,*)'t�tulo da lista e dados dos n�s'
	 write(8,201)
       write(8,*) 'DESL y',ipt
	 do ibet=1,nnos
	  write(8,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),0.,p(glgl(ibet,2))
     #-p0(glgl(ibet,2))
	 enddo
	 write(8,*)
	 write(8,*)'t�tulo da lista e dados dos n�s'
	 write(8,201)
       write(8,*) 'DESL z',ipt
	 do ibet=1,nnos
	  write(8,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),0.,p(glgl(ibet,3))
     #-p0(glgl(ibet,3))
	 enddo
       else
	endif

 20   format(18i6)
      
      return
      end

c      
c      *******************************************
c      Subrotina que calcula as funcoes de forma
c      e suas derivadas para qualquer ordem de aproxima��o
c      novidade em MEF
c      *******************************************   
      subroutine formaederi0
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

c     o valor do ponto de interesse xsi1 vem de fora da subrotina
c     pode ser um ponto de Gauss (no procedimento de integra��o)
c     Ou uma coordenada qualquer (usualmente nodal)
c     para c�lculo de esfor�o interno por exemplo   
c$$$
c$$$	drl1=0.
c$$$      do ks=1,notl(j)
c$$$	 rl1(ks)=1.
c$$$	 do is=1,notl(j)
c$$$	  drl1e=1.
c$$$	  if(ks.ne.is) then
c$$$	   rl1(ks)=rl1(ks)*(xsi1-xs1(is))/(xs1(ks)-xs1(is)) !fun��o
c$$$	   do js=1,notl(j)
c$$$	    if ((is.ne.js).and.(ks.ne.js)) then
c$$$	     drl1e=drl1e*(xsi1-xs1(js)) !derivada
c$$$           else
c$$$	    endif
c$$$	   enddo
c$$$         drl1(ks)=drl1(ks)+drl1e !derivada
c$$$	   else
c$$$	  endif
c$$$	 enddo
c$$$	enddo
c$$$
c$$$      do is=1,notl(j)
c$$$	 do js=1,notl(j)
c$$$	  if(is.ne.js) then
c$$$	   drl1(is)=drl1(is)/(xs1(is)-xs1(js)) !derivada
c$$$	   else
c$$$	  endif
c$$$	 enddo
c$$$	enddo
c$$$
c$$$	do is=1,notl(j)
c$$$	 fi(is)=rl1(is) !fun��o
c$$$	enddo
c$$$
c$$$	do is=1,notl(j)
c$$$	 dfi(is,1)=drl1(is) !derivada
c$$$	enddo

c     ////////////////////////////////////////////////////////////
c     CALCULATING B-SPLINES
c     ////////////////////////////////////////////////////////////
     
         x3 = xsi1*xsi1*xsi1
         x2 = xsi1*xsi1
         x1 = xsi1
                          
         if (itipoel(j).eq.1) then  !funcoes de forma para o elemento do centro (tipo 01)
      
            fi(1)=-1.d0/48.d0*x3+1.d0/16.d0*x2-1.d0/16.d0*x1+1.d0/48.d0
            
            dfi(1,1)=-3.d0/48.d0*x2+2.d0/16.d0*x1-1.d0/16.d0
              
            fi(2)=(1.d0/16.d0*x3-1.d0/16.d0*x2-5.d0/16.d0*x1+23.d0
     $           /48.d0)
            
            dfi(2,1)=(3.d0/16.d0*x2-2.d0/16.d0*x1-5.d0/16.d0)

            
            fi(3)=(-1.d0/16.d0*x3-1.d0/16.d0*x2+5.d0/16.d0*x1+23.d0
     $           /48.d0)
            
            dfi(3,1)=(-3.d0/16.d0*x2-2.d0/16.d0*x1+5.d0/16.d0)
                          
            fi(4)=1.d0/48.d0*x3+1.d0/16.d0*x2+1.d0/16.d0*x1+1.d0/48.d0
            
            dfi(4,1)=3.d0/48.d0*x2+2.d0/16.d0*x1+1.d0/16.d0

c     write(*,*)x1,fi(1),dfi(1,1),itipoel(j)
         endif

         if (itipoel(j).eq.2)then !funcoes de forma para o elemento da ponta a esquerda (tipo 02)
            fi(1)=(-1.d0/8.d0*x3)+(3.d0/8.d0*x2) -(3.d0/8.d0*x1)+(1.d0
     $           /8.d0)
            
            dfi(1,1)=(-3.d0/8.d0*x2)+(6.d0/8.d0*x1) -(3.d0/8.d0)
           
            fi(2)=(7.d0/32.d0*x3)-(15.d0/32.d0*x2)-(3.d0/32.d0*x1)
     $           +(19.d0/32.d0)
            dfi(2,1)=(21.d0/32.d0*x2)-(30.d0/32.d0*x1)-(3.d0/32.d0) 
            
            fi(3)=(-22.d0/192.d0*x3)+(1.d0/32.d0*x2)+(13.d0/32.d0*x1)
     $           +(25.d0/96.d0)
            dfi(3,1)=(-66.d0/192.d0*x2) +(2.d0/32*x1)+(13.d0/32.d0)
            
            fi(4)= (1.d0/48.d0*x3)+(1.d0/16.d0*x2)+(1.d0/16.d0*x1)+(1.d0
     $           /48.d0)
            dfi(4,1)=(3.d0/48*x2)+(2.d0/16*x1)+(1.d0/16)
            
c           write(*,*)x1,fi(1),dfi(1,1),itipoel(j)
         endif
      
               
                 
         if (itipoel(j).eq.3)then !funcoes de forma para o elemento da ponta a direita (tipo 03)
           
            fi(1)= (-1.d0/48.d0*x3)+ (1.d0/16.d0*x2)-(1.d0/16.d0*x1)+
     $           (1.d0/48.d0)
            dfi(1,1)=(-3.d0/48.d0*x2)+(2.d0/16*x1)-(1.d0/16)

            fi(2)=(22.d0/192.d0*x3) +( 1.d0/32.d0*x2) -(13.d0/32.d0*x1)
     $           +(25.d0/96.d0)
            dfi(2,1)=(66.d0/192.d0*x2) +(2.d0/32*x1) - (13.d0/32.d0)

    
            fi(3)=-(7.d0/32.d0*x3) -(15.d0/32.d0*x2) + (3.d0/32.d0*x1)
     $           +(19.d0/32.d0)
            dfi(3,1)=-(21.d0/32.d0*x2) -(30.d0/32.d0*x1) +(3.d0/32.d0) 

           
            fi(4)=(1.d0/8.d0*x3)+(3.d0/8.d0*x2) + (3.d0/8.d0*x1)+(1.d0
     $           /8.d0)
            dfi(4,1)=(3.d0/8.d0*x2)+(6.d0/8.d0*x1) +(3.d0/8.d0)         
c         write(*,*)x1,fi(1),dfi(1,1),itipoel(j)
         endif
 
          if (itipoel(j).eq.4)then  !funcoes de forma para o segundo elemento da ponta a esquerda (tipo 04)  
            
             fi(1)=(-2.d0/64.d0*x3)+(6.d0/64.d0*x2) - (6.d0/64.d0*x1)
     $            +(2.d0/64.d0)
             dfi(1,1)=(-6.d0/64.d0*x2)+(12.d0/64.d0*x1) - (6.d0/64.d0)

             fi(2)=(14.d0/192.d0*x3)-(3.d0/32.d0*x2)-(9.d0/32.d0*x1)
     $            +(15.d0/32.d0)
             dfi(2,1)=(42.d0/192.d0*x2) -(6.d0/32.d0*x1) -(9.d0/32.d0)

                                                 
             fi(3)=(-1.d0/16.d0*x3)-(1.d0/16.d0*x2)+(5.d0/16.d0*x1)
     $            +(23.d0/48.d0)
             dfi(3,1)=(-3.d0/16.d0*x2)-(2.d0/16.d0*x1) +(5.d0/16.d0)
                           
             fi(4)= (1.d0/48.d0*x3)+ (3.d0/48.d0*x2)+(3.d0/48.d0*x1)+
     $            (1.d0/48.d0)
             dfi(4,1)=(3.d0/48.d0*x2)+(6.d0/48*x1)+(3.d0/48.d0)
            
c             write(*,*)x1,fi(1),dfi(1,1),itipoel(j)
         endif

         if (itipoel(j).eq.5)then !funcoes de forma para o segundo elemento da ponta a direita (tipo 05)
            fi(1)= (-1.d0/48.d0*x3)+ (3.d0/48.d0*x2)-(3.d0/48.d0*x1)+
     $           (1.d0/48.d0)
            dfi(1,1)=-(3.d0/48.d0*x2)+(6.d0/48.d0*x1)-(3.d0/48.d0)


            fi(2)=(1.d0/16.d0*x3) - (1.d0/16.d0*x2) -(5.d0/16.d0*x1)
     $           +(23.d0/48.d0)
            dfi(2,1)=(3.d0/16.d0*x2) - (2.d0/16.d0*x1) -(5.d0/16.d0)

            
            fi(3)=(-14.d0/192.d0*x3) -(3.d0/32.d0*x2) +(9.d0/32.d0*x1)
     $           +(15.d0/32.d0)
            dfi(3,1)=(-42.d0/192.d0*x2) -(6.d0/32.d0*x1) +(9.d0/32.d0)
            
            
           
            fi(4)=(2.d0/64.d0*x3)+(6.d0/64.d0*x2)+(6.d0/64.d0*x1)+(2.d0
     $           /64.d0)
            dfi(4,1)=(6.d0/64.d0*x2)+(12.d0/64.d0*x1) +(6.d0/64.d0)

c            write(*,*)x1,fi(1),dfi(1,1),itipoel(j)
         endif

         sumwbspl=0.d0
         sumdwbspldxsi=0.d0

         do inurbs=1,4
            sumwbspl=sumwbspl+fi(inurbs)*pesoNURBS(ic(j,inurbs))
            sumdwbspldxsi=sumdwbspldxsi+dfi(inurbs,1)*pesoNURBS(ic(j
     $           ,inurbs))

         enddo
         do inurbs=1,4
            dfi(inurbs,1)=(dfi(inurbs,1)*pesoNURBS(ic(j,inurbs))
     $           *sumwbspl-pesoNURBS(ic(j,inurbs))*fi(inurbs)
     $           *sumdwbspldxsi) /(sumwbspl*sumwbspl)
         enddo
c$$$
c$$$         
         do inurbs=1,4
            fi(inurbs)=fi(inurbs)*pesoNURBS(ic(j,inurbs))/sumwbspl
         enddo


         
c$$$
c$$$
c$$$       x3=((xsi1+1.d0)/2.d0)**3.
c$$$       x2=((xsi1+1.d0)/2.d0)**2.
c$$$       x1=((xsi1+1.d0)/2.d0)
c$$$
c$$$      if (itipoel(j).eq.1)then
c$$$
c$$$         fi(1)=-1./6.*x3+1./2.*x2-1./2.*x1+1./6.
c$$$         dfi(1,1)=(-1./2.*x2+x1-1./2.)/2.d0
c$$$         fi(2)=(1./2.*x3-x2+2./3.)
c$$$         dfi(2,1)=(3./2.*x2-2.*x1)/2.
c$$$         fi(3)=(-1./2.*x3+1./2.*x2+1./2.*x1+1./6.)
c$$$         dfi(3,1)=(-3./2.*x2+x1+1./2.)/2.
c$$$         fi(4)=1./6.*x3
c$$$         dfi(4,1)=(1./2.*x2)/2.
c$$$c         write(*,*)x1,fi(1),dfi(1,1)
c$$$      endif
c$$$
c$$$      if (itipoel(j).eq.2)then        
c$$$         fi(1)=(1./6.*x3-x1+1.)
c$$$         dfi(1,1)=((1./2.*x2-1.)/2.)
c$$$         fi(2)=(-1./3.*x3+x1)
c$$$         dfi(2,1)=(-x2+1.)/2.
c$$$         fi(3)=1./6.*x3
c$$$         dfi(3,1)=(1./2.*x2)/2.
c$$$      endif
c$$$
c$$$      if (itipoel(j).eq.3)then
c$$$         fi(1)=-1./6.*x3+1./2.*x2-1./2.*x1+1./6.
c$$$         dfi(1,1)=(-1./2.*x2+x1-1./2.)/2.d0
c$$$         fi(2)=1./3.*x3-x2+2./3.
c$$$         dfi(2,1)=(x2-2.*x1)/2.
c$$$         fi(3)=-1./6.*x3+1./2.*x2+1./2.*x1+1./6.
c$$$         dfi(3,1)=(-1./2.*x2+x1+1./2.)/2.
c$$$
c$$$      endif
      return
      end


c     ************************************************
c     subrotina que define a lei constitutiva para cada elemento
c     ****************************************************
      Subroutine leielast
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'


c     Constantes el�sticas anisotropica, 
c     aqui serao isotropicas por enquanto
c     ex,ey,ez,rxy,rxz,ryz,gxy,gxz,gyz
c     estamos usando o modulo de elasticidade
c     do elemento 1 depois generalizaremos

      ex=re(j)
      ey=ex
      ez=ex
      rxy=ri(j)
      rxz=ri(j)
      ryz=ri(j)
      ryx=rxy*ey/ex
      rzx=rxz*ez/ex
      rzy=ryz*ez/ey
      de=0.d0
      de(1,1)=1.d0/ex
      de(2,2)=1.d0/ey
      de(3,3)=1.d0/ez      
      de(1,2)=-de(1,1)*rxy
      de(1,3)=-de(1,1)*rxz
      de(2,1)=de(1,2)
      de(3,1)=de(1,3)
      de(2,3)=-de(2,2)*ryz
      de(3,2)=de(2,3)
      gxy=dsqrt(ex*ey)/(2.d0*dsqrt((1.d0+rxy)*(1.d0+ryx)))
      gxz=dsqrt(ex*ez)/(2.d0*dsqrt((1.d0+rxz)*(1.d0+rzx)))
      gyz=dsqrt(ey*ez)/(2.d0*dsqrt((1.d0+ryz)*(1.d0+rzy)))
      de(4,4)=1.d0/gxy
      de(5,5)=1.d0/gxz
      de(6,6)=1.d0/gyz

c     a inversa da de � a ce que na realidade � a hook
c      call dlinrg(6,de,6,ceret,6)
      call inverse(6,de,ceret)

c     Constantes Pl�sticas ex,ey,ez,rxy,rxz,ryz,gxy,gxz,gyz
c     aqui o retorno ser� na mesma dire��o da tentativa el�stica
c     no futuro teremos que criar a pl�stica diferente da elastica
c     e ler os valores dos m�dulos pl�sticos (que s� influenciam na dire��o do fluxo)
c       read(3,*)
c       read(3,*) ex,ey,ez,rxy,rxz,ryz

      ryx=rxy*ey/ex
      rzx=rxz*ez/ex
      rzy=ryz*ez/ey
      dpret=0.d0
      dpret(1,1)=1.d0/ex
      dpret(2,2)=1.d0/ey
      dpret(3,3)=1.d0/ez      
      dpret(1,2)=-dpret(1,1)*rxy
      dpret(1,3)=-dpret(1,1)*rxz
      dpret(2,1)=dpret(1,2)
      dpret(3,1)=dpret(1,3)
      dpret(2,3)=-dpret(2,2)*ryz
      dpret(3,2)=dpret(2,3)
      gxy=dsqrt(ex*ey)/(2.d0*dsqrt((1.d0+rxy)*(1.d0+ryx)))
      gxz=dsqrt(ex*ez)/(2.d0*dsqrt((1.d0+rxz)*(1.d0+rzx)))
      gyz=dsqrt(ey*ez)/(2.d0*dsqrt((1.d0+ryz)*(1.d0+rzy)))
      dpret(4,4)=1.d0/gxy
      dpret(5,5)=1.d0/gxz
      dpret(6,6)=1.d0/gyz

c      call dlinrg(6,dpret,6,cp,6)
      call inverse(6,de,ceret)

      return
      end


c     ************************************************
c     subrotina que define a lei constitutiva
c     dever� ser unificada em nota��o com a anterior
c     ****************************************************
      Subroutine definelei
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      rhook=0.d0
      do j=1,nel
         rrk=re(j)/(1.d0+ri(j))/(1.d0-2.d0*ri(j))
         rri=ri(j)
         rrg=rg(j)
         
         rhook(j,1,1)=rrk*(1.d0-rri)
         rhook(j,1,2)=rrk*rri
         rhook(j,2,1)=rrk*rri
         rhook(j,2,2)=rrk*(1.d0-rri)
         rhook(j,3,3)=rrg
         
c     rhook � o mesmo que a cret, dever� ser unificada a nota��o n programa
c     em breve
         
      enddo
      return
      end


c     ************************************************
c     Subrotina que define a lei constitutiva estado plano de tens�es
c     desativado
c     ****************************************************
      Subroutine defineleiEPT
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      rhook=0.d0
      do j=1,nel
         rrk=re(j)/(1.d0+ri(j))
         rri=ri(j)
         rrg=rg(j)
         
         rhook(j,1,1)=rrk*1.d0/(1.d0-rri)
         rhook(j,1,2)=rrk*rri*1.d0/(1.d0-rri)
         rhook(j,2,1)=rrk*rri*1.d0/(1.d0-rri)
         rhook(j,2,2)=rrk*1.d0/(1.d0-rri)
         rhook(j,3,3)=rrg
         
      enddo
      return
      end


c     ************************************************
c     Subrotina que calcula segundas derivadas das deformacoes
c     em funcao das posicoes nodais (gls) locais
c     para cada no li
c     novidade semelhante � deridef
c     ****************************************************
      Subroutine deri2def
      USE BIG_STUFF	 
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     precisamos registrar qual e o no do elemento que estamos considerando
c     para usar na segunda derivadaS
      d2rex=0.d0
      d2rey=0.d0
      d2rez=0.d0
      d2rexy=0.d0
      d2rexz=0.d0
      d2reyz=0.d0

      do k1=1,3
         do k2=1,3
            
c     derivada da deformacao ex
            
            aup=0.d0
            do i1=1,2
               do j1=1,2
                  aup(i1)=aup(i1)+d2ch(i1,j1,k1,k2)*up(j1)
               enddo
            enddo
            do i1=1,2
               d2rex(k1,k2)=d2rex(k1,k2)+aup(i1)*up(i1)/rlex/2.d0
            enddo
            d2rex(k1,k2)=d2rex(k1,k2)-drex(k1,li)*drex(k2,lj)/rlex
            
c     derivada da deformacao ey
            
            aup=0.d0
            do i1=1,2
               do j1=1,2
                  aup(i1)=aup(i1)+d2ch(i1,j1,k1,k2)*vp(j1)
               enddo
            enddo
            do i1=1,2
               d2rey(k1,k2)=d2rey(k1,k2)+aup(i1)*vp(i1)/rley/2.d0
            enddo
            d2rey(k1,k2)=d2rey(k1,k2)-drey(k1,li)*drey(k2,lj)/rley

            goto 10             !2d

c     derivada da deformacao ez
            
            aup=0.d0
            do i1=1,2
               do j1=1,2
                  aup(i1)=aup(i1)+d2ch(i1,j1,k1,k2)*wp(j1)
               enddo
            enddo
            do i1=1,2
               d2rez(k1,k2)=d2rez(k1,k2)+aup(i1)*wp(i1)/rlez/2.d0
            enddo
            d2rez(k1,k2)=d2rez(k1,k2)-drez(k1,li)*drez(k2,lj)/rlez
 10         continue
            
            goto 40 
c     derivada da deformacao exy
            
            aup=0.d0
            aup1=0.d0
            aup2=0.d0
            aup12=0.d0
            do i1=1,2
               do j1=1,2
                  aup(i1)=aup(i1)+dch(i1,j1,k1,li)*vp(j1)
                  aup2(i1)=aup(i1)+dch(i1,j1,k2,lj)*vp(j1)
                  aup1(i1)=aup1(i1)+ch(i1,j1)*vp(j1)
                  aup12(i1)=aup12(i1)+d2ch(i1,j1,k1,k2)*vp(j1)
               enddo
            enddo
            
            aur=0.d0
            aur1=0.d0
            aur2=0.d0        
            aur12=0.d0

            do i1=1,2
               aur=aur+aup(i1)*up(i1)
               aur2=aur2+aup2(i1)*up(i1)
               aur1=aur1+aup1(i1)*up(i1)
               aur12=aur12+aup12(i1)*up(i1)
            enddo

            rdent=dabs((rlex*rley)**2-aur1**2)
            roxy=drex(k2,lj)*rley+drey(k2,lj)*rlex
        
            aa=1.d0/(2.d0*rlex*rley*dsqrt(dabs(rdent)))

            aa1=-aa**2.d0*(2.d0*roxy*rdent**.5d0+rlex*rley*rdent**(
     $           -.5d0)*(2.d0*rlex*rley*roxy-2.d0*aur1*aur2))

            d2rexy(k1,k2)=aa1*(aur*rlex*rley -aur1*(rley*drex(k1,li)
     $           +rlex*drey(k1,li)))+ aa*(aur12*rley*rlex+aur*(rley
     $           *drex(k2,lj)+rlex*drey(k2,lj))- aur2*(rley*drex(k1,li)
     $           +rlex*drey(k1,li)) -aur1*(drey(k2,lj)*drex(k1,li) +rley
     $           *d2rex(k1,k2)+drex(k2,lj)*drey(k1,li) +rlex*d2rey(k1
     $           ,k2)))

            d2rexy(k1,k2)=d2rexy(k1,k2)*2.d0 !gama
            
            goto 20             !2d

c      derivada da deformacao exz

            aup=0.d0
            aup1=0.d0
            aup2=0.d0
            aup12=0.d0
            do i1=1,2
               do j1=1,2
                  aup(i1)=aup(i1)+dch(i1,j1,k1,li)*wp(j1)
                  aup2(i1)=aup(i1)+dch(i1,j1,k2,lj)*wp(j1)
                  aup1(i1)=aup1(i1)+ch(i1,j1)*wp(j1)
                  aup12(i1)=aup12(i1)+d2ch(i1,j1,k1,k2)*wp(j1)
               enddo
            enddo

            aur=0.d0
            aur1=0.d0
            aur2=0.d0
            aur12=0.d0

            do i1=1,2
               aur=aur+aup(i1)*up(i1)
               aur2=aur2+aup2(i1)*up(i1)
               aur1=aur1+aup1(i1)*up(i1)
               aur12=aur12+aup12(i1)*up(i1)
            enddo

            rdent=dabs((rlex*rlez)**2.d0-aur1**2.d0)
            roxz=drex(k2,lj)*rlez+drez(k2,lj)*rlex

            aa=1.d0/(2.d0*rlex*rlez*dsqrt(dabs(rdent)))
            
            aa1=-aa**2.d0*(2.d0*roxz*rdent**.5d0+rlex*rlez*rdent**(
     $           -.5d0)*(2.d0*rlex*rlez*roxz-2.d0*aur1*aur2))

            d2rexz(k1,k2)=aa1*(aur*rlex*rlez -aur1*(rlez*drex(k1,li)
     $           +rlex*drez(k1,li)))+ aa*(aur12*rlez*rlex+aur*(rlez
     $           *drex(k2,lj)+rlex*drez(k2,lj))- aur2*(rlez*drex(k1,li)
     $           +rlex*drez(k1,li)) -aur1*(drez(k2,lj)*drex(k1,li) +rlez
     $           *d2rex(k1,k2)+drex(k2,lj)*drez(k1,li) +rlex*d2rez(k1
     $           ,k2)))

            d2rexz(k1,k2)=d2rexz(k1,k2)*2.d0 !gama

c      derivada da deformacao eyz

            aup=0.d0
            aup1=0.d0
            aup2=0.d0
            aup12=0.d0

            do i1=1,2
               do j1=1,2
                  aup(i1)=aup(i1)+dch(i1,j1,k1,li)*wp(j1)
                  aup2(i1)=aup(i1)+dch(i1,j1,k2,lj)*wp(j1)
                  aup1(i1)=aup1(i1)+ch(i1,j1)*wp(j1)
                  aup12(i1)=aup12(i1)+d2ch(i1,j1,k1,k2)*wp(j1)
               enddo
            enddo

            aur=0.d0
            aur1=0.d0
            aur2=0.d0
            aur12=0.d0
            do i1=1,2
               aur=aur+aup(i1)*vp(i1)
               aur2=aur2+aup2(i1)*vp(i1)
               aur1=aur1+aup1(i1)*vp(i1)
               aur12=aur12+aup12(i1)*vp(i1)
            enddo

            rdent=dabs((rley*rlez)*(rley*rlez)-aur1*aur1)
            royz=drey(k2,lj)*rlez+drez(k2,lj)*rley

            aa=1.d0/(2.d0*rley*rlez*dsqrt(dabs(rdent)))

            aa1=-aa*aa*(2.d0*royz*rdent**.5d0+rley*rlez*rdent**(-.5d0)*
     $           (2.d0*rley*rlez*royz-2.d0*aur1*aur2))

            d2reyz(k1,k2)=aa1*(aur*rley*rlez -aur1*(rlez*drey(k1,li)
     $           +rley*drez(k1,li)))+ aa*(aur12*rlez*rley+aur*(rlez
     $           *drey(k2,lj)+rley*drez(k2,lj))- aur2*(rlez*drey(k1,li)
     $           +rley*drez(k1,li)) -aur1*(drez(k2,lj)*drey(k1,li) +rlez
     $           *d2rey(k1,k2)+drey(k2,lj)*drez(k1,li) +rley*d2rez(k1
     $           ,k2)))

            d2reyz(k1,k2)=d2reyz(k1,k2)*2.d0 !gama

 20         continue
 40         continue

         enddo
      enddo


      return
      end

c     ************************************************
c     subrotina que calcula a segunda derivadas do tensor de cauchy-green
c     em relacao a cada gl 1,2,3 para um determind no
c     ****************************************************
      Subroutine d2cauchy
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      d2af=0.d0
      
      do k1=1,3
         do k2=1,3
            aaux=0.d0
            do ii=1,2
               do jj=1,2
                  do kk=1,2
                     aaux(ii,jj)=aaux(ii,jj)+daf(ii,kk,lj,k2)*a0inv(kk
     $                    ,jj)
                  enddo
               enddo
            enddo
            aaux1=0.d0
            do ii=1,2
               do jj=1,2
                  do kk=1,2
                     aaux1(ii,jj)=aaux1(ii,jj)+daf(ii,kk,li,k1)*a0inv(kk
     $                    ,jj)
                  enddo
               enddo
            enddo

c     nao zero para acumular
            
            do ii=1,2
               do jj=1,2
                  do kk=1,2
                     d2ch(ii,jj,k1,k2)=d2ch(ii,jj,k1,k2)+aaux1(kk,ii)
     $                    *aaux(kk,jj)+aaux(kk,ii)*aaux1(kk,jj)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      return
      end

c      
c     *******************************************
c     Gerencia o calculo de tensoes
c     para as fibras superiores de uma barra
c     intuito apenas de prepara��o
c     para p�s-orocessamento
c     *******************************************   
      subroutine tensaoacima
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     nel � o numero de elementos
      
      do j=1,nel
         call eletensaoacima
      enddo
      
      return
      end
c     
c     *******************************************
c     calcula tensos nos elementos
c     *******************************************   
      subroutine eletensaoacima
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     atualizar configura��o nos elementos
      
      call leielast
      
      ng=notl(j) !notl(j)
c     aquisi��o das coordenadas e pesos de gaus j� fiz no in�cio do programa
      call gauss(ng,qsi,w)
      call gauss(ng3,qsi3,w3)
c     
c$$$  call dgqrul(ng,1,0.d0,0.d0,0,qsi,qsi,w)
c$$$  call dgqrul(ng3,1,0.d0,0.d0,0,qsi3,qsi3,w3)
      call geometriaatual 
      
      do ig1=1,ng
         
         up=0.
         vp=0.
         
c     xsi1=(2.*ig1-(1.+notl(j)))/(notl(j)-1.)
         
         xsi1=qsi(ig1)
c     ***********************
c     variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c         call suportefuncoesdeforma
         call formaederi0      
         
         igf=3                  ! a terceira fibra � a de cima
      ig3=ng3
      
      
c     xsi 3 para face inferior 
      xsi3=qsi3(ig3)
      
      call calc_A0              !matriz gradiente do adimensional para inicial      
      
      up(1)=cos(salfa-pi/2.)
      up(2)=sin(salfa-pi/2.)
      
      vp(1)=cos(salfa)
      vp(2)=sin(salfa)
      
      
      call calc_Af              !matriz do adimensonal para o atual
      call cauchy               !gradiente da transformacao total
      
      
      call deformacoes          !deformacoes longitudinais, distorcoes e alongamentos
      call energia
      call derivenertensao      !vetor de forcas internas
      
      enddo
      
      
      return
      end




c      
c     *******************************************
c     Gerencia o c�lculo os esfor�os solicitantes internos
c     prepara��o de resultados para p�s-processamento
c     de engenharia
c     *******************************************   
      subroutine momento
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     nel � o numero de elementos
      
      rmz=0
      rvy=0.
      rnx=0.
      
      do j=1,nel
         call elemom
      enddo
      
      return
      end
c     
c     *******************************************
c     Calcula os esfor�os
c     *******************************************   
      subroutine elemom
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      call leielast
      
c        atualizar configura��o nos elementos

      
      ng=notl(j) !notl(j)
c     aquisi��o das coordenadas e pesos de gaus j� fiz no in�cio do programa
      call gauss(ng,qsi,w)
      call gauss(ng3,qsi3,w3)
c$$$c
c$$$      call dgqrul(ng,1,0.d0,0.d0,0,qsi,qsi,w)
c$$$      call dgqrul(ng3,1,0.d0,0.d0,0,qsi3,qsi3,w3)


      call geometriaatual 

      do ig1=1,ng
         up=0.
         vp=0.

c      xsi1=(2.*ig1-(1.+notl(j)))/(notl(j)-1.)

         xsi1=qsi(ig1)
c       ***********************
c       variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c         call suportefuncoesdeforma
         call formaederi0      

       do igf=1,3 ! a terceira fibra � a de cima
       do ig3=1,ng3
          up=0.
          vp=0.

c      xsi 3 para face inferior 
         xsi3=qsi3(ig3)

         call calc_A0           !matriz gradiente do adimensional para inicial      

         up(1)=cos(salfa-pi/2.)
         up(2)=sin(salfa-pi/2.)
         
         vp(1)=cos(salfa)
         vp(2)=sin(salfa)
         
         
         call calc_Af           !matriz do adimensonal para o atual
         call cauchy            !gradiente da transformacao total
         
         
         call deformacoes       !deformacoes longitudinais, distorcoes e alongamentos
         call energia
         call derivenermom      !vetor de forcas internas
         
      enddo
      enddo
      enddo
      
      return
      end

c      
c     *******************************************
c     Gerencia o calculo de tensoes
c     para as fibras inferiores de uma barra
c     intuito apenas de prepara��o
c     para p�s-orocessamento
c     *******************************************   
      subroutine tensaoabaixo
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     nel � o numero de elementos
      
      do j=1,nel
         call eletensaoabaixo
      enddo
      
      return
      end
      
c     
c     *******************************************
c     Calcula as tensoea para as fibras inferiores
c     *******************************************   
      subroutine eletensaoabaixo
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
c     atualizar configura��o nos elementos
      
      call leielast
      
      ng=notl(j) !notl(j)
c     aquisi��o das coordenadas e pesos de gaus j� fiz no in�cio do programa
      call gauss(ng,qsi,w)
      call gauss(ng3,qsi3,w3)     
c$$$  call dgqrul(ng,1,0.d0,0.d0,0,qsi,qsi,w)
c$$$  call dgqrul(ng3,1,0.d0,0.d0,0,qsi3,qsi3,w3)
      call geometriaatual 
      
      do ig1=1,ng
         up=0.
         vp=0.
         
c     xsi1=(2.*ig1-(1.+notl(j)))/(notl(j)-1.)
         
         xsi1=qsi(ig1)
c     ***********************
c     variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c         call suportefuncoesdeforma
         call formaederi0      
         
         igf=1                  ! a terceira fibra � a de cima
         ig3=1
         
c     xsi 3 para face inferior 
         xsi3=qsi3(ig3)
         
         call calc_A0           !matriz gradiente do adimensional para inicial      
         
         up(1)=cos(salfa-pi/2.)
         up(2)=sin(salfa-pi/2.)
         
         vp(1)=cos(salfa)
         vp(2)=sin(salfa)
         
         
         call calc_Af           !matriz do adimensonal para o atual
         call cauchy            !gradiente da transformacao total
         
         
         call deformacoes       !deformacoes longitudinais, distorcoes e alongamentos
         call energia
         call derivenertensao   !vetor de forcas internas
         
      enddo
      
      return
      end


c     ************************************************
c     Subrotina que calcula as tens�es de engenharia
c     para p�s-processamento semelhante � derivener
c     ****************************************************
      Subroutine derivenertensao
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

c      duex, duey, duez etc sao as tensoes x, y, z etc el�sticas fict�cias
c      retirando destas as tens�es pl�sticas acumuladas teremos a tens�o tentativa
c      entrando na rotina de retorno, teremos como resposta o n�vel de tens�o real
c      que dever� ocupar a posi��o dos duex, duey, duez etc no c�culo da 
c      forca interna.


c     c�lculo do valor da temperatura para o ponto de integra��o
      
      tpos=(a1*xsi3+fd(j,igf))*(temc(j)-temb(j))/2.+
     #(temc(j)+temb(j))/2.

      et(1)=rex-tpos*calt(j)
      et(2)=rey-tpos*0.
      et(3)=rez
      et(4)=rexy
      et(5)=rexz
      et(6)=reyz

c     CHAMADA DA SUBROTINA DE RETORNO na derivener
c     vamos corrigir o valor das tens�es calculadas

      do iep=1,6
         ee(iep)=et(iep)-ep(j,ig1,ig3,igf,iep)
      enddo
      
      s=0.
      do i=1,6
         do jr=1,6
            s(i)=s(i)+ceret(i,jr)*ee(jr)
         enddo
      enddo
      
      duex=s(1)
      duey=s(2)
      duez=s(3)
      duexy=s(4)
      duexz=s(5)
      dueyz=s(6)


      sx(j,ig1)=duex
      sy(j,ig1)=duey
      sxy(j,ig1)=duexy

      return
      end

c     ************************************************
c     Subrotina que calcula as tens�es de engenharia
c     para integra��o dos esfor�os solicitantes
c     ****************************************************
      Subroutine derivenermom
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

c      duex, duey, duez etc sao as tensoes x, y, z etc el�sticas fict�cias
c      retirando destas as tens�es pl�sticas acumuladas teremos a tens�o tentativa
c      entrando na rotina de retorno, teremos como resposta o n�vel de tens�o real
c      que dever� ocupar a posi��o dos duex, duey, duez etc no c�culo da 
c      forca interna.

c     c�lculo do valor da temperatura para o ponto de integra��o
      
      tpos=(a1*xsi3+fd(j,igf))*(temc(j)-temb(j))/2.+
     #(temc(j)+temb(j))/2.

      et(1)=rex-tpos*calt(j)
      et(2)=rey-tpos*0.
      et(3)=rez
      et(4)=rexy
      et(5)=rexz
      et(6)=reyz


c     CHAMADA DA SUBROTINA DE RETORNO na derivener
c     vamos corrigir o valor das tens�es calculadas

      do iep=1,6
         ee(iep)=et(iep)-ep(j,ig1,ig3,igf,iep)
      enddo
      
      s=0.
      do i=1,6
         do jr=1,6
            s(i)=s(i)+ceret(i,jr)*ee(jr)
         enddo
      enddo
      
      duex=s(1)
      duey=s(2)
      duez=s(3)
      duexy=s(4)
      duexz=s(5)
      dueyz=s(6)

      rmz(j,ig1)=rmz(j,ig1)+duex*rjac2*w3(ig3)*(a1*xsi3+fd(j,igf))
      rvy(j,ig1)=rvy(j,ig1)+duexy*rjac2*w3(ig3)
      rnx(j,ig1)=rnx(j,ig1)+duex*rjac2*w3(ig3)
      
      return
      end


c      *******************************************
c      Saida para pos processamento tens�es nas fibras superiores
c      ******************************************* 
      subroutine saidaposnovtensaoacima
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      if(ipt.eq.1) then
       write(9,200)
	 write(9,*)
	 write(9,*)'teste-geral-padrao - txt'
	 write(9,*)
	 write(9,*)'ites de linhas em branco ou com textos,'
	 write(9,*)'a diferente de #'
	 write(9,*)
       write(9,*)'nnos nelem nlistas'
	 write(9,201)
	 write(9,*) nnos,nel,3*(npt/nprint)
       write(9,*)
	 write(9,*)'coordx coordy coordz deslx delsy deslz'
	 write(9,201)
	 do i=1,nnos
	  write(9,111) p0(glgl(i,1)),p0(glgl(i,2)),p0(glgl(i,3)),
     #  p(glgl(i,1))-p0(glgl(i,1)),p(glgl(i,2))
     #-p0(glgl(i,2)),p(glgl(i,3))-p0(glgl(i,3))
	 enddo
	 write(9,*)
	 write(9,*)'tipoetc     *grauaprox n�1 n�2...non'
	 write(9,201)
       do j=1,nel
        write(9,20) 1,(notl(j)-1),(ic(j,k),k=1,notl(j))      
	 enddo
       else
	endif

 120  format(18i6)
 111  format (6(1x,e17.6))
 112  format (4(1x,e17.6))
 200  format(20hpegar comentarios em)
 201  format(1h#)


      if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     #.lt.1.0d-6)) then

c     inicio das listas    
c     preparando o tensao x
         valor=0.
         do j=1,nel
            into=notl(j)            
            call prepara
            
            do ino=1,notl(j)	   
               do ih=1,ng
                  valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sx(j,ih)
	   enddo
      enddo
      deALLOCATE (rmfi)
      
      enddo
      
c     valor=sx
      
      fno=0.
      do j=1,nel
         DO ino=1,notl(j)       !monta o vetor de valores nodais
            fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
         ENDDO
      enddo
      do j=1,nnos
         fno(j)=fno(j)/(kn(j)*1.)
      enddo
      write(9,*)
      write(9,*)'t�tulo da lista e dados dos n�s'
      write(9,201)
      write(9,*) 'sx',ipt
      do ibet=1,nnos
         write(9,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet)
      enddo
      
c     preparando tensao y
      valor=0.
      do j=1,nel
         ng=notl(j) !notl(j)
            call prepara
	  do ino=1,notl(j)
	   do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sy(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)
	 enddo
       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
 	 write(9,*)
	 write(9,*)'t�tulo da lista e dados dos n�s'
	 write(9,201)
       write(9,*) 'sy',ipt
	 do ibet=1,nnos
	  write(9,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet)
       enddo

c     preparando tensaoxy
       valor=0.
       do j=1,nel
        ng=notl(j) !notl(j)
            call prepara
	  do ino=1,notl(j)
	   do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sxy(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)
	 enddo
       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
   	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
 	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
	 write(9,*)
	 write(9,*)'t�tulo da lista e dados dos n�s'
	 write(9,201)
       write(9,*) 'sxy',ipt
	 do ibet=1,nnos
	  write(9,112) p(glgl(ibet,1))-p0(glgl(ibet,1))
     #,p(glgl(ibet,2))-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet) 
       enddo
       else
	endif

 20   format(18i6)

      return
      end

c      *******************************************
c      Saida para pos processamento de esfor�os solicitantes
c       por n�
c      ******************************************* 
      subroutine saidaposmomento
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      if(ipt.eq.1) then
       write(90,200)
	 write(90,*)
	 write(90,*)'teste-geral-padrao - txt'
	 write(90,*)
	 write(90,*)'ites de linhas em branco ou com textos,'
	 write(90,*)'a diferente de #'
	 write(90,*)
       write(90,*)'nnos nelem nlistas'
	 write(90,201)
	 write(90,*) nnos,nel,3*(npt/nprint)
       write(90,*)
	 write(90,*)'coordx coordy coordz deslx delsy deslz'
	 write(90,201)
	 do i=1,nnos
	  write(90,111) p0(glgl(i,1)),p0(glgl(i,2)),p0(glgl(i,3)),
     #  p(glgl(i,1))-p0(glgl(i,1)),p(glgl(i,2))
     #-p0(glgl(i,2)),p(glgl(i,3))-p0(glgl(i,3))
	 enddo
	 write(90,*)
	 write(90,*)'tipoetc     *grauaprox n�1 n�2...non'
	 write(90,201)
       do j=1,nel
        write(90,20) 1,(notl(j)-1),(ic(j,k),k=1,notl(j))      
	 enddo
       else
	endif

 120  format(18i6)
 111  format (6(1x,e17.6))
 112  format (4(1x,e17.6))
 200  format(20hpegar comentarios em)
 201  format(1h#)


      if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     #.lt.1.0d-6)) then

c     inicio das listas    
c     preparando a forca normal
       valor=0.
       do j=1,nel
	      ng=notl(j) !notl(j)

            call prepara

	  do ino=1,notl(j)	   
	  do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*rnx(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)

	 enddo

       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
	 write(90,*)
	 write(90,*)'t�tulo da lista e dados dos n�s'
	 write(90,201)
       write(90,*) 'N',ipt
	 do ibet=1,nnos
	  write(90,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet)
       enddo

c     preparando forca cortant
       valor=0.
       do j=1,nel
	      ng=notl(j) !notl(j)
            call prepara
	  do ino=1,notl(j)
	   do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*rvy(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)
	 enddo
       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
 	 write(90,*)
	 write(90,*)'t�tulo da lista e dados dos n�s'
	 write(90,201)
       write(90,*) 'V',ipt
	 do ibet=1,nnos
	  write(90,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet)
       enddo

c     preparando momento fletor
       valor=0.
       do j=1,nel
        ng=notl(j) !notl(j)
            call prepara
	  do ino=1,notl(j)
	   do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*rmz(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)
	 enddo
       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
   	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
 	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
	 write(90,*)
	 write(90,*)'t�tulo da lista e dados dos n�s'
	 write(90,201)
       write(90,*) 'M',ipt
	 do ibet=1,nnos
	  write(90,112) p(glgl(ibet,1))-p0(glgl(ibet,1))
     #,p(glgl(ibet,2))-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet) 
       enddo
       else
	endif

 20   format(18i6)

      return
      end

c      *******************************************
c      Saida para pos processamento esfor;os solicitantes por elemento
c      ******************************************* 
      subroutine momentonovo
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
       k=0
	 do ig=1,nel
	  do jl=1,notl(ig)
      	k=k+1
        enddo
	 enddo

      if(ipt.eq.1) then
       write(90,200)
	 write(90,*)
	 write(90,*)'teste-geral-padrao - txt'
	 write(90,*)
	 write(90,*)'ites de linhas em branco ou com textos,'
	 write(90,*)'a diferente de #'
	 write(90,*)
       write(90,*)'nnos nelem nlistas'
	 write(90,201)
	 write(90,*) k,nel,3*(npt/nprint)
       write(90,*)
	 write(90,*)'coordx coordy coordz deslx delsy deslz'
	 write(90,201)


	 do ig=1,nel
	  do jl=1,notl(ig)
        i=ic(ig,jl)
	  write(90,111) p0(glgl(i,1)),p0(glgl(i,2)),p0(glgl(i,3)),
     #  p(glgl(i,1))-p0(glgl(i,1)),p(glgl(i,2))
     #-p0(glgl(i,2)),p(glgl(i,3))-p0(glgl(i,3))
	  enddo
	 enddo


	 write(90,*)
	 write(90,*)'tipoetc     *grauaprox n�1 n�2...non'
	 write(90,201)
        kl=0
       do j=1,nel
           write(90,20) 1,(notl(j)-1),(kl+k,k=1,notl(j))      
	  kl=kl+notl(j)
	 enddo
       else
	endif

 120  format(18i6)
 111  format (6(1x,e17.6))
 112  format (4(1x,e17.6))
 122  format (1i6,4(1x,e17.6))
 200  format(20hpegar comentarios em)
 201  format(1h#)


      if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     #.lt.1.0d-6)) then

c     inicio das listas    
c     preparando a forca normal nos nos do elemento
       valor=0.
       do j=1,nel

	      ng=notl(j) !notl(j) !numero total de nos do elemento = ao numero de pontos de gauss

	  do ino=1,notl(j)	   
	    valor(j,ino)=rnx(j,ino)
	  enddo
      enddo


	 write(90,*)
	 write(90,*)'t�tulo da lista e dados dos n�s'
	 write(90,201)
       write(90,*) 'N',ipt
       write(5,*) 'N',ipt
	 do ig=1,nel
	  do il=1,notl(ig)
	   ibet=ic(ig,il)

	  write(90,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),valor(ig,il) ! imprimindo a lista

	  write(5,122) ig,p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),valor(ig,il) ! imprimindo a lista

       enddo
	enddo


c     preparando forca cortant no elemento
       valor=0.
       do j=1,nel

	      ng=notl(j) !notl(j) !numero total de nos do elemento = ao numero de pontos de gauss

	  do ino=1,notl(j)	   
	    valor(j,ino)=rvy(j,ino)
	  enddo

	 enddo


	 write(90,*)
	 write(90,*)'t�tulo da lista e dados dos n�s'
	 write(90,201)
       write(90,*) 'V',ipt
       write(5,*) 'V',ipt
	 do ig=1,nel
	  do il=1,notl(ig)
	   ibet=ic(ig,il)
	  write(90,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),valor(ig,il) ! imprimindo a lista

	  write(5,122) ig,p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),valor(ig,il) ! imprimindo a lista

       enddo
	enddo

c     preparando momento fletor no elemento
       valor=0.
       do j=1,nel
        ng=notl(j) !notl(j)
	  do ino=1,notl(j)
	    valor(j,ino)=rmz(j,ino)
	  enddo
	 enddo


	 write(90,*)
	 write(90,*)'t�tulo da lista e dados dos n�s'
	 write(90,201)
       write(90,*) 'M',ipt
       write(5,*) 'M',ipt
	 do ig=1,nel
	  do il=1,notl(ig)
	   ibet=ic(ig,il)

	  write(90,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),valor(ig,il) ! imprimindo a lista

	  write(5,122) ig,p(glgl(ibet,1))-p0(glgl(ibet,1)),p(glgl(ibet,2))
     #-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),valor(ig,il) ! imprimindo a lista


       enddo
	enddo

       else
	endif

 20   format(18i6)

      return
      end

c      *******************************************
c      Saida para pos processamento tens�es nas fibras inferiores
c      ******************************************* 
      subroutine saidaposnovtensaoabaixo
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      if(ipt.eq.1) then
       write(19,200) 
	 write(19,*)
	 write(19,*)'teste-geral-padrao - txt'
	 write(19,*)
	 write(19,*)'ites de linhas em branco ou com textos,'
	 write(19,*)'a diferente de #'
	 write(19,*)
       write(19,*)'nnos nelem nlistas'
	 write(19,201)
	 write(19,*) nnos,nel,2*(npt/nprint)
       write(19,*)
	 write(19,*)'coordx coordy coordz deslx delsy deslz'
	 write(19,201)
	 do i=1,nnos
	  write(19,111) p0(glgl(i,1)),p0(glgl(i,2)),p0(glgl(i,3)),
     #  p(glgl(i,1))-p0(glgl(i,1)),p(glgl(i,2))-p0(glgl(i,2)),
     #p(glgl(i,3))-p0(glgl(i,3))
	 enddo
	 write(19,*)
	 write(19,*)'tipoetc     *grauaprox n�1 n�2...non'
	 write(19,201)
       do j=1,nel
        write(19,20) 1,(notl(j)-1),(ic(j,k),k=1,notl(j))      
	 enddo
       else
	endif

 120   format(18i6)
 111  format (6(1x,e17.6))
 112  format (4(1x,e17.6))
 200  format(20hpegar comentarios em)
 201  format(1h#)

      if ((ipt.eq.1).or.((ipt/(nprint*1.)-int(ipt/(nprint*1.)))
     #.lt.1.0d-6)) then

c      inicio das listas    
c      preparando o tensao principal1
       valor=0.
       do j=1,nel
	      ng=notl(j) !notl(j)
            call prepara
	  do ino=1,notl(j)
	   do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sx(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)
	 enddo

      valor=sx      
       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
	 write(19,*)
	 write(19,*)'t�tulo da lista e dados dos n�s'
	 write(19,201)
       write(19,*) 's1',ipt
	 do ibet=1,nnos
	  write(19,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),
     #p(glgl(ibet,2))-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet)
       enddo

c      preparando tensao principal 2
       valor=0.
       do j=1,nel
	      ng=notl(j) !notl(j)
            call prepara
	  do ino=1,notl(j)
	   do ih=1,ng
	    valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sy(j,ih)
	   enddo
	  enddo
	      deALLOCATE (rmfi)
	 enddo
       fno=0.
       do j=1,nel
        DO ino=1,notl(j)  !monta o vetor de valores nodais
	   fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
	  ENDDO
	 enddo
       do j=1,nnos
	  fno(j)=fno(j)/(kn(j)*1.)
       enddo
	 write(19,*)
	 write(19,*)'t�tulo da lista e dados dos n�s'
	 write(19,201)
       write(19,*) 's2',ipt
	 do ibet=1,nnos
	  write(19,112) p(glgl(ibet,1))-p0(glgl(ibet,1)),
     #p(glgl(ibet,2))-p0(glgl(ibet,2)),
     #  p(glgl(ibet,3))-p0(glgl(ibet,3)),fno(ibet)
       enddo

      else
	endif
 20   format(18i6)

      return
      end
  
c     *******************************************
c     Cria��o de matriz tranforma��o de valores nos pontos de Gauss
c     para valores nos n�s do elemento
c     novidade para mef
c     *******************************************   
      subroutine prepara
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      ALLOCATE (rmfi(ng,ng))
      rmfi=0.d0
      
c     Aquisi��o das coordenadas e pesos de Gauss
      call gauss(ng,qsi,w)
c     call dgqrul(ng,1,0.d0,0.d0,0,qsi,qsi,w)
      
      do ig1=1,ng
         
c     xsi1=qsi(ig1)
        
!         xsi1=(2.*ig1-(1.+notl(j)))/(notl(j)-1.) !neste caso as funcoes de forma sao em numero igual ao pontos de gaus
         
c     ***********************
c     variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c     call suportefuncoesdeforma
         call formaederi0
         
         do ifunc=1,ng
            rmfi(ig1,ifunc)=fi(ifunc)
         enddo
      enddo
      
c     CALL DLINRG(ng,rmfi,ng,rmfi,ng) !problema com a dimensao
      call inverse(ng,rmfi,rmfi)
      return
      end


      

c      *******************************************
c      Subrotina que calcula a matriz de massa (lagrangeana, constante)
c      e tamb�m a matriz que transforma carregamentos 
c      distribu�dos em concentrados equivalentes.
c      *******************************************   
      subroutine massa
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      ng=notl(j)                !notl(j)
c     Aquisi��o das coordenadas e pesos de Gauss
      call gauss(ng,qsi,w)
      call gauss(ng3,qsi3,w3)
c     
c     
c      call dgqrul(ng,1,0.d0,0.d0,0,qsi,qsi,w)
c      call dgqrul(ng3,1,0.d0,0.d0,0,qsi3,qsi3,w3)

      do ig1=1,ng
         xsi1=qsi(ig1)
c     ***********************
c     Variaveis necess�rias para funcao de forma e derivadas de ordem qualquer
c     call suportefuncoesdeforma
         call formaederi0      
         do igf=1,3             !Para as faixas (fibras)
            do ig3=1,ng3        ! para altura (de cada fibra)
               xsi3=qsi3(ig3)
               call calc_A0     !Matriz gradiente do adimensional para inicial      
               do ir=1,notl(j)
                  do jm=1,notl(j)
                     rmp(ir,jm,j)=rmp(ir,jm,j)+fi(ir)*fi(jm)*w(ig1)
     $                    *w3(ig3)*rjac0   !     massa
                                                                                  !consistente
                                                                                  !pequena
c$$$                       write(*,*)' rkl',j
c$$$            do iaaa=1,notl(j)
c$$$               do ibbb=1,notl(j)
c$$$                  write(*,*)iaaa,ibbb,rmp(iaaa,ibbb,j)
c$$$               enddo
c$$$            enddo
                  enddo
               enddo
            enddo               ! faixas
         enddo                  !altura
      enddo                     !comprimento
       
c     Expandindo para massa grande (a memos de in�rcia rotacional (n�o incluida)
      do ir=1,notl(j)
         ib=(ir-1)*3
         do ik=1,2
            igb=ib+ik
            do jm=1,notl(j)
               jb=(jm-1)*3+ik
               rmg(igb,jb,j)=rmp(ir,jm,j)*ro(j)
            enddo
         enddo
      enddo
      
c     Calculando rmp para usar na carga (sem densidade)
      do ir=1,notl(j)
         do jm=1,notl(j)
            rmp(ir,jm,j)=0.d0
         enddo
      enddo

      do ig1=1,ng               !A integral � apenas na longitudinal
         xsi1=qsi(ig1)
c     call suportefuncoesdeforma
         call formaederi0  
         igf=1                  !so para usar calc_a0    
         ig3=1                  !so para usar calc_a0    
         call calc_A0           !matriz gradiente do adimensional para inicial      
         do ir=1,notl(j)
            do jm=1,notl(j)
               rmp(ir,jm,j)=rmp(ir,jm,j)+fi(ir)*fi(jm)*w(ig1)*rjac1
            enddo
         enddo
      enddo                     !comprimento
      
      return
      end


c      *******************************************
c      subrotina que introduz renumeracao para rotulas
c     *******************************************   
      subroutine rotular
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     primeira etapa: alocar algumas vari�veis para calcular outras
      ALLOCATE (ngl(nnos)) 
      ALLOCATE (iglv(nnos,nel))
      ALLOCATE (glgl(nnos,17))
      ALLOCATE (itr(nel,30))            
      ALLOCATE (kn(nnos))
      ALLOCATE (im(nnos,17))
      ALLOCATE (ic(nel,30))
      
      glgl=0
      iglv=0
      ngl=0
      itr=0
      
      kn=0
      im=0
      ic=0
      
      read (3,*)
      read (3,*)	  
c     incidenca dos elementos
c     elemento, no1, no2,...no6
      write(5,*) 'incidenca dos elementos'
      write(5,*) 'elemento, no1, no2,...no6'

      do j=1,nel !La�o para todos os elementos
       read(3,*) k,(ic(k,i),i=1,notl(k))
       write(5,*) k,(ic(k,i),i=1,notl(k))

c     Novidade incid�ncia inversa para montagem r�pida ma Matriz-vetor MA27
c     A quais elementos pertence o no?
       do i=1,notl(k)           ! notl(j)=numero de nos do elemento
c       kn(nnos) informa quantos elementos estao conectados no n�      
        kn(ic(k,i))=kn(ic(k,i))+1 !k elemento i � o no local e ic e a incidencia
c       incidencia inversa im(no-global,numero do elento)
        im(ic(k,i),kn(ic(k,i)))=k
      enddo
      enddo                     !fim do la�o para todos os elementos


c     Informamos as r�tulas atrav�s dos n�s de extremidade dos elementos
      read (3,*)
      read (3,*)	  
c     Tipo do no de extremidade
c     Elemento, noi, tipo, nof, tipo ! 0=restrito 1=livre (rotulado)
      write(5,*) 'tipo do no de extremidade'
      write(5,*) 'elemento, tipo do noi, tipo do nof' ! 0=restrito 1=livre'
      do j=1,nel
       read(3,*) k,itr(k,1),itr(k,notl(k)) !itr(elemento,no-local)
       write(5,*) k,itr(k,1),itr(k,notl(k))
      enddo

      do ib=1,nnos
         do it=1,kn(ib)
            if(ic(im(ib,it),1).eq.ib) then
               marca(it)=itr(im(ib,it),1)
            else
               marca(it)=itr(im(ib,it),notl(im(ib,it)))
            endif
         enddo
         
         iicont=0
         do it=1,kn(ib)
            if(marca(it).eq.0) then
               iicont=iicont+1
               nordem(iicont)=im(ib,it)
            endif
         enddo
         
         do it=1,kn(ib)
            if(marca(it).eq.1) then
               iicont=iicont+1
               nordem(iicont)=im(ib,it)
            endif
         enddo
         
         do it=1,kn(ib)
            im(ib,it)=nordem(it)
         enddo
         
      enddo
      
c     Informando a incid�ncia inversa
      write(5,*) 'no,numero de elementos a que pertence'
      do i=1,nnos
         write(5,*) i,kn(i)
      enddo
      write(5,*) 'no,elemeto 1 , elemento 2 ...'
      do i=1,nnos
         write(5,*) i,(im(i,j),j=1,kn(i))
      enddo
      


c     In�cio do processo de defini��o da conectividade entre n�s 
c     e graus de liberdade.
c     Um n� que pertence a mais de um emento possuir� duas transla��es
c     Por�m poder� asumir mais de um giro
c     Assim, se pertence apenas a um elemento possui tr�s graus de liberdade
c     Se pertence a dois elementos e no segundo
c     elemento � rotulado ter� quatro graus de liberdade
c     Duas transla��es e dois giros
c     Caso Perten�a a tr�s elementos e for rotulado apenas em um possui 
c     quatro gls. Caso perten�a a tr�s elementos e for rotulado em dois
c      ter� 5 gls. E assim sucessivamente
 
      do ib=1,nnos              !La�o sobre os n�s
         iloc=2                 !duas transla��es impostas
         do jj=1,kn(ib)         !La�o sobre os elementos conectados ao n� jj = elemento local
            j=im(ib,jj)         !j � o elemento global
            do ie=1,notl(j)     !La�o sobre os n�s do elemento em quest�o
               if (ic(j,ie).eq.ib) then !Descobrindo o n� local correspondente ao n� ib investigado
                  i=ie          !no local 
                  goto 100
               else
               endif
            enddo               !Fim do la�o sobre os n�s do elemento em quest�o      
 100        continue
            iloc=iloc+itr(j,i)  !Somando o grau de liberdade de giro para barras rotuladas
            if(jj.eq.1) iloc=3  !se � aprimeira barra esse grau � o 3 sempre
            iglv(ib,j)=iloc     ! Registrando o n�mero local do grau de liberdade de giro para o n� global e elemento global
         enddo                  !Fim do la�o sobre os elementos conectados ao n� jj = elemento local
         ngl(ib)=iloc           !registrando o n�mero de graus de liberdade para cada n�
      enddo                     !Fim do la�o sobre o n�s
      
      gl0=0.d0                    !contador base de graus de liberdade
      do ib=1,nnos              !La�o sobre os n�s da estrutura
         glgl(ib,1)=gl0+1       !registrando transla��o 1
         glgl(ib,2)=gl0+2       !registrando transla��o 2
         do jj=3,ngl(ib)        !La�o sobre o n�mero giros
            glgl(ib,jj)=glgl(ib,jj-1)+1 !registrando giro jj
         enddo
         gl0=glgl(ib,ngl(ib))   !atualizando o contador base para cada n�
      enddo                     !fim do la�o sobre os n�s da estrutura
      n=gl0                     !registrando o n�mero de graus de liberdade total do problema
      return
      end
      
c      *******************************************
c      Subrotina que ajeita as secoes transversais
c      O usu�rio de engenharia b�sica entrar� com
c      dados de �rea, momento de in�rcia, posicao 
c      do cg e altura da secao e o programa 
c      transforma esses dados em tr�s faixas gen�ricas.
c      em vers�o posterior um usu�rio mais avan�ado
c      poder� informar caracter�sticas para cada fibra
c      de sua se��o transversal. Esta op��o estar� dispon�vel
c      em uma segunda vers�o do software.
c      *******************************************   
      subroutine secoes
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      do j=1,nel
         hh=resp(j)             !altura total
         aa=ra(j)               !�rea total
         rii=rib(j)             !Momento de in�rcia total
         yyb=yb(j)              !posi��o do cg (da se��o total)
         rde=dsqrt(aa*(432.d0*rii+aa*(-47.d0*hh*hh+108.d0*hh*yyb-108.d0
     $        *yyb*yyb)))
         hd1=-(aa*(hh-18.d0*yyb)+rde)/12.d0/aa
         hd2=-(aa*(hh-18.d0*yyb)-rde)/12.d0/aa
         if((hd1.gt.0.).and.(hd1.lt.hh)) then
            hb(j,1)=hd1         !altura da fibra 1
         else
            hb(j,1)=hd2
            if((hd2.lt.0.d0).or.(hd2.gt.hh)) then
               write(*,*) 'secao com problema para o elemento',j
            endif
         endif
         hb(j,2)=3.d0*yyb-hh/2.d0-2.d0*hb(j,1) !altura da fibra 2
         hb(j,3)=hh-hb(j,1)-hb(j,2) !altura da fibra 3
         brd(j,1)=aa/hb(j,1)/3.d0 !base da fibra 1
         brd(j,2)=aa/hb(j,2)/3.d0 !base da fibra 2
         brd(j,3)=aa/hb(j,3)/3.d0 !base da fibra 3
         y1b=hb(j,1)/2.d0 
         y2b=hb(j,1)+hb(j,2)/2.d0
         y3b=hb(j,1)+hb(j,2)+hb(j,3)/2.d0
         fd(j,1)=y1b-yyb        !posi��o do cg da fibra 1
         fd(j,2)=y2b-yyb        !posi��o do cg da fibra 2
         fd(j,3)=y3b-yyb        !posi��o do cg da fibra 1
      enddo
      return
      end
      
c      *******************************************
c      Subrotina que transforma as cargas distribuidas conservativas
c      em concentradas equivalentes
c     *******************************************   
      subroutine alocadistr
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      do j=1,nel
         do i=1,notl(j)
            igl=ic(j,i)
            do k=1,notl(j)
               df(glgl(igl,1))=df(glgl(igl,1))+rmp(i,k,j)*qx(j,k)
               df(glgl(igl,2))=df(glgl(igl,2))+rmp(i,k,j)*qy(j,k)
            enddo
         enddo
      enddo
      
      return
      end


c      *******************************************
c      subrotina que transforma as cargas distribuidas
c      n�o conservativas em nodais equivalentes
c      *******************************************   
      subroutine alocanco
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
c     O comportamento temporal da carga ainda segue o mesmo padr�o
c     da carga consevativa, pode ser melhorado em vers�o futura
      raux= (ca+cb*ipc*dt+cc*(ipc*dt)*(ipc*dt)+cd*dsin(ce*ipc*dt)+cf
     $     *dcos(cg*ipc*dt) +crh*exp(cri*ipc*dt)+cj*exp(ck*ipc*dt))
      
      do j=1,nel
         do i=1,notl(j)
            igl=ic(j,i)
            do k=1,notl(j)
               if((ipc.eq.1).and.(ia.eq.0))then
                  f(glgl(igl,1))=f(glgl(igl,1))+rmp(i,k,j)
     #*(qt(j,k)*sin(pl0(j,3,k))+qn(j,k)*cos(pl0(j,3,k)))*raux
                  f(glgl(igl,2))=f(glgl(igl,2))+rmp(i,k,j)
     #*(qn(j,k)*sin(pl0(j,3,k))-qt(j,k)*cos(pl0(j,3,k)))*raux
               else
                  f(glgl(igl,1))=f(glgl(igl,1))+rmp(i,k,j)
     #*(qt(j,k)*sin(pl(j,3,k))+qn(j,k)*cos(pl(j,3,k)))*raux
                  f(glgl(igl,2))=f(glgl(igl,2))+rmp(i,k,j)
     #*(qn(j,k)*sin(pl(j,3,k))-qt(j,k)*cos(pl(j,3,k)))*raux
               endif
            enddo
         enddo
      enddo
      
      return
      end

c     *****************************************************
c     Subrotina que faz aleitura de dados para plasticidade
c     *****************************************************
      subroutine lerplast
      USE BIG_STUFF
	implicit real*8(a-h),integer*4(i-n),real*8(o-z)
	include 'bloclagatbig.for'

c     Leitura dos par�metros de resist�ncia
c     nos moldes do modelo constitutivo de Drucke-Prager
c     com possibilidade de caracter�sticas anisotr�picas 
c     Cone n�o sim�trico
c     Esta subrotina foi deixada para ser usada
c     na vers�o de casca a ser implementada

      read(4,*) ! matriz dos parametros de resistencia
      read(4,*) rpret(1,1),rpret(1,2),rpret(1,3),rpret(1,4)
     #,rpret(1,5),rpret(1,6)
      read(4,*) rpret(2,1),rpret(2,2),rpret(2,3),rpret(2,4)
     #,rpret(2,5),rpret(2,6)
      read(4,*) rpret(3,1),rpret(3,2),rpret(3,3),rpret(3,4)
     #,rpret(3,5),rpret(3,6)
      read(4,*) rpret(4,1),rpret(4,2),rpret(4,3),rpret(4,4)
     #,rpret(4,5),rpret(4,6)
      read(4,*) rpret(5,1),rpret(5,2),rpret(5,3),rpret(5,4)
     #,rpret(5,5),rpret(5,6)
      read(4,*) rpret(6,1),rpret(6,2),rpret(6,3),rpret(6,4)
     #,rpret(6,5),rpret(6,6)

      rpret=rpret/2.

      read(4,*) ! matriz dos par�metros de abertura do cone
      read(4,*) rl(1),rl(2),rl(3),rl(4),rl(5),rl(6)
      read(4,*) !numero de trechs  adicionais nL
c               (o primeiro trecho segue comportamento elastopl�stico comum)
	read(4,*) nt
      
c      Aloca��o din�mica de termos de plasticidade
c      observar que como pertence ao big stuff est� presente em todo o 
c      programa como vari�vel global
	ALLOCATE (defl(2,nt+2))
	ALLOCATE (rhdg(2,nt+2))
	ALLOCATE (rhvg(2,nt+2))
      defl=0.
	rhdg=0.
	rhvg=0.
c     Novidade (Endurecimento volum�trico e desviador)
      read(4,*) ! Resist�ncia, hardening vol, hard dev (tracao e compressao)
	read(4,*) st,sc,rhvg(1,1),rhdg(1,1),rhvg(2,1),rhdg(2,1)

	do it=1,nt !n�mero de trechos n�o lineares al�m do primeiro
       read(4,*) ! hardeningv hardeningd def (tracao e compressao)
	 read(4,*) rhvg(1,it+1),rhdg(1,it+1),defl(1,it)
     #,rhvg(2,it+1),rhdg(2,it+1),defl(2,it)
      enddo

      defl(1,nt+1)=defl(1,nt)
      defl(2,nt+1)=defl(2,nt)
      rrmi=(sc-st)/sqrt(3.)/(sc+st)
      rl=rl*rrmi
	s0=st*(1./sqrt(3.)+rrmi)
	defl=defl*(1./sqrt(3.)+rrmi)

      return 
	end

c     **************************************************************
c     Subrotina de retorno para calculo da tens�o real (plasticidade)
c     *************************************************************
      subroutine retorno
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      rtest=et(1)+et(2)+et(3)   !vari�vel para testar trecho do cone
c     (tra��o ou compress�o)
      imret=1
      
      if (rtest.lt.0.) imret=2 
      
      
      if (ipc.eq.1) then        !Inicializando o endurecimento
         rhv(1)=rhvg(imret,1)
         rhv(2)=rhvg(imret,2)
         rhd(1)=rhdg(imret,1)
         rhd(2)=rhdg(imret,2)
      endif
      
c     Calculo da deforma��o total equivalente
c     usando  norma do modelo constitutivo anisotr�pico
      
      rpe=0.
      do i=1,6
         do jr=1,6
            rpe(i)=rpe(i)+rpret(i,jr)*et(jr)
         enddo
      enddo
      
      rde=0.
      do i=1,6
         rde=rde+rpe(i)*et(i)
      enddo
      
      rle=0.
      do i=1,6
         rle=rle+rl(i)*et(i)
      enddo
      
      eeqa=eeq(j,ig1,ig3,igf)   !deforma��o equivalente anterior
      eeq(j,ig1,ig3,igf)=dsqrt(rde)+rle !deforma��o equivalente atual
      eeq1=eeq(j,ig1,ig3,igf) 
      
      auxesc=1.
      na=1

      do itret=1,nt 
C      estabelecendo mudan�a de hardening nos cantos
c      mudan�a de trecho de hardening lineares
       if((eeq1.gt.defl(imret,itret)).and.
     #(eeqa.lt.defl(imret,itret))) then
          rhv(1)=rhvg(imret,itret+1*0)
          rhv(2)=rhvg(imret,itret+1)
          rhd(1)=rhdg(imret,itret+1*0)
          rhd(2)=rhdg(imret,itret+1)
          na=2
          eeq2=eeq1
          auxesc=defl(imret,itret)/eeq1	 
          eeq1=auxesc*eeq1
          et=auxesc*et
       else
       endif
       
       if((eeq1.gt.defl(imret,itret)).and.
     #(eeqa.gt.defl(imret,itret))) then
          na=1
          rhv(1)=rhvg(imret,itret+1) 
          rhv(2)=rhvg(imret,itret+1)
          rhd(1)=rhdg(imret,itret+1)
          rhd(2)=rhdg(imret,itret+1)
       else
       endif
      enddo

      
      do iaret=1,na             !caso n�o haja mudan�a de hardening na=1
         dlam=0.                !zerando multiplicador pl�stico
         if(iaret.eq.2) et=et/auxesc
         prsi=0.                !zerando for�a interna (conjugado da deforma��o equivalete)
         do iep=1,6
            ee(iep)=et(iep)-ep(j,ig1,ig3,igf,iep)
         enddo
         s=0.
         do i=1,6
            do jr=1,6
               s(i)=s(i)+ceret(i,jr)*ee(jr)
            enddo
         enddo
         rls=0.
         do i=1,6
            rls=rls+rl(i)*s(i)
         enddo
         f0=rls-s0-prsi0(j,ig1,ig3,igf)
         if(f0.gt.0.) then      !teste de degenera��o da superf�cie pl�stica
            f0=0.
            rhv(1)=0. 
            rhv(2)=0.
            rhd(1)=0.
            rhd(2)=0.
         else
         endif
         rhret=0.
         do i=1,6
            do jr=1,6
               rhret(i)=rhret(i)+dpret(i,jr)*s(jr)
            enddo
         enddo
         rt=0.
         do i=1,6
            do jr=1,6
               rt(i)=rt(i)+ceret(i,jr)*rhret(jr)
            enddo
         enddo
         rpt=0.
         do i=1,6
            do jr=1,6
               rpt(i)=rpt(i)+rpret(i,jr)*rt(jr)
            enddo
         enddo
         reret=0.
         do i=1,6
            reret=reret+rpt(i)*s(i)
         enddo
         rph=0.
         do i=1,6
            do jr=1,6
               rph(i)=rph(i)+rpret(i,jr)*rhret(jr)
            enddo
         enddo
         RD1=0.
         do i=1,6
            RD1=RD1+rph(i)*rhret(i)
         enddo
         rps=0.
         do i=1,6
            do jr=1,6
               rps(i)=rps(i)+rpret(i,jr)*s(jr)
            enddo
         enddo
         RD=0.
         do i=1,6
            RD=RD+rps(i)*s(i)
         enddo
         rfret=0.
         do i=1,6
            rfret=rfret+rps(i)*rt(i)
         enddo
         rg1=0.
         do i=1,6
            rg1=rg1+rpt(i)*rt(i)
         enddo
         rlh=0.
         do i=1,6
            rlh=rlh+rl(i)*rhret(i)
         enddo
         rbret=rhd(iaret)*dsqrt(rd1)+rhv(iaret)*rlh
         raret=0.
         do i=1,6
            raret=raret+rl(i)*rt(i)
         enddo
         rcret=raret+rbret
         seq=dsqrt(rd)+rls
         fret=f0+dsqrt(rd)      !c�lculo do valor da fun��o objetiva
         alfret=0.
         betret=0.
         gamret=0.
         delt=0.
         if(fret.gt.1.0d-10) then !violou o crit�rio de plastifica��o?
            alfret=rg1-rcret**2  
            betret=2.*rcret*f0-reret-rfret
            gamret=rd-f0**2.
            delt=betret**2.-4.*alfret*gamret
            if(dabs(alfret).gt.1.d-6) then
               dl1=(-betret+dsqrt(delt))/2./alfret
               dl2=(-betret-dsqrt(delt))/2./alfret
c     caso delta seja nulo o bascara resulta em uma unica raiz dada por (cuidado com sinal alf e gamret)
               if (dabs(delt).lt.1.0d-11) then
                  dl1=dsqrt(gamret)/dsqrt(alfret)
                  dl2=dl1
               endif
            else
               dl1=-gamret/betret
               dl2=dl1
            endif
            dlam=0.
            if ((dl1.gt.0.).and.(dl1.lt.dl2)) then
               dlam=dl1         !multiplicador pl�stico
            else
               dlam=dl2         !multiplicador pl�stico
            endif
            if((dlam).lt.1.d-15) then
               dlam=0.
            else
               prsi=dlam*rbret  !varia��o da for�a interna
               do iep=1,6
                  ep(j,ig1,ig3,igf,iep)=ep(j,ig1,ig3,igf,iep)+dlam
     $                 *rhret(iep)
!     deforma��o pl�stica
               enddo
            endif
            prsi0(j,ig1,ig3,igf)=prsi0(j,ig1,ig3,igf)+prsi !atualiza��o da for�a interna
            do iep=1,6
               ee(iep)=et(iep)-ep(j,ig1,ig3,igf,iep) !deforma��o el�stica
            enddo
c     goto 10 !fim das itera��es se o programa 
c     de retorno fosse desacoplado de um problema mec�nicao
         else
         endif
         if(iaret.eq.2) eeq1=eeq2 
      enddo
      
      return
      end
      


c      *************************************
c       Resolve o sistema subrotina de conhecimento p�blico
c      *************************************
      subroutine ma27
      USE BIG_STUFF	 
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      iw=0
      iw1=0
      icntl=0
      cntl=0
      ikeep=0
      info=0
      iflag=0
      nsteps=0
      maxfrt=0
      w27=0
      
      rhs=v
c      write(*,*)'rhs'
c      write(*,*)rhs
C     
      
C     INITIALIZE MA27
      CALL MA27ID(ICNTL,CNTL)
      
C     
C     STORE ARRAY LENGTHS
      
      
c     nzexp1=Int((2*nz+3*n+1)*1.2)
      
      LIW=nzexp1
      LA=nzexp
C     
      
C     ASK FOR FULL PRINTING FROM MA27 PACKAGE
      ICNTL(3) = 0
C     
C     SET IFLAG TO INDICATE PIVOT SEQUENCE IS TO BE FOUND BY MA27AD
      IFLAG = 0
C     
      
      
      
      
C     READ MATRIX AND RIGHT-HAND SIDE
c     READ (5,*) N,NZ
c     READ (5,*) (IRN(I),ICN(I),A(I),I=1,NZ)
c     READ (5,*) (RHS(I),I=1,N)
C     
C     ANALYSE SPARSITY PATTERN
      CALL MA27AD(N27,NZ,IRN,ICN,IW,LIW,IKEEP,IW1,NSTEPS,IFLAG,
     +     ICNTL,CNTL,INFO,OPS)
      
      
      
      
      
      
C     FACTORIZE MATRIX
      CALL MA27BD(N27,NZ,IRN,ICN,A,LA,IW,LIW,IKEEP,NSTEPS,MAXFRT,
     *     IW1,ICNTL,CNTL,INFO)
      
      
C     
C     SOLVE THE EQUATIONS
      CALL MA27CD(N27,A,LA,IW,LIW,W27,MAXFRT,RHS,IW1,NSTEPS,ICNTL,
     +     INFO)
      v=rhs
c     write(*,*)'solution'
c      write(*,*)rhs
c      read(*,*)
      return
      END

      
* *******************************************************************
* COPYRIGHT (c) 1982 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 20 September 2001
C  September 2001: threadsafe version of MA27
C  19/3/03. Array ICNTL in MA27GD made assumed size.

      SUBROUTINE MA27ID(ICNTL,CNTL)
      INTEGER ICNTL(30)
      DOUBLE PRECISION CNTL(5)
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
      ICNTL(1) = 66
      ICNTL(2) = 66
      ICNTL(3) = 0
      ICNTL(4) = 2139062143
      ICNTL(5) = 1
      ICNTL(IFRLVL+1)  = 32639
      ICNTL(IFRLVL+2)  = 32639
      ICNTL(IFRLVL+3)  = 32639
      ICNTL(IFRLVL+4)  = 32639
      ICNTL(IFRLVL+5)  = 14
      ICNTL(IFRLVL+6)  = 9
      ICNTL(IFRLVL+7)  = 8
      ICNTL(IFRLVL+8)  = 8
      ICNTL(IFRLVL+9)  = 9
      ICNTL(IFRLVL+10) = 10
      ICNTL(IFRLVL+11) = 32639
      ICNTL(IFRLVL+12) = 32639
      ICNTL(IFRLVL+13) = 32639
      ICNTL(IFRLVL+14) = 32689
      ICNTL(IFRLVL+15) = 24
      ICNTL(IFRLVL+16) = 11
      ICNTL(IFRLVL+17) = 9
      ICNTL(IFRLVL+18) = 8
      ICNTL(IFRLVL+19) = 9
      ICNTL(IFRLVL+20) = 10
      ICNTL(26) = 0
      ICNTL(27) = 0
      ICNTL(28) = 0
      ICNTL(29) = 0
      ICNTL(30) = 0
      CNTL(1) = 0.1D0
      CNTL(2) = 1.0D0
      CNTL(3) = 0.0D0
      CNTL(4) = 0.0
      CNTL(5) = 0.0
      RETURN
      END
      SUBROUTINE MA27AD(N,NZ,IRN,ICN,IW,LIW,IKEEP,IW1,NSTEPS,IFLAG,
     +                 ICNTL,CNTL,INFO,OPS)
      INTEGER IFLAG,LIW,N,NSTEPS,NZ
      DOUBLE PRECISION CNTL(5),OPS
      INTEGER ICNTL(30),INFO(20)
      INTEGER ICN(*),IKEEP(N,3),IRN(*),IW(LIW),IW1(N,2)
      INTEGER I,IWFR,K,L1,L2,LLIW
      EXTERNAL MA27GD,MA27HD,MA27JD,MA27KD,MA27LD,MA27MD,MA27UD
      INTRINSIC MIN
      DO 5 I = 1,15
        INFO(I) = 0
    5 CONTINUE
      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 40
      WRITE (ICNTL(2),FMT=10) N,NZ,LIW,IFLAG
   10 FORMAT(/,/,' ENTERING MA27AD WITH      N     NZ      LIW  IFLAG',
     +       /,21X,I7,I7,I9,I7)
      NSTEPS = 0
      K = MIN(8,NZ)
      IF (ICNTL(3).GT.1) K = NZ
      IF (K.GT.0) WRITE (ICNTL(2),FMT=20) (IRN(I),ICN(I),I=1,K)
   20 FORMAT (' MATRIX NON-ZEROS',/,4 (I9,I6),/,
     +       (I9,I6,I9,I6,I9,I6,I9,I6))
      K = MIN(10,N)
      IF (ICNTL(3).GT.1) K = N
      IF (IFLAG.EQ.1 .AND. K.GT.0) THEN
        WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
      END IF
   30 FORMAT (' IKEEP(.,1)=',10I6,/, (12X,10I6))
   40 IF (N.LT.1 .OR. N.GT.ICNTL(4)) GO TO 70
      IF (NZ.LT.0) GO TO 100
      LLIW = LIW - 2*N
      L1 = LLIW + 1
      L2 = L1 + N
      IF (IFLAG.EQ.1) GO TO 50
      IF (LIW.LT.2*NZ+3*N+1) GO TO 130
      CALL MA27GD(N,NZ,IRN,ICN,IW,LLIW,IW1,IW1(1,2),IW(L1),IWFR,
     +           ICNTL,INFO)
      CALL MA27HD(N,IW1,IW,LLIW,IWFR,IW(L1),IW(L2),IKEEP(1,2),
     +           IKEEP(1,3),IKEEP,ICNTL(4),INFO(11),CNTL(2))
      GO TO 60
   50 IF (LIW.LT.NZ+3*N+1) GO TO 120
      CALL MA27JD(N,NZ,IRN,ICN,IKEEP,IW,LLIW,IW1,IW1(1,2),IW(L1),IWFR,
     +           ICNTL,INFO)
      CALL MA27KD(N,IW1,IW,LLIW,IWFR,IKEEP,IKEEP(1,2),IW(L1),IW(L2),
     +           INFO(11))
   60 CALL MA27LD(N,IW1,IW(L1),IKEEP,IKEEP(1,2),IKEEP(1,3),IW(L2),
     +           NSTEPS,ICNTL(5))
      IF(NZ.GE.1) IW(1) = IRN(1) + 1
      CALL MA27MD(N,NZ,IRN,ICN,IKEEP,IKEEP(1,3),IKEEP(1,2),IW(L2),
     +           NSTEPS,IW1,IW1(1,2),IW,INFO,OPS)
      GO TO 160
   70 INFO(1) = -1
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
   80 FORMAT (' **** ERROR RETURN FROM MA27AD **** INFO(1)=',I3)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=90) N
   90 FORMAT (' VALUE OF N OUT OF RANGE ... =',I10)
      GO TO 160
  100 INFO(1) = -2
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=110) NZ
  110 FORMAT (' VALUE OF NZ OUT OF RANGE .. =',I10)
      GO TO 160
  120 INFO(2) = NZ + 3*N + 1
      GO TO 140
  130 INFO(2) = 2*NZ + 3*N + 1
  140 INFO(1) = -3
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=150) LIW,INFO(2)
  150 FORMAT (' LIW TOO SMALL, MUST BE INCREASED FROM',I10,
     +       ' TO AT LEAST',I10)
  160 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0 .OR. INFO(1).LT.0) GO TO 200
      WRITE (ICNTL(2),FMT=170) NSTEPS,INFO(1),OPS,INFO(2),INFO(3),
     +  INFO(4),INFO(5),INFO(6),INFO(7),INFO(8),INFO(11)
  170 FORMAT (/,' LEAVING MA27AD WITH NSTEPS  INFO(1)    OPS IERROR',
     +          ' NRLTOT NIRTOT',
     +        /,20X,2I7,F7.0,3I7,
     +        /,20X,' NRLNEC NIRNEC NRLADU NIRADU  NCMPA',
     +        /,20X,6I7)
      K = MIN(9,N)
      IF (ICNTL(3).GT.1) K = N
      IF (K.GT.0) WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
      K = MIN(K,NSTEPS)
      IF (K.GT.0) WRITE (ICNTL(2),FMT=180) (IKEEP(I,2),I=1,K)
  180 FORMAT (' IKEEP(.,2)=',10I6,/, (12X,10I6))
      IF (K.GT.0) WRITE (ICNTL(2),FMT=190) (IKEEP(I,3),I=1,K)
  190 FORMAT (' IKEEP(.,3)=',10I6,/, (12X,10I6))
  200 CONTINUE
      RETURN
      END


      SUBROUTINE MA27BD(N,NZ,IRN,ICN,A,LA,IW,LIW,IKEEP,NSTEPS,MAXFRT,
     +                 IW1,ICNTL,CNTL,INFO)
      INTEGER LA,LIW,MAXFRT,N,NSTEPS,NZ
      DOUBLE PRECISION A(LA),CNTL(5)
      INTEGER ICN(*),IKEEP(N,3),IRN(*),IW(LIW),IW1(N)
      INTEGER ICNTL(30),INFO(20)
      INTEGER I,IAPOS,IBLK,IPOS,IROWS,J1,J2,JJ,K,KBLK,KZ,LEN,NCOLS,
     +        NROWS,NZ1
      EXTERNAL MA27ND,MA27OD
      INTRINSIC ABS,MIN

      INFO(1) = 0
      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 60
      WRITE (ICNTL(2),FMT=10) N,NZ,LA,LIW,NSTEPS,CNTL(1)
   10 FORMAT (/,/,
     + ' ENTERING MA27BD WITH      N     NZ       LA      LIW',
     +       ' NSTEPS      U',/,21X,I7,I7,I9,I9,I7,1PD10.2)
      KZ = MIN(6,NZ)

      IF (ICNTL(3).GT.1) KZ = NZ

      write(*,*) 'aqui',icntl(2),kz,la
      IF (NZ.GT.0) WRITE (ICNTL(2),FMT=20) (A(K),IRN(K),ICN(K),K=1,KZ)
   20 FORMAT (' MATRIX NON-ZEROS',/,1X,2 (1P,D16.3,2I6),/,
     +       (1X,1P,D16.3,2I6,1P,D16.3,2I6))


      K = MIN(9,N)
      IF (ICNTL(3).GT.1) K = N

      IF (K.GT.0) WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
   30 FORMAT (' IKEEP(.,1)=',10I6,/, (12X,10I6))
      K = MIN(K,NSTEPS)
      IF (K.GT.0) WRITE (ICNTL(2),FMT=40) (IKEEP(I,2),I=1,K)
   40 FORMAT (' IKEEP(.,2)=',10I6,/, (12X,10I6))
      IF (K.GT.0) WRITE (ICNTL(2),FMT=50) (IKEEP(I,3),I=1,K)
   50 FORMAT (' IKEEP(.,3)=',10I6,/, (12X,10I6))


   60 IF (N.LT.1 .OR. N.GT.ICNTL(4)) GO TO 70

      IF (NZ.LT.0) GO TO 100
      IF (LIW.LT.NZ) GO TO 120
      IF (LA.LT.NZ+N) GO TO 150
      IF (NSTEPS.LT.1 .OR. NSTEPS.GT.N) GO TO 175
      CALL MA27ND(N,NZ,NZ1,A,LA,IRN,ICN,IW,LIW,IKEEP,IW1,ICNTL,INFO)
      IF (INFO(1).EQ.-3) GO TO 130
      IF (INFO(1).EQ.-4) GO TO 160
      CALL MA27OD(N,NZ1,A,LA,IW,LIW,IKEEP,IKEEP(1,3),NSTEPS,MAXFRT,
     +           IKEEP(1,2),IW1,ICNTL,CNTL,INFO)
      IF (INFO(1).EQ.-3) GO TO 130
      IF (INFO(1).EQ.-4) GO TO 160
      IF (INFO(1).EQ.-5) GO TO 180
      IF (INFO(1).EQ.-6) GO TO 200
      IF (INFO(1).EQ.3 .AND. ICNTL(2).GT.0) THEN
        WRITE (ICNTL(2),FMT=65) INFO(1),INFO(2)
      END IF

   65 FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +        '  *** INFO(1) =',I2,
     +        /,5X,'MATRIX IS SINGULAR. RANK=',I5)
      GO TO 220
   70 INFO(1) = -1

      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
   80 FORMAT (' **** ERROR RETURN FROM MA27BD **** INFO(1)=',I3)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=90) N
   90 FORMAT (' VALUE OF N OUT OF RANGE ... =',I10)
      GO TO 220
  100 INFO(1) = -2
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=110) NZ
  110 FORMAT (' VALUE OF NZ OUT OF RANGE .. =',I10)
      GO TO 220
  120 INFO(1) = -3
      INFO(2) = NZ
  130 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=140) LIW,INFO(2)
  140 FORMAT (' LIW TOO SMALL, MUST BE INCREASED FROM',I10,' TO',
     +       ' AT LEAST',I10)
      GO TO 220
  150 INFO(1) = -4
      INFO(2) = NZ + N
  160 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=170) LA,INFO(2)
  170 FORMAT (' LA TOO SMALL, MUST BE INCREASED FROM ',I10,' TO',
     +       ' AT LEAST',I10)
      GO TO 220
  175 INFO(1) = -7
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) THEN
        WRITE (ICNTL(1),FMT='(A)') ' NSTEPS is out of range'
      END IF
      GO TO 220
  180 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=190) INFO(2)
  190 FORMAT (' ZERO PIVOT AT STAGE',I10,
     +        ' WHEN INPUT MATRIX DECLARED DEFINITE')
      GO TO 220
  200 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=210)
  210 FORMAT (' CHANGE IN SIGN OF PIVOT ENCOUNTERED',
     +        ' WHEN FACTORING ALLEGEDLY DEFINITE MATRIX')
  220 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0 .OR. INFO(1).LT.0) GO TO 310
      WRITE (ICNTL(2),FMT=230) MAXFRT,INFO(1),INFO(9),INFO(10),INFO(12),
     +  INFO(13),INFO(14),INFO(2)
  230 FORMAT (/,' LEAVING MA27BD WITH',
     +        /,10X,'  MAXFRT  INFO(1) NRLBDU NIRBDU NCMPBR',
     +         ' NCMPBI   NTWO IERROR',
     +        /,11X,8I7)
      IF (INFO(1).LT.0) GO TO 310
      KBLK = ABS(IW(1)+0)
      IF (KBLK.EQ.0) GO TO 310
      IF (ICNTL(3).EQ.1) KBLK = 1
      IPOS = 2
      IAPOS = 1
      DO 300 IBLK = 1,KBLK
        NCOLS = IW(IPOS)
        NROWS = IW(IPOS+1)
        J1 = IPOS + 2
        IF (NCOLS.GT.0) GO TO 240
        NCOLS = -NCOLS
        NROWS = 1
        J1 = J1 - 1
  240   WRITE (ICNTL(2),FMT=250) IBLK,NROWS,NCOLS
  250   FORMAT (' BLOCK PIVOT =',I8,' NROWS =',I8,' NCOLS =',I8)
        J2 = J1 + NCOLS - 1
        IPOS = J2 + 1
        WRITE (ICNTL(2),FMT=260) (IW(JJ),JJ=J1,J2)
  260   FORMAT (' COLUMN INDICES =',10I6,/, (17X,10I6))
        WRITE (ICNTL(2),FMT=270)
  270   FORMAT (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')
        LEN = NCOLS
        DO 290 IROWS = 1,NROWS
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (ICNTL(2),FMT=280) (A(JJ),JJ=J1,J2)
  280     FORMAT (1P,5D13.3)
          LEN = LEN - 1
          IAPOS = J2 + 1
  290   CONTINUE
  300 CONTINUE
  310 RETURN
      END
      SUBROUTINE MA27CD(N,A,LA,IW,LIW,W,MAXFRT,RHS,IW1,NSTEPS,
     + ICNTL,INFO)
      INTEGER LA,LIW,MAXFRT,N,NSTEPS
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFRT)
      INTEGER IW(LIW),IW1(NSTEPS),ICNTL(30),INFO(20)
      INTEGER I,IAPOS,IBLK,IPOS,IROWS,J1,J2,JJ,K,KBLK,LATOP,LEN,NBLK,
     +        NCOLS,NROWS
      EXTERNAL MA27QD,MA27RD
      INTRINSIC ABS,MIN
      INFO(1) = 0
      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 110
      WRITE (ICNTL(2),FMT=10) N,LA,LIW,MAXFRT,NSTEPS
   10 FORMAT (/,/,' ENTERING MA27CD WITH      N     LA    LIW MAXFRT',
     +       '  NSTEPS',/,21X,5I7)
      KBLK = ABS(IW(1)+0)
      IF (KBLK.EQ.0) GO TO 90
      IF (ICNTL(3).EQ.1) KBLK = 1
      IPOS = 2
      IAPOS = 1
      DO 80 IBLK = 1,KBLK
        NCOLS = IW(IPOS)
        NROWS = IW(IPOS+1)
        J1 = IPOS + 2
        IF (NCOLS.GT.0) GO TO 20
        NCOLS = -NCOLS
        NROWS = 1
        J1 = J1 - 1
   20   WRITE (ICNTL(2),FMT=30) IBLK,NROWS,NCOLS
   30   FORMAT (' BLOCK PIVOT =',I8,' NROWS =',I8,' NCOLS =',I8)
        J2 = J1 + NCOLS - 1
        IPOS = J2 + 1
        WRITE (ICNTL(2),FMT=40) (IW(JJ),JJ=J1,J2)
   40   FORMAT (' COLUMN INDICES =',10I6,/, (17X,10I6))
        WRITE (ICNTL(2),FMT=50)
   50   FORMAT (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')
        LEN = NCOLS
        DO 70 IROWS = 1,NROWS
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (ICNTL(2),FMT=60) (A(JJ),JJ=J1,J2)
   60     FORMAT (1P,5D13.3)
          LEN = LEN - 1
          IAPOS = J2 + 1
   70   CONTINUE
   80 CONTINUE
   90 K = MIN(10,N)
      IF (ICNTL(3).GT.1) K = N
      IF (N.GT.0) WRITE (ICNTL(2),FMT=100) (RHS(I),I=1,K)
  100 FORMAT (' RHS',1P,5D13.3,/, (4X,1P,5D13.3))
  110 IF (IW(1).LT.0) GO TO 130
      NBLK = IW(1)
      IF (NBLK.GT.0) GO TO 140
      DO 120 I = 1,N
        RHS(I) = 0.0D0
  120 CONTINUE
      GO TO 150
  130 NBLK = -IW(1)
  140 CALL MA27QD(N,A,LA,IW(2),LIW-1,W,MAXFRT,RHS,IW1,NBLK,LATOP,ICNTL)
      CALL MA27RD(N,A,LA,IW(2),LIW-1,W,MAXFRT,RHS,IW1,NBLK,LATOP,ICNTL)
  150 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 170
      WRITE (ICNTL(2),FMT=160)
  160 FORMAT (/,/,' LEAVING MA27CD WITH')
      IF (N.GT.0) WRITE (ICNTL(2),FMT=100) (RHS(I),I=1,K)
  170 CONTINUE
      RETURN
      END
      SUBROUTINE MA27GD(N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
      INTEGER IWFR,LW,N,NZ
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
      INTEGER ICNTL(*),INFO(20)
      INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LR,N1,NDUP
      INFO(2) = 0
      DO 10 I = 1,N
        IPE(I) = 0
   10 CONTINUE
      LR = NZ
      IF (NZ.EQ.0) GO TO 120
      DO 110 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 90
        ELSE IF (I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 90
        ELSE
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF
   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD',
     +          '  *** INFO(1) =',I2)
        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF
   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')
   80   I = 0
        J = 0
        GO TO 100
   90   IPE(I) = IPE(I) + 1
        IPE(J) = IPE(J) + 1
  100   IW(K) = J
        LR = LR + 1
        IW(LR) = I
  110 CONTINUE
  120 IQ(1) = 1
      N1 = N - 1
      IF (N1.LE.0) GO TO 140
      DO 130 I = 1,N1
        FLAG(I) = 0
        IF (IPE(I).EQ.0) IPE(I) = -1
        IQ(I+1) = IPE(I) + IQ(I) + 1
        IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
      FLAG(N) = 0
      IF (LR.GE.LAST) GO TO 160
      K1 = LR + 1
      DO 150 K = K1,LAST
        IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
      IWFR = LAST + 1
      IF (NZ.EQ.0) GO TO 230
      DO 220 K = 1,NZ
        J = IW(K)
        IF (J.LE.0) GO TO 220
        L = K
        IW(K) = 0
        DO 210 ID = 1,NZ
          IF (L.GT.NZ) GO TO 170
          L = L + NZ
          GO TO 180
  170     L = L - NZ
  180     I = IW(L)
          IW(L) = 0
          IF (I.LT.J) GO TO 190
          L = IQ(J) + 1
          IQ(J) = L
          JN = IW(L)
          IW(L) = -I
          GO TO 200
  190     L = IQ(I) + 1
          IQ(I) = L
          JN = IW(L)
          IW(L) = -J
  200     J = JN
          IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
  230 NDUP = 0
      DO 280 I = 1,N
        K1 = IPE(I) + 1
        K2 = IQ(I)
        IF (K1.LE.K2) GO TO 240
        IPE(I) = 0
        IQ(I) = 0
        GO TO 280
  240   DO 260 K = K1,K2
          J = -IW(K)
          IF (J.LE.0) GO TO 270
          L = IQ(J) + 1
          IQ(J) = L
          IW(L) = I
          IW(K) = J
          IF (FLAG(J).NE.I) GO TO 250
          NDUP = NDUP + 1
          IW(L) = 0
          IW(K) = 0
  250     FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
        IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
      IF (NDUP.EQ.0) GO TO 310
      IWFR = 1
      DO 300 I = 1,N
        K1 = IPE(I) + 1
        IF (K1.EQ.1) GO TO 300
        K2 = IQ(I) + IPE(I)
        L = IWFR
        IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 290 K = K1,K2
          IF (IW(K).EQ.0) GO TO 290
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
  290   CONTINUE
        IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN
      END
      SUBROUTINE MA27HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
     +                 NCMPA,FRATIO)
      DOUBLE PRECISION FRATIO
      INTEGER IWFR,LW,N,IOVFLO,NCMPA
      INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
      INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
     +        KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
     +        NP0,NS,NVPIV,NVROOT,ROOT
      EXTERNAL MA27UD
      INTRINSIC ABS,MIN
      DO 10 I = 1,N
        IPD(I) = 0
        NV(I) = 1
        FLAG(I) = IOVFLO
   10 CONTINUE
      MD = 1
      NCMPA = 0
      NFLG = IOVFLO
      NEL = 0
      ROOT = N+1
      NVROOT = 0
      DO 30 IS = 1,N
        K = IPE(IS)
        IF (K.LE.0) GO TO 20
        ID = IW(K) + 1
        NS = IPD(ID)
        IF (NS.GT.0) LST(NS) = IS
        NXT(IS) = NS
        IPD(ID) = IS
        LST(IS) = -ID
        GO TO 30
   20   NEL = NEL + 1
        FLAG(IS) = -1
        NXT(IS) = 0
        LST(IS) = 0
   30 CONTINUE
      DO 340 ML = 1,N
        IF (NEL+NVROOT+1.GE.N) GO TO 350
        DO 40 ID = MD,N
          MS = IPD(ID)
          IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
        NVPIV = NV(MS)
        NS = NXT(MS)
        NXT(MS) = 0
        LST(MS) = 0
        IF (NS.GT.0) LST(NS) = -ID
        IPD(ID) = NS
        ME = MS
        NEL = NEL + NVPIV
        IDN = 0
        KP = IPE(ME)
        FLAG(MS) = -1
        IP = IWFR
        LEN = IW(KP)
        DO 140 KP1 = 1,LEN
          KP = KP + 1
          KE = IW(KP)
          IF (FLAG(KE).LE.-2) GO TO 60
          IF (FLAG(KE).LE.0) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 140
             KE = ROOT
             IF (FLAG(KE).LE.0) GO TO 140
          END IF
          JP = KP - 1
          LN = LEN - KP1 + 1
          IE = MS
          GO TO 70
   60     IE = KE
          JP = IPE(IE)
          LN = IW(JP)
   70     DO 130 JP1 = 1,LN
            JP = JP + 1
            IS = IW(JP)
            IF (FLAG(IS).LE.0) THEN
               IF (IPE(IS).EQ.-ROOT) THEN
                  IS = ROOT
                  IW(JP) = ROOT
                  IF (FLAG(IS).LE.0) GO TO 130
               ELSE
                  GO TO 130
               END IF
            END IF
            FLAG(IS) = 0
            IF (IWFR.LT.LW) GO TO 100
            IPE(MS) = KP
            IW(KP) = LEN - KP1
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA27UD(N,IPE,IW,IP-1,LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 90
            DO 80 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   80       CONTINUE
   90       IP = LWFR
            JP = IPE(IE)
            KP = IPE(ME)
  100       IW(IWFR) = IS
            IDN = IDN + NV(IS)
            IWFR = IWFR + 1
            LS = LST(IS)
            LST(IS) = 0
            NS = NXT(IS)
            NXT(IS) = 0
            IF (NS.GT.0) LST(NS) = LS
            IF (LS.LT.0) THEN
              LS = -LS
              IPD(LS) = NS
            ELSE IF (LS.GT.0) THEN
              NXT(LS) = NS
            END IF
  130     CONTINUE
          IF (IE.EQ.MS) GO TO 150
          IPE(IE) = -ME
          FLAG(IE) = -1
  140   CONTINUE
  150   NV(MS) = IDN + NVPIV
        IF (IWFR.EQ.IP) GO TO 330
        K1 = IP
        K2 = IWFR - 1
        LIMIT = NINT(FRATIO*(N-NEL))
        DO 310 K = K1,K2
          IS = IW(K)
          IF (IS.EQ.ROOT) GO TO 310
          IF (NFLG.GT.2) GO TO 170
          DO 160 I = 1,N
            IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
            IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160     CONTINUE
          NFLG = IOVFLO
  170     NFLG = NFLG - 1
          ID = IDN
          KP1 = IPE(IS) + 1
          NP = KP1
          KP2 = IW(KP1-1) + KP1 - 1
          DO 220 KP = KP1,KP2
            KE = IW(KP)
          IF (FLAG(KE).EQ.-1) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 220
             KE = ROOT
             IW(KP) = ROOT
             IF (FLAG(KE).EQ.-1) GO TO 220
          END IF
          IF (FLAG(KE).GE.0) GO TO 230
            JP1 = IPE(KE) + 1
            JP2 = IW(JP1-1) + JP1 - 1
            IDL = ID
            DO 190 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).LE.NFLG) GO TO 190
              ID = ID + NV(JS)
              FLAG(JS) = NFLG
  190       CONTINUE
            IF (ID.GT.IDL) GO TO 210
            DO 200 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).NE.0) GO TO 210
  200       CONTINUE
            IPE(KE) = -ME
            FLAG(KE) = -1
            GO TO 220
  210       IW(NP) = KE
            FLAG(KE) = -NFLG
            NP = NP + 1
  220     CONTINUE
          NP0 = NP
          GO TO 250
  230     KP0 = KP
          NP0 = NP
          DO 240 KP = KP0,KP2
            KS = IW(KP)
            IF (FLAG(KS).LE.NFLG) THEN
               IF (IPE(KS).EQ.-ROOT) THEN
                  KS = ROOT
                  IW(KP) = ROOT
                  IF (FLAG(KS).LE.NFLG) GO TO 240
               ELSE
                  GO TO 240
               END IF
            END IF
            ID = ID + NV(KS)
            FLAG(KS) = NFLG
            IW(NP) = KS
            NP = NP + 1
  240     CONTINUE
  250     IF (ID.GE.LIMIT) GO TO 295
          IW(NP) = IW(NP0)
          IW(NP0) = IW(KP1)
          IW(KP1) = ME
          IW(KP1-1) = NP - KP1 + 1
          JS = IPD(ID)
          DO 280 L = 1,N
            IF (JS.LE.0) GO TO 300
            KP1 = IPE(JS) + 1
            IF (IW(KP1).NE.ME) GO TO 300
            KP2 = KP1 - 1 + IW(KP1-1)
            DO 260 KP = KP1,KP2
              IE = IW(KP)
              IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260       CONTINUE
            GO TO 290
  270       JS = NXT(JS)
  280     CONTINUE
  290     IPE(JS) = -IS
          NV(IS) = NV(IS) + NV(JS)
          NV(JS) = 0
          FLAG(JS) = -1
          NS = NXT(JS)
          LS = LST(JS)
          IF (NS.GT.0) LST(NS) = IS
          IF (LS.GT.0) NXT(LS) = IS
          LST(IS) = LS
          NXT(IS) = NS
          LST(JS) = 0
          NXT(JS) = 0
          IF (IPD(ID).EQ.JS) IPD(ID) = IS
          GO TO 310
  295     IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IW(K) = ROOT
            IPE(IS) = -ROOT
            NV(ROOT) = NV(ROOT) + NV(IS)
            NV(IS) = 0
            FLAG(IS) = -1
          END IF
          NVROOT = NV(ROOT)
          GO TO 310
  300     NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
          MD = MIN(MD,ID)
  310   CONTINUE
        DO 320 K = K1,K2
          IS = IW(K)
          IF (NV(IS).EQ.0) GO TO 320
          FLAG(IS) = NFLG
          IW(IP) = IS
          IP = IP + 1
  320   CONTINUE
        IWFR = K1
        FLAG(ME) = -NFLG
        IW(IP) = IW(K1)
        IW(K1) = IP - K1
        IPE(ME) = K1
        IWFR = IP + 1
        GO TO 335
  330   IPE(ME) = 0
  335   CONTINUE
  340 CONTINUE
  350 DO 360 IS = 1,N
        IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
          IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IPE(IS) = -ROOT
          END IF
          NVROOT = NVROOT + NV(IS)
          NV(IS) = 0
         END IF
  360 CONTINUE
      DO 370 IE = 1,N
        IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
      IF(NVROOT.GT.0)NV(ROOT)=NVROOT
      END
      SUBROUTINE MA27UD(N,IPE,IW,LW,IWFR,NCMPA)
      INTEGER IWFR,LW,N,NCMPA
      INTEGER IPE(N),IW(LW)
      INTEGER I,IR,K,K1,K2,LWFR
      NCMPA = NCMPA + 1
      DO 10 I = 1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR = 1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K = LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K = K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
      SUBROUTINE MA27JD(N,NZ,IRN,ICN,PERM,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
      INTEGER IWFR,LW,N,NZ
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW),PERM(N)
      INTEGER ICNTL(30),INFO(20)
      INTEGER I,ID,IN,J,JDUMMY,K,K1,K2,L,LBIG,LEN
      INTRINSIC MAX
      INFO(1) = 0
      INFO(2) = 0
      DO 10 I = 1,N
        IQ(I) = 0
   10 CONTINUE
      IF (NZ.EQ.0) GO TO 110
      DO 100 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF(I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 80
        ELSE IF(I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 80
        ELSE
          IW(K) = 0
          IF (I.GE.1 .AND. I.LE.N) GO TO 100
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IW(K) = 0
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF
   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD',
     +          '  *** INFO(1) =',I2)
        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF
   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')
        GO TO 100
   80   IF (PERM(J).GT.PERM(I)) GO TO 90
        IQ(J) = IQ(J) + 1
        GO TO 100
   90   IQ(I) = IQ(I) + 1
  100 CONTINUE
  110 IWFR = 1
      LBIG = 0
      DO 120 I = 1,N
        L = IQ(I)
        LBIG = MAX(L,LBIG)
        IWFR = IWFR + L
        IPE(I) = IWFR - 1
  120 CONTINUE
      IF (NZ.EQ.0) GO TO 250
      DO 160 K = 1,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 160
        L = K
        IW(K) = 0
        DO 150 ID = 1,NZ
          J = ICN(L)
          IF (PERM(I).LT.PERM(J)) GO TO 130
          L = IPE(J)
          IPE(J) = L - 1
          IN = IW(L)
          IW(L) = I
          GO TO 140
  130     L = IPE(I)
          IPE(I) = L - 1
          IN = IW(L)
          IW(L) = J
  140     I = -IN
          IF (I.LE.0) GO TO 160
  150   CONTINUE
  160 CONTINUE
      K = IWFR - 1
      L = K + N
      IWFR = L + 1
      DO 190 I = 1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 180
        DO 170 JDUMMY = 1,LEN
          IW(L) = IW(K)
          K = K - 1
          L = L - 1
  170   CONTINUE
  180   IPE(J) = L
        L = L - 1
  190 CONTINUE
      IF (LBIG.GE.ICNTL(4)) GO TO 210
      DO 200 I = 1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0
  200 CONTINUE
      GO TO 250
  210 IWFR = 1
      DO 240 I = 1,N
        K1 = IPE(I) + 1
        K2 = IPE(I) + IQ(I)
        IF (K1.LE.K2) GO TO 220
        IPE(I) = 0
        GO TO 240
  220   IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 230 K = K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 230
          IW(IWFR) = J
          IWFR = IWFR + 1
          FLAG(J) = I
  230   CONTINUE
        K = IPE(I)
        IW(K) = IWFR - K - 1
  240 CONTINUE
  250 RETURN
      END
      SUBROUTINE MA27KD(N,IPE,IW,LW,IWFR,IPS,IPV,NV,FLAG,NCMPA)
      INTEGER IWFR,LW,N,NCMPA
      INTEGER FLAG(N),IPE(N),IPS(N),IPV(N),IW(LW),NV(N)
      INTEGER I,IE,IP,J,JE,JP,JP1,JP2,JS,KDUMMY,LN,LWFR,ME,MINJS,ML,MS
      EXTERNAL MA27UD
      INTRINSIC MIN
      DO 10 I = 1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML = 1,N
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY = 1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
          DO 50 JP1 = 1,LN
            JP = JP + 1
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA27UD(N,IPE,IW,IP-1,LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = MIN(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
   60     IPE(IE) = -ME
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      RETURN
      END
      SUBROUTINE MA27LD(N,IPE,NV,IPS,NE,NA,ND,NSTEPS,NEMIN)
      INTEGER N,NSTEPS,NEMIN
      INTEGER IPE(N),IPS(N),NA(N),ND(N),NE(N),NV(N)
      INTEGER I,IB,IF,IL,IS,ISON,K,L,NR
      DO 10 I = 1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
      DO 20 I = 1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
      NR = N + 1
      DO 50 I = 1,N
        IF (NV(I).LE.0) GO TO 50
        IF = -IPE(I)
        IF (IF.EQ.0) GO TO 40
        IS = -IPS(IF)
        IF (IS.LE.0) GO TO 30
        IPE(I) = IS
   30   IPS(IF) = -I
        GO TO 50
   40   NR = NR - 1
        NE(NR) = I
   50 CONTINUE
      IS = 1
      I = 0
      DO 160 K = 1,N
        IF (I.GT.0) GO TO 60
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   DO 70 L = 1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
        IF (NV(I).LE.0) GO TO 120
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        IF (NA(IS).NE.1) GO TO 90
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
   90   IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NA(IS).EQ.0) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
          I = -IB
          IL = IL + 1
        END IF
  160 CONTINUE
      NSTEPS = IS - 1
      RETURN
      END
      SUBROUTINE MA27MD(N,NZ,IRN,ICN,PERM,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
     +                 IW,INFO,OPS)
      DOUBLE PRECISION OPS
      INTEGER N,NSTEPS,NZ
      INTEGER ICN(*),IRN(*),IW(*),LSTKI(N),LSTKR(N),NA(NSTEPS),
     +        ND(NSTEPS),NE(NSTEPS),PERM(N),INFO(20)
      INTEGER I,INEW,IOLD,IORG,IROW,ISTKI,ISTKR,ITOP,ITREE,JOLD,JORG,K,
     +        LSTK,NASSR,NELIM,NFR,NSTK,NUMORG,NZ1,NZ2
      DOUBLE PRECISION DELIM
      INTEGER NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC
      INTRINSIC MAX,MIN
      IF (NZ.EQ.0) GO TO 20
      IF (IRN(1).NE.IW(1)) GO TO 20
      IRN(1) = IW(1) - 1
      NZ2 = 0
      DO 10 IOLD = 1,N
        INEW = PERM(IOLD)
        LSTKI(INEW) = LSTKR(IOLD) + 1
        NZ2 = NZ2 + LSTKR(IOLD)
   10 CONTINUE
      NZ1 = NZ2/2 + N
      NZ2 = NZ2 + N
      GO TO 60
   20 DO 30 I = 1,N
        LSTKI(I) = 1
   30 CONTINUE
      NZ1 = N
      IF (NZ.EQ.0) GO TO 50
      DO 40 I = 1,NZ
        IOLD = IRN(I)
        JOLD = ICN(I)
        IF (IOLD.LT.1 .OR. IOLD.GT.N) GO TO 40
        IF (JOLD.LT.1 .OR. JOLD.GT.N) GO TO 40
        IF (IOLD.EQ.JOLD) GO TO 40
        NZ1 = NZ1 + 1
        IROW = MIN(PERM(IOLD)+0,PERM(JOLD)+0)
        LSTKI(IROW) = LSTKI(IROW) + 1
   40 CONTINUE
   50 NZ2 = NZ1
   60 ISTKI = 0
      ISTKR = 0
      OPS = 0.0D0
      NRLADU = 0
      NIRADU = 1
      NIRTOT = NZ1
      NRLTOT = NZ1
      NIRNEC = NZ2
      NRLNEC = NZ2
      NUMORG = 0
      ITOP = 0
      DO 100 ITREE = 1,NSTEPS
        NELIM = NE(ITREE)
        DELIM = NELIM
        NFR = ND(ITREE)
        NSTK = NA(ITREE)
        NASSR = NFR* (NFR+1)/2
        IF (NSTK.NE.0) NASSR = NASSR - LSTKR(ITOP) + 1
        NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
        NIRTOT = MAX(NIRTOT,NIRADU+NFR+2+ISTKI+NZ1)
        NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)
        NIRNEC = MAX(NIRNEC,NIRADU+NFR+2+ISTKI+NZ2)
        DO 70 IORG = 1,NELIM
          JORG = NUMORG + IORG
          NZ2 = NZ2 - LSTKI(JORG)
   70   CONTINUE
        NUMORG = NUMORG + NELIM
        IF (NSTK.LE.0) GO TO 90
        DO 80 K = 1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          LSTK = LSTKI(ITOP)
          ISTKI = ISTKI - LSTK
          ITOP = ITOP - 1
   80   CONTINUE
   90   NRLADU = NRLADU + (NELIM* (2*NFR-NELIM+1))/2
        NIRADU = NIRADU + 2 + NFR
        IF (NELIM.EQ.1) NIRADU = NIRADU - 1
        OPS = OPS + ((NFR*DELIM*(NFR+1)-(2*NFR+1)*DELIM*(DELIM+1)/2+
     +        DELIM* (DELIM+1)* (2*DELIM+1)/6)/2)
        IF (ITREE.EQ.NSTEPS) GO TO 100
        IF (NFR.EQ.NELIM) GO TO 100
        ITOP = ITOP + 1
        LSTKR(ITOP) = (NFR-NELIM)* (NFR-NELIM+1)/2
        LSTKI(ITOP) = NFR - NELIM + 1
        ISTKI = ISTKI + LSTKI(ITOP)
        ISTKR = ISTKR + LSTKR(ITOP)
        NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
        NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
  100 CONTINUE
      NRLNEC = MAX(NRLNEC,N+MAX(NZ,NZ1))
      NRLTOT = MAX(NRLTOT,N+MAX(NZ,NZ1))
      NRLNEC = MIN(NRLNEC,NRLTOT)
      NIRNEC = MAX(NZ,NIRNEC)
      NIRTOT = MAX(NZ,NIRTOT)
      NIRNEC = MIN(NIRNEC,NIRTOT)
      INFO(3) = NRLTOT
      INFO(4) = NIRTOT
      INFO(5) = NRLNEC
      INFO(6) = NIRNEC
      INFO(7) = NRLADU
      INFO(8) = NIRADU
      RETURN
      END
      SUBROUTINE MA27ND(N,NZ,NZ1,A,LA,IRN,ICN,IW,LIW,PERM,IW2,ICNTL,
     +                 INFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LA,LIW,N,NZ,NZ1
      DOUBLE PRECISION A(LA)
      INTEGER ICN(*),IRN(*),IW(LIW),IW2(N),PERM(N),ICNTL(30),INFO(20)
      DOUBLE PRECISION ANEXT,ANOW
      INTEGER I,IA,ICH,II,IIW,INEW,IOLD,IPOS,J1,J2,JJ,JNEW,JOLD,JPOS,K
      INTRINSIC MIN
      INFO(1) = 0
      IA = LA
      DO 10 IOLD = 1,N
        IW2(IOLD) = 1
        A(IA) = ZERO
        IA = IA - 1
   10 CONTINUE
      INFO(2) = 0
      NZ1 = N
      IF (NZ.EQ.0) GO TO 80
      DO 70 K = 1,NZ
        IOLD = IRN(K)
        IF (IOLD.GT.N .OR. IOLD.LE.0) GO TO 30
        JOLD = ICN(K)
        IF (JOLD.GT.N .OR. JOLD.LE.0) GO TO 30
        INEW = PERM(IOLD)
        JNEW = PERM(JOLD)
        IF (INEW.NE.JNEW) GO TO 20
        IA = LA - N + IOLD
        A(IA) = A(IA) + A(K)
        GO TO 60
   20   INEW = MIN(INEW,JNEW)
        IW2(INEW) = IW2(INEW) + 1
        IW(K) = -IOLD
        NZ1 = NZ1 + 1
        GO TO 70
   30   INFO(1) = 1
        INFO(2) = INFO(2) + 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=40) INFO(1)
        ENDIF
   40   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +          '  *** INFO(1) =',I2)
        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=50) K,IRN(K),ICN(K)
        END IF
   50   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')
   60   IW(K) = 0
   70 CONTINUE
   80 IF (NZ.LT.NZ1 .AND. NZ1.NE.N) GO TO 100
      K = 1
      DO 90 I = 1,N
        K = K + IW2(I)
        IW2(I) = K
   90 CONTINUE
      GO TO 120
  100 K = 1
      DO 110 I = 1,N
        K = K + IW2(I) - 1
        IW2(I) = K
  110 CONTINUE
  120 IF (NZ1.GT.LIW) GO TO 210
      IF (NZ1+N.GT.LA) GO TO 220
      IF (NZ1.EQ.N) GO TO 180
      DO 140 K = 1,NZ
        IOLD = -IW(K)
        IF (IOLD.LE.0) GO TO 140
        JOLD = ICN(K)
        ANOW = A(K)
        IW(K) = 0
        DO 130 ICH = 1,NZ
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          INEW = MIN(INEW,JNEW)
          IF (INEW.EQ.PERM(JOLD)) JOLD = IOLD
          JPOS = IW2(INEW) - 1
          IOLD = -IW(JPOS)
          ANEXT = A(JPOS)
          A(JPOS) = ANOW
          IW(JPOS) = JOLD
          IW2(INEW) = JPOS
          IF (IOLD.EQ.0) GO TO 140
          ANOW = ANEXT
          JOLD = ICN(JPOS)
  130   CONTINUE
  140 CONTINUE
      IF (NZ.GE.NZ1) GO TO 180
      IPOS = NZ1
      JPOS = NZ1 - N
      DO 170 II = 1,N
        I = N - II + 1
        J1 = IW2(I)
        J2 = JPOS
        IF (J1.GT.JPOS) GO TO 160
        DO 150 JJ = J1,J2
          IW(IPOS) = IW(JPOS)
          A(IPOS) = A(JPOS)
          IPOS = IPOS - 1
          JPOS = JPOS - 1
  150   CONTINUE
  160   IW2(I) = IPOS + 1
        IPOS = IPOS - 1
  170 CONTINUE
  180 DO 190 IOLD = 1,N
        INEW = PERM(IOLD)
        JPOS = IW2(INEW) - 1
        IA = LA - N + IOLD
        A(JPOS) = A(IA)
        IW(JPOS) = -IOLD
  190 CONTINUE
      IPOS = NZ1
      IA = LA
      IIW = LIW
      DO 200 I = 1,NZ1
        A(IA) = A(IPOS)
        IW(IIW) = IW(IPOS)
        IPOS = IPOS - 1
        IA = IA - 1
        IIW = IIW - 1
  200 CONTINUE
      GO TO 230
  210 INFO(1) = -3
      INFO(2) = NZ1
      GO TO 230
  220 INFO(1) = -4
      INFO(2) = NZ1 + N
  230 RETURN
      END
      SUBROUTINE MA27OD(N,NZ,A,LA,IW,LIW,PERM,NSTK,NSTEPS,MAXFRT,NELIM,
     +                 IW2,ICNTL,CNTL,INFO)
      DOUBLE PRECISION ZERO,HALF,ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
      INTEGER LA,LIW,MAXFRT,N,NSTEPS,NZ
      DOUBLE PRECISION A(LA),CNTL(5)
      INTEGER IW(LIW),IW2(N),NELIM(NSTEPS),NSTK(NSTEPS),PERM(N)
      INTEGER ICNTL(30),INFO(20)
      DOUBLE PRECISION AMAX,AMULT,AMULT1,AMULT2,DETPIV,RMAX,SWOP,
     +        THRESH,TMAX,UU
      INTEGER AINPUT,APOS,APOS1,APOS2,APOS3,ASTK,ASTK2,AZERO,I,IASS,
     +        IBEG,IDUMMY,IELL,IEND,IEXCH,IFR,IINPUT,IOLDPS,IORG,IPIV,
     +        IPMNP,IPOS,IROW,ISNPIV,ISTK,ISTK2,ISWOP,IWPOS,IX,IY,J,J1,
     +        J2,JCOL,JDUMMY,JFIRST,JJ,JJ1,JJJ,JLAST,JMAX,JMXMIP,JNEW,
     +        JNEXT,JPIV,JPOS,K,K1,K2,KDUMMY,KK,KMAX,KROW,LAELL,LAPOS2,
     +        LIELL,LNASS,LNPIV,LT,LTOPST,NASS,NBLK,NEWEL,NFRONT,NPIV,
     +        NPIVP1,NTOTPV,NUMASS,NUMORG,NUMSTK,PIVSIZ,POSFAC,POSPV1,
     +        POSPV2
      INTEGER NTWO,NEIG,NCMPBI,NCMPBR,NRLBDU,NIRBDU
      EXTERNAL MA27PD
      INTRINSIC ABS,MAX,MIN
      INTEGER IDIAG
      IDIAG(IX,IY) = ((IY-1)* (2*IX-IY+2))/2
      NBLK = 0
      NTWO = 0
      NEIG = 0
      NCMPBI = 0
      NCMPBR = 0
      MAXFRT = 0
      NRLBDU = 0
      NIRBDU = 0
      UU = MIN(CNTL(1),HALF)
      UU = MAX(UU,-HALF)
      DO 10 I = 1,N
        IW2(I) = 0
   10 CONTINUE
      IWPOS = 2
      POSFAC = 1
      ISTK = LIW - NZ + 1
      ISTK2 = ISTK - 1
      ASTK = LA - NZ + 1
      ASTK2 = ASTK - 1
      IINPUT = ISTK
      AINPUT = ASTK
      AZERO = 0
      NTOTPV = 0
      NUMASS = 0
      DO 760 IASS = 1,NSTEPS
        NASS = NELIM(IASS)
        NEWEL = IWPOS + 1
        JFIRST = N + 1
        NFRONT = 0
        NUMSTK = NSTK(IASS)
        LTOPST = 1
        LNASS = 0
        IF (NUMSTK.EQ.0) GO TO 80
        J2 = ISTK - 1
        LNASS = NASS
        LTOPST = ((IW(ISTK)+1)*IW(ISTK))/2
        DO 70 IELL = 1,NUMSTK
          JNEXT = JFIRST
          JLAST = N + 1
          J1 = J2 + 2
          J2 = J1 - 1 + IW(J1-1)
          DO 60 JJ = J1,J2
            J = IW(JJ)
            IF (IW2(J).GT.0) GO TO 60
            JNEW = PERM(J)
            IF (JNEW.LE.NUMASS) NASS = NASS + 1
            DO 20 IDUMMY = 1,N
              IF (JNEXT.EQ.N+1) GO TO 30
              IF (PERM(JNEXT).GT.JNEW) GO TO 30
              JLAST = JNEXT
              JNEXT = IW2(JLAST)
   20       CONTINUE
   30       IF (JLAST.NE.N+1) GO TO 40
            JFIRST = J
            GO TO 50
   40       IW2(JLAST) = J
   50       IW2(J) = JNEXT
            JLAST = J
            NFRONT = NFRONT + 1
   60     CONTINUE
   70   CONTINUE
        LNASS = NASS - LNASS
   80   NUMORG = NELIM(IASS)
        J1 = IINPUT
        DO 150 IORG = 1,NUMORG
          J = -IW(J1)
          DO 140 IDUMMY = 1,LIW
            JNEW = PERM(J)
            IF (IW2(J).GT.0) GO TO 130
            JLAST = N + 1
            JNEXT = JFIRST
            DO 90 JDUMMY = 1,N
              IF (JNEXT.EQ.N+1) GO TO 100
              IF (PERM(JNEXT).GT.JNEW) GO TO 100
              JLAST = JNEXT
              JNEXT = IW2(JLAST)
   90       CONTINUE
  100       IF (JLAST.NE.N+1) GO TO 110
            JFIRST = J
            GO TO 120
  110       IW2(JLAST) = J
  120       IW2(J) = JNEXT
            NFRONT = NFRONT + 1
  130       J1 = J1 + 1
            IF (J1.GT.LIW) GO TO 150
            J = IW(J1)
            IF (J.LT.0) GO TO 150
  140     CONTINUE
  150   CONTINUE
        IF (NEWEL+NFRONT.LT.ISTK) GO TO 160
        CALL MA27PD(A,IW,ISTK,ISTK2,IINPUT,2,NCMPBR,NCMPBI)
        IF (NEWEL+NFRONT.LT.ISTK) GO TO 160
        INFO(2) = LIW + 1 + NEWEL + NFRONT - ISTK
        GO TO 770
  160   J = JFIRST
        DO 170 IFR = 1,NFRONT
          NEWEL = NEWEL + 1
          IW(NEWEL) = J
          JNEXT = IW2(J)
          IW2(J) = NEWEL - (IWPOS+1)
          J = JNEXT
  170   CONTINUE
        MAXFRT = MAX(MAXFRT,NFRONT)
        IW(IWPOS) = NFRONT
        LAELL = ((NFRONT+1)*NFRONT)/2
        APOS2 = POSFAC + LAELL - 1
        IF (NUMSTK.NE.0) LNASS = LNASS* (2*NFRONT-LNASS+1)/2
        IF (POSFAC+LNASS-1.GE.ASTK) GO TO 180
        IF (APOS2.LT.ASTK+LTOPST-1) GO TO 190
  180   CALL MA27PD(A,IW,ASTK,ASTK2,AINPUT,1,NCMPBR,NCMPBI)
        IF (POSFAC+LNASS-1.GE.ASTK) GO TO 780
        IF (APOS2.GE.ASTK+LTOPST-1) GO TO 780
  190   IF (APOS2.LE.AZERO) GO TO 220
        APOS = AZERO + 1
        LAPOS2 = MIN(APOS2,ASTK-1)
        IF (LAPOS2.LT.APOS) GO TO 210
        DO 200 K = APOS,LAPOS2
          A(K) = ZERO
  200   CONTINUE
  210   AZERO = APOS2
  220   IF (NUMSTK.EQ.0) GO TO 260
        DO 250 IELL = 1,NUMSTK
          J1 = ISTK + 1
          J2 = ISTK + IW(ISTK)
          DO 240 JJ = J1,J2
            IROW = IW(JJ)
            IROW = IW2(IROW)
            APOS = POSFAC + IDIAG(NFRONT,IROW)
            DO 230 JJJ = JJ,J2
              J = IW(JJJ)
              APOS2 = APOS + IW2(J) - IROW
              A(APOS2) = A(APOS2) + A(ASTK)
              A(ASTK) = ZERO
              ASTK = ASTK + 1
  230       CONTINUE
  240     CONTINUE
          ISTK = J2 + 1
  250   CONTINUE
  260   DO 280 IORG = 1,NUMORG
          J = -IW(IINPUT)
          IROW = IW2(J)
          APOS = POSFAC + IDIAG(NFRONT,IROW)
          DO 270 IDUMMY = 1,NZ
            APOS2 = APOS + IW2(J) - IROW
            A(APOS2) = A(APOS2) + A(AINPUT)
            AINPUT = AINPUT + 1
            IINPUT = IINPUT + 1
            IF (IINPUT.GT.LIW) GO TO 280
            J = IW(IINPUT)
            IF (J.LT.0) GO TO 280
  270     CONTINUE
  280   CONTINUE
        NUMASS = NUMASS + NUMORG
        J1 = IWPOS + 2
        J2 = IWPOS + NFRONT + 1
        DO 290 K = J1,J2
          J = IW(K)
          IW2(J) = 0
  290   CONTINUE
        LNPIV = -1
        NPIV = 0
        DO 650 KDUMMY = 1,NASS
          IF (NPIV.EQ.NASS) GO TO 660
          IF (NPIV.EQ.LNPIV) GO TO 660
          LNPIV = NPIV
          NPIVP1 = NPIV + 1
          JPIV = 1
          DO 640 IPIV = NPIVP1,NASS
            JPIV = JPIV - 1
            IF (JPIV.EQ.1) GO TO 640
            APOS = POSFAC + IDIAG(NFRONT-NPIV,IPIV-NPIV)
            IF (UU.GT.ZERO) GO TO 320
            IF (ABS(A(APOS)).LE.CNTL(3)) GO TO 790
            IF (NTOTPV.GT.0) GO TO 300
            IF (A(APOS).GT.ZERO) ISNPIV = 1
            IF (A(APOS).LT.ZERO) ISNPIV = -1
  300       IF (A(APOS).GT.ZERO .AND. ISNPIV.EQ.1) GO TO 560
            IF (A(APOS).LT.ZERO .AND. ISNPIV.EQ.-1) GO TO 560
            IF (INFO(1).NE.2) INFO(2) = 0
            INFO(2) = INFO(2) + 1
            INFO(1) = 2
            I = NTOTPV + 1
            IF (ICNTL(2).GT.0 .AND. INFO(2).LE.10) THEN
              WRITE (ICNTL(2),FMT=310) INFO(1),I
            END IF
  310       FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +              '  *** INFO(1) =',I2,/,' PIVOT',I6,
     +             ' HAS DIFFERENT SIGN FROM THE PREVIOUS ONE')
            ISNPIV = -ISNPIV
            IF (UU.EQ.ZERO) GO TO 560
            GO TO 800
  320       AMAX = ZERO
            TMAX = AMAX
            J1 = APOS + 1
            J2 = APOS + NASS - IPIV
            IF (J2.LT.J1) GO TO 340
            DO 330 JJ = J1,J2
              IF (ABS(A(JJ)).LE.AMAX) GO TO 330
              JMAX = IPIV + JJ - J1 + 1
              AMAX = ABS(A(JJ))
  330       CONTINUE
  340       J1 = J2 + 1
            J2 = APOS + NFRONT - IPIV
            IF (J2.LT.J1) GO TO 360
            DO 350 JJ = J1,J2
              TMAX = MAX(ABS(A(JJ)),TMAX)
  350       CONTINUE
  360       RMAX = MAX(TMAX,AMAX)
            APOS1 = APOS
            KK = NFRONT - IPIV
            LT = IPIV - (NPIV+1)
            IF (LT.EQ.0) GO TO 380
            DO 370 K = 1,LT
              KK = KK + 1
              APOS1 = APOS1 - KK
              RMAX = MAX(RMAX,ABS(A(APOS1)))
  370       CONTINUE
  380       IF (ABS(A(APOS)).GT.MAX(CNTL(3),UU*RMAX)) GO TO 450
            IF (ABS(AMAX).LE.CNTL(3)) GO TO 640
            APOS2 = POSFAC + IDIAG(NFRONT-NPIV,JMAX-NPIV)
            DETPIV = A(APOS)*A(APOS2) - AMAX*AMAX
            THRESH = ABS(DETPIV)
            THRESH = THRESH/ (UU*MAX(ABS(A(APOS))+AMAX,
     +               ABS(A(APOS2))+AMAX))
            IF (THRESH.LE.RMAX) GO TO 640
            RMAX = ZERO
            J1 = APOS2 + 1
            J2 = APOS2 + NFRONT - JMAX
            IF (J2.LT.J1) GO TO 400
            DO 390 JJ = J1,J2
              RMAX = MAX(RMAX,ABS(A(JJ)))
  390       CONTINUE
  400       KK = NFRONT - JMAX + 1
            APOS3 = APOS2
            JMXMIP = JMAX - IPIV - 1
            IF (JMXMIP.EQ.0) GO TO 420
            DO 410 K = 1,JMXMIP
              APOS2 = APOS2 - KK
              KK = KK + 1
              RMAX = MAX(RMAX,ABS(A(APOS2)))
  410       CONTINUE
  420       IPMNP = IPIV - NPIV - 1
            IF (IPMNP.EQ.0) GO TO 440
            APOS2 = APOS2 - KK
            KK = KK + 1
            DO 430 K = 1,IPMNP
              APOS2 = APOS2 - KK
              KK = KK + 1
              RMAX = MAX(RMAX,ABS(A(APOS2)))
  430       CONTINUE
  440       IF (THRESH.LE.RMAX) GO TO 640
            PIVSIZ = 2
            GO TO 460
  450       PIVSIZ = 1
  460       IROW = IPIV - NPIV
            DO 550 KROW = 1,PIVSIZ
              IF (IROW.EQ.1) GO TO 530
              J1 = POSFAC + IROW
              J2 = POSFAC + NFRONT - (NPIV+1)
              IF (J2.LT.J1) GO TO 480
              APOS2 = APOS + 1
              DO 470 JJ = J1,J2
                SWOP = A(APOS2)
                A(APOS2) = A(JJ)
                A(JJ) = SWOP
                APOS2 = APOS2 + 1
  470         CONTINUE
  480         J1 = POSFAC + 1
              J2 = POSFAC + IROW - 2
              APOS2 = APOS
              KK = NFRONT - (IROW+NPIV)
              IF (J2.LT.J1) GO TO 500
              DO 490 JJJ = J1,J2
                JJ = J2 - JJJ + J1
                KK = KK + 1
                APOS2 = APOS2 - KK
                SWOP = A(APOS2)
                A(APOS2) = A(JJ)
                A(JJ) = SWOP
  490         CONTINUE
  500         IF (NPIV.EQ.0) GO TO 520
              APOS1 = POSFAC
              KK = KK + 1
              APOS2 = APOS2 - KK
              DO 510 JJ = 1,NPIV
                KK = KK + 1
                APOS1 = APOS1 - KK
                APOS2 = APOS2 - KK
                SWOP = A(APOS2)
                A(APOS2) = A(APOS1)
                A(APOS1) = SWOP
  510         CONTINUE
  520         SWOP = A(APOS)
              A(APOS) = A(POSFAC)
              A(POSFAC) = SWOP
              IPOS = IWPOS + NPIV + 2
              IEXCH = IWPOS + IROW + NPIV + 1
              ISWOP = IW(IPOS)
              IW(IPOS) = IW(IEXCH)
              IW(IEXCH) = ISWOP
  530         IF (PIVSIZ.EQ.1) GO TO 550
              IF (KROW.EQ.2) GO TO 540
              IROW = JMAX - (NPIV+1)
              JPOS = POSFAC
              POSFAC = POSFAC + NFRONT - NPIV
              NPIV = NPIV + 1
              APOS = APOS3
              GO TO 550
  540         NPIV = NPIV - 1
              POSFAC = JPOS
  550       CONTINUE
            IF (PIVSIZ.EQ.2) GO TO 600
  560       A(POSFAC) = ONE/A(POSFAC)
            IF (A(POSFAC).LT.ZERO) NEIG = NEIG + 1
            J1 = POSFAC + 1
            J2 = POSFAC + NFRONT - (NPIV+1)
            IF (J2.LT.J1) GO TO 590
            IBEG = J2 + 1
            DO 580 JJ = J1,J2
              AMULT = -A(JJ)*A(POSFAC)
              IEND = IBEG + NFRONT - (NPIV+JJ-J1+2)
CDIR$ IVDEP
              DO 570 IROW = IBEG,IEND
                JCOL = JJ + IROW - IBEG
                A(IROW) = A(IROW) + AMULT*A(JCOL)
  570         CONTINUE
              IBEG = IEND + 1
              A(JJ) = AMULT
  580       CONTINUE
  590       NPIV = NPIV + 1
            NTOTPV = NTOTPV + 1
            JPIV = 1
            POSFAC = POSFAC + NFRONT - NPIV + 1
            GO TO 640
  600       IPOS = IWPOS + NPIV + 2
            NTWO = NTWO + 1
            IW(IPOS) = -IW(IPOS)
            POSPV1 = POSFAC
            POSPV2 = POSFAC + NFRONT - NPIV
            SWOP = A(POSPV2)
            IF (DETPIV.LT.ZERO) NEIG = NEIG + 1
            IF (DETPIV.GT.ZERO .AND. SWOP.LT.ZERO) NEIG = NEIG + 2
            A(POSPV2) = A(POSPV1)/DETPIV
            A(POSPV1) = SWOP/DETPIV
            A(POSPV1+1) = -A(POSPV1+1)/DETPIV
            J1 = POSPV1 + 2
            J2 = POSPV1 + NFRONT - (NPIV+1)
            IF (J2.LT.J1) GO TO 630
            JJ1 = POSPV2
            IBEG = POSPV2 + NFRONT - (NPIV+1)
            DO 620 JJ = J1,J2
              JJ1 = JJ1 + 1
              AMULT1 = - (A(POSPV1)*A(JJ)+A(POSPV1+1)*A(JJ1))
              AMULT2 = - (A(POSPV1+1)*A(JJ)+A(POSPV2)*A(JJ1))
              IEND = IBEG + NFRONT - (NPIV+JJ-J1+3)
CDIR$ IVDEP
              DO 610 IROW = IBEG,IEND
                K1 = JJ + IROW - IBEG
                K2 = JJ1 + IROW - IBEG
                A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
  610         CONTINUE
              IBEG = IEND + 1
              A(JJ) = AMULT1
              A(JJ1) = AMULT2
  620       CONTINUE
  630       NPIV = NPIV + 2
            NTOTPV = NTOTPV + 2
            JPIV = 2
            POSFAC = POSPV2 + NFRONT - NPIV + 1
  640     CONTINUE
  650   CONTINUE
  660   IF (NPIV.NE.0) NBLK = NBLK + 1
        IOLDPS = IWPOS
        IWPOS = IWPOS + NFRONT + 2
        IF (NPIV.EQ.0) GO TO 690
        IF (NPIV.GT.1) GO TO 680
        IW(IOLDPS) = -IW(IOLDPS)
        DO 670 K = 1,NFRONT
          J1 = IOLDPS + K
          IW(J1) = IW(J1+1)
  670   CONTINUE
        IWPOS = IWPOS - 1
        GO TO 690
  680   IW(IOLDPS+1) = NPIV
  690   LIELL = NFRONT - NPIV
        IF (LIELL.EQ.0 .OR. IASS.EQ.NSTEPS) GO TO 750
        IF (IWPOS+LIELL.LT.ISTK) GO TO 700
        CALL MA27PD(A,IW,ISTK,ISTK2,IINPUT,2,NCMPBR,NCMPBI)
  700   ISTK = ISTK - LIELL - 1
        IW(ISTK) = LIELL
        J1 = ISTK
        KK = IWPOS - LIELL - 1
CDIR$ IVDEP
        DO 710 K = 1,LIELL
          J1 = J1 + 1
          KK = KK + 1
          IW(J1) = IW(KK)
  710   CONTINUE
        LAELL = ((LIELL+1)*LIELL)/2
        KK = POSFAC + LAELL
        IF (KK.NE.ASTK) GO TO 720
        ASTK = ASTK - LAELL
        GO TO 740
  720   KMAX = KK - 1
CDIR$ IVDEP
        DO 730 K = 1,LAELL
          KK = KK - 1
          ASTK = ASTK - 1
          A(ASTK) = A(KK)
  730   CONTINUE
        KMAX = MIN(KMAX,ASTK-1)
        DO 735 K = KK,KMAX
          A(K) = ZERO
  735   CONTINUE
  740   AZERO = MIN(AZERO,ASTK-1)
  750   IF (NPIV.EQ.0) IWPOS = IOLDPS
  760 CONTINUE
      IW(1) = NBLK
      IF (NTWO.GT.0) IW(1) = -NBLK
      NRLBDU = POSFAC - 1
      NIRBDU = IWPOS - 1
      IF (NTOTPV.EQ.N) GO TO 810
      INFO(1) = 3
      INFO(2) = NTOTPV
      GO TO 810
  770 INFO(1) = -3
      GO TO 810
  780 INFO(1) = -4
      INFO(2) = LA + MAX(POSFAC+LNASS,APOS2-LTOPST+2) - ASTK
      GO TO 810
  790 INFO(1) = -5
      INFO(2) = NTOTPV + 1
      GO TO 810
  800 INFO(1) = -6
      INFO(2) = NTOTPV + 1
  810 CONTINUE
      INFO(9) = NRLBDU
      INFO(10) = NIRBDU
      INFO(12) = NCMPBR
      INFO(13) = NCMPBI
      INFO(14) = NTWO
      INFO(15) = NEIG
      RETURN
      END
      SUBROUTINE MA27PD(A,IW,J1,J2,ITOP,IREAL,NCMPBR,NCMPBI)
      INTEGER IREAL,ITOP,J1,J2,NCMPBR,NCMPBI
      DOUBLE PRECISION A(*)
      INTEGER IW(*)
      INTEGER IPOS,JJ,JJJ
      IPOS = ITOP - 1
      IF (J2.EQ.IPOS) GO TO 50
      IF (IREAL.EQ.2) GO TO 20
      NCMPBR = NCMPBR + 1
      IF (J1.GT.J2) GO TO 40
      DO 10 JJJ = J1,J2
        JJ = J2 - JJJ + J1
        A(IPOS) = A(JJ)
        IPOS = IPOS - 1
   10 CONTINUE
      GO TO 40
   20 NCMPBI = NCMPBI + 1
      IF (J1.GT.J2) GO TO 40
      DO 30 JJJ = J1,J2
        JJ = J2 - JJJ + J1
        IW(IPOS) = IW(JJ)
        IPOS = IPOS - 1
   30 CONTINUE
   40 J2 = ITOP - 1
      J1 = IPOS + 1
   50 RETURN
      END
      SUBROUTINE MA27QD(N,A,LA,IW,LIW,W,MAXFNT,RHS,IW2,NBLK,LATOP,ICNTL)
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
      INTEGER LA,LATOP,LIW,MAXFNT,N,NBLK
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFNT)
      INTEGER IW(LIW),IW2(NBLK),ICNTL(30)
      DOUBLE PRECISION W1,W2
      INTEGER APOS,IBLK,IFR,ILVL,IPIV,IPOS,IRHS,IROW,IST,J,J1,J2,J3,JJ,
     +        JPIV,K,K1,K2,K3,LIELL,NPIV
      INTRINSIC ABS,MIN
      APOS = 1
      IPOS = 1
      J2 = 0
      IBLK = 0
      NPIV = 0
      DO 140 IROW = 1,N
        IF (NPIV.GT.0) GO TO 90
        IBLK = IBLK + 1
        IF (IBLK.GT.NBLK) GO TO 150
        IPOS = J2 + 1
        IW2(IBLK) = IPOS
        LIELL = -IW(IPOS)
        NPIV = 1
        IF (LIELL.GT.0) GO TO 10
        LIELL = -LIELL
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
   10   J1 = IPOS + 1
        J2 = IPOS + LIELL
        ILVL = MIN(NPIV,10)
        IF (LIELL.LT.ICNTL(IFRLVL+ILVL)) GO TO 90
        IFR = 0
        DO 20 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   20   CONTINUE
        JPIV = 1
        J3 = J1
        DO 70 IPIV = 1,NPIV
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 70
          IF (IW(J3).LT.0) GO TO 40
          JPIV = 1
          J3 = J3 + 1
          APOS = APOS + 1
          IST = IPIV + 1
          IF (LIELL.LT.IST) GO TO 70
          W1 = W(IPIV)
          K = APOS
          DO 30 J = IST,LIELL
            W(J) = W(J) + A(K)*W1
            K = K + 1
   30     CONTINUE
          APOS = APOS + LIELL - IST + 1
          GO TO 70
   40     JPIV = 2
          J3 = J3 + 2
          APOS = APOS + 2
          IST = IPIV + 2
          IF (LIELL.LT.IST) GO TO 60
          W1 = W(IPIV)
          W2 = W(IPIV+1)
          K1 = APOS
          K2 = APOS + LIELL - IPIV
          DO 50 J = IST,LIELL
            W(J) = W(J) + W1*A(K1) + W2*A(K2)
            K1 = K1 + 1
            K2 = K2 + 1
   50     CONTINUE
   60     APOS = APOS + 2* (LIELL-IST+1) + 1
   70   CONTINUE
        IFR = 0
        DO 80 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          RHS(J) = W(IFR)
   80   CONTINUE
        NPIV = 0
        GO TO 140
   90   IF (IW(J1).LT.0) GO TO 110
        NPIV = NPIV - 1
        APOS = APOS + 1
        J1 = J1 + 1
        IF (J1.GT.J2) GO TO 140
        IRHS = IW(J1-1)
        W1 = RHS(IRHS)
        K = APOS
        DO 100 J = J1,J2
          IRHS = ABS(IW(J)+0)
          RHS(IRHS) = RHS(IRHS) + A(K)*W1
          K = K + 1
  100   CONTINUE
        APOS = APOS + J2 - J1 + 1
        GO TO 140
  110   NPIV = NPIV - 2
        J1 = J1 + 2
        APOS = APOS + 2
        IF (J1.GT.J2) GO TO 130
        IRHS = -IW(J1-2)
        W1 = RHS(IRHS)
        IRHS = IW(J1-1)
        W2 = RHS(IRHS)
        K1 = APOS
        K3 = APOS + J2 - J1 + 2
        DO 120 J = J1,J2
          IRHS = ABS(IW(J)+0)
          RHS(IRHS) = RHS(IRHS) + W1*A(K1) + W2*A(K3)
          K1 = K1 + 1
          K3 = K3 + 1
  120   CONTINUE
  130   APOS = APOS + 2* (J2-J1+1) + 1
  140 CONTINUE
  150 LATOP = APOS - 1
      RETURN
      END
      SUBROUTINE MA27RD(N,A,LA,IW,LIW,W,MAXFNT,RHS,IW2,NBLK,LATOP,ICNTL)
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
      INTEGER LA,LATOP,LIW,MAXFNT,N,NBLK
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFNT)
      INTEGER IW(LIW),IW2(NBLK),ICNTL(30)
      DOUBLE PRECISION W1,W2
      INTEGER APOS,APOS2,I1RHS,I2RHS,IBLK,IFR,IIPIV,IIRHS,ILVL,IPIV,
     +        IPOS,IRHS,IST,J,J1,J2,JJ,JJ1,JJ2,JPIV,JPOS,K,LIELL,LOOP,
     +        NPIV
      INTRINSIC ABS,MIN
      APOS = LATOP + 1
      NPIV = 0
      IBLK = NBLK + 1
      DO 180 LOOP = 1,N
        IF (NPIV.GT.0) GO TO 110
        IBLK = IBLK - 1
        IF (IBLK.LT.1) GO TO 190
        IPOS = IW2(IBLK)
        LIELL = -IW(IPOS)
        NPIV = 1
        IF (LIELL.GT.0) GO TO 10
        LIELL = -LIELL
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
   10   JPOS = IPOS + NPIV
        J2 = IPOS + LIELL
        ILVL = MIN(10,NPIV) + 10
        IF (LIELL.LT.ICNTL(IFRLVL+ILVL)) GO TO 110
        J1 = IPOS + 1
        IFR = 0
        DO 20 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   20   CONTINUE
        JPIV = 1
        DO 90 IIPIV = 1,NPIV
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 90
          IPIV = NPIV - IIPIV + 1
          IF (IPIV.EQ.1) GO TO 30
          IF (IW(JPOS-1).LT.0) GO TO 60
   30     JPIV = 1
          APOS = APOS - (LIELL+1-IPIV)
          IST = IPIV + 1
          W1 = W(IPIV)*A(APOS)
          IF (LIELL.LT.IST) GO TO 50
          JJ1 = APOS + 1
          DO 40 J = IST,LIELL
            W1 = W1 + A(JJ1)*W(J)
            JJ1 = JJ1 + 1
   40     CONTINUE
   50     W(IPIV) = W1
          JPOS = JPOS - 1
          GO TO 90
   60     JPIV = 2
          APOS2 = APOS - (LIELL+1-IPIV)
          APOS = APOS2 - (LIELL+2-IPIV)
          IST = IPIV + 1
          W1 = W(IPIV-1)*A(APOS) + W(IPIV)*A(APOS+1)
          W2 = W(IPIV-1)*A(APOS+1) + W(IPIV)*A(APOS2)
          IF (LIELL.LT.IST) GO TO 80
          JJ1 = APOS + 2
          JJ2 = APOS2 + 1
          DO 70 J = IST,LIELL
            W1 = W1 + W(J)*A(JJ1)
            W2 = W2 + W(J)*A(JJ2)
            JJ1 = JJ1 + 1
            JJ2 = JJ2 + 1
   70     CONTINUE
   80     W(IPIV-1) = W1
          W(IPIV) = W2
          JPOS = JPOS - 2
   90   CONTINUE
        IFR = 0
        DO 100 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          RHS(J) = W(IFR)
  100   CONTINUE
        NPIV = 0
        GO TO 180
  110   IF (NPIV.EQ.1) GO TO 120
        IF (IW(JPOS-1).LT.0) GO TO 150
  120   NPIV = NPIV - 1
        APOS = APOS - (J2-JPOS+1)
        IIRHS = IW(JPOS)
        W1 = RHS(IIRHS)*A(APOS)
        J1 = JPOS + 1
        IF (J1.GT.J2) GO TO 140
        K = APOS + 1
        DO 130 J = J1,J2
          IRHS = ABS(IW(J)+0)
          W1 = W1 + A(K)*RHS(IRHS)
          K = K + 1
  130   CONTINUE
  140   RHS(IIRHS) = W1
        JPOS = JPOS - 1
        GO TO 180
  150   NPIV = NPIV - 2
        APOS2 = APOS - (J2-JPOS+1)
        APOS = APOS2 - (J2-JPOS+2)
        I1RHS = -IW(JPOS-1)
        I2RHS = IW(JPOS)
        W1 = RHS(I1RHS)*A(APOS) + RHS(I2RHS)*A(APOS+1)
        W2 = RHS(I1RHS)*A(APOS+1) + RHS(I2RHS)*A(APOS2)
        J1 = JPOS + 1
        IF (J1.GT.J2) GO TO 170
        JJ1 = APOS + 2
        JJ2 = APOS2 + 1
        DO 160 J = J1,J2
          IRHS = ABS(IW(J)+0)
          W1 = W1 + RHS(IRHS)*A(JJ1)
          W2 = W2 + RHS(IRHS)*A(JJ2)
          JJ1 = JJ1 + 1
          JJ2 = JJ2 + 1
  160   CONTINUE
  170   RHS(I1RHS) = W1
        RHS(I2RHS) = W2
        JPOS = JPOS - 2
  180 CONTINUE
  190 RETURN
      END



* *******************************************************************
* COPYRIGHT (c) 1967 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 4 Oct 1992
C       Toolpack tool decs employed.
C       SAVE statement for COMMON FA01ED added.
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C
C
      DOUBLE PRECISION FUNCTION FA01AD(I)
      INTEGER I
      DOUBLE PRECISION R,S
      INTRINSIC DINT,MOD
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      EXTERNAL FA01FD
      SAVE /FA01ED/
      R = GR*9228907D0/65536D0
      S = DINT(R)
      GL = MOD(S+GL*9228907D0,65536D0)
      GR = R - S
      IF (I.GE.0) FA01AD = (GL+GR)/65536D0
      IF (I.LT.0) FA01AD = (GL+GR)/32768D0 - 1.D0
      GR = GR*65536D0
      RETURN
      END
      SUBROUTINE FA01BD(MAX,NRAND)
      INTEGER MAX,NRAND
      DOUBLE PRECISION FA01AD
      EXTERNAL FA01AD
      INTRINSIC DBLE,INT
      NRAND = INT(FA01AD(1)*DBLE(MAX)) + 1
      RETURN
      END
      SUBROUTINE FA01CD(IL,IR)
      INTEGER IL,IR
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      IL = GL
      IR = GR
      RETURN
      END
      SUBROUTINE FA01DD(IL,IR)
      INTEGER IL,IR
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      GL = IL
      GR = IR
      RETURN
      END
      BLOCK DATA FA01FD
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      DATA GL/21845D0/
      DATA GR/21845D0/
      END

c     *******************************************
c     Pontos e pesos de gauss
c     ******************************************* 
      subroutine gauss(nga,qsiga,wga)

      integer*4 nga,iga,jga,mga
      real*8 x1ga,x2ga,qsiga(20),wga(20),p1ga,p2ga,p3ga,ppga
     #,xlga,xmga,zga,z1ga,pi
      pi=2.d0*datan2(1.d0,0.d0) !definido pi radianos
      x1ga=-1.d0
      x2ga=1.d0
      mga=(nga+1)/2
      xmga=0.d0
      xlga=1.d0
      do iga=1,mga
        zga=cos(pi*(iga-.25d0)/(nga+.5d0))
 1      continue !Newton rapson para refinar raiz
        p1ga=1.d0
        p2ga=0.d0        
        do jga=1,nga
            p3ga=p2ga
            p2ga=p1ga
            p1ga=((2.d0*jga-1.d0)*zga*p2ga-(jga-1.d0)*p3ga)/(jga*1.d0)
        enddo
        ppga=nga*(zga*p1ga-p2ga)/(zga*zga-1.d0)
        z1ga=zga
        zga=z1ga-p1ga/ppga
        if(dabs(zga-z1ga).gt.1.d-15) goto 1 !Newton Rapson
        qsiga(iga)=xmga-xlga*zga
        qsiga(nga+1-iga)=xmga+xlga*zga        
        wga(iga)=2.d0*xlga/((1.d0-zga*zga)*ppga*ppga)
        wga(nga+1-iga)=wga(iga)
      enddo

c      do iga=1,nga
c      write(*,*) iga,qsiga(iga),wga(iga)
c      enddo

      return
      end


!****************************************************************************
!
!  ROTINA:     INVERSE
!
!  PROP�SITO:  Inverter uma matriz
!
!
!****************************************************************************
      subroutine inverse(iEq, dbMin,dbMout)

      implicit none
      integer ii, ij, ik, iEq
      real*8 dbMin(iEq,iEq), dbMout(iEq,iEq)
      real*8 Id(iEq,iEq), dbAux


      ! iEq            - numero de equa��es
      ! dbMin            - Matriz a ser resolvida 

      ! Gera a matriz identidade (onde estar� a resposta - dbMout)
      ! Copia dbMin para n�o alter�-la - Id no final ser� a matriz identidade

      do ii = 1,iEq
            do ij = 1,iEq
                  if (ii .eq. ij) then 
                        dbMout(ii,ij) = 1.d0
                  else
                        dbMout(ii,ij) = 0.d0
                  end if

c                  Id(ii,ij) = dbMin(ii,ij)
            end do
      end do

      do ii=1, iEq
       if (dabs(dbmin(ii,ii)) .lt. 1.d-10) then
         do ij=1, iEq
             if (dabs(dbmin(ij,ii)).gt.1.d-10) then
            do ik=1,iEq
             dbmin(ii,ik)=dbmin(ii,ik)+dbmin(ij,ik)
             dbmout(ii,ik)=dbmout(ii,ik)+dbmout(ij,ik)
            enddo
          endif
        enddo
       else
       endif
      enddo


      id=dbmin

      do ii = 1,iEq
            dbAux = Id(ii,ii)
            do ij=1,iEq
                  Id(ii,ij) = Id(ii,ij) / dbAux
                  dbMout(ii,ij) = dbMout(ii,ij) / dbAux
            end do

            do ik = ii+1,iEq
                  if (id(ik,ii) .ne. 0.d0) then
                        dbAux = Id(ik,ii)


                        do ij = 1,iEq
                              Id(ik,ij) = Id(ik,ij) - Id(ii,ij) * dbAux
                              dbMout(ik,ij) = dbMout(ik,ij) - 
     #                                            dbMout(ii,ij) * dbAux
                        end do
                  end if
            end do
      end do

      do ii = iEq,1,-1
            do ik = ii-1,1,-1
                  if (id(ik,ii) .ne. 0.d0) then
                        dbAux = Id(ik,ii)

                        do ij = iEq, 1, -1
                              Id(ik,ij) = Id(ik,ij) - Id(ii,ij) * dbAux
                              dbMout(ik,ij) = dbMout(ik,ij) - 
     #                                            dbMout(ii,ij) * dbAux
                        end do
                  end if
            end do
      end do

      return 
      end

c     *******************************************
c     Saida para pos processamento
c     ******************************************* 
      subroutine saidaposparaview
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      character*10 iprintc
      write(iprintc,'(I10)') ipt+1000000000
      
      open(unit=8,file='desnov'//iprintc//'.vtu',status='unknown')
      
      write(8,500)
      write(8,*)'<VTKFile type="UnstructuredGrid">'
      write(8,*)'  <UnstructuredGrid>'
      write(8,*)'  <Piece NumberOfPoints="', nnos,'
     $"  NumberOfCells= "',nel,'">'
      
      
c     imprime coordenadas dos nos
      write(8,*) '    <Points>'
      write(8,*) '      <DataArray type="Float64"
     $ NumberOfComponents ="3" format="ascii">'
      do i=1,nnos
c     p(glgl(i,1)),' ',p(glgl(i,2)),' ',0!
         write(8,*)     solint(glgl(i,1)),' ',solint(glgl(i,2)),' ',0! coordint(1,i),' ',coordint(2,i),' ',0
      enddo
      write(8,*)'      </DataArray>'
      write(8,*)'    </Points>'
      
c     write elemente connectivity
      write(8,*)'    <Cells>'
      write(8,*)'      <DataArray type="Int32" Name= "connectivity"
     $ format="ascii">'
      do k=1,nel
         if(itipoel(k).ne.2)then
            if(itipoel(k).ne.3)then
               write(8,*) (ic(k,j)-1,j=2,3)
         endif
      endif
      if(itipoel(k).eq.2)then
         write(8,*) ic(k,1)-1,ic(k,3)-1
      endif
      if(itipoel(k).eq.3)then
         write(8,*) ic(k,2)-1,ic(k,4)-1
      endif
c$$$          if(itipoel(k).eq.4)then
c$$$            write(8,*) (ic(k,j)-1,j=2,3)
c$$$         endif
c$$$         if(itipoel(k).eq.3)then
c$$$            write(8,*) (ic(k,j)-1,j=2,3)
c$$$         endif
      enddo
      write(8,*) '      </DataArray>'
      
c     write offsets in the data array
      write(8,*) '      <DataArray type="Int32" Name="offsets"
     $ format="ascii">'
      nant=0
      do jel=1,nel
         nant=nant+2
         write(8,*) nant
      enddo
      write(8,*)'      </DataArray>'
      
c     write element types
      write(8,*) '      <DataArray type="UInt8" Name="types"
     $ format ="ascii">'
      do jel=1,nel
         write(8,*) 4
      enddo
      write(8,*)'      </DataArray>'
      write(8,*)'    </Cells>'
      
c     write deslocamentos
      write(8,*)'    <PointData>'
      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
     $=" 3" Name="desloc"           format="ascii">'
      do i=1,nnos
c$$$         write(8,*) p(glgl(i,1))-p0(glgl(i,1)),' ',p(glgl(i
c$$$     $        ,2))-p0(glgl(i,2)),' ',0
         
         write(8,*) solint(glgl(i,1))- confini(1,i),' ',solint(glgl(i
     $        ,2))-confini(2,i),' ',0

c     write(8,*) solint(glgl(i,1)),' ',solint(glgl(i,2)),' ',0

c        write(8,*) confini(1,i),' ', confini(2,i),' ',0
         
      enddo
      write(8,*)'      </DataArray> '

c     write deslocamentos do pontos de controle
      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
     $=" 3" Name="deslpc"           format="ascii">'
      do i=1,nnos
          write(8,*) p(glgl(i,1))-p0(glgl(i,1)),' ',p(glgl(i
     $        ,2))-p0(glgl(i,2)),' ',0

      enddo
      write(8,*)'      </DataArray> '


c     write tensoes
      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
     $=" 3" Name="confin" format="ascii">'
      do i=1,nnos
         write(8,*) confini(1,i),' ',confini(2,i),' ',0
         
      enddo
      write(8,*)'      </DataArray> '

c     write tensoes
      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
     $=" 3" Name="tensao" format="ascii">'
      do i=1,nnos
        do ii=1,3
         if (sigmaacima(i,ii).lt.(10**(-30.)))sigmaacima(i,ii)=0.d0
        enddo
         write(8,*) sigmaacima(i,1),' ',sigmaacima(i,2),' ',sigmaacima(i
     $        ,3)
      enddo
      write(8,*)'      </DataArray> '
      
      
c$$$c     write tensao abaixo
c$$$      call tensaoabaixo         !Cacula tensoe principais
c$$$      
c$$$      
c     inicio das listas    
c     preparando o tensao x
c      valor=0.
c      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         
c$$$         call prepara
c$$$         
c$$$         do ino=2,3 !notl(j)	   
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sx(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$         
c$$$      enddo
c$$$      
c$$$c     valor=sx
c$$$      
c$$$      fno=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno(j)=fno(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$      
c$$$      
c$$$      
c$$$      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
c$$$     $=" 3" Name="tensaoabaixo" format="ascii">'
c$$$      do i=1,nnos
c$$$         write(8,*) fno(i),' ',fno1(i),' ',fno2(i)
c$$$      enddo
c$$$      write(8,*)'      </DataArray> '
c$$$      
c$$$      
c$$$c     write tensao acima
c$$$      call tensaoacima          !Cacula tensoe principais
c$$$      
c$$$      
c$$$c     inicio das listas    
c$$$c     preparando o tensao x
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         
c$$$         call prepara
c$$$         
c$$$         do ino=1,notl(j)	   
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sx(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$         
c$$$      enddo
c$$$      
c$$$c     valor=sx
c$$$      
c$$$      fno=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno(j)=fno(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$c     preparando tensao y
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         call prepara
c$$$         do ino=1,notl(j)
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sy(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$      enddo
c$$$      fno1=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno1(ic(j,ino))=fno1(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno1(j)=fno1(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$      
c$$$c     preparando tensaoxy
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         call prepara
c$$$         do ino=1,notl(j)
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sxy(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$      enddo
c$$$      fno2=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno2(ic(j,ino))=fno2(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno2(j)=fno2(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$      
c$$$      
c$$$      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
c$$$     $=" 3" Name="tensaoacima" format="ascii">'
c$$$      do i=1,nnos
c$$$         write(8,*) fno(i),' ',fno1(i),' ',fno2(i)
c$$$      enddo
c$$$      write(8,*)'      </DataArray> '
      write(8,*)'    </PointData>'
      
c     finalise
      write(8,*)'  </Piece>'
      write(8,*)'  </UnstructuredGrid>'
      write(8,*)'</VTKFile>'
      
 500  format(21h<?xml version="1.0"?>)
      return
      end

c     *******************************************
c     Saida para pos processamento
c     ******************************************* 
      subroutine saidaposparaviewo
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      character*10 iprintc
      write(iprintc,'(I10)') ipt+1000000000
      
      open(unit=8,file='desnov'//iprintc//'.vtu',status='unknown')
      
      write(8,500)
      write(8,*)'<VTKFile type="UnstructuredGrid">'
      write(8,*)'  <UnstructuredGrid>'
      write(8,*)'  <Piece NumberOfPoints="', nnos,'
     $"  NumberOfCells= "',nel,'">'
      
      
c     imprime coordenadas dos nos
      write(8,*) '    <Points>'
      write(8,*) '      <DataArray type="Float64" 
     $NumberOfComponents ="3" format="ascii">'
      do i=1,nnos
         write(8,*) p0(glgl(i,1)),' ',p0(glgl(i,2)),' ',0
      enddo
      write(8,*)'      </DataArray>'
      write(8,*)'    </Points>'
      
c     write elemente connectivity
      write(8,*)'    <Cells>'
      write(8,*)'      <DataArray type="Int32" Name= "connectivity" 
     $format="ascii">'
      do k=1,nel
         if(itipoel(k).eq.1)then
            write(8,*) (ic(k,j)-1,j=2,3)
         endif
         if(itipoel(k).eq.2)then
            write(8,*) (ic(k,j)-1,j=1,2)
         endif
         if(itipoel(k).eq.3)then
            write(8,*) (ic(k,j)-1,j=2,3)
         endif
      enddo
      write(8,*) '      </DataArray>'
      
c     write offsets in the data array
      write(8,*) '      <DataArray type="Int32" Name="offsets" 
     $format="ascii">'
      nant=0
      do jel=1,nel
         nant=nant+2
         write(8,*) nant
      enddo
      write(8,*)'      </DataArray>'
      
c     write element types
      write(8,*) '      <DataArray type="UInt8" Name="types" 
     $format ="ascii">'
      do jel=1,nel
         write(8,*) 4
      enddo
      write(8,*)'      </DataArray>'
      write(8,*)'    </Cells>'
      
c     write deslocamentos
      write(8,*)'    <PointData>'
      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
     $=" 3" Name="desloc" format="ascii">' 
      do i=1,nnos
         write(8,*) p(glgl(i,1))-p0(glgl(i,1)),' ',p(glgl(i
     $        ,2))-p0(glgl(i,2)),' ',0
         
      enddo
      write(8,*)'      </DataArray> '

c     write tensoes
      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
     $=" 3" Name="tensao" format="ascii">'
      do i=1,nnos
         write(8,*) sigmaacima(i,1),' ',sigmaacima(i,2),' ',sigmaacima(i
     $        ,3)
         
      enddo
      write(8,*)'      </DataArray> '
      
      
c$$$c     write tensao abaixo
c$$$      call tensaoabaixo         !Cacula tensoe principais
c$$$      
c$$$      
c$$$c     inicio das listas    
c$$$c     preparando o tensao x
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         
c$$$         call prepara
c$$$         
c$$$         do ino=1,notl(j)	   
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sx(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$         
c$$$      enddo
c$$$      
c$$$c     valor=sx
c$$$      
c$$$      fno=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno(j)=fno(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$      
c$$$      
c$$$      
c$$$      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
c$$$     $=" 3" Name="tensaoabaixo" format="ascii">'
c$$$      do i=1,nnos
c$$$         write(8,*) fno(i),' ',fno1(i),' ',fno2(i)
c$$$      enddo
c$$$      write(8,*)'      </DataArray> '
c$$$      
c$$$      
c$$$c     write tensao acima
c$$$      call tensaoacima          !Cacula tensoe principais
c$$$      
c$$$      
c$$$c     inicio das listas    
c$$$c     preparando o tensao x
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         
c$$$         call prepara
c$$$         
c$$$         do ino=1,notl(j)	   
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sx(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$         
c$$$      enddo
c$$$      
c$$$c     valor=sx
c$$$      
c$$$      fno=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno(ic(j,ino))=fno(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno(j)=fno(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$c     preparando tensao y
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         call prepara
c$$$         do ino=1,notl(j)
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sy(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$      enddo
c$$$      fno1=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno1(ic(j,ino))=fno1(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno1(j)=fno1(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$      
c$$$c     preparando tensaoxy
c$$$      valor=0.
c$$$      do j=1,nel
c$$$         ng=notl(j)                   !notl(j)
c$$$         call prepara
c$$$         do ino=1,notl(j)
c$$$            do ih=1,ng
c$$$               valor(j,ino)=valor(j,ino)+rmfi(ino,ih)*sxy(j,ih)
c$$$            enddo
c$$$         enddo
c$$$         deALLOCATE (rmfi)
c$$$      enddo
c$$$      fno2=0.
c$$$      do j=1,nel
c$$$         DO ino=1,notl(j)       !monta o vetor de valores nodais
c$$$            fno2(ic(j,ino))=fno2(ic(j,ino))+valor(j,ino)
c$$$         ENDDO
c$$$      enddo
c$$$      do j=1,nnos
c$$$         fno2(j)=fno2(j)/(kn(j)*1.)
c$$$      enddo
c$$$      
c$$$      
c$$$      
c$$$      write(8,*)'      <DataArray type="Float64" NumberOfComponents 
c$$$     $=" 3" Name="tensaoacima" format="ascii">'
c$$$      do i=1,nnos
c$$$         write(8,*) fno(i),' ',fno1(i),' ',fno2(i)
c$$$      enddo
c$$$      write(8,*)'      </DataArray> '
      write(8,*)'    </PointData>'
      
c     finalise
      write(8,*)'  </Piece>'
      write(8,*)'  </UnstructuredGrid>'
      write(8,*)'</VTKFile>'
      
 500  format(21h<?xml version="1.0"?>)
      return
      end
      

c     *******************************************
c     subrotina que imprime derivadas primeira e segunda
c     *******************************************   
      subroutine imprege
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      write(*,*) 'x', '    y', '       dy_dx'
      do j=1,nel
         do i=1,ng

            xsi1=qsi(i)
            call formaederi0
            call calc_A0 

            
            x=0.d0
            y=0.d0
            dx_dxsi=0.d0
            dy_dxsi=0.d0
            do k=1,notl(j)
               nog=ic(j,k)!encontra numero global do n�x
               x=p0(glgl(nog,1))*fi(k)+x
               y=p0(glgl(nog,2))*fi(k)+y

c              derivada de x em rela�ao a xsi
               dx_dxsi=dfi(k,1)*p0(glgl(nog,1))+dx_dxsi
c              derivada de y em rela�ao a xsi
               dy_dxsi=dfi(k,1)*p0(glgl(nog,2))+dy_dxsi


           
            enddo


c           derivada de y em rela��o a x
            dy_dx=(dy_dxsi)*(1.d0/(dx_dxsi))


            write(*,*) x,y,dy_dx
         enddo
      enddo
      end



c     ***************************************************
c     Subrotina que calcula normal no ponto
c     ***************************************************
      subroutine vetnorm
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

      write(*,*) 'x','y', 'nx_', 'ny_'

      Anx(1)=-dsin(1.487779048189210000d0)
      Any(1)=dcos(1.487779048189210000d0)
      Anx(2)=-dsin(1.461364298564610000d0)
      Any(2)=dcos(1.461364298564610000d0)
      Anx(3)=-dsin(1.341508350980550000d0)
      Any(3)=dcos(1.341508350980550000d0)
      Anx(4)=-dsin(0.985309944313188000d0)
      Any(4)=dcos(0.985309944313188000d0)
      Anx(5)=-dsin(0.446785799812043000d0)
      Any(5)=dcos(0.446785799812043000d0)
      Anx(6)=-dsin(0.156127161102117000d0)
      Any(6)=dcos(0.156127161102117000d0)
      Anx(7)=-dsin(0.057964534640109600d0)
      Any(7)=dcos(0.057964534640109600d0)
      Anx(8)=-dsin(0.025049497627388500d0)
      Any(8)=dcos(0.025049497627388500d0)
      Anx(9)=-dsin(0.012927537191682300d0)
      Any(9)=dcos(0.012927537191682300d0)
      Anx(10)=-dsin(0.008039942090901970d0)
      Any(10)=dcos(0.008039942090901970d0)

      do j=1,nel
         do i=1,ng

            call formaederi0
            call calc_A0 

            x=0
            y=0
            anx_=0
            any_=0
            do k=1,notl(j)
               nog=ic(j,k)!encontra numero global do n�x

               x=p0(glgl(nog,1))*fi(k)+x
               y=p0(glgl(nog,2))*fi(k)+y

c              normal x
               anx_=anx(nog)*fi(k)+anx_
c              normal y
               any_=any(nog)*fi(k)+any_

            enddo

            write(*,*) x,y,anx_,any_
         enddo
      enddo
      end

c     ********************************************************
c     Subrotina imprime tensao pontos de gauss
c     ******************************************************** 
      subroutine impreten
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      
      write(*,*) 'x', '    y', '       sx'

      do j=1,nel
c     atualizar configura��o nos elementos
      call gauss(ng3,qsi3,w3)

      call leielast

      call geometriaatual 

         do ig1=1,ng
   
               up=0.d0
               vp=0.d0
               x=0.d0
               y=0.d0
               
               xsi1=qsi(ig1)
               call formaederi0
               call calc_A0 
               

               do k=1,notl(j)
                  nog=ic(j,k)   !encontra numero global do n� k
                  x=p(glgl(nog,1))*fi(k)/p(glgl(nnos,1))+x
                  y=p(glgl(nog,2))*fi(k)+y
               enddo
               igf=1            ! a terceira fibra � a de cima
               ig3=1!ng3
               
c     xsi 3 para face inferior 
               xsi3=0.d0    !qsi3(ig3) -1.d0 = abaixo 1.d0=acima
               
               call calc_A0     !matriz gradiente do adimensional para inicial  
               
               up(1)=cos(salfa-pi/2.)
               up(2)=sin(salfa-pi/2.)
               
               vp(1)=cos(salfa)
               vp(2)=sin(salfa)
               
               call calc_Af     !matriz do adimensonal para o atual
               call cauchy      !gradiente da transformacao total
               
               call deformacoes !deformacoes longitudinais, distorcoes e
!alongamentos
               call energia
               call derivenertensaobspl !vetor de forcas internas
               !sx,sy,sxy(elemento,ptgauss)

           write(*,*) x,y,sigmax
         enddo
      enddo
    
      
      return
      end


c     ********************************************************
c     Subrotina calcula tensao n�
c     ******************************************************** 
      subroutine tennodal
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'
      

      do j=1,nel
c     atualizar configura��o nos elementos
         call gauss(ng3,qsi3,w3)
         
         call leielast
         
         call geometriaatual 
         
         do ig1=1,2
            
            up=0.d0
            vp=0.d0
            x=0.d0
            y=0.d0
            
            xsi1=2.d0*ig1-3.d0
            call formaederi0
            call calc_A0 
c     acima
c$$$            
c$$$            igf=3               ! a terceira fibra � a de cima
c$$$            ig3=3
c$$$            
c$$$            xsi3=1.d0
c     abaixo

            igf=3               ! a terceira fibra � a de cima
            ig3=3
            
            xsi3=1.d0

            call calc_A0        !matriz gradiente do adimensional para inicial  
            
            up(1)=cos(salfa-pi/2.)
            up(2)=sin(salfa-pi/2.)
            
            vp(1)=cos(salfa)
            vp(2)=sin(salfa)
            
            call calc_Af        !matriz do adimensonal para o atual
            call cauchy         !gradiente da transformacao total
            
            call deformacoes    !deformacoes longitudinais, distorcoes e
!     alongamentos
            call energia
            call derivenertensaobspl !vetor de forcas internas
!sx,sy,sxy(elemento,ptgauss)
c            if(itipoel(j).eq.2)then


               nog=ic(j,ig1+1)
               if (itipoel(j).eq.2)nog=ic(j,2*ig1-1)
               if (itipoel(j).eq.3)nog=ic(j,2*ig1)
             
c            else
c               nog=ic(j,ig1+1)
c            endif
            
            sigmaacima(nog,1)=sigmax
            sigmaacima(nog,2)=sigmay
            sigmaacima(nog,3)=sigmaxy
            
         enddo
      enddo
      
      
      return
      end



c     ************************************************
c     Subrotina que calcula as tens�es de engenharia
c     para p�s-processamento semelhante � derivener
c     ****************************************************
      Subroutine derivenertensaobspl
      USE BIG_STUFF
      implicit real*8(a-h),integer*4(i-n),real*8(o-z)
      include 'bloclagatbig.for'

c      duex, duey, duez etc sao as tensoes x, y, z etc el�sticas
c      fict�cias retirando destas as tens�es pl�sticas acumuladas
c      teremos a tens�o tentativa entrando na rotina de retorno, teremos
c      como resposta o n�vel de tens�o real que dever� ocupar a posi��o
c      dos duex, duey, duez etc no c�culo da forca interna.

c     c�lculo do valor da temperatura para o ponto de integra��o
      
      tpos=(a1*xsi3+fd(j,igf))*(temc(j)-temb(j))/2.+
     $(temc(j)+temb(j))/2.
      
      et(1)=rex-tpos*calt(j)
      et(2)=rey-tpos*0.
      et(3)=rez
      et(4)=rexy
      et(5)=rexz
      et(6)=reyz

c     CHAMADA DA SUBROTINA DE RETORNO na derivener
c     vamos corrigir o valor das tens�es calculadas

      do iep=1,6
         ee(iep)=et(iep)-ep(j,ig1,ig3,igf,iep)
      enddo
      
      s=0.
      do i=1,6
         do jr=1,6
            s(i)=s(i)+ceret(i,jr)*ee(jr)
         enddo
      enddo

      duex=s(1)
      duey=s(2)
      duez=s(3)
      duexy=s(4)
      duexz=s(5)
      dueyz=s(6)
     
      sigmax=duex
      sigmay=duey
      tauxy=duexy
      
      return
      end

