      common /bloco0/ nnos,nel,nnr,nnc,n,ig1,ig2,nfc,ipt
     #,k2,li,lj,ig3,neltr,kkk,nzexp,nzexp1,kkka,nprint
      common /bloco17/ nelcs,ngs,jlado,idv,ig,ia,it,ih1
     #,igf,nelcar,nelnco,iposic,nnoel

      common /bloco1/ npt,ngf,ipc,j,noel,ng,ng3,ngf1,not
     	common /bloco2/ rnorma,tol

      common /bloco3/ t,t1,t2,ba,bb,bc,bd,be,bf,bg,bh,bri,bj,bk
     #,ca,cb,cc,cd,ce,cf,cg,crh,cri,cj,ck
     	common /bloco4/ dxg1,dxg2,dxp1,dxp2,dyg1,dyg2
     	common /bloco5/ dyp1,dyp2,rjac,rjac0,salfa,calfa,tf1,tf2,tfa,tfb
      common /bloco9/ pi,teta,rjac2,a1
      
      common /bloco11/ re1,re2,ren,rex,rey,rexy
     #,rez,rexz,reyz

	common /bloco14/ rr,alfa,bea,gam,alfa1,alfa2,beta1,beta2
     #,gama1,gama2,dgaxydx,dgaxydy,dgaxzdx,dgaxzdz,dgayzdy,
     #dgayzdz,d2gxydxdy,d2gxydydx,d2gxydxdx,d2gxydydy,
     #d2gxzdxdz,d2gxzdzdx,d2gxzdxdx,d2gxydzdz,
     #d2gyzdydz,d2gyzdzdy,d2gyzdydy,d2gyzdzdz


	common /bloco19/ velo,rbn,rgn,dt,xc,yc,xp,yp,rimp


      common /bloco20/ dpl(3,2),xsi1,xsi2,xsi3,qi3
     

      common /bloco21/ pxsi(3),peta(3),pdel(3),a0(2,2),af(2,2),
     #a0inv(2,2),aaux(3,3),ch(2,2),up(3),vp(3),wp(3),aup(3),aaux1(3,3)
      common /bloco22/ rco1,rco2,aup2(3),rco3
     #,drco1x,drco2x,drco1y,drco2y,d2ue(6,6),aval(2),acp(2,2),dcp(3,3)
     #,drco3x,drco3y,rjac1,marca(17),nordem(17)

c      common /bloco23/ due(6),daf1(3,3),d2af(3,3),
c     #aup1(3),rlex,rley,rlez,rlexy,aup12(3)
c     #,rlexz,rleyz,ham(7,3),wh(7)

      common /bloco23/ due(6),daf1(3,3),d2ch(3,3,6,6),d2af(3,3), aup1(3)
     $     ,rlex,rley,rlez,rlexy,aup12(3),ham(12,3),wh(12) ,rlexz,rleyz
     $     ,pp1(3,6,6),pp2(3,6,6),pp3(3,6,6),d2rex(6,6), d2rey(6,6)
     $     ,d2rez(6,6),d2rexy(6,6),d2rexz(6,6),d2reyz(6,6),sigmax,sigmay
     $     ,tauxy



      common /bloco24/ICNTL(30),INFO(20),CNTL(5),OPS,
     #LIW,LA,IFLAG,I,NSTEPS,MAXFRT,kac
      
	common/bloco25/n27,nz

 
c     blocos para a plasticidade e temperatura
      common /bloco26/  np,ip,iaret,nt,ntem
      common /bloco27/ s(6),rl(6),rhret(6),rt(6),ceret(6,6),cp(6,6)
     #,dpret(6,6),rpret(6,6),s0,tpos
     #,ds(6),rps(6),rph(6),rpt(6),ee(6),de(6,6),et(6)
     #,rpe(6),et1(6),rhd(2),rhv(2)
