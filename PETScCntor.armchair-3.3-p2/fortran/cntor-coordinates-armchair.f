	program cntor
c	use numerical_libraries
	implicit real*8 (a-h,o-z)
	
	parameter(imax=44480)	
	parameter(idim=40)
	
	parameter(ndiag=imax/idim)
	
	parameter(imaxp1=imax+1)
	 		
      	parameter(idimp1=idim+1)
      	parameter(idimm1=idim-1)
      	parameter(idim2=2*idim)
      	parameter(idimd2=idim/2)

  	parameter(imax1=imax-idim+1)
	parameter(imax2=imax-idim+2)
	parameter(imaxd21=imax/2+1)
	parameter(imaxd22=imax/2+2)
	parameter(imaxd21i=imax/2-idim+1)
	parameter(imaxd22i=imax/2-idim+2)	

      	parameter(itor=idim/2)
      		 
	parameter(ig=4)
	parameter(ntofe=8001)	

      	complex*16 hmod(imax,imax)
      	complex*16 gleft(imax,imax)
	complex*16 work(imaxp1),ipiv(imax)
      	complex*16 gl(ig,ig),gr(ig,ig),html(ig,ig),htmr(ig,ig)
	complex*16 gaml(ig,ig),gamrr(ig,ig),gtl(ig,ig),gtr(ig,ig)
      	complex*16 unitz,bphs,test,tofe,trans
      	complex*16 tb,tbc

      	dimension ic(itor),imaxc(itor)
      	dimension yl(ig),zl(ig),yr(ig),zr(ig),isl(ig),isr(ig),r(imax,3)
        dimension rnew(imax,3), index(imax)
        dimension thetaphi(imax,2)
	dimension yl180(ig),zl180(ig),yr180(ig)
	dimension zr180(ig),isl180(ig),isr180(ig)
 	dimension venergy(ntofe),vtofe(ntofe),currtofe(ntofe)		
 	     	        
	data zip,one,two /0.0d0,1.0d0,2.0d0/
	data info /1/

c---------------------------------------------------------------------------------------
C 'imax' carbon atoms, 'idim' C atoms in small ring, 
C 'nrng' C atoms in large rings, (m,m) armchair torus:

c 180 degrees (layers 51,50):
      	data isl180/1,2,imax1,imax2/
      	data isr180/imaxd21,imaxd22,imaxd21i,imaxd22i/
		
c 180 degrees (layers 51,50):
      	data isl/1,2,imax1,imax2/
      	data isr/imaxd21,imaxd22,imaxd21i,imaxd22i/
c--------------------------------------------------------------------------------------- 
 	
c---------------------------------------------------------------------------------------
C 3600 carbon atoms, (3,3) armchair torus:

c 180 degrees (layers 151,150):
c      	data isl180/1,2,3589,3590/
c      	data isr180/1801,1802,1789,1790/

c Relative position of metallic leads (angle between leads; left lead always attached 
c at isl/1,2,3589,3590/)
c Example: 3600 sites with 12 C atoms per layer/ring -> 300 layers
		
c 180 degrees (layers 151,150):
c      	data isl/1,2,3589,3590/
c      	data isr/1801,1802,1789,1790/

c 90 degrees (layers 76,75):
c     	data isl/1,2,3589,3590/
c    	data isr/901,902,889,890/
c---------------------------------------------------------------------------------------

c---------------------------------------------------------------------------------------
C 1800 carbon atoms, (3,3) armchair torus:

c 180 degrees (layers 76,75):
c      	data isl180/1,2,1789,1790/
c      	data isr180/901,902,889,890/

c Relative position of metallic leads (angle between leads; left lead always attached 
c at isl/1,2,1789,1790/)
c Example: 1800 sites with 12 C atoms per layer/ring -> 150 layers

c 180 degrees (layers 76,75):
c      	data isl/1,2,1789,1790/
c        data isr/901,902,889,890/

c 90 degrees (layers 38,37):
c     	data isl/1,2,1789,1790/
c      	data isr/457,458,445,446/
c---------------------------------------------------------------------------------------

c---------------------------------------------------------------------------------------
C 4000 carbon atoms, (5,5) armchair torus:

c 180 degrees (layers 1,100):
c      	data isl180/1,2,3981,3982/
c      	data isr180/2001,2002,1981,1982/

c Relative position of metallic leads (angle between leads; left lead always attached 
c at isl/1,2,3981,3982/)
c Example: 4000 sites with 20 C atoms per layer/ring -> 200 layers
		
c 180 degrees (layers 1,100):
c      	data isl/1,2,3981,3982/
c      	data isr/2001,2002,1981,1982/

c 90 degrees (layers 51,50):
c     	data isl/1,2,3981,3982/
c    	data isr/1001,1002,981,982/
c--------------------------------------------------------------------------------------- 
    
c---------------------------------------------------------------------------------------
C 2000 carbon atoms, (5,5) armchair torus:

c 180 degrees (layers 51,50):
c      	data isl180/1,2,1981,1982/
c      	data isr180/1001,1002,981,982/

c Relative position of metallic leads (angle between leads; left lead always attached 
c at isl/1,2,1981,1982/)
c Example: 2000 sites with 20 C atoms per layer/ring -> 100 layers
		
c 180 degrees (layers 51,50):
c      	data isl/1,2,1981,1982/
c      	data isr/1001,1002,981,982/

c 90 degrees (layers 26,25):
c     	data isl/1,2,1981,1982/
c    	data isr/501,502,481,482/
c--------------------------------------------------------------------------------------- 

c---------------------------------------------------------------------------------------
C 12000 carbon atoms, (10,10) armchair torus:

c 180 degrees (layers 151,150):
c      	data isl180/1,2,11961,11962/
c      	data isr180/6001,6002,5961,5962/
		
c 180 degrees (layers 151,150):
c      	data isl/1,2,11961,11962/
c      	data isr/6001,6002,5961,5962/

c 90 degrees (layers 76,75):
c     	data isl/1,2,11961,11962/
c    	data isr/3001,3002,2961,2962/
c---------------------------------------------------------------------------------------    
    
c---------------------------------------------------------------------------------------
C 4000 carbon atoms, (10,10) armchair torus:

c 180 degrees (layers 1,100):
c      	data isl180/1,2,3961,3962/
c      	data isr180/2001,2002,1961,1962/
		
c 180 degrees (layers 51,50):
c      	data isl/1,2,3961,3962/
c      	data isr/2001,2002,1961,1962/

c 90 degrees (layers 26,25):
c     	data isl/1,2,3961,3962/
c    	data isr/1001,1002,961,962/
c---------------------------------------------------------------------------------------

c---------------------------------------------------------------------------------------
C 2000 carbon atoms, (10,10) armchair torus:

c 180 degrees (layers 26,25):
c      	data isl180/1,2,1961,1962/
c      	data isr180/1001,1002,961,962/

c Relative position of metallic leads (angle between leads; left lead always attached 
c at isl/1,2,1961,1962/)
c Example: 2000 sites with 40 C atoms per layer/ring -> 100 layers
		
c 180 degrees (layers 26,25):
c      	data isl/1,2,1961,1962/
c      	data isr/1001,1002,961,962/

c 90 degrees (layers 13,12):
c     	data isl/1,2,1961,1962/
c    	data isr/501,502,461,462/
c---------------------------------------------------------------------------------------
 
	common/hbloc/hmod
	common/gbloc/gleft
      	common/rblk/bigr,smallr,r
      	common/rnewblk/rnew
      	common/thphi/thetaphi
      	common/nblk/ndm1,ndm2
	common/emubloc/emu1,emu2
	common/vblk/vtofe,venergy
	common/cblk/corrsd
	common/thblk/therm	
		
c       read(5,*)a,mtop,ntop
c       read(5,*)thop 
c	read(5,*)energy,deng,energymax
c       read(5,*)bigr,smallr
c       read(5,*)bfield
c       read(5,*)eta,dparam

	a = 2.0d0
	
	mtop = 100
	ntop = 100
	
	thop = -0.25d0

c 22400 atoms, 560 atoms per major ring, 40 atoms per minor ring, (10,10) armchair torus:
ccc outer radius R = bigr+smallr = 116.2395
ccc	bigr   = 280*1.42*sqrt(3)/2/pi = 109.474
ccc	smallr = 20*1.42*3/2/2/pi = 6.765
c	bigr = 109.474d0
c	smallr = 6.76481d0

c 44400 atoms, 1100 atoms per major ring, 40 atoms per minor ring, (10,10) armchair torus:
ccc outer radius R = bigr+smallr = 114.409
ccc	bigr   = 1100/4*1.42*sqrt(3)/2/pi = 107.643
ccc	smallr = 40/2*(1.42*3/2)/2/pi = 6.76481
	bigr = 107.643d0
	smallr = 6.76481d0
	
c 12000 atoms, 300 rings, 40 atoms per ring, (10,10) armchair torus:	
c	bigr = 117.4d0
c	smallr = 6.765d0

c 4000 atoms, 100 rings, 40 atoms per ring, (10,10) armchair torus:	
c	bigr = 39.15d0
c	smallr = 6.765d0

c 2000 atoms, 50 rings, 40 atoms per ring, (10,10) armchair torus:	
c	bigr = 19.58d0
c	smallr = 3.383d0
		
c 2000 atoms, 100 rings, 20 atoms per ring, (5,5) armchair torus:
c	bigr = 39.15d0
c	smallr = 3.383d0

c 3600 atoms, 300 rings, 12 atoms per ring, (3,3) armchair torus:		
c	bigr = 117.4d0
c	smallr = 2.0d0

c 1800 atoms, 150 rings, 12 atoms per ring, (3,3) armchair torus:		
c	bigr = 53.7d0
c	smallr = 2d0

c Important (M. Jack): Original code had been written with bigr = 116 
c but value was adjusted to bigr = 117.4 to accurately reproduce 
c spacing between rings of 2.46 A at 1.42 A C-C bond length like 
c in M. Leamy's finite-element phonon modeling code; 07/23/11.

	bfield = 0.0d0
		
c	eta = 5.0d-2
c	eta = 5.0d-3
c	eta = 5.0d-4
c	eta = 5.0d-5
c       eta = 5.0d-6
c       eta = 5.0d-7
c	eta = 5.0d-8
	eta = 5.0d-9
c       eta = 5.0d-10
c       eta = 5.0d-11
 	
	therm=0.03d0	
	        
	vsdmax=0.1d0
	vsdmin=0d0
	nvsd=100
		
	energymin = -0.2d0
	energymax =  0.2d0
	deng = 0.00005d0

	energy = energymin	

	dparam=0d0
	
	vhop0=-3.1d0
	vhop0t=vhop0
c	vhop0t=0d0
c	vhop=-3.1d0
	vphon=5.3d0
c	vphon=0d0
		   
	e0=0.0d0
	dos=0.0d0
	sum=0.0d0
	corrsd=0.0d0
	
	ndm1=ndiag-1
	ndm2=ndiag-2

        print 12,bigr,smallr
 12     format('bigr = ',d10.4,1x,'A','   smallr = ',d10.4,1x,'A'/)

        print 13,vhop0,thop
 13     format('vhop0 = ',d10.4,1x,'eV','   thop = ',d10.4,1x,'eV'/)

        print 14,vsdmax
 14     format('vsd = ',d10.4,1x,'V'/)

        print 15,imax,idim,ndiag
 15     format('N_torus = ',i6,1x,'   n_ring = ',i4,1x,
     &  '   n_diag = ',i4,1x/)

	print 16,energymax,energymin
 16     format('e_max = ',d10.4,1x,'eV','   e_min = ',d10.4,1x,'eV'/)

	print 17,deng
 17     format('del_e = ',d10.4,1x,'eV'/)

c build the array of nanotube positions and y,z planar coordinates for the wire/tube
c interface. this is y,z,coordinates are only approximate; the curvature of the tube
c is neglected.

c       sum2=0d0
c        
c       n=10000
c	
c       do 10 l=1,3
c       n=10*n		
c	call bodesrule2(energymin,energymax,n,sum2)
c	print 20,energymin,energymax,sum2
c 20	format(3x,d10.4,3x,d10.4,3x,d10.4)
c 10 	continue

	call tpos(bigr,smallr,r,thetaphi)
	
	OPEN(1335, FILE='ARRAY_thetaphi.dat', FORM='formatted')	
	OPEN(1336, FILE='ARRAY_posr.dat', FORM='formatted')
c	REWIND(1335)
c	REWIND(1336)
	
	do 22 i=1,imax
		write(unit=1335, fmt='(3x,f12.6,3x,f12.6,3x,i6)') 
     &          thetaphi(i,1),thetaphi(i,2),i		
		write(unit=1336, fmt='(3x,f12.6,3x,f12.6,3x,f12.6,3x,i6)') 
     &  	r(i,1),r(i,2),r(i,3),i
 22	continue
	CLOSE(1335)	
	CLOSE(1336)

c	OPEN(1337, FILE='ARRAY_posr_ph.dat', FORM='formatted')
c	REWIND(1337)
	
	OPEN(1337, FILE='ARRAY_posr.dat', FORM='formatted')
c	REWIND(1337)
					
 	do 24 i=1,imax  
    		read(1337, 25) rnew(i,1),rnew(i,2),rnew(i,3),index(i)	
 25 		format(3x,f12.6,3x,f12.6,3x,f12.6,3x,i6)		
 24	continue
	CLOSE(1337)
	
	OPEN(1338, FILE='ARRAY_posr_new.dat', FORM='formatted')
c	REWIND(1338)
				
	do 28 i=1,imax	
		write(unit=1338, fmt='(3x,f12.6,3x,f12.6,3x,f12.6,3x,i6)') 
     &          rnew(i,1),rnew(i,2),rnew(i,3),index(i)
 28	continue	
	CLOSE(1338)	


	CLOSE(1339)
 	CLOSE(1340)

	end
	
	
*******************************************************************************************************

	subroutine tpos(bigr,smallr,r,thetaphi)
	implicit real*8 (a-h,o-z)
	parameter (imax=44480)
	parameter (idim=40)
	dimension r(imax,3)
	dimension thetaphi(imax,2)
	dimension temp1(idim),temp2(idim)
	dimension theta(imax),phi(imax)
	common/nblk/ndm1,ndm2

	xmax=1.0d0*imax
	xdim=1.0d0*idim
	xdim15=1.5d0/2.0d0*xdim
	
	pi=dacos(-1.0d0)
	
	utheta=pi/xdim15
	uphi=pi/(xmax/xdim)
	
	idimd2=idim/2
	idimd4=(idim-4)/4
	
	do 100 k=1,idimd2
	   temp1(2*k-1)=3*(k-1)*utheta
 	   temp1(2*k)=(3*k-2)*utheta
 100 	continue
 
	temp2(1)=0.0d0
	temp2(2)=uphi
	temp2(3)=uphi
	temp2(4)=0.0d0
	   
	do 150 k=1,idimd4
	      do 150 i=1,4
		 temp2(4*k+i)=temp2(i)
 150	continue
	      
	      do 300 k=0,ndm1
		 do 300 i=1,idim
		    theta(idim*k+i)=temp1(i)
 300    	continue
		
        	do 350 k=0,ndm1
		   do 350 i=1,idim
		      phi(idim*k+i)=2*k*uphi+temp2(i)
 350		   continue
		   
		   do 400 i=1,imax
		      r(i,1)=(bigr+smallr*dcos(theta(i)))*dcos(phi(i))
		      r(i,2)=(bigr+smallr*dcos(theta(i)))*dsin(phi(i))
		      r(i,3)=smallr*dsin(theta(i))
		      thetaphi(i,1)=theta(i)
		      thetaphi(i,2)=phi(i)
 400		   continue
			
	
c	do 405 i=1,imax
c	print 410,i,j,r(i,1),r(i,2),r(i,3)
c       410	format(i4,2x,i4,3x,3(d10.4,3x))
c       405   continue
		      
		      return
		      end

********************************************************************************

	function delpos(i,j)
	implicit real*8 (a-h,o-z)
	parameter(imax=44480)
        dimension r(imax,3),rnew(imax,3)
	
	data acc/1.42d0/
	
	common/rblk/bigr,smallr,r
      	common/rnewblk/rnew
      	
      	del=(r(i,1)-r(j,1))**2d0 + (r(i,2)-r(j,2))**2d0 
     &  + (r(i,3)-r(j,3))**2d0
	del=sqrt(del)
      	
      	delnew=(rnew(i,1)-rnew(j,1))**2d0 + (rnew(i,2)-rnew(j,2))**2d0 
     &  + (rnew(i,3)-rnew(j,3))**2d0
	delnew=sqrt(delnew)
	
	delpos=delnew-del
 
	return
	end


********************************************************************************

	subroutine posyz(yl,zl,yr,zr,isl,isr)
	implicit real*8 (a-h,o-z)
	parameter(imax=44480)
	parameter(ig=4)
	dimension yl(ig),zl(ig),yr(ig),zr(ig),isl(ig),isr(ig)
	data acc/1.42d0/
	common/rblk/bigr,smallr,r(imax,3)
	
	do 10 i=1,ig
	   yl(i)=r(isl(i),2)
	   zl(i)=r(isl(i),3)
	   yr(i)=r(isr(i),2)
	   zr(i)=r(isr(i),3)
 
 10	continue
	
	return
	end

********************************************************************************

	subroutine posyz180(yl180,zl180,yr180,zr180,isl180,isr180)
	implicit real*8 (a-h,o-z)
	parameter(imax=44480)
	parameter(ig=4)
	dimension yl180(ig),zl180(ig),yr180(ig) 
	dimension zr180(ig),isl180(ig),isr180(ig)
	data acc/1.42d0/
	common/rblk/bigr,smallr,r(imax,3)
	
	do 10 i=1,ig
	   yl180(i)=r(isl180(i),2)
	   zl180(i)=r(isl180(i),3)
	   yr180(i)=r(isr180(i),2)
	   zr180(i)=r(isr180(i),3)
 
 10	continue
	
	return
	end

********************************************************************************

	function bphs(i,j,bfield)
	implicit real*8 (a-h,o-z)
	complex*16 bphs
	parameter(imax=44480)
	dimension r(imax,3)
	common/rblk/bigr,smallr,r
	
	arg=earg(r(i,1),r(i,2),r(i,3),
     1              r(j,1),r(j,2),r(j,3),bfield)
	bphs=dcmplx(dcos(arg),dsin(arg))
	return
	end
	
********************************************************************************
	
	function earg(xi,yi,zi,xj,yj,zj,bfield)
	implicit real*8 (a-h,o-z)
	dimension phihat(2)
	dimension rij(3)
	earg=0.0d0
	rhoij=(xi+xj)*(xi+xj)+(yi+yj)*(yi+yj)
	rhoij=.5d0*dsqrt(rhoij)
	if(rhoij.lt.0.00001)then
	   earg=0.0d0
	else
	   phihat(1)=-(yi+yj)/rhoij/2.0d0
	   phihat(2)= (xi+xj)/rhoij/2.0d0
	   rij(1)=xi-xj
	   rij(2)=yi-yj
	   do 100 i=1,2
	      earg=phihat(i)*rij(i)+earg
 100	   continue
	endif
	earg=.5d0*bfield*rhoij*earg
	earg=(1.52d-05)*earg
	return
	end

********************************************************************************
	function dweight(i,j,dparam)
	implicit real*8 (a-h,o-z)
	parameter(imax=44480)
	common/rblk/bigr,smallr,r(imax,3)
	dist=0.0d0
	do 200 k=1,3
 200	   dist=(r(i,k)-r(j,k))**2.0d0+dist
	   dweight=1.0d0+dparam*dsqrt(dist)
	   return
	   end

*********************************************************************************

	subroutine invert(a,b)
	implicit real*8 (a-h,o-z)
	parameter(idim=40)
c	dimension ipiv(idim)
	complex*16 work(idim)
	complex*16 a(idim,idim),atmp(idim,idim),b(idim,idim)
	integer ipiv(idim)
	integer info=1
	call copmat(a,atmp)
cc      call zgetrf(idim,idim,atmp,idim,ipiv,info)
cc      call zgetri(idim,atmp,idim,ipiv,work,idim,info)
cc	call dlincg(idim, atmp, idim, atmp, idim)

ccccccccccccccccccccccccccc
	
	call zgetrf(idim, idim, atmp, idim, ipiv, info);
	call zgetri( idim, atmp, idim, ipiv, work, idim, info);
	
ccccccccccccccccccccccccccc
	call copmat(atmp,b)
	return
	end

*****************************************************************************

	subroutine copmat(a,b)
	parameter(idim=40)
	complex*16 a(idim,idim),b(idim,idim)
	do 10 i=1,idim
	do 10 j=1,idim
	b(i,j)=a(i,j)
 10	continue
	return
	end

******************************************************************************

	subroutine mult5(z1,z2,z3,z4,z5,zres)
	parameter (idim=40)
	complex*16 z1(idim,idim),z2(idim,idim),z3(idim,idim),
     1             z4(idim,idim),z5(idim,idim)
	complex*16 prodx(idim,idim),prody(idim,idim),zres(idim,idim)

c start the matrix multiplications.

	call mult3(z1,z2,z3,prodx)
	call mult2(z4,z5,prody)
	call mult2(prodx,prody,zres)

	return
	end

******************************************************************************

	subroutine mult4(z1,z2,z3,z4,zres)
	parameter (idim=40)
	complex*16 z1(idim,idim),z2(idim,idim)
	complex*16 z3(idim,idim),z4(idim,idim)
	complex*16 prodl(idim,idim),prodr(idim,idim),zres(idim,idim)

c start the matrix multiplications.

	do 75 i=1,idim
	   do 75 j=1,idim
	      prodl(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,idim
 76		 prodl(i,j)=z1(i,m)*z2(m,j)+prodl(i,j)
 75	      continue

        	do 85 i=1,idim
        	do 85 j=1,idim
        	prodr(i,j)=dcmplx(0.0d0,0.0d0)
        	do 86  m=1,idim
 86		   prodr(i,j)=z3(i,m)*z4(m,j)+prodr(i,j)
 85		continue
		
		do 105 i=1,idim
		   do 105 j=1,idim
		      zres(i,j)=dcmplx(0.0d0,0.0d0)
		      do 106 m=1,idim
 106			 zres(i,j)=prodl(i,m)*prodr(m,j)+zres(i,j)
 105		      continue
		      return
		      end

******************************************************************************

	subroutine mult3(t1,t2,t3,tres)
	parameter (idim=40)
	complex*16 t1(idim,idim),t2(idim,idim),t3(idim,idim)
	complex*16 prod12(idim,idim),tres(idim,idim)
	
	do 75 i=1,idim
	   do 75 j=1,idim
	      prod12(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,idim
 76		 prod12(i,j)=t1(i,m)*t2(m,j)+prod12(i,j)
 75	      continue
	      
	      do 105 i=1,idim
		 do 105 j=1,idim
		    tres(i,j)=dcmplx(0.0d0,0.0d0)
		    do 106 m=1,idim
 106		       tres(i,j)=prod12(i,m)*t3(m,j)+tres(i,j)
 105		    continue
		    return
		    end

**********************************************************************************

	subroutine mult2(za,zb,zout)
	parameter (idim=40)
	complex*16 za(idim,idim),zb(idim,idim),zout(idim,idim)
	do 75 i=1,idim
	   do 75 j=1,idim
	      zout(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,idim
 76		 zout(i,j)=za(i,m)*zb(m,j)+zout(i,j)
 75	      continue
	      return
	      end

******************************************************************************

	subroutine tmult4(z1,z2,z3,z4,tofe)
	parameter (ig=4)
	complex*16 z1(ig,ig),z2(ig,ig),z3(ig,ig),z4(ig,ig)
	complex*16 prodl(ig,ig),prodr(ig,ig),zres(ig,ig),tofe

	tofe=dcmplx(0.0d0,0.0d0)
	
	do 75 i=1,ig
	   do 75 j=1,ig
	      prodl(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,ig
 76		 prodl(i,j)=z1(i,m)*z2(m,j)+prodl(i,j)
 75	      continue

	      do 85 i=1,ig
		 do 85 j=1,ig
		    prodr(i,j)=dcmplx(0.0d0,0.0d0)
		    do 86  m=1,ig
 86		       prodr(i,j)=z3(i,m)*z4(m,j)+prodr(i,j)
 85		    continue
		    
		    do 105 i=1,ig
		       do 105 j=1,ig
			  zres(i,j)=dcmplx(0.0d0,0.0d0)
			  do 106 m=1,ig
 106			     zres(i,j)=prodl(i,m)*prodr(m,j)+zres(i,j)
 105			  continue
	
			  do 115 i=1,ig
 115			     tofe=tofe+zres(i,i)
			     
			     return
			     end
	
******************************************************************************

	subroutine add5(z1,z2,z3,z4,z5,zres)
	parameter (idim=40)
	complex*16 z1(idim,idim),z2(idim,idim),z3(idim,idim),
     1             z4(idim,idim),z5(idim,idim),zres(idim,idim)
	do 100 i=1,idim
	do 100 j=1,idim
 100	zres(i,j)=z1(i,j)+z2(i,j)+z3(i,j)+z4(i,j)+z5(i,j)
	return
	end

******************************************************************************

	subroutine add2(z1,z2,zres)
	parameter (idim=40)
	complex*16 z1(idim,idim),z2(idim,idim),zres(idim,idim)
	do 100 i=1,idim
	do 100 j=1,idim
 100	zres(i,j)=z1(i,j)+z2(i,j)
	return
	end

******************************************************************************


	subroutine prmute(a,b)
	parameter(idim=40)
	complex*16 a(idim,idim),b(idim,idim)
	do 100 i=1,idim
	do 100 j=1,idim
 100	b(j,i)=a(i,j)
	return
	end

*********************************************************************************

	subroutine hinit(hmat)
	implicit real*8 (a-h,o-z)
	parameter(idim=40)
	complex*16 hmat(idim,idim)
	do 55 i=1,idim
	   do 55 j=1,idim
        	hmat(i,j)=dcmplx(0.0d0,0.0)
 55	     continue
	     return
	     end

**********************************************************************

	subroutine inipv(ipivg)
	implicit real*8 (a-h,o-z)
	parameter(idim=40)	
      	parameter(idimp1=idim+1)
      	parameter(idimm1=idim-1)
      	parameter(idim2=2*idim)
	dimension ipivg(idimp1)
	do 167 j=1,idim+1
 167	   ipivg(j)=0
	   return
	   end

**********************************************************************************

	subroutine iniwg(workg)
	parameter(idim=40)	
      	parameter(idimp1=idim+1)
      	parameter(idimm1=idim-1)
      	parameter(idim2=2*idim)
	complex*16 workg(idim2)
	do 233 j=1,idim2
 233	   workg(j)=dcmplx(0.d0,0.d0)
	   return
	   end


********************************print blocks for testing*************************** 
c       
c       do 540 i=1,idim
c       do 540 j=1,idim
c       print 541,i,j,gl(i,j),gr(i,j)
c       541	format(2x,i3,2x,i3,4x,4(d10.4))
c       540	continue
c
c       itest=1
c       if(itest.eq.0)then
c       do 997 i=1,imax
c       do 997 j=1,imax
c       testv=hmod(i,j)*dconjg(hmod(i,j))
c       if(testv.gt.zip)then
c       print 998,i,j,hmod(i,j)
c       998   format(2x,i3,2x,i3,3x,2(d10.3,3x))
c       endif
c       997   continue
c       endif
c
c       do 40 i=1,imax
c       do 40 j=1,imax
c       print 41,i,j,dreal(hmod(i,j))
c       41     format(2x,i3,2x,i3,2x,1(d10.4))
c       40     continue
c
c       do 64,i=1,idim
c       do 64,j=1,idim
c       print 74,i,j,dreal(gkk(i,j)),dreal(gk0(i,j)),
c       1               dreal(g0k(i,j)),dreal(g00(i,j))
c       74    	format(2x,i3,2x,i3,3x,4(d10.4,3x))
c       64    	continue
c       do 6 i=25,36
c       do 6 j=25,36
c       testv=gleft(i,j)*dconjg(gleft(i,j))
c       if(testv.ge.0.00000001)then
c       print 7,i,j,gleft(i,j)
c       7      	format(2x,i3,2x,i3,3x,2(d10.4,3x))
c       endif
c       6      	continue
c       
c       do 1050 i=1,12
c       do 1050 j=1,12
c       print 1055,i,j,dreal(gjj1(i,j)),dreal(gjj4(i,j))
c       1055	format(i3,3x,i3,3x,2(d10.4,3x))
c       1050	continue
c       
c       do 505 i=1,12
c       do 505 j=1,12
c       print 515, i,j,grnj(i,j)
c       515	format(i3,3x,i3,3x,2(d10.4,2x))
c       505	continue
c       
c        	do 500 i=1,imax
c        	do 500 j=1,imax
c        	print 510, i,j,dreal(biggr(i,j)),dreal(biggr(j,i))
c       510	format(i3,3x,i3,3x,2(d10.4,2x))
c       500 	continue
c       
c       do 300 i=1,12
c       do 300 j=1,12
c       print 44,i,j,gl(i,j),gr(i,j)
c       44    	format(i3,2x,i3,2x,4(d11.5,3x))
c       300   	continue
c       
c       do 500 i=1,imax
c       print 550, i, r(i,1),r(i,2),r(i,3)
c       550	format(i3,3x,3(d10.4,3x))
c       500 	continue


**********************************************************************************
	
	function fermi(e,emu)
	implicit real*8 (a-h,o-z)
 	common/thblk/therm	

 	fermi = 1.0d0/(exp((e-emu)/therm)+1.0d0)
	
	end

**********************************************************************************
	
	subroutine bodesrule(a,b)
	implicit real*8 (a-h,o-z)
	parameter(ntofe=8001)
	dimension venergy(ntofe),vtofe(ntofe),currtofe(ntofe)	  
	common/emubloc/emu1,emu2
	common/vblk/vtofe,venergy
	common/cblk/corrsd	

	do 45 k=1,ntofe
	      currtofe(k)=0d0
 45	continue
		      	
      	c1=14./45.

      	c2=64./45.

      	c3=24./45.

      	c4=c2

      	c5=c1

      	corrsd=0.0d0
      	
      	ntofem=ntofe-1
        ntofem4=ntofe-4

      	h=(b-a)/ntofem
      	
      	do 50 i=1,ntofe
      	        e=venergy(i)
      		fermi1=fermi(e,emu1)
      		fermi2=fermi(e,emu2)
c		currtofe(i)=e*e
		currtofe(i)=2d0*vtofe(i)*(fermi1-fermi2)
		if (abs(e).lt.h) currtofe(i)=0;				
c		print*,e
c		print 60,e,currtofe(i)
c 60		format(3x,d10.4,3x,d10.4)
 50     continue

      	do 100 j=1,ntofem4,4

      	corrsd=corrsd+h*(c1*currtofe(j)+c2*currtofe(j+1)+c3*currtofe(j+2)
     &+c4*currtofe(j+3)+c5*currtofe(j+4))
    		
c    	print 70,corrsd
c 70	format(3x,d10.4,3x)
   
 100  	continue
       

      	end


**********************************************************************************
	
	subroutine bodesrule2(a,b,n,sum2)
	implicit double precision (a-h,o-z)

      	c1=14./45.

     	c2=64./45.	

      	c3=24./45.

     	c4=c2

      	c5=c1

      	sum2=0.0d0

      	h=(b-a)/n

      	ns=n/4

      	x=a

      	do 100 j=1,ns

      	sum2=sum2+h*(c1*f(x)+c2*f(x+h)+c3*f(x+2.0*h)

     *         + c4*f(x+3.0*h)+c5*f(x+4.0*h))

      	x=x+4.0*h

 100  	continue

      	print*,sum2

      	end


      	function f(x)

      	implicit double precision (a-h,o-z)

      	f=x*x

     	return

      	end
	
	



















