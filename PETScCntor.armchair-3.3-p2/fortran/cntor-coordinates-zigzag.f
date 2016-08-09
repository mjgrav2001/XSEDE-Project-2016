	program cntor coordinates
c	use numerical_libraries
	implicit real*8 (a-h,o-z)
	
	parameter(imax=24000)	
	parameter(idim=40)	
	parameter(nrng=imax/idim)
	
	parameter(imaxp1=imax+1)

      	parameter(idimp1=idim+1)
      	parameter(idimm1=idim-1)
      	parameter(idim2=2*idim)
      	parameter(idimd2=idim/2)
      		 		
      	parameter(nrngp1=nrng+1)
      	parameter(nrngm1=nrng-1)
      	parameter(nrng2=2*nrng)
 
 	parameter(imax1=imax-nrng+1)
 	parameter(imax2=imax-nrng+2)		
	parameter(nrngd1=nrng/2+1)
	parameter(nrngd2=nrng/2+2)
	parameter(imaxn1=imax-nrng/2+1)
	parameter(imaxn2=imax-nrng/2+2)
 
      	parameter(itor=nrng/2)
		 
	parameter(ig=4)
	parameter(ntofe=8001)	

      	complex*16 hmod(imax,imax)
      	complex*16 gleft(imax,imax)
	complex*16 work(imaxp1),ipiv(imax)
      	complex*16 gl(ig,ig),gr(ig,ig),html(ig,ig),htmr(ig,ig)
	complex*16 gaml(ig,ig),gamrr(ig,ig),gtl(ig,ig),gtr(ig,ig)
      	complex*16 unitz,bphs,test,tofe,trans
      	complex*16 phsarr(nrngd2,nrngd2),tb,tbc

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
C 'imax' carbon atoms, large rings with 'nrng' atoms, small rings 
C with 'idim' atoms, (m,0) zigzag torus:

c 180 degrees:
      	data isl180/1,2,imax1,imax2/
      	data isr180/nrngd2,nrngd1,imaxn2,imaxn1/
      			
c 180 degrees:
      	data isl/1,2,imax1,imax2/
      	data isr/nrngd2,nrngd1,imaxn2,imaxn1/
      	
c 180 degrees:
c      	data isl180/2000,1901,100,1/
c     	data isr180/1951,1950,51,50/

c      	data isl180/1901,1,2000,100/
c     	data isr180/1950,50,1951,51/
      			
c 180 degrees:
c      	data isl/2000,1901,100,1/
c      	data isr/1951,1950,51,50/	

c      	data isl/1901,1,2000,100/
c     	data isr/1950,50,1951,51/
c--------------------------------------------------------------------------------------- 
      		
c---------------------------------------------------------------------------------------
C 2000 carbon atoms, large rings with 100 atoms, small rings 
C with 20 atoms, (m,0) zigzag torus:

c 180 degrees (positions 50,51):
c      	data isl180/100,1,2000,1901/
c      	data isr180/50,51,1950,1951/

c Relative position of metallic leads (angle between leads; left lead always attached 
c at isl/100,1,2000,1901/)
c Example: 2000 sites with 20 C atoms per small ring
c -> 100 C atoms per large ring
		
c 180 degrees (positions 50,51):
c      	data isl/100,1,2000,1901/
c      	data isr/50,51,1950,1951/

c 90 degrees (layers 25,26):
c      	data isl/100,1,2000,1901/
c      	data isr/25,26,1925,1926/
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

c 16960 atoms, 40 atoms per small rings (width), 424 atoms per large ring, 
c zigzag torus, outer radius = bigr+smallr = 79.6876
cc	bigr   = 212/2*(1.42*3/2)/2/pi = 71.8680
cc	smallr = 40/2*1.42*sqrt(3)/2/pi = 7.8196
c	bigr   = 71.8680
c	smallr = 7.8196
	
c 24000 atoms, 40 atoms per small rings (width), 212 atoms per large ring, 
c zigzag torus, outer radius = bigr+smallr = 105.606
cc	bigr   = 600/2*(1.42*3/2)/2/pi = 101.700
cc	smallr = 40/4*1.42*sqrt(3)/2/pi = 3.90566
	bigr   = 101.700d0
	smallr = 3.90566d0
	
c example:

c Important (M. Jack): Original code had been written with bigr = 116 
c but value was adjusted to bigr = 117.4 to accurately reproduce 
c spacing between rings of 2.46 A at 1.42 A C-C bond length like 
c in M. Leamy's finite-element phonon modeling code; 07/23/11.

	bfield = 0.0d0
	 	
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
	vphon=0d0
		   
	e0=0.0d0
	dos=0.0d0
	sum=0.0d0
	corrsd=0.0d0
	
	ndm1=idim-1
	ndm2=idim-2

        print 12,bigr,smallr
 12     format('bigr = ',d10.4,1x,'A','   smallr = ',d10.4,1x,'A'/)

        print 13,vhop0,thop
 13     format('vhop0 = ',d10.4,1x,'eV','   thop = ',d10.4,1x,'eV'/)

        print 14,vsdmax
 14     format('vsd = ',d10.4,1x,'V'/)

        print 15,imax,nrng,idim
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


*******************************************************************************

	subroutine tpos(bigr,smallr,r,thetaphi)
	implicit real*8 (a-h,o-z)
	parameter (imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)	
	dimension r(imax,3)
	dimension thetaphi(imax,2)
	dimension temp1(nrng),temp2(nrng)
	dimension theta(imax),phi(imax)
      	common/nblk/ndm1,ndm2
      	
	xmax=1.0d0*imax
	xrng=1.0d0*nrng
	xrng15=3.0d0/2.0d0*xrng
	
	pi=dacos(-1.0d0)
	
	uphi=2.0d0*pi/xrng15
	utheta=pi/(xmax/xrng)
	
	nrngd2=nrng/2
	nrngd4=(nrng-4)/4
	
	do 100 k=1,nrngd2
	   temp1(2*k-1)=3*(k-1)*uphi
 	   temp1(2*k)=(3*k-2)*uphi
 100 	continue
 
	temp2(1)=0.0d0
	temp2(2)=utheta
	temp2(3)=utheta
	temp2(4)=0.0d0
	   
	do 150 k=1,nrngd4
	      do 150 i=1,4
		 temp2(4*k+i)=temp2(i)
 150	continue
	      
	      do 300 k=0,ndm1
		 do 300 i=1,nrng
		    phi(nrng*k+i) = -uphi/2.d0+temp1(i)
 300    	continue
		
        	do 350 k=0,ndm1
		   do 350 i=1,nrng
		      theta(nrng*k+i) = utheta/2.d0+k*2*utheta+temp2(i)
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
	parameter(imax=24000)
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
	parameter(imax=24000)
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
	parameter(imax=24000)
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
	parameter(imax=24000)
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
	parameter(imax=24000)
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
	parameter (imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
c	dimension ipiv(nrng)
	complex*16 work(nrng)
	complex*16 a(nrng,nrng),atmp(nrng,nrng),b(nrng,nrng)
	integer ipiv(nrng)
	integer info=1
	call copmat(a,atmp)

ccccccccccccccccccccccccccc
	
	call zgetrf(nrng, nrng, atmp, nrng, ipiv, info);
	call zgetri(nrng, atmp, nrng, ipiv, work, nrng, info);
	
ccccccccccccccccccccccccccc
	call copmat(atmp,b)
	return
	end

*****************************************************************************

	subroutine copmat(a,b)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 a(nrng,nrng),b(nrng,nrng)
	do 10 i=1,nrng
	do 10 j=1,nrng
	b(i,j)=a(i,j)
 10	continue
	return
	end

******************************************************************************

	subroutine mult5(z1,z2,z3,z4,z5,zres)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 z1(nrng,nrng),z2(nrng,nrng),z3(nrng,nrng),
     1             z4(nrng,nrng),z5(nrng,nrng)
	complex*16 prodx(nrng,nrng),prody(nrng,nrng),zres(nrng,nrng)

c start the matrix multiplications.

	call mult3(z1,z2,z3,prodx)
	call mult2(z4,z5,prody)
	call mult2(prodx,prody,zres)

	return
	end

******************************************************************************

	subroutine mult4(z1,z2,z3,z4,zres)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 z1(nrng,nrng),z2(nrng,nrng)
	complex*16 z3(nrng,nrng),z4(nrng,nrng)
	complex*16 prodl(nrng,nrng),prodr(nrng,nrng),zres(nrng,nrng)

c start the matrix multiplications.

	do 75 i=1,nrng
	   do 75 j=1,nrng
	      prodl(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,nrng
 76		 prodl(i,j)=z1(i,m)*z2(m,j)+prodl(i,j)
 75	      continue

        	do 85 i=1,nrng
        	do 85 j=1,nrng
        	prodr(i,j)=dcmplx(0.0d0,0.0d0)
        	do 86  m=1,nrng
 86		   prodr(i,j)=z3(i,m)*z4(m,j)+prodr(i,j)
 85		continue
		
		do 105 i=1,nrng
		   do 105 j=1,nrng
		      zres(i,j)=dcmplx(0.0d0,0.0d0)
		      do 106 m=1,nrng
 106			 zres(i,j)=prodl(i,m)*prodr(m,j)+zres(i,j)
 105		      continue
		      return
		      end

******************************************************************************

	subroutine mult3(t1,t2,t3,tres)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 t1(nrng,nrng),t2(nrng,nrng),t3(nrng,nrng)
	complex*16 prod12(nrng,nrng),tres(nrng,nrng)
	
	do 75 i=1,nrng
	   do 75 j=1,nrng
	      prod12(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,nrng
 76		 prod12(i,j)=t1(i,m)*t2(m,j)+prod12(i,j)
 75	      continue
	      
	      do 105 i=1,nrng
		 do 105 j=1,nrng
		    tres(i,j)=dcmplx(0.0d0,0.0d0)
		    do 106 m=1,nrng
 106		       tres(i,j)=prod12(i,m)*t3(m,j)+tres(i,j)
 105		    continue
		    return
		    end

**********************************************************************************

	subroutine mult2(za,zb,zout)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 za(nrng,nrng),zb(nrng,nrng),zout(nrng,nrng)
	do 75 i=1,nrng
	   do 75 j=1,nrng
	      zout(i,j)=dcmplx(0.0d0,0.0d0)
	      do 76  m=1,nrng
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
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 z1(nrng,nrng),z2(nrng,nrng),z3(nrng,nrng),
     1             z4(nrng,nrng),z5(nrng,nrng),zres(nrng,nrng)
	do 100 i=1,nrng
	do 100 j=1,nrng
 100	zres(i,j)=z1(i,j)+z2(i,j)+z3(i,j)+z4(i,j)+z5(i,j)
	return
	end

******************************************************************************

	subroutine add2(z1,z2,zres)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 z1(nrng,nrng),z2(nrng,nrng),zres(nrng,nrng)
	do 100 i=1,nrng
	do 100 j=1,nrng
 100	zres(i,j)=z1(i,j)+z2(i,j)
	return
	end

******************************************************************************

	subroutine prmute(a,b)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 a(nrng,nrng),b(nrng,nrng)
	do 100 i=1,nrng
	do 100 j=1,nrng
 100	b(j,i)=a(i,j)
	return
	end

*********************************************************************************

	subroutine hinit(hmat)
	implicit real*8 (a-h,o-z)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)
	complex*16 hmat(nrng,nrng)
	do 55 i=1,nrng
	   do 55 j=1,nrng
        	hmat(i,j)=dcmplx(0.0d0,0.0)
 55	     continue
	     return
	     end

**********************************************************************

	subroutine inipv(ipivg)
	implicit real*8 (a-h,o-z)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)	
      	parameter(nrngp1=nrng+1)
      	parameter(nrngm1=nrng-1)
      	parameter(nrng2=2*nrng)
	dimension ipivg(nrngp1)
	do 167 j=1,nrng+1
 167	   ipivg(j)=0
	   return
	   end

**********************************************************************************

	subroutine iniwg(workg)
	parameter(imax=24000)
	parameter(idim=40)
	parameter(nrng=imax/idim)	
      	parameter(nrngp1=nrng+1)
      	parameter(nrngm1=nrng-1)
      	parameter(nrng2=2*nrng)
	complex*16 workg(nrng2)
	do 233 j=1,nrng2
 233	   workg(j)=dcmplx(0.d0,0.d0)
	   return
	   end


********************************print blocks for testing*************************** 
c       
c       do 540 i=1,nrng
c       do 540 j=1,nrng
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
c       do 64,i=1,nrng
c       do 64,j=1,nrng
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
     & +c4*currtofe(j+3)+c5*currtofe(j+4))
    		
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
	
	















