!***** Program: Selection the points in three dirence geometris in no Euler space  (comovil space)
!Purpose: Identifies groups and filaments of galaxies in a comovil space.  

!Inputs:  
!- initial_parameters.F90  ;
!- RAN_DES_n7mill_nside1024.dat;  randoms points in the comovil space
!- read_filaments_read: geometry kind filament (cylinder)
!Outputs:  
!- A file listing the geometry (circles, rings and cone) and their associated points.  

!Method:  
!- Implements the Gitano-Friends-of-Friends clustering algorithm with a specified linking length.  

!Requirements:  
!- Fortran 90 compiler (e.g., gfortran).  

!Usage:  
!1. Compile the program using:  
!   gfortran -o geo_selec_points.x geo_selec_points.f90  
!2. Run the program:  
!   ./geo_selec_points.x 

!Author: [Salerno Juan Manuel]  
!Date: June 2022 

!Notes:  
!- Future versions will include an adaptive linking length feature.



!-----------------------------------------------------------
    include 'programas/count_lines.f90'
    include 'programas/DLOCATE.FOR'
    include 'initial_parameters.F90'
    include 'programas/distanciaytheta.f90'
!-----------------------------------------------------------
    module ran
     integer ::  N_ran
     integer,allocatable,dimension(:) :: id_ran
     real*4,allocatable,dimension(:) :: ra_ran,dec_ran
     integer,allocatable,dimension(:) :: flag
     character ::color*4
     end module ran
!*******************************************************************************
!*******************************************************************************
      module fil_grid
        real :: xg_min,xg_max,yg_min,yg_max
        real,parameter :: lgrid=5. !Mpc/h comoving
        integer :: ngrid,nx,ny
        integer,allocatable,dimension(:) :: grp_1st,grp_nxt! <- 6df grp_1st,grp_nxt
        integer,allocatable,dimension(:) :: gal_1st,gal_nxt ! second nodo
        integer :: ngal
      end module fil_grid
!**************************
!-----------------------------------------------------------
    module fil_parametres
     integer :: N_fil
     integer, allocatable, dimension(:) ::  ID_fil
     integer, allocatable, dimension(:,:) :: Id_clu
     real*4,allocatable,dimension(:,:) :: ra_fil, dec_fil,Rvir_fil,z_fil,Factor
     character, allocatable,dimension (:) :: chori_fil*269
     real*4, allocatable,dimension (:) :: rp, Dc
     real*4, allocatable,dimension (:,:) :: x_nod
     integer, allocatable,dimension(:) :: flag_fil
     integer, allocatable, dimension(:,:) :: counts_ran
      character ::  salida*8
    end module fil_parametres

!*******************************************************************************
!*******************************************************************************
    program campact_filaments
     call read_ran
     call reader_fil
     call grid_limits
     call grid_ran
     call grupos
     call filaments
     call Infall
    end program campact_filaments
!***********************************************************
    subroutine read_ran
    use ran
    use condiciones, only : Num_ran
    implicit none
    real*4 dc1,da1,z5
    real*4,allocatable,dimension(:):: aa,bb
    integer :: i,j,k, ll, count_gal,ss,ss1,N

    
    open(10,file='RAN_DES_n7mill_nside1024.dat',status='old')
    call count_lines (10, N)
    
    allocate( aa(N),bb(N))
    N_ran=0 
    do i=1,Num_ran
     read(10,*)  aa(i), bb(i)
     if(bb(i) .LE. 5.995) N_ran=N_ran+1
    enddo

    allocate( ra_ran(N_ran),dec_ran(N_ran))
    N_ran=0

    do i=1,Num_ran
     if(bb(i) .LE. 5.995) then
     N_ran=N_ran+1
      ra_ran(N_ran)=aa(i)
     dec_ran(N_ran)=bb(i)
     endif
   enddo

    close(10)
    print*, N_ran
    end subroutine read_ran
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
   subroutine reader_fil
    use fil_parametres
    use condiciones
    implicit none
    integer, allocatable, dimension(:,:) :: counts
    real*4 dc1,da1,z5
    real*4 daa(200001),ez(200001),ddz,suma,zzz,sumav,zz
    integer d,dd,git
    real*4:: dh=c/h
    integer :: i,j,k, N
    real*4,allocatable,dimension(:,:) :: ra, dec,Rvir
    character, allocatable,dimension (:) :: chori*281
    real*4,allocatable,dimension(:) :: rp1,ov,z,largor,zz1
    real*4, allocatable,dimension (:,:) :: x_nod1
    real*4  :: da

    
   if(ang_den .gE. 19) then 
    salida='out_all/'
   elseif(ang_den .lT. 19) then
    salida='out_cut/'
   endif


    open(10,file='read_filaments_read.dat')
    call count_lines (10, N)
    allocate(chori(N),ra(2,N),dec(2,N))
    allocate(Rvir(2,N))
    allocate(rp1(N),ov(N),largor(N))
    allocate(x_nod1(6,N),z(N),zz1(N))
   allocate( counts(2,N) )
    do i=1,N 
     read(10,'(a281)') chori(i)
     read(chori(i),'(19x,f10.6,1x,f10.6)') ra(1,i),dec(1,i)
     read(chori(i),'(140x,f10.6,1x,f10.6)') ra(2,i),dec(2,i)
     read(chori(i),'(86x,f5.3)') Rvir(1,i)
     read(chori(i),'(207x,f5.3)') Rvir(2,i)
     read(chori(i),'(249x,f5.2)') rp1(i)
     read(chori(i),'(98x,3(f8.3,2x))') x_nod1(1:3, i)
     read(chori(i),'(220x,3(f8.3,2x))') x_nod1(4:6, i)
     read(chori(i),'(262x,f6.3)') ov(i)
     read(chori(i),'(63x,f7.5)') z(i)
     read(chori(i),'(184x,f7.5)') zz1(i)
     read(chori(i),'(272x,I3,3x,I3)')  counts(1:2,i)


    N_fil=0
    do i=1,N
     if(ov(i) .ge. 2.0 .and. (ra(1,i) < 150 .or. ra(1,i) >290) .and. z(i) < zmax ) N_fil=N_fil+1
    enddo
    print*, 'Numero de filamentos empleados', N_fil
    allocate(chori_fil(N_fil),ra_fil(2,N_fil),dec_fil(2,N_fil))
    allocate(Rvir_fil(2,N_fil),Id_fil(N_fil), Id_clu(2,N_fil))
    allocate(rp(N_fil))
    allocate(x_nod(6,N_fil),z_fil(2,N_fil),Factor(2,N_fil),Dc(N_fil))
   allocate( counts_ran(2,N_fil) )
    N_fil=0

    do i=1,N
     if(ov(i) .ge. 2.0 .and.  (ra(1,i) < 150 .or. ra(1,i) >290) .and. z(i) < zmax ) then
      N_fil=N_fil+1
      chori_fil(N_fil)= chori(i)
      ra_fil(1,N_fil)= ra(1,i); dec_fil(1,N_fil)= dec(1,i)
      ra_fil(2,N_fil)= ra(2,i); dec_fil(2,N_fil)= dec(2,i)
      Rvir_fil(1,N_fil)= Rvir(1,i) ;  Rvir_fil(2,N_fil)=  Rvir(2,i)
      rp(N_fil)= rp1(i)	
      z_fil(1,N_fil)=z(i)
      z_fil(2,N_fil)=zz1(i)
      counts_ran(1,N_fil)= counts(1,i)
      counts_ran(2,N_fil)= counts(2,i)
      read(chori_fil(N_fil),'(3x,I4,2x,I4,117x,I4)') Id_fil(N_fil), Id_clu(1,N_fil), Id_clu(2,N_fil)
     end if
    enddo

    do i=1,N_fil
       if(ra_fil(1,i) > 150) ra_fil(1,i) = ra_fil(1,i) - 360
       if(ra_fil(2,i) > 150) ra_fil(2,i) = ra_fil(2,i) - 360
    enddo

    do i=1,N_fil
     do j=1,2
        z5=z_fil(j,i)
        ddz=0.00001  !no modificar
         zzz=0.
         suma=0
         do d=0,200000
          zzz=zzz+ddz
          if(zzz.gt.z5) then
           dd=d-1
           goto 34
          end if 
          ez(j)=(om*(1+zzz)**3.+ol)**.5
          suma=suma+ddz/ez(j)
          dc1=dh*suma
          daa(j)=dc1/(1+zzz)
         end do       

 34      da1=dh*suma/(1.+z5)
          dc1=dh*suma
         Dc(i)=dc1
         
         da=da1
         Factor(j,i)= (da/206.264806)*3.6
       enddo
       
         x_nod(1,i)=  Dc(i)*cos(dec_fil(1,i)*(pi/180.))*cos(ra_fil(1,i)*(pi/180.))
         x_nod(2,i) = Dc(i)*cos(dec_fil(1,i)*(pi/180.))*sin(ra_fil(1,i)*(pi/180.))
         x_nod(3,i) = Dc(i)*sin(dec_fil(1,i)*(pi/180.))
           
         x_nod(4,i)=  Dc(i)*cos(dec_fil(2,i)*(pi/180.))*cos(ra_fil(2,i)*(pi/180.))
         x_nod(5,i) = Dc(i)*cos(dec_fil(2,i)*(pi/180.))*sin(ra_fil(2,i)*(pi/180.))
         x_nod(6,i) = Dc(i)*sin(dec_fil(2,i)*(pi/180.))
       
       
    enddo

    
    close(10)
    end subroutine reader_fil


!*******************************************************************************
!***********************************************************
!*******************************************************************************
    subroutine grid_limits
    use fil_grid
    use fil_parametres, only:  ra_fil, dec_fil,N_fil
    use ran
    implicit none
    xg_min=min(minval(ra_ran), minval(ra_fil(1,1:N_fil)))      
    xg_max=max(maxval(ra_ran), maxval(ra_fil(1,1:N_fil)))
    write(*,*) xg_min,'< X <',xg_max
    yg_min=min(minval(dec_ran), minval(dec_fil(1,1:N_fil)) )
    yg_max=max(maxval(dec_ran), maxval(dec_fil(2,1:N_fil)) )
    write(*,*) yg_min,'< Y <',yg_max
    nx=int((xg_max-xg_min)/lgrid)+1
    ny=int((yg_max-yg_min)/lgrid)+1

    ngrid=nx*ny
    write(*,*) '# Grids X = ',nx
    write(*,*) '# Grids Y = ',ny
    write(*,*) '# Grids = ',ngrid
    end subroutine grid_limits

!*******************************************************************************

!*******************************************************************************
!*******************************************************************************
      subroutine grid_ran
      use ran
      use fil_grid, only:gal_1st,gal_nxt,ngrid
      implicit none
      integer :: i,j
      integer,external :: igrid
      real*4,dimension(2,N_ran) :: x_nod2

       allocate(gal_1st(ngrid),gal_nxt(N_ran))
      gal_nxt=0
      gal_1st=0
      
      do i=1,N_ran
      x_nod2(1,i) = ra_ran(i) 
      x_nod2(2,i) = dec_ran(i)
       j=igrid(x_nod2(1:2,i))
       gal_nxt(i)=gal_1st(j)
       gal_1st(j) = i           
      end do
      write(*,*) 'Galaxy grid ready'

     end subroutine grid_ran
!...............................................................................
!*******************************************************************************

!*******************************************************************************

!*******************************************************************************
   subroutine grupos
    use fil_parametres, only:  ra_fil, dec_fil,N_fil,flag_fil,Factor,rp,Rvir_fil,&
    Id_clu, counts_ran,salida
    use ran
    use condiciones
    use fil_grid
    implicit none
    integer :: i,j,k, region,inode,mm,inc,ig,ix,iy,ixmin,ixmax,iymin,iymax,ix0,iy0
    integer,dimension(N_fil) :: flag2
    real*4 ::  cx0,cx1
    integer :: ng_search, nn
    real*4,dimension(2) :: x0,dx,dxg,iv,jv,kv
    real*4,dimension(N_ran,N_fil ) :: cdist,cdist2
    real*4 :: x,y,z,xg,yg,zg,xmin,xmax,x1,x2
    real*4, dimension(2) :: x_n, y_n,Nu_n  
    real,external :: z_length, dot
    real*4, dimension(N_ran):: y_g,x_g,nu_g
    real*4 :: R_nor, dist_dentro, dist_dentro2
    integer, allocatable, dimension(:,:) :: counts
    integer, allocatable, dimension(:) :: counts2
    real*8 ::alf1,del1,alf2,del2,tita,theta,z1,z2

!primero buscamos las galaxias que caen dentro del radio virial de cada
!nodo

   if(corte_col .eq. 1) color="rojo" 
   if(corte_col .eq. 0) color="azul"
   if(corte_col .eq. -1) color="tota"

   open(12,file=''//salida//'ran_LSBG_groups_information_'//color//'.dat')

   allocate(flag_fil(N_fil), counts(2,N_fil), counts2(N_fil) )
   allocate(flag(N_ran))
   flag(1:N_ran)=500
   flag_fil(1:N_fil)=1000
   counts(1,1:N_fil)=0
   counts(2,1:N_fil)=0
   ng_search = 2
   z1=1.
   z2=1.0

   do j=1,N_fil
    if( rp(j)/ Factor(1,j) .le. ang_den ) then
     x0(1)= ra_fil(1,j)
     x0(2)=dec_fil(1,j)  
     ix0=int((x0(1)-xg_min)/lgrid)+1
     iy0=int((x0(2)-yg_min)/lgrid)+1
     ixmin=max(1,ix0-ng_search) 
     ixmax=min(nx,ix0+ng_search) 
     iymin=max(1,iy0-ng_search) 
     iymax=min(ny,iy0+ng_search)     

     do ix=ixmin,ixmax
       do iy=iymin,iymax
         ig = ix+(iy-1)*nx
          k = gal_1st(ig)
         do while(k/=0)      
           cx0=ra_fil(1,j); cx1=dec_fil(1,j)
           alf1= cx0;    alf2= ra_ran(k)
           del1= cx1 ;  del2= dec_ran(k)
!           call distanciaytheta(alf1,del1,z1,alf2,del2,z2,theta)       
           cdist(k,j)=sqrt((cx0-ra_ran(k))**2.+ (cx1-dec_ran(k))**2.)!theta

            if( flag(k) .eq. 500 .and.  counts_ran(1,j) .GT. 0 ) then
             if (cdist(k,j) .LE.  Rvir_fil(1,j)/Factor(1,j) ) then
              flag(k)=100
              write(12, 200 ) ra_ran(k),dec_ran(k),Factor(1,j),Id_clu(1,j)
              flag_fil(j)=2000
              counts(1,j) = counts(1,j) + 1          
            endif
          endif
          k=gal_nxt(k)
        enddo
       enddo
      enddo
      
      x0(1)= ra_fil(2,j)
      x0(2)=dec_fil(2,j)  
      ix0=int((x0(1)-xg_min)/lgrid)+1
      iy0=int((x0(2)-yg_min)/lgrid)+1
      ixmin=max(1,ix0-ng_search) 
      ixmax=min(nx,ix0+ng_search) 
      iymin=max(1,iy0-ng_search) 
      iymax=min(ny,iy0+ng_search)
      do ix=ixmin,ixmax
        do iy=iymin,iymax
         ig=ix+(iy-1)*nx
         k=gal_1st(ig)
         do while(k/=0)      
          cx0=ra_fil(2,j); cx1=dec_fil(2,j)
           alf1= cx0;    alf2= ra_ran(k)
           del1= cx1 ;  del2= dec_ran(k)
!           call distanciaytheta(alf1,del1,z1,alf2,del2,z2,theta)       
           cdist2(k,j)=sqrt((cx0-ra_ran(k))**2.+ (cx1-dec_ran(k))**2.)!theta
           if (cdist2(k,j) .LE.  Rvir_fil(2,j)/Factor(2,j) .and. flag(k) .eq. 500 .and.  counts_ran(2,j) .GT. 0 ) then
            flag(k)=100
            write(12, 200 )  ra_ran(k),dec_ran(k),Factor(2,j),Id_clu(2,j)
            flag_fil(j)=2000
            counts(2,j) = counts(2,j) + 1
         endif
       k=gal_nxt(k)
       enddo
      enddo  
     enddo 
     endif 
    enddo

  

    close(12)

    
  200  format(2(F9.3,1x),1x, 1x, f7.3,1x,I5)

    end subroutine grupos
!*******************************************************************************
!*******************************************************************************
    subroutine filaments
    use fil_parametres
    use ran
    use condiciones
    implicit none
    integer :: i,j,k, region,inode,mm,inc
    real*4, dimension(N_ran):: y_g,x_g,nu_g
    real*4 ::  cx0,cx1,xlength,zlength
    real,external :: z_length, dot
    real*4,dimension(3) :: x0,dx,dxg,iv,jv,kv
    integer, allocatable, dimension(:,:) :: counts
    integer, allocatable, dimension(:) :: counts2
    integer,dimension(N_fil) :: flag2
    real*4 :: x,y,z,xg,yg,zg,xmin,xmax,x1,x2

    allocate( counts(2,N_fil), counts2(N_fil) )

    counts2=0
    flag2(1:N_fil)=70
    open(10,file='ran_filaments_data.dat')

    do j=1,N_fil
      if( flag_fil(j) .eq. 2000 ) then
        x0(1:3)=0.5*(x_nod(1:3,j)+ x_nod(4:6,j))
        call versors(x_nod(1:3,j),x_nod(4:6,j),iv,jv,kv)
        zlength=z_length(x_nod(1:6,j))
        xmax=rp(j)/2.
        xmin=-xmax
        x1= Rvir_fil(1,j)+xmin
        x2=xmax- Rvir_fil(2,j)
        xlength=rp(j)

        do i=1, N_ran 
         if(  flag(i) .eq. 500  .or. flag(i) .eq. 400  ) then
		 x_g(i)  = Dc(j)*cos(dec_ran(i)*(pi/180.))*cos(ra_ran(i)*(pi/180.))
		 y_g(i)  = Dc(j)*cos(dec_ran(i)*(pi/180.))*sin(ra_ran(i)*(pi/180.))
		 nu_g(i) = Dc(j)*sin(dec_ran(i)*(pi/180.))
		    dx(1)= x_g(i)-x0(1); dx(2)=y_g(i)-x0(2);  dx(3)=nu_g(i)-x0(3)       
		     z=dot(dx,kv)
		        if(abs(z)<=zlength) then
		          y=dot(dx,jv)
		           if(abs(y)<=ylength) then
		            x=dot(dx,iv)
		            if(abs(x)<=xlength) then
		            region=2
		            if(x<=xmin.or.x>=xmax) region=1 
		            if(x>x1.and.x<x2) region=3
		             if(x<=0) then 
		                 inode=1
		             else
		                  inode=2
		             endif
		            if(abs(y)>=filwidth) region=1
		             if(region .GE. 2 ) then
		                flag2(j)=10
		                 counts2(j)=counts2(j)+1
		                 flag(i) = 400
      		                write(10,*) ra_ran(i),dec_ran(i),flag(i)
		              endif
		          end if
		        endif  
                end if
             endif
         enddo             	
       endif
    enddo



    close(13)
    close(10)
    
    
    end subroutine filaments

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
   subroutine Infall
   use fil_parametres
   use ran
   use condiciones
   use fil_grid
   implicit none
   integer :: i,j,k, region,inode,ng_search,N_ff,N_ff2
   integer :: ig,ix,iy,ixmin,ixmax,iymin,iymax,ix0,iy0
   integer,dimension(N_fil) :: flag2
   real*4 ::  cx0,cx1,xlength,zlength
   real*4,dimension(3) :: x0,dx,dxg,iv,jv,kv
   real*4,dimension(N_ran,N_fil ) :: cdist,cdist2
   real*4 ::  Rnor,Rnor2
   real*4 :: x,y,z,xg,yg,zg,xmin,xmax,x1,x2
   real*4, dimension(2) :: x_n, y_n,Nu_n  
   real,external :: z_length, dot
   real*4, dimension(N_ran):: y_g,x_g,nu_g
   real*4 :: dist_dentro, dist_dentro2
   integer, allocatable, dimension(:,:) :: counts
   integer, allocatable, dimension(:) :: counts2
   real*8 ::alf1,del1,alf2,del2,tita,theta,z1,z2
   z1=1.
   z2=1.0
    print*, N_ran, N_fil

 
    open(13,file='ran_infall.dat')
    allocate(counts(2,N_fil), counts2(N_fil) )
   
    counts(1,1:N_fil)=0
    counts(2,1:N_fil)=0
  
   ng_search = 7
   print*, R_infall
   do j=1,N_fil
    if( rp(j)/ Factor(1,j) .le. ang_den  .and. flag_fil(j) .eq. 2000 ) then 
     x0(1)= ra_fil(1,j)
     x0(2)=dec_fil(1,j)  
     ix0=int((x0(1)-xg_min)/lgrid)+1
     iy0=int((x0(2)-yg_min)/lgrid)+1
     ixmin=max(1,ix0-ng_search) 
     ixmax=min(nx,ix0+ng_search) 
     iymin=max(1,iy0-ng_search) 
     iymax=min(ny,iy0+ng_search)     
     do ix=ixmin,ixmax
      do iy=iymin,iymax
        ig=ix+(iy-1)*nx
        i=gal_1st(ig)
        do while(i/=0)     
         if(flag(i) .eq. 500 .or. flag(i) .eq. 300) then
          cx0=ra_fil(1,j); cx1=dec_fil(1,j)
          alf1= cx0;    alf2= ra_ran(i)
          del1= cx1 ;  del2= dec_ran(i)
 !         call distanciaytheta(alf1,del1,z1,alf2,del2,z2,theta)
          cdist(i,j)=sqrt((cx0-ra_ran(i))**2.+ (cx1-dec_ran(i))**2.)!theta
          if (cdist(i,j) .LE.  (R_infall *Rvir_fil(1,j) )/Factor(1,j) .and.  counts_ran(1,j) .GT. 0 ) then
           flag(i)=300
           write(13, 200 )  ra_ran(i),dec_ran(i),R_infall * Rvir_fil(1,j),Id_clu(1,j)
           counts(1,j) = counts(1,j) + 1
         endif
       endif
      i=gal_nxt(i)
      enddo
     enddo
    enddo 
   
  close(13)


    open(15,file=''//salida//'ran_infall_area_LSBG.dat')
 
    do j=1,N_fil
     if(counts(1,j) > 0 ) write(15,300) Id_clu(1,j), counts(1,j), Rvir_fil(1,j)/ Factor(1,j), ra_fil(1,j), dec_fil(1,j)
     if(counts(2,j) > 0 ) write(15,300) Id_clu(2,j), counts(2,j), Rvir_fil(2,j)/ Factor(2,j), ra_fil(2,j), dec_fil(2,j)
    enddo

      print*, '# infall', sum( counts(1,1:N_fil))


   


 200  format(2(F9.3,1x),1x,f5.2,2x,I4)
 300 Format(I4,1x, I6, f7.3, 1x, 2(1x,F7.3)  ) 
    end subroutine Infall
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************



!*******************************************************************************
!     Dot product between the vectors V1 and V2
      function dot(v1,v2)
!...............................................................................
      implicit none
!...............................................................................
      real :: dot
      real,dimension(3) :: v1,v2
!...............................................................................
      dot=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
!...............................................................................
      end function dot
!*******************************************************************************


!*******************************************************************************
!     Computes the projected distance at the baricentre
      function projected(v1,v2)
!...............................................................................
      implicit none
      real :: projected
      real,dimension(3) :: v1,v2,dv,iv,jv,kv
      real,external :: dot
!...............................................................................
      call versors(v1,v2,iv,jv,kv)
      dv=v2-v1
      projected=abs(dot(dv,iv))
!...............................................................................
      end function projected
!*******************************************************************************
!*******************************************************************************
!     Computes the i,j,k versors centred in the baricentre with
!     z-axis in the l.o.s. and x-axis alongside the projected separation
!     in the sky
      subroutine versors(v1,v2,i,j,k)
!...............................................................................
      implicit none
      real,dimension(3) :: v1,v2,i,j,k,xm,dx
      real,external :: dot
!...............................................................................
      dx=v2-v1
      xm=0.5*(v1+v2)
      call versor(xm,k)
      xm=dx-(dot(dx,k))*k
      call versor(xm,i)
      call cross(k,i,j)
!...............................................................................
      end subroutine versors
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!     Gives the versor V2 in the direction of V1
      subroutine versor(v1,v2)
!...............................................................................
      implicit none
      real,dimension(3) :: v1,v2
      real :: vv
!...............................................................................
      vv=sqrt(sum(v1**2))
      v2=v1/vv
!...............................................................................
      end subroutine versor

!*******************************************************************************
!     Vectorial product between vectors V1 and V2
      subroutine cross(v1,v2,v3)
!...............................................................................
      implicit none
      real,dimension(3) :: v1,v2,v3
!...............................................................................
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
!...............................................................................
      end subroutine cross
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!     Computes the maximum and minimum redshift (in Mpc/h) of the box 
!     within which we search for galaxies.
      function z_length(x)
!...............................................................................
      implicit none
      real :: z_length
      integer :: i,j
      real,dimension(6) :: x
      real :: z1,z2
!...............................................................................
      z1=sqrt(x(1)**2+x(2)**2+x(3)**2)
      z2=sqrt(x(4)**2+x(5)**2+x(6)**2)
      z_length=0.5*abs(z1-z2)+20.
!...............................................................................
      end function z_length
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

!     Computes grid number to a position vector X
      function igrid(x)

      use fil_grid, only: xg_min,yg_min, lgrid,nx,ny

      implicit none
      integer :: ix,iy,igrid
      real,dimension(2) :: x

        ix=int((x(1)-xg_min)/lgrid)+1 
        iy=int((x(2)-yg_min)/lgrid)+1 

        igrid=ix+(iy-1)*nx
!...............................................................................
      end function igrid
!*******************************************************************************
