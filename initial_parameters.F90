      module condiciones
!      real*8,parameter :: ra_min=149.47, ra_max=150.77
!      real*8,parameter :: dec_min=1.62,  dec_max=2.83
      real*4,parameter  :: om=0.3,ol=0.7, h=100, c=299792.458
      real*4, parameter ::  zmin=0.009,  zmax=0.065
      real*4, parameter ::  R_infall = 2.5
      real*4, parameter ::  ylength=7, filwidth=1.3!/0.7
      real*4, parameter :: pi=4.0*atan(1.0)
      real*4, parameter ::  z_deep = 40., ang_den= 18
!      integer,parameter ::  corte_col = -1!or 0 or 1: 1  -1 for all, 1 for rojas, 0 por azules

      integer, parameter:: Num_ran= 7000000
      end module condiciones
 





