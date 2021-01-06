!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_pk02

CONTAINS

!------------------------------------------------------------------
   SUBROUTINE PK02RAD(RTHRATEN,p8w,p_phy,pi_phy,dz8w,t_phy,          &
                     t8w, rho_phy, z_phy,R_d,G,CP,dt,xlat,degrad,  &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     )

!------------------------------------------------------------------
   IMPLICIT NONE
!------------------------------------------------------------------
   INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte  

   REAL, INTENT(IN    )      ::        DEGRAD

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                              RTHRATEN

   REAL, INTENT(IN   )   ::                   R_d,CP,G,dt

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                     dz8w, &
                                                             p8w, &
                                                           p_phy, &
                                                          pi_phy, &
                                                           t_phy, &
                                                             t8w, &
                                                         rho_phy, &
                                                           z_phy  
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN ) ::                                     xlat

   INTEGER :: i,j,K,NK
   real :: delty,delthez,p0,sec_p_d,sigb,kka,kks,kkf,rcp
   real :: ttmp,teq,sig,sigterm,kkt,t_tend
   real :: eps, p_strat, t_strat, w_lat, phi0, del_phi, t_us, t_pv, sigma, delta, gam, lapse

!------------------------------------------------------------------
! Newtonian relaxation scheme from Held and Suarez, Bull. Amer. Met.
! Soc., Vol. 75, No. 10., p1825-1830, 1994.  (box on page 1826)
! CEN and MIR  31-JUL-04
! +
! Add Stratospheric Structrue followed by Polvani and Kushner (2002, GRL, Vol. 29, No. 7., 10.1029/2001GL014284)

   delty   = 60.0

   !eps     = -10.0  !for southern-hemisphere winter
   eps     = 10.0   !for northern-hemisphere winter

   delthez = 10.0
   p0      = 100000.0
   sec_p_d = 86400.
   sigb    = 0.7
   kka     = 1.0/40.0   ! units of per day
   kks     = 0.25
   kkf     = 1.0
   rcp     = R_d/CP
 
   p_strat = 10000.0
   t_strat = 216.649 ! T(p_strat)=216.649 -> US Standard Atmos.
   phi0    = -50.
   del_phi = 10.
   !gam     = 0.  !no-Vortex
   gam     = 2. 
   lapse   = -gam*0.001 ! K/km

   !j_loop: DO J=jts,MIN(jte,jde-1)
   !k_loop: DO K=kts,MIN(kte,kde-1)
   !i_loop: DO I=its,MIN(ite,ide-1)
   j_loop: DO J=jts,jte
   k_loop: DO K=kts,kte
   i_loop: DO I=its,ite

      if (p_phy(i,k,j) >= p_strat) then 
          !for troposphere
          ttmp = 315.0 - delty*(sin(xlat(i,j)*degrad))**2.0- &
                   eps*(sin(xlat(i,j)*degrad))- &
                   delthez*alog(p_phy(i,k,j)/p0)*(cos(xlat(i,j)*degrad))**2.0
                   
          teq=max(t_strat,ttmp*(p_phy(i,k,j)/p0)**rcp) 
      else 
          !for stratosphere
          !w_lat = 0.5*(1.-tanh((xlat(i,j)-phi0)/del_phi))  !for southern hemisphere winter
          w_lat = 0.5*(1.-tanh((-xlat(i,j)-phi0)/del_phi))  !for northern hemisphere winter
     
          t_pv = t_strat*(p_phy(i,k,j)/p_strat)**(-R_d*lapse/G)

          call US_atmosphere(z_phy(i,k,j), sigma, delta, t_us) !t_us is ratio with 288.15          

          teq=(1.-w_lat)*(t_us*288.15) + w_lat*t_pv 

          !write (0,*) 'SHP-HS1:', w_lat, p_phy(i,k,j), p_strat
          !write (0,*) 'SHP-HS2:', t_strat, (p_phy(i,k,j)/p_strat), R_d, gam, G, (-R_d*gam/G), t_pv
          !write (0,*) 'SHP-HS3:', z_phy(i,k,j), sigma, delta, t_us
      end if


      sig=p_phy(i,k,j)/p8w(i,1,j)
      sigterm=max(0.0,(sig-sigb)/(1.0-sigb))

      kkt=kka+(kks-kka)*sigterm*(cos(xlat(i,j)*degrad))**4.0

      t_tend=-kkt*(t_phy(i,k,j)-teq)/sec_p_d  ! t_tend in kinetic K/s

      RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)

   ENDDO i_loop
   ENDDO k_loop
   ENDDO j_loop                                          

   END SUBROUTINE PK02RAD

!====================================================================
!+
SUBROUTINE US_atmosphere(alt, sigma, delta, theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
! AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
! NOTE - If alt > 86, the values returned will not be correct, but they will
!   not be too far removed from the correct values for density.
!   The reference document does not use the terms pressure and temperature
!   above 86 km.
  IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
  REAL,INTENT(IN)::  alt        ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma      ! density/sea-level standard density
  REAL,INTENT(OUT):: delta      ! pressure/sea-level standard pressure
  REAL,INTENT(OUT):: theta      ! temperature/sea-level standard temperature
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
  REAL,PARAMETER:: REARTH = 6369.0                 ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                     ! hydrostatic constant
  INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
  INTEGER:: i,j,k                                                  ! counters
  REAL:: h                                       ! geopotential altitude (km)
  REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
  REAL:: tlocal                                           ! local temperature
  REAL:: deltah                             ! height above base of this layer
!============================================================================
!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
!============================================================================
  REAL,DIMENSION(NTAB),PARAMETER:: htab= &
                          (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
  REAL,DIMENSION(NTAB),PARAMETER:: ttab= &
          (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
  REAL,DIMENSION(NTAB),PARAMETER:: ptab= &
               (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
                                     6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
  REAL,DIMENSION(NTAB),PARAMETER:: gtab= &
                                (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  i=1
  j=NTAB                                       ! setting up for binary search
  DO
    k=(i+j)/2                                              ! integer division
    IF (h < htab(k)) THEN
      j=k
    ELSE
      i=k
    END IF
    IF (j <= i+1) EXIT
  END DO

  tgrad=gtab(i)                                     ! i will be in 1...NTAB-1
  tbase=ttab(i)
  deltah=h-htab(i)
  tlocal=tbase+tgrad*deltah
  theta=tlocal/ttab(1)                                    ! temperature ratio

  IF (tgrad == 0.0) THEN                                     ! pressure ratio
    delta=ptab(i)*EXP(-GMR*deltah/tbase)
  ELSE
    delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
  END IF

  sigma=delta/theta                                           ! density ratio
  RETURN
END Subroutine US_atmosphere   ! -----------------------------------------------

END MODULE module_ra_pk02
