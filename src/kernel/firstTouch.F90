  SUBROUTINE firstTouch(t, q                                      &   
                   ,qci, qrs, den, p, delz                        &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                                                                  )
!  USE module_mp_wsm6
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!  Mimic memory access patterns of wsm62D() while setting physics arrays 
!  to zero.  
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL, DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL, DIMENSION( its:ite , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci
  REAL, DIMENSION( its:ite , kts:kte, 3 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qrs
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
! LOCAL VAR
  INTEGER :: i, k
!
      do k = kts, kte
        do i = its, ite
          t(i,k) = 0.
          q(i,k) = 0.
          qci(i,k,1) = 0.
          qci(i,k,2) = 0.
          qrs(i,k,1) = 0.
          qrs(i,k,2) = 0.
          qrs(i,k,3) = 0.
          den(i,k) = 0.
          p(i,k) = 0.
          delz(i,k) = 0.
        enddo
      enddo
      do i = its, ite
        rain(i) = 0.
        rainncv(i) = 0.
        sr(i) = 0.
      enddo
      if (PRESENT(snow)) then
        do i = its, ite
          snow(i,lat) = 0.
        enddo
      endif
      if (PRESENT(snowncv)) then
        do i = its, ite
          snowncv(i,lat) = 0.
        enddo
      endif
      if (PRESENT(graupel)) then
        do i = its, ite
          graupel(i,lat) = 0.
        enddo
      endif
      if (PRESENT(graupelncv)) then
        do i = its, ite
          graupelncv(i,lat) = 0.
        enddo
      endif
  END SUBROUTINE firstTouch

