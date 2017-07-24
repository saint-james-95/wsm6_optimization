  PROGRAM wsm6
!-------------------------------------------------------------------
  USE module_mp_wsm6
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------

  INTERFACE
    subroutine readarray2(arr,arrname,unitno,ips,ipe)
      REAL,             INTENT(OUT) :: arr(:,:)
      CHARACTER(LEN=*), INTENT(IN)  :: arrname
      INTEGER,          INTENT(IN)  :: unitno
      INTEGER,          INTENT(IN)  :: ips,ipe
      REAL, ALLOCATABLE :: tmparr(:)
      INTEGER :: i,j,ij,jsize,CHUNK,ipn
    end subroutine readarray2

    subroutine readarray3(arr,arrname,unitno,ips,ipe)
      REAL,             INTENT(OUT) :: arr(:,:,:)
      CHARACTER(LEN=*), INTENT(IN)  :: arrname
      INTEGER,          INTENT(IN)  :: unitno
      INTEGER,          INTENT(IN)  :: ips,ipe
      REAL, ALLOCATABLE :: tmparr(:,:)
      INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn
    end subroutine readarray3

    subroutine readarray4(arr,arrname,unitno,ips,ipe)
      REAL,             INTENT(OUT) :: arr(:,:,:,:)
      CHARACTER(LEN=*), INTENT(IN)  :: arrname
      INTEGER,          INTENT(IN)  :: unitno
      INTEGER,          INTENT(IN)  :: ips,ipe
      REAL, ALLOCATABLE :: tmparr(:,:,:)
      INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn,m,msize
    end subroutine readarray4

    SUBROUTINE writearray2(arr,arrname,unitno,ips,ipe)
      REAL,             INTENT(IN) :: arr(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: arrname
      INTEGER,          INTENT(IN) :: unitno
      INTEGER,          INTENT(IN) :: ips,ipe
      REAL, ALLOCATABLE :: tmparr(:)
      INTEGER :: i,j,ij,jsize,CHUNK,ipn
    end subroutine writearray2

    SUBROUTINE writearray3(arr,arrname,unitno,ips,ipe)
      REAL,             INTENT(IN) :: arr(:,:,:)
      CHARACTER(LEN=*), INTENT(IN) :: arrname
      INTEGER,          INTENT(IN) :: unitno
      INTEGER,          INTENT(IN) :: ips,ipe
      REAL, ALLOCATABLE :: tmparr(:,:)
      INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn
    end subroutine writearray3

    SUBROUTINE writearray4(arr,arrname,unitno,ips,ipe)
      REAL,             INTENT(IN) :: arr(:,:,:,:)
      CHARACTER(LEN=*), INTENT(IN) :: arrname
      INTEGER,          INTENT(IN) :: unitno
      INTEGER,          INTENT(IN) :: ips,ipe
      REAL, ALLOCATABLE :: tmparr(:,:,:)
      INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn,m,msize
    end subroutine writearray4  
  end INTERFACE

  INTEGER ::   ids, ide,  jds, jde, kds,kde , &
               ims, ime,  jms, jme, kms,kme , &
               its, ite,  jts, jte, kts,kte
  ! "chunk" indices
  INTEGER ::  iids,iide, jjds,jjde,           &
              iims,iime, jjms,jjme,           &
              iits,iite, jjts,jjte
! timers
#include <gptl.inc>
  INTEGER :: ret
  REAL*8  :: totaltime
  INTEGER :: handle = 0
! Contains definition of _CHUNK_
#include <chunk.h>
! Contains definition of _NZ_
#include <literalk.h>
  INTEGER :: CHUNK, NZ_BUILD
  INTEGER :: num_tiles_C
  INTEGER ::  itimestep
  REAL, ALLOCATABLE :: t(:,:,:)
  REAL, ALLOCATABLE :: qci(:,:,:,:)
  REAL, ALLOCATABLE :: qrs(:,:,:,:)
  REAL, ALLOCATABLE :: q(:,:,:)
  REAL, ALLOCATABLE :: den(:,:,:)
  REAL, ALLOCATABLE :: p(:,:,:)
  REAL, ALLOCATABLE :: delz(:,:,:)
  REAL :: delt, g, rd, rv, t0c, den0, cpd, cpv, ep1, ep2, qmin, &
          XLS, XLV0, XLF0, cliq, cice, psat, denr
  REAL, ALLOCATABLE :: rain(:,:)
  REAL, ALLOCATABLE :: rainncv(:,:)
  REAL, ALLOCATABLE :: sr(:,:)
  REAL, ALLOCATABLE :: snow(:,:)
  REAL, ALLOCATABLE :: snowncv(:,:)
  REAL, ALLOCATABLE :: graupel(:,:)
  REAL, ALLOCATABLE :: graupelncv(:,:)

!+---+-----------------------------------------------------------------+
! LOCAL VAR
  INTEGER ::               i,j,k
  INTEGER :: ios, unitno
  CHARACTER(LEN=64) :: fn  ! file name

!+---+-----------------------------------------------------------------+
! EXTERNAL FUNCTIONS
#ifdef _OPENMP
  integer, external :: omp_get_max_threads
#endif
!+---+-----------------------------------------------------------------+

  call gptlprocess_namelist ('GPTLnamelist', 77, ret)
  ret = gptlinitialize ()
  ret = gptlstart('Total')

  ! read constants
  PRINT *,'wsm6init():  read constants'
  fn = 'wsm6_constants.dat'
  unitno=31
  open (unitno,file=trim(fn),form="unformatted",action='read', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open constant file ',trim(fn), &
               ' . stopping'
    stop
  endif
  read(unitno) pi, xlv1
  read(unitno) qc0, qck1
  read(unitno) bvtr1, bvtr2, bvtr3, bvtr4, bvtr6, g1pbr, g3pbr, &
               g4pbr, g6pbr, g5pbro2, pvtr, eacrr, pacrr, &
               precr1, precr2, roqimax
  read(unitno) bvts1, bvts2, bvts3, bvts4, g1pbs, g3pbs, g4pbs,  &
               g5pbso2, pvts, pacrs, precs1, precs2, pidn0r, pidn0s
  read(unitno) pacrc
  read(unitno) bvtg1, bvtg2, bvtg3, bvtg4, g1pbg, g3pbg, g4pbg,  &
               pacrg, g5pbgo2, pvtg, precg1, precg2, pidn0g
  read(unitno) rslopermax, rslopesmax, rslopegmax, rsloperbmax,  &
               rslopesbmax, rslopegbmax, rsloper2max, rslopes2max,  &
               rslopeg2max, rsloper3max, rslopes3max, rslopeg3max
  close(unitno)

  ! read input data
  PRINT *,'wsm62D():  read input state'
  fn = 'wsm6_input.dat'
  unitno=31
  open (unitno,file=trim(fn),form="unformatted",action='read', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open input file ',trim(fn), &
               ' . stopping'
    stop
  endif
  read(unitno) itimestep
  PRINT *,'wsm62D():  itimestep == ',itimestep
  ! read serial versions of indices
  read(unitno) ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte

  ! assert constraints on indices
  ! NOTE that [ikj]d[se] are ignored
  if ((ims/=its).or.(ime/=ite).or.(jms/=jts).or.(jme/=jte)) then
    print *,'ERROR:  index mismatch reading file ',trim(fn), &
            ' . stopping'
    stop
  endif
  if ((ims/=1).or.(jms/=1)) then
    print *,'ERROR:  incorrect start index reading file ',trim(fn), &
            ' . stopping'
    stop
  endif

  ! set default values of "chunk" indices
  iids = ids
  iide = ide
  iims = ims
  iime = ime
  iits = its
  iite = ite
  jjds = jds
  jjde = jde
  jjms = jms
  jjme = jme
  jjts = jts
  jjte = jte

  ! set up optional "i" chunking and optional fixed vertical extent
#ifdef _CHUNK_
   CHUNK = _CHUNK_
#else
   CHUNK = iite-iits+1
#endif
#ifdef _NZ_
   NZ_BUILD = _NZ_
   ! if specified, NZ_BUILD must match namelist setting for nvl
   if (NZ_BUILD/=kte-kts+1) then
     print *, 'ERROR:  Build-time-specified NZ must equal namelist nz, values are: ',NZ_BUILD,kte-kts+1
     call flush(6)
     stop
   endif
#else
   NZ_BUILD = kte-kts+1
#endif
   num_tiles_C = (iite-iits+1) / CHUNK
   if (mod((iite-iits+1),CHUNK) > 0) then
     num_tiles_C = num_tiles_C + 1
   endif
   iime = CHUNK
   iite = CHUNK
   jjme = num_tiles_C
   jjte = num_tiles_C
   PRINT *,'ims,ime,iims,iime',ims,ime,iims,iime
#ifdef _OPENMP
   PRINT "('omp_get_max_threads() returned:  ',I9)",omp_get_max_threads()
#endif
   PRINT *,'CHUNK =', CHUNK
   PRINT *,'NUMBER OF CHUNKS =', num_tiles_C
#ifdef _NZ_
   PRINT *,'NZ_BUILD =', NZ_BUILD
#endif

  ! allocate arrays
  ALLOCATE(t(iits:iite,kts:kte,jjts:jjte))
  ALLOCATE(qci(iits:iite,kts:kte,2,jjts:jjte))
  ALLOCATE(qrs(iits:iite,kts:kte,3,jjts:jjte))
  ALLOCATE(q(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(den(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(p(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(delz(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(rain(iims:iime,jjms:jjme))
  ALLOCATE(rainncv(iims:iime,jjms:jjme))
  ALLOCATE(sr(iims:iime,jjms:jjme))
  ALLOCATE(snow(iims:iime,jjms:jjme))
  ALLOCATE(snowncv(iims:iime,jjms:jjme))
  ALLOCATE(graupel(iims:iime,jjms:jjme))
  ALLOCATE(graupelncv(iims:iime,jjms:jjme))

!$OMP PARALLEL DO &
!$OMP PRIVATE ( j ) &
!$OMP SCHEDULE(runtime)
  do j = jjts,jjte
    CALL firstTouch(t(iits,kts,j), q(iims,kms,j)              &
               ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
               ,den(iims,kms,j)                               &
               ,p(iims,kms,j), delz(iims,kms,j)               &
               ,j                                             &
               ,rain(iims,j),rainncv(iims,j)                  &
               ,sr(iims,j)                                    &
               ,iids,iide, jjds,jjde, kds,kde                 &
               ,iims,iime, jjms,jjme, kms,kme                 &
               ,iits,iite, jjts,jjte, kts,kte                 &
               ,snow,snowncv                                  &
               ,graupel,graupelncv                            &
                                                              )
  enddo  ! j loop
!$OMP END PARALLEL DO

  ! read remaining input data
  call readarray3(t,'t',unitno,its,ite)
  call readarray4(qci,'qci',unitno,its,ite)
  call readarray4(qrs,'qrs',unitno,its,ite)
  call readarray3(q,'q',unitno,its,ite)
  call readarray3(den,'den',unitno,its,ite)
  call readarray3(p,'p',unitno,its,ite)
  call readarray3(delz,'delz',unitno,its,ite)
  read(unitno) delt
  print *,' delt = ',delt
  read(unitno) g
  print *,' g = ',g
  read(unitno) cpd
  print *,' cpd = ',cpd
  read(unitno) cpv
  print *,' cpv = ',cpv
  read(unitno) t0c
  print *,' t0c = ',t0c
  read(unitno) den0
  print *,' den0 = ',den0
  read(unitno) rd
  print *,' rd = ',rd
  read(unitno) rv
  print *,' rv = ',rv
  read(unitno) ep1
  print *,' ep1 = ',ep1
  read(unitno) ep2
  print *,' ep2 = ',ep2
  read(unitno) qmin
  print *,' qmin = ',qmin
  read(unitno) XLS
  print *,' XLS = ',XLS
  read(unitno) XLV0
  print *,' XLV0 = ',XLV0
  read(unitno) XLF0
  print *,' XLF0 = ',XLF0
  read(unitno) cliq
  print *,' cliq = ',cliq
  read(unitno) cice
  print *,' cice = ',cice
  read(unitno) psat
  print *,' psat = ',psat
  read(unitno) denr
  print *,' denr = ',denr
  call readarray2(rain,'rain',unitno,its,ite)
  call readarray2(rainncv,'rainncv',unitno,its,ite)
  call readarray2(sr,'sr',unitno,its,ite)
  call readarray2(snow,'snow',unitno,its,ite)
  call readarray2(snowncv,'snowncv',unitno,its,ite)
  call readarray2(graupel,'graupel',unitno,its,ite)
  call readarray2(graupelncv,'graupelncv',unitno,its,ite)
  close(unitno)

  ! minimize timer overhead inside OpenMP loop
  ret = gptlinit_handle ('WSM62D', handle)
  ! call WSM6
  ret = gptlstart('WSM62D+OpenMP')
!$OMP PARALLEL DO &
!$OMP PRIVATE ( j,ret ) &
!$OMP SCHEDULE(runtime)
  do j = jjts,jjte
    ret = gptlstart_handle ('WSM62D', handle)
    CALL wsm62D(t(iits,kts,j), q(iims,kms,j)                  &
               ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
               ,den(iims,kms,j)                               &
               ,p(iims,kms,j), delz(iims,kms,j)               &
               ,delt,g, cpd, cpv, rd, rv, t0c                 &
               ,ep1, ep2, qmin                                &
               ,XLS, XLV0, XLF0, den0, denr                   &
               ,cliq,cice,psat                                &
               ,j                                             &
               ,rain(iims,j),rainncv(iims,j)                  &
               ,sr(iims,j)                                    &
               ,iids,iide, jjds,jjde, kds,kde                 &
               ,iims,iime, jjms,jjme, kms,kme                 &
               ,iits,iite, jjts,jjte, kts,kte                 &
               ,snow,snowncv                                  &
               ,graupel,graupelncv                            &
                                                              )
    ret = gptlstop_handle ('WSM62D', handle)
  enddo  ! j loop
!$OMP END PARALLEL DO
   ret = gptlstop('WSM62D+OpenMP')

  ! write output data
  PRINT *,'wsm62D():  itimestep == ',itimestep,', write output state'
  fn = 'wsm6_output.dat'
  unitno=31
  open (unitno,file=trim(fn),form="unformatted",action='write', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open output file ',trim(fn), &
               ' . stopping'
    stop
  endif
  call writearray3(t,'t',unitno,its,ite)
  call writearray4(qci,'qci',unitno,its,ite)
  call writearray4(qrs,'qrs',unitno,its,ite)
  call writearray3(q,'q',unitno,its,ite)
  call writearray2(rain,'rain',unitno,its,ite)
  call writearray2(rainncv,'rainncv',unitno,its,ite)
  call writearray2(sr,'sr',unitno,its,ite)
  call writearray2(snow,'snow',unitno,its,ite)
  call writearray2(snowncv,'snowncv',unitno,its,ite)
  call writearray2(graupel,'graupel',unitno,its,ite)
  call writearray2(graupelncv,'graupelncv',unitno,its,ite)
  close(unitno)

  ! deallocate arrays
  DEALLOCATE(t)
  DEALLOCATE(qci)
  DEALLOCATE(qrs)
  DEALLOCATE(q)
  DEALLOCATE(den)
  DEALLOCATE(p)
  DEALLOCATE(delz)
  DEALLOCATE(rain)
  DEALLOCATE(rainncv)
  DEALLOCATE(sr)
  DEALLOCATE(snow)
  DEALLOCATE(snowncv)
  DEALLOCATE(graupel)
  DEALLOCATE(graupelncv)

  ret = gptlstop('Total')
  ret = gptlget_wallclock ('Total', 0, totaltime)  ! The "0" is thread number
  print*,''
  print*,'Total time =' , totaltime

!+---+-----------------------------------------------------------------+

  ! print timing info
  ret = gptlpr (0)
  ret = gptlpr_summary (0)

  END PROGRAM wsm6

