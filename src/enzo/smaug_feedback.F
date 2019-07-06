#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c

      subroutine star_maker_smaug(nx, ny, nz,
     &     d, dm, temp, u, v, w,
     &     dt, r, metal, dx, t, z, procnum,
     &     d1, x1, v1, t1,
     &     nmax, xstart, ystart, zstart, ibuff,
     &     imetal, imethod, tindsf,
     &     odthresh, useodthresh, masseff, smthresh, level, np,
     &     xp, yp, zp, up, vp, wp,
     &     mp, tdp, tcp, metalf,
     &     imetalSNIa, metalSNIa, metalfSNIa,
     &     imetalSNII, metalSNII, metalfSNII,
     &     mfcell)

c
c  CREATES STAR PARTICLES
c
c  Created 2014, adapted from star_maker2
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    tindsf   - if true, remove dt/t_dyn term in starfraction
c    odthresh - overdensity threshold (some number * avg. density)
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    level - current level of refinement
c    procnum - processor number (for output)
c    imetalSNIa - SN Ia metallicity flag (0 - none, 1 - yes)
c    mfcell - the maximum baryonic mass of a cell we allow to become a star particle.
c
c  OUTPUTS:
c
c    np   - number of particles created
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    metalfSNIa - metallicity fraction of particle (from SN Ia)
c    nmax     - particle array size specified by calling routine
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC imetalSNIa, procnum, tindsf, useodthresh
      INTG_PREC imetalSNII
      R_PREC    metalSNII(nx,ny,nz)
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    metalSNIa(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax),
     &     metalfSNIa(nmax), metalfSNII(nmax)
      R_PREC    odthresh, masseff, smthresh, mfcell
c
      R_PREC   sformsum
      save   sformsum
      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii
      INTG_PREC, save :: ran1_init = 0
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      R_PREC   random, pstar
      real*8 msolar
      parameter (pi=pi_val, G=GravConst,
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass)
c
      ii = np

c     Initialize RNG the first time this routine is called
c     Use MPI rank to seed the RNG, ensuring different streams
c     on each processor.
      if(ran1_init .eq. 0) then
         call init_random_seed(procnum)
         ran1_init = 1
      endif

c
c  for each zone, : "star" particle is created if the density is
c  greater than a critical density
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              1) is density greater than threshold (if we are using it)?
c
               if (useodthresh .eq. 1) then
                  if (d(i,j,k) .lt. odthresh) then
                     goto 10
                  endif
               else
c              Make sure we never use more than half the mass in the cell to
c              create a particle.  One can think of this as an implicit
c              resolution-dependent density threshold.
                  if (mfcell*bmass .lt. smthresh) then
                     goto 10
                  endif
               endif

c
c              2) Calculate useful intermediate values.
c

               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar
               isosndsp2 = sndspdC * temp(i,j,k)
               jeanmass = pi/(6._RKIND*sqrt(d(i,j,k)*dble(d1))) *
     &           dble(pi * isosndsp2 / G)**1.5_RKIND / msolar

c               if (jeanmass > bmass) then
c                  goto 10
c               endif

c
c              3) Calculate fraction of cell mass to convert into
c                 star particles.
c

               if (tindsf .eq. 1) then
                  starfraction = masseff
               else
                  starfraction = min(masseff*dt/tdyn, 1.0_RKIND)
               endif

c
c              4) If the star particle mass is below the minimum mass,
c                 form stars stochastically such that the average SFR
c                 is correct.
c
               if ((starfraction*bmass .lt. smthresh) .or.
     &             (useodthresh .ne. 1)) then
                  call random_number(random)
                  pstar = starfraction*bmass/smthresh
                  if (random .gt. pstar) then
                     goto 10
                  else
                     starfraction = smthresh/bmass
                  endif
               endif

               starfraction = min(starfraction, mfcell)

c
c              5) Create a star particle
c
               ii = ii + 1
               mp(ii)  = starfraction*d(i,j,k)
               tcp(ii) = t
               tdp(ii) = tdyn
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
               else
                  metalf(ii) = 0._RKIND
               endif

c              Metallicity from Type Ia SNe

               if (imetalSNIa .eq. 1) then
                  metalfSNIa(ii) = metalSNIa(i,j,k)    ! in here metal is a fraction
               endif
               if (imetalSNII .eq. 1) then
                   metalfSNII(ii) = metalSNII(i,j,k)
               endif
c
c              Remove mass from grid
c
               d(i,j,k) = (1.0_RKIND - starfraction)*d(i,j,k)

c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20

10          continue

            enddo
         enddo
      enddo
 20   continue
c
      if (ii .ge. nmax) then
         write(6,*) 'star_maker_ssn: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii
c
c      if (np .ne. 0) then
c         write(6,*) 'Stars created: number,time,level: ', np, t, level
c      endif
cc
      return
      end


c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine smaug_feedback(nx, ny, nz,
     &                       d, dm, te, ge, u, v, w,
     &                       metal, zfield1, zfield2,
     &                       idual, imetal, imulti_metals, imethod, 
     &                       dt, r, dx, t, z,
     &                       d1, x1, v1, t1, sn_param, m_eject, yield,
     &                       distrad, diststep, distcells, 
     &                       npart, xstart, ystart, zstart, ibuff,
     &                       xp, yp, zp, up, vp, wp,
     &                       mp, tdp, tcp, metalf, type, justburn)
c
c  RELEASES "STAR" PARTICLE ENERGY, MASS AND METALS
c
c  written by: Brian O,Shea
c  date:       13 November 2002
c  	this is a copy of star_maker2.src, and is designated
c	as of this writing to be for unigrid runs with star 
c	formation, and for experimentation.  See star_maker2.src
c	for the changelog prior to 13 November 2002.
c
c  modified1:  Britton Smith
c  date:       August 2010
c              Added distributed feedback where gas, metals, and 
c              energy are injected evenly over multiple cells.
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metallicity density field
c    r     - refinement field (0 if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imulti_metals - flag to use multi metals zfield 1 and 2
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    npart    - particle array size specified by calling routine
c    sn_param - fraction of stellar rest mass that goes to feedback
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c    justburn     - time-weighted mass of star formation (code units)
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod,
     &      distrad, diststep, distcells
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    metal(nx,ny,nz)
      R_PREC    zfield1(nx,ny,nz), zfield2(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart)
c
c  Locals
c    (msolar_e51 is one solar rest mass energy divided by 10^51 erg)
c
      INTG_PREC i, j, k, n, ic, jc, kc, stepk, stepj, cellstep
      R_PREC mform, tfactor, clight, energy, sn_param, msolar_e51,
     &     m_eject, yield, minitial, xv1, xv2, dratio,
     &     distmass
      parameter (clight = c_light, msolar_e51 = 1800._RKIND)
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
c      write(6,*) 'star_feedback3: start'
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2
     &           .and. ((t-tcp(n)) .gt. (10.0*3.155E13/t1)) ) then

c
c          change to negative
           tcp(n) = -1.0d0 * tcp(n)

c        The star particle creation algorithm partnered with this 
c          feedback algorithm creates a star particle instantaneously.
c          However, we do feedback as if the star particles are created 
c          over a long period of time (in code units), so the particle
c          actually loses mass over time in an exponentially decaying
c          way.

c        Determine how much of a given star particle would have been 
c          turned into stars during this timestep.  Then calculate the mass
c          which should have formed during this timestel dt using the integral
c          form of the Cen & Ostriker formula.

c
c         Compute index of the cell that the star particle
c           resides in.
c 
            i = int((xp(n) - xstart)/dx) + 1
            j = int((yp(n) - ystart)/dx) + 1
            k = int((zp(n) - zstart)/dx) + 1
c
c         check bounds - if star particle is outside of this grid
c         then exit and give a warning.
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'warning: star particle out of grid',i,j,k
               goto 100
            endif
c
c          skip if very little mass is formed.
c
c          calculate mass added to each cell
c
            mform    = mp(n)
            distmass = (mform * m_eject) / distcells
c
c          if using distributed feedback, check if particle is
c          too close to the boundary
c
            if (distrad .gt. 0) then
               i = max((1 + ibuff + distrad), 
     &              min((nx - ibuff - distrad), i))
               j = max((1 + ibuff + distrad), 
     &              min((ny - ibuff - distrad), j))
               k = max((1 + ibuff + distrad), 
     &              min((nz - ibuff - distrad), k))
            endif
c
c           subtract ejected mass from particle (ejection due
c           to winds, supernovae)
c
            mp(n) = mp(n) - mform * m_eject
c
c           Record amount of star formation in this grid.
c
            justburn = justburn + mform * dt * dx**3
c
c           Calculate how much of the star formation in this
c           timestep would have gone into supernova energy.
c
c            energy = sn_param * mform * (clight/v1)**2 / distcells


c           Needs to be an energy density in code units

            energy = 1.0d51 / (dx**3) / (d1*x1**3)/v1/v1 / distcells

c
c             (aug 7/00 addition: share ejected energy between gas
c              being ejected and that remaining in star particle)
c
#define NO_SHARE_ENERGY
#ifdef SHARE_ENERGY
            energy = energy*mform*m_eject/
     &              (mform*m_eject + minitial*exp(-xv2)*(1.0 + xv2))
#endif /* SHARE_ENERGY */
c
c           Add energy to energy field
c
            do kc = k-distrad,k+distrad
               stepk = abs(kc-k)
               do jc = j-distrad,j+distrad
                  stepj = stepk + abs(jc-j)
                  do ic = i-distrad,i+distrad
                     cellstep = stepj + abs(ic-i)
                     if (cellstep .le. diststep) then
                        dratio = 1.0/(d(ic,jc,kc) + distmass)
                        te(ic,jc,kc) = ((te(ic,jc,kc)*d(ic,jc,kc)) +
     &                       energy) * dratio
                        if (idual .eq. 1)
     &                       ge(ic,jc,kc) = 
     &                       ((ge(ic,jc,kc)*d(ic,jc,kc)) +
     &                       energy) * dratio
c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
                        if (imetal .eq. 1) then
c
c           "Cen method".  This takes into account gas recycling.
c
                           metal(ic,jc,kc) = 
     &                          (metal(ic,jc,kc)*d(ic,jc,kc) + 
     &                          (mform / distcells) * 
     &                          (yield * (1.0-metalf(n)) + 
     &                          m_eject * metalf(n))) * dratio
            !metal is a fraction
c
                           if (imulti_metals .eq. 1) then
                              zfield1(ic,jc,kc) = metal(ic,jc,kc)
                              zfield2(ic,jc,kc) = metal(ic,jc,kc)
                           endif
c
                        endif
c
c           Mass and momentum feedback
c
                        u(ic,jc,kc) = u(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * up(n)
                        v(ic,jc,kc) = v(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * vp(n)
                        w(ic,jc,kc) = w(ic,jc,kc)*d(ic,jc,kc) +
     &                       distmass * wp(n)
                        d(ic,jc,kc) = d(ic,jc,kc) + distmass
                        u(ic,jc,kc) = u(ic,jc,kc)/d(ic,jc,kc)
                        v(ic,jc,kc) = v(ic,jc,kc)/d(ic,jc,kc)
                        w(ic,jc,kc) = w(ic,jc,kc)/d(ic,jc,kc)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
                        if (imethod .ne. 2 .and. idual .eq. 1) then
                           te(ic,jc,kc) = 0.5*(u(ic,jc,kc)**2 + 
     &                          v(ic,jc,kc)**2 + w(ic,jc,kc)**2) +
     &                          ge(ic,jc,kc)
                        endif
                     endif
                  enddo
               enddo
            enddo
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
c      write(6,*) 'star_feedback3: end'
      return
      end
