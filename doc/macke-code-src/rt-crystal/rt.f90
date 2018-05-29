
! Short Description: 
! 'rt' calculates the reflection and transmission processes of a bundle of
! parallel and equidistant incoming rays at a polyedric shaped particle:
! A ray is defined by direction of propagation, direction of electric field
! vectors, and a Stokes vector, which defines intensity and polarization
! state. 
! Additionally, the diffraction at the crystals geometrical cross section
! is calculated for each crystal orientation. 

! Andreas Macke, Feb. '96

! Special thanks to Jouni Peltoniemi whose rotation routines for the 
! Stokes matrix have been adopted, and to Karri Muinonen for an important
! hint regarding the random-rotations at direct forward and backscattering.

! Diffraction now along the entire range of scattering angles by using 
! formular for inclination factor. 

! Variable angular bin size for ray tracing AND diffraction 
! sideproduct: diffraction at exact angels: ifn + '.diff_discrete'
! NOTE 1: This can be very time consuming for high angular resolution!!
! Therefore, go_every_n controls frequency of calls to polydiff
! NOTE 2: *.diff.discrete can differ from .diff if angular resolution is to small

!!$Reference: 
!!$@Article{Macke96,
!!$  author = 	 "Macke, A. and Mueller, J. and Raschke, E.",
!!$  title = 	 "Single scattering properties of atmospheric ice crystals",
!!$  journal =	 "J. Atmos. Sci."
!!$  year =	 1996,
!!$  volume =	 53,
!!$  pages =	 "2813-2825"
!!$}

! Input file (logical name 'ifn') must contain:
!   afn : name of crystal coordinates file
!   This crystal file must contain (see PDAL):
!     um : number of surfaces 
!     vm(1..,i,..,um) : number of points in surface i
!     p1(i,j), p2(i,j), p3(i,j) : edge coordinates of point j in surface i
!   la, rbi, ibi : wavelength (\mu m), real & imag. part of refractive index
!   rays : Number of incident rays per crystal orientation
!   orientmax : Number of crystal orientations
!   recdeepmax : Number of ray recursions
!   totreflmax : # of total internal reflections 
!   distortion : degree of crystal distortion, 0 .lt. distortion .lt. 1 
!   rtddim :  # of angular bins - 1
!   diff_seperate = 1 if diffraction NOT added to ray tracing part, 0 else
!   no_diff = 1 if diffraction is omitted, 0 else

! Output files:
! ifn + '.status' : runtime status each 10000 rays, optical scalars like
!                   cross sections, single scattering albedo and 
!                   asymmetry parameter.
! ifn + '.save'   : saves intermediate results (binary!). Read with 'get.save.f'
! ifn + '.angles' : scattering angles == center angles of angular bins
! 
! ifn + '.1dr'    : phase function excluding diffraction
! ifn + '.diff', ifn + 'diff_discrete'   : diffraction pattern 
! ifn + '.p(r)11'    : total phase function 
! ifn + '.p(r)12'    : degree of linear polarization
!        .p(r)22/p(r)33/p(r)34/p(r)44 phase matrix elements
! 'r' instead of 'p' if diffraction is excluded from final output.

! Questions, comments, ... -> 
!!$  Andreas Macke, Institut fuer Meereskunde, Abt. Maritime Meteorologie
!!$  Duesternbrooker Weg 20, 24105 Kiel, Germany
!!$  phone: +49-(0)431-597-3875, fax: +49-(0)431-565876
!!$  email: amacke@ifm.uni-kiel.de
!!$  anonymous ftp: ftp.ifm.uni-kiel.de/pub/macke
!!$  www: http://www.ifm.uni-kiel.de/me/research/Projekte/RemSens/RemSens-overview.html

program rt
implicit none

  integer(8), parameter ::   rtddim_max = 1800  ! max # of scatt. angle bins   
  integer(8), parameter ::   un = 30              ! max. # of plane surfaces 
  integer(8), parameter ::   vn = 10               ! max. # of points per surface 
  integer(8), parameter ::   raytreeim = 1000 ! max. # of subrays per incom. ray  
  integer(8), parameter ::   angle_max = 10000! -> max # of angles -> Fresnel & Snell discretization 
  integer(8), parameter ::   theta_diff_max = 900 ! # of zenith angles for diffr.
  integer(8), parameter ::   n_phi_diff = 90  ! # of azimuth angles for diffr.
  integer(8), parameter ::   nstokes = 4      !  # of Stokes vectors per ray
  integer(8), parameter ::   go_every_n = 1      !  controlls polydiff-calls

  type raytype                                       ! ray - pointer      
    real(8) in_stokes(1:4,1:4)
    real(8) es1, es2, es3, ep1, ep2, ep3, ek1, ek2, ek3 
    integer(8) utrold
    real(8) xold, yold, zold
    character(1) aorb
    integer(8) recdeep
    integer(8) totrefl
  end type raytype 
  type(raytype), dimension(raytreeim) :: raytree 

  real(8) fn1(1:un), fn2(1:un), fn3(1:un)               ! surface normal 
  real(8) lq(1:un), nk1(1:un), nk2(1:un), nk3(1:un)       ! distances 
  real(8) p1(1:un, 1:vn), p2(1:un, 1:vn), p3(1:un, 1:vn)  ! crystal points 
  real(8) p4(1:un, 1:vn), p5(1:un, 1:vn), p6(1:un, 1:vn)
  real(8) p7(1:un, 1:vn), p8(1:un, 1:vn), p9(1:un, 1:vn) 
  real(8) kano1(1:un, 1:vn), kano2(1:un, 1:vn), kano3(1:un, 1:vn) ! segment norm.
  real(8) rtd(0:rtddim_max,1:nstokes, 1:4) ! scattered Stokes vectors 
  real(8) pij(0:rtddim_max,1:4, 1:4) ! phase matrix
  real(8) diff_total(0:rtddim_max)    ! diffraction function 
  integer(8) umtr(1:un)                   ! lists hitted surfaces 
  integer(8) vm(1:un)                ! # of points per surface 
  real(8) xs1(1:un), xs2(1:un), xs3(1:un) ! intersections with all surfaces 
  real(8) refrwa(0:angle_max), refrwb(0:angle_max)      ! angle of refraction  
  real(8) ra(1:4,1:4,0:angle_max), ta(1:4,1:4,0:angle_max)
  real(8) rb(1:4,1:4,0:angle_max), tb(1:4,1:4,0:angle_max)  ! Fresnel coeff.
  real(8) rot_stokes(1:4,1:4), trans_r(1:4,1:4), trans_t(1:4,1:4)
  real(8) id_mat(1:4,1:4), rotate_out_stokes(1:4,1:4), rotate_in_stokes(1:4,1:4)
! transf. matrices for rotation, reflection and transmission of Stokes v. 
  real(8) istep_stat(0:100), ref_stat(0:100), utr_stat(0:100)
! for statistics: to analyse origin of scattering features at a given angle 
  real(8) in_stokes(1:4,1:4), new_in_stokes(1:4,1:4), refl_stokes(1:4,1:4), trans_stokes(1:4,1:4), out_stokes(1:4,1:4) ! stokes matrices
  real(8) wfa(0:angle_max), wfb(0:angle_max)                       !s. PINIT 
  real(8) cos_theta_t_a(0:angle_max), cos_theta_t_b(0:angle_max)    !s. PINIT 
  real(8) poly_x(1:1000), poly_y(1:1000) ! polygon of crystal projection   
  real(8) phi_diff(1:n_phi_diff), sqrint(1:n_phi_diff) ! az. int. in polydiff
  real(8) simpson_fak(0:2), simpson
  real(8) hvec1, hvec2, hvec3
  real(8) cross_1, cross_2, cross_3 ! cross product
  real(8) ac, bc, cc                            ! new hit point - old hit point 
  real(8)  co11, co12, co13, co14, co21, co22  ! matrix coefficients for determ. 
  real(8)  co23, co24, co31, co32, co33, co34    ! ray/surface intersection 
  real(8)  d1, d2, d3, d4                                    ! -> PGS 
  real(8)  ev11, ev12, ev13, ev21, ev22, ev23  ! electric field- and...  
  real(8)  ev31, ev32, ev33       ! ... propagation -vector for incident...
  real(8)  vk71, vk72, vk73, vk81, vk82, vk83
  real(8)  vk91, vk92, vk93   ! ...reflected, and...
  real(8)  vk101, vk102, vk103, vk111, vk112, vk113
  real(8)  vk121, vk122, vk123  ! ... refracted rays 
  real(8)  vk51, vk52, vk53                   ! local surface normal 
  real(8)  gwu    ! math. epsilon to avoid numerical errors 
  real(8)  nf                         ! normalizes direction vectors 
  real(8)  snk11, snk12, snk13, snk21, snk22, snk23  ! parameter for surfaces 
  real(8)  sp1, sp2, sp3                 ! crystals center of gravity 
  real(8)  vf             ! indicates wether Cramer gives solution or not
  real(8)  w                               ! dcos(angle of incidence) 
  real(8)  wd                           ! scattering angle in degree 
  real(8)  wf, wvz                                     ! -> PEG* 
  real(8)  xmax, xmin, ymax, ymin, zmax, zmin  ! extreme coordinates of crystal 
  real(8)  x1, x2, x3       ! solution of 3x3 lin eq. system via Cramer 
  real(8)  x_bound, y_bound, xy_bound  ! area overlapped by the inc. rays 
  real(8)  grenzw, cosgrenzw     ! angle for total reflection , dcos( :. ) 
  real(8)  mabs, absfak        ! absorption coeff., exp(mabs*pathlength) 
  real(8)  la, rbi, ibi ! wavelength, real(8)- and imag. part of refractive index 
  real(8)  waveno                         ! 2 wavenumber: pi/wavelength 
  real(8)  rayconst1                  ! 2/pi*angle_max: for binning 
  real(8)  diffconst1, diffconst2  ! 1/(/k^2), k^2/(4pi^2) -> diffraction 
  real(8)  diffconst3              ! rad*360/n_phi_diff -> diffraction 
  real(8)  poly_sum, gcs_mean_poly         ! projection area summation 
  real(8)  alpha_euler, beta_euler, gamma_euler ! Euler angles for crystal rot. 
  real(8)  xold, yold, zold           ! 'old' point of intersection 
  real(8)  gcs, gcs_mean, gcs_sum     ! 'scanning' geometrical cross section 
  real(8)  area                        ! crystal projection polygon surface 
  real(8)  direction ! 
  real(8)  before_sys, after_sys, before_user, after_user ! time ###
  real(8)  pi,twopi, pitwo                    ! pi, 2*pi, pi/2 
  real(8)  rad, deg                                 ! 180/pi, pi/180 
  real(8)  step_width                  ! distance between incoming rays 
  real(8)  help1, help2, help3, help4, help5, help6, help7, help8, test1
  real(8)  sum, sum1, sum2
  real(8)  wegl                            ! path length inside crystal 
  real(8)  myfunc_arg                  ! single argument for functions 
  real(8)  lq_max                  ! max length between crystal points 
  real(8)  energy_in, energy_out, energy_abs ! total inc,scatt,abs energies 
  real(8)  energy_violation ! adds up numerical errors in Fresnel-eq.s
  real(8)  distortion                              ! surface distortion 
  real(8)  resol_fak                    ! angular resolution for sampling 
  real(8)  w_0                           ! single scattering albedo 
  real(8)  g, g_rt, g_diff, g_diff_discrete      ! asymmetry parameter 
  real(8) lost_tr, lost_rm ! lost due to finite # of total refl, subrays 
  real(8) sin_rot, cos_rot ! rotation matrix for stokes vector
  real(8) w1, w2 ! random cos and sin
  real(8) e_0, e_pi, epsilon ! delta peak energy and epsilon
  real(8) la2 ! = la*la
  real(8) rad2resol_fak ! rad/resol_fak
  real(8) sintheta, costheta, angle, d_o1
  real(8) rand
  integer(8)  rtddim ! # of angular rt bins
  integer(8)  theta_diff ! = rtddim/2, diffraction bins in [0,90]
  integer(8)  i_glob, j_glob,  mainu, i_ray, k_glob, l_glob, m_glob, &
              m, i0 , h          ! global loop indices 
  integer(8)  i_sub, j_sub ! lokal 
  integer(8)  nt, ntm_sum, nterror, nerror          ! temporary ray counter  
  integer(8)  ntr, ntotal ! # of hits per orient., total # of hits
  integer(8)  nwd                         ! address of scattering angle in rtd 
  integer(8)  u, um                         ! Number of surface out of 1:um 
  integer(8)  utr, utrold                       ! number of hitted surface 
  integer(8)  vmu                            ! a surface under consideration 
  integer(8)  ntm, orient ! revised number of inc. rays per orientation 'orient' 
  integer(8)  orientmax, rays    ! max # of orientations, rays per orientation 
  integer(8)  rayhaben, raysoll  ! actual and total # of subrays per incom. ray 
  integer(8)  istep, recdeepMax, totreflmax  ! actual & maximum # of recursions and tital reflections for a subray 
  integer(8)  dummy !
  integer(8)  angle_i                            ! address for incident angle 
  integer(8)  angle_g                     !     "       angle of total reflection 
  integer(8)  x_range, y_range               ! determ. loop over incoming rays 
  integer(8)  n_poly           ! # of polygon points for projected crystal area 
  integer(8)  totrefl               ! number of total reflections of a subray
  integer(8)  pure_str ! get rid of trailing blanks  
  integer(8)  start_time, end_time
  integer(8)  poly_orient ! # of "good" orientations in polydiff 
  integer(8)  iangle
  integer(8)  diff_seperate, no_diff

  character(1) aorb
  character(100)  afn, ifn
  real(4) random_arg

  read(5,*) ifn ! input filename 
  h = 1
  do while (ifn(h:h) .ne. ' ')
    h = h + 1
  end do
  pure_str = h - 1 ! store for all other concatenations
  open(unit = 10, file = ifn(1:pure_str), status = 'old', form = 'formatted')

  call PDAL      ! Input of crystal geometry 
  call SDATIN    ! Input of wave characteristic 
  call PINIT ! Initialize variables, constants, Fresnel-coefficients, ...

  do orient  = 1, orientmax  ! crystal orientation loop 
    do i_glob  = 1, um
      do j_glob  = 1, vm(i_glob)
        p1(i_glob, j_glob) =  p7(i_glob, j_glob) ! get original orientation...
        p2(i_glob, j_glob) =  p8(i_glob, j_glob)
        p3(i_glob, j_glob) =  p9(i_glob, j_glob)
      end do
    end do
    call PSP ! ... and center of gravity...
    call PROT ! ...before rotating 
    call PEX ! extreme coordinates: xmax, xmin, ymax, ymin
    call PSP
    do mainu  = 1, um
      call PNK ! crystal facets coordinates
      call PFN ! surface normals fn1/2/3 (surface)
      call PPD ! kano1/2/3(line,surface)
    end do

!    call PDAS  ! write rotated crystal coordinates to 'rot.crystal'

    x_bound = xmax - xmin 
    y_bound = ymax - ymin
    xy_bound = x_bound*y_bound ! area, which is scanned by the rays 
    ntm_sum = ntm_sum + rays ! total # of incoming rays 
    ntr = 0 ! # of hits per orientation
 
    do i_ray = 1, rays  ! ray - loop, what else
      call random_number(rand)
      xold = xmin + rand*x_bound ! random point for initial ray
      call random_number(rand)
      yold = ymin + rand*y_bound
      zold = 0.0
    
        nt = nt + 1 ! temporary ray counter  

        raysoll = 0 ! actual position in raytree
        rayhaben = 1 ! number of entries in raytree 

! INCOMING RAY : 

        raytree(1)%in_stokes = id_mat
        raytree(1)%es1 = 1. ! e-perp
        raytree(1)%es2 = 0.
        raytree(1)%es3 = 0.
        raytree(1)%ep1 = 0. ! e-par
        raytree(1)%ep2 = 1.
        raytree(1)%ep3 = 0.
        raytree(1)%ek1 = 0. ! direction of propagation 
        raytree(1)%ek2 = 0.
        raytree(1)%ek3 = -1.
        raytree(1)%utrold = um + 1 ! last hitted surface 
        raytree(1)%xold = xold ! last hitted point 
        raytree(1)%yold = yold
        raytree(1)%zold = zold
        raytree(1)%aorb = 'A' ! A: ray is outside-, B: inside the crystal 
        raytree(1)%recdeep = 0 ! # of so far interactions do THIS ray 
        raytree(1)%totrefl = 0 ! # of total reflections    

        do while (raysoll .lt. rayhaben) ! raytree loop
          raysoll = raysoll + 1
          in_stokes = raytree(raysoll)%in_stokes
          ev11 = raytree(raysoll)%es1
          ev12 = raytree(raysoll)%es2
          ev13 = raytree(raysoll)%es3
          ev21 = raytree(raysoll)%ep1
          ev22 = raytree(raysoll)%ep2
          ev23 = raytree(raysoll)%ep3
          ev31 = raytree(raysoll)%ek1
          ev32 = raytree(raysoll)%ek2
          ev33 = raytree(raysoll)%ek3
          utrold = raytree(raysoll)%utrold 
          xold = raytree(raysoll)%xold
          yold = raytree(raysoll)%yold
          zold = raytree(raysoll)%zold
          aorb = raytree(raysoll)%aorb  
          istep = raytree(raysoll)%recdeep
          totrefl = raytree(raysoll)%totrefl

          if (istep .le. recdeepMax .and. totrefl .lt. totreflmax .and. in_stokes(1,1) .gt. 1.0e-10) then
            if (aorb .eq. 'A') then             ! external ray 
              call PTXA ! -> utr 
              if (utr .eq. 0 .and. istep .eq. 0) then ! missed crystals mark 
                rayhaben = raysoll             ! -> skip this ray 
              end if
              if (utr .gt. 0) then          ! hit 
                if (istep .eq. 0) then ! if incoming ray = initial ray then 
                  energy_in = energy_in + in_stokes(1,1)
                  ntr = ntr + 1 ! # of hits per orientation
                end if
                call PRTA ! Ray tracing for external rays 
              end if ! hit
              if (utr .eq. 0 .and. istep .gt. 0) then ! no more hits: 

!SAMPLING: 
                wd = dacos(-ev33) ! cos(theta) = 0*ev31 + 0*ev32 + -1*ev33 
                wd = wd * rtddim/pi
                nwd = idnint(wd) ! index for sampling array rtd(0:180)

! Get energy for delta forward and backward scattering
if (abs(ev33 + 1.0d+0) .lt. epsilon) e_0 = e_0 + out_stokes(1,1)
if (abs(ev33 - 1.0d+0) .lt. epsilon) e_pi = e_pi + out_stokes(1,1)

                if (nwd .eq. 0 .or. nwd .eq. rtddim) then 
                call random_rotation
                  hvec1 = cos_rot ! perp to incoming direction 
                  hvec2 =-sin_rot ! with random azimuth
                  hvec3 = 0.
                else
                  hvec1 = -ev32 
                  hvec2 =  ev31
                  hvec3 = 0.0  ! nur der Vollstaendigkeit halber
                  nf = sqrt(hvec1*hvec1 + hvec2*hvec2) 
                  hvec1 = hvec1/nf 
                  hvec2 = hvec2/nf
                end if

                call rotation_matrix(1d0, 0d0, 0d0, hvec1, hvec2, hvec3, 0d0, 0d0, -1d0)
! Rotate initial Stokes matrix into scattering plane

                rotate_in_stokes = rot_stokes

                call rotation_matrix(hvec1, hvec2, hvec3, ev11, ev12, ev13, ev31, ev32, ev33)
! 2nd rotation:
                rotate_out_stokes = rot_stokes

                out_stokes = matmul(rotate_out_stokes, matmul(in_stokes, rotate_in_stokes))

                do m_glob  = 1, nstokes  ! stokes vectors 
                  do k_glob  = 1, 4  ! stokes vect. components 
                    rtd(nwd, m_glob, k_glob) = rtd(nwd, m_glob, k_glob) + &
                                               out_stokes(m_glob, k_glob)
                  end do
                end do
                energy_out = energy_out + out_stokes(1,1)
              end if! no more hits 
            end if ! external ray 
            if (istep .le. recdeepMax .and. aorb .eq. 'B') then
              call PTXB                                    ! -> utr
              if (utr .ne. 0) call PRTB ! ### 
            end if
          end if ! recdeep loop

          if (istep .gt. recdeepmax)  lost_rm = lost_rm + in_stokes(1,1)
          if (totrefl .eq. totreflmax) lost_tr = lost_tr + in_stokes(1,1)

        end do ! do while raysoll < rayhaben 
    end do ! ray loop

    ntotal = ntotal + ntr ! total # of hits 
    gcs = xy_bound*ntr/rays ! actual 'scanning' geometrical cross section 
    gcs_sum = gcs_sum + gcs
    gcs_mean = gcs_sum/dble(orient) ! ### leave this here !!  
    if (no_diff .eq. 0) then 
      call PGCS     ! geometrical cross section area 
      if (mod(orient,go_every_n) .eq. 0) call polydiff  
! diffraction at geometrical cross section every go_every_n
    end if
    if (nt .gt. 10000) then 
      call PZS
    endif

  end do ! orientation-loop
  call PZS
  call PRTNORM
  call PZS

  contains ! for "association" of global variables
           ! in all subroutines:  

  subroutine GLS            ! Solve 3x3 linear equation system via Cramer 
   
    d1 = co11 * co22 * co33 + co21 * co32 * co13 + co31 * co12 * co23 - &
         co31 * co22 * co13 - co11 * co32 * co23 - co21 * co12 * co33
    if (d1 .eq. 0) then
      vf = 2
    else
      d2 = co14 * co22 * co33 + co24 * co32 * co13 + co34 * co12 * co23 - &
           co34 * co22 * co13 - co14 * co32 * co23 - co24 * co12 * co33
      d3 = co11 * co24 * co33 + co21 * co34 * co13 + co31 * co14 * co23 - &
           co31 * co24 * co13 - co11 * co34 * co23 - co21 * co14 * co33
      d4 = co11 * co22 * co34 + co21 * co32 * co14 + co31 * co12 * co24 - &
           co31 * co22 * co14 - co11 * co32 * co24 - co21 * co12 * co34
      x1 = d2 / d1
      x2 = d3 / d1
      x3 = d4 / d1
      vf = 0
    end if
  end subroutine GLS

  subroutine PXS ! Points of intersection ray/crystal surface 
  integer(8) sm, sl
  real(8) lqh
    co11 = snk11
    co12 = snk12
    co13 = snk13
    co14 = 1
    co21 = snk21
    co22 = snk22
    co23 = snk23
    co24 = 1
    co31 = nk1(u)
    co32 = nk2(u)
    co33 = nk3(u)
    co34 = 1
    call GLS
    ac = x1 - xold
    bc = x2 - yold
    cc = x3 - zold
    direction = ev31*ac + ev32*bc + ev33*cc 
    if (direction .gt. 0 .and. vf .eq. 0 .and. u .ne. utrold) then
!intersection 'in front of' ray & intersection exists & not old plane 
      xs1(u) = x1 ! sort with respect to distance 
      xs2(u) = x2
      xs3(u) = x3
      lqh = ac * ac + bc * bc + cc * cc
      h = h + 1
      lq(h) = lqh
      umtr(h) = u
      do sm = 1, h - 1
        if (lq(h) .lt. lq(sm)) then
          do sl = h, (sm + 1), -1
            lq(sl) = lq(sl - 1)
            umtr(sl) = umtr(sl - 1) ! order of hitted surfaces 
          end do
          lq(sm) = lqh
          umtr(sm) = u
        end if
      end do
    end if  
  end subroutine PXS

  subroutine PWS
! determines wether xs is within the boundary of the surface 
  integer(8)  v, forlim, um1
  real(8) scalpr
  logical inside

    u = umtr(m)
    um1 = u
    inside = .true.
    vmu = vm(um1)
    v = 0
    do while (inside .and. v .ne. vmu)
      v = v + 1
      scalpr = kano1(um1,v)*(xs1(um1) - p1(um1,v)) + &
               kano2(um1,v)*(xs2(um1) - p2(um1,v)) + &
               kano3(um1,v)*(xs3(um1) - p3(um1,v))
      if (scalpr .gt. 0.0) inside = .false.
    end do
    if (inside) utr = u
    m = m + 1
  end subroutine pws

  subroutine PSG            ! Ray equations 
! a ray is determined by the intersection of two plane equations, 
!  where plane I is defined by direction of propagation & 1st electric 
!  field vector and plane  II by            "           & 2nd electric 
!  field vector
  
    co11 = xold
    co12 = yold
    co13 = zold
    co14 = 1
    co21 = xold + ev21
    co22 = yold + ev22
    co23 = zold + ev23
    co24 = 1
    co31 = xold + ev31
    co32 = yold + ev32
    co33 = zold + ev33
    co34 = 1
    call GLS
    snk11 = x1
    snk12 = x2
    snk13 = x3
    co11 = xold + ev11
    co12 = yold + ev12
    co13 = zold + ev13
    co14 = 1
    co21 = xold
    co22 = yold
    co23 = zold
    co24 = 1
    co31 = xold + ev31
    co32 = yold + ev32
    co33 = zold + ev33
    co34 = 1
    call GLS
    snk21 = x1
    snk22 = x2
    snk23 = x3
  end subroutine PSG

  subroutine PTXA
! hitted surface 'utr' and point of intersection 'xs(utr)' 
  integer(8) su
    call PSG            ! ray equation 
    do su = 1, um 
      lq(su) = 0  ! initialization 
      umtr(su) = 0
    end do
    h = 0
    do su = 1, um ! point of intersection with unbounded surface 
      u = su
      if (u .ne. utrold) then 
        call PXS ! -> lq(u), xs(u), u <> utrold 
      end if
    end do
    m = 1
    utr = 0
    do while (m .le. h .and. utr .ne. u .and. lq(m) .lt. lq_max) 
      call PWS ! determ. wether xs(m) lies within the bounded surface 
    end do
  end subroutine PTXA

  subroutine PTXB ! (almost) the same as PTXA but for internal rays 
  integer(8) su
    call PSG
    h = 0
    do su = 1, um
       umtr(su) = 0
    end do
    do su = 1, um
      u = su
      call PXS             ! --> lq(:),umtr(:)   
    end do
    m = 1
    utr = 0
    do while (m .le. h .and. utr .ne. u) 
      call PWS
    end do
    wegl = sqrt(lq(m-1))
    help3 = mabs*wegl
    if (help3 .lt. - 50) then ! ###
      absfak = 0.0
    else
      absfak = exp(help3)     ! Absorption 
    end if


  end subroutine PTXB

  subroutine cross_product(a1, a2, a3, b1, b2, b3)
  real(8) a1, a2, a3, b1, b2, b3
    cross_1 = a2*b3 - a3*b2
    cross_2 = a3*b1 - a1*b3
    cross_3 = a1*b2 - a2*b1
    nf = sqrt(cross_1*cross_1 + cross_2*cross_2 + cross_3*cross_3)
    cross_1 = cross_1/nf
    cross_2 = cross_2/nf
    cross_3 = cross_3/nf
  end subroutine cross_product

  subroutine rotation_matrix(en1, en2, en3, eo1, eo2, eo3, ko1, ko2, ko3)
  real(8) eo1, eo2, eo3 ! rotation_matrix: perp. to reference plane
  real(8) en1, en2, en3 !        "       : perp. e-field 
  real(8) ko1, ko2, ko3 !        '       : direction

    myfunc_arg = en1*eo1 + en2*eo2 + en3*eo3 ! perp_new * perp_old
    if (abs(myfunc_arg) .gt. 1.0) myfunc_arg = sign(1.0d0, myfunc_arg)
    ev21 = ko2*en3 - ko3*en2 
    ev22 = ko3*en1 - ko1*en3
    ev23 = ko1*en2 - ko2*en1
    test1 = myfunc_arg*(eo1*ev21 + eo2*ev22 + eo3*ev23) ! perp_new * par_old 
    help1 = myfunc_arg*myfunc_arg
    cos_rot = 2.*help1 - 1.
    help2 = 4.*help1*(1. - help1) ! sin^2(2*psi_rot) 
    help3 = sqrt(help2)
    sin_rot = sign(help3, test1)
    rot_stokes(2,2) = cos_rot
    rot_stokes(2,3) =-sin_rot
    rot_stokes(3,2) = sin_rot
    rot_stokes(3,3) = cos_rot
  end subroutine rotation_matrix

  subroutine tilt
  integer counter
  real(8) kh1, kh2, kh3
    counter = 0 ! random distortion of the rays propagation vector
    w = 1.0 
    do while (w .ge. 0 .and. counter .lt. 10)  ! equivalent to displacemet
      counter = counter + 1  ! of crystal coordinates
      call random_number(rand)
      kh1 = ev31 + distortion*(2*rand-1)
      call random_number(rand)
      kh2 = ev32 + distortion*(2*rand-1)
      call random_number(rand)
      kh3 = ev33 + distortion*(2*rand-1)
      test1 = sqrt(kh1*kh1 + kh2*kh2 + kh3*kh3)
      kh1 = kh1/test1 
      kh2 = kh2/test1 
      kh3 = kh3/test1
      w = kh1*fn1(utr) + kh2*fn2(utr) + kh3*fn3(utr) !-dcos(theta_i
    end do

    if (counter .lt. 10) then
      ev31 = kh1 
      ev32 = kh2 
      ev33 = kh3
    else ! w = negative cosine of angle of incidence 
      w = ev31 * fn1(utr) + ev32 * fn2(utr) + ev33 * fn3(utr)
    end if

  end subroutine tilt

  subroutine reflection_vectors
    help1 = -2.*w ! Attention: if w not < 0 in PRTA then @#$
    vk91 = ev31 + help1*fn1(utr) ! direction vector for reflection 
    vk92 = ev32 + help1*fn2(utr)
    vk93 = ev33 + help1*fn3(utr)   
    nf = SQRT(vk91*vk91 + vk92*vk92 + vk93*vk93)
    vk91 = vk91 / nf 
    vk92 = vk92 / nf 
    vk93 = vk93 / nf

    vk71 = fn2(utr)*ev33 - fn3(utr)*ev32 ! vector perp. to plane of incidence
    vk72 = fn3(utr)*ev31 - fn1(utr)*ev33 ! perp = f x prop_in
    vk73 = fn1(utr)*ev32 - fn2(utr)*ev31
    nf = SQRT(vk71*vk71 + vk72*vk72 + vk73*vk73)
    help4 = -abs(w)/w ! sim. f pointing in(out)ward for ex(in)ternal reflection
    vk71 = help4*vk71/nf 
    vk72 = help4*vk72/nf 
    vk73 = help4*vk73/nf
  end subroutine reflection_vectors

  subroutine refraction_vectors
    wvz = w / abs(w)
      !  ---- direction of refraction ----  
    vk121 = wf * (ev31 -w * fn1(utr)) + wvz * help2 * fn1(utr)
    vk122 = wf * (ev32 -w * fn2(utr)) + wvz * help2 * fn2(utr)
    vk123 = wf * (ev33 -w * fn3(utr)) + wvz * help2 * fn3(utr)
    nf = SQRT(vk121 * vk121 + vk122 * vk122 + vk123 * vk123)
    vk121 = vk121 / nf 
    vk122 = vk122 / nf 
    vk123 = vk123 / nf
  end subroutine refraction_vectors

  subroutine random_rotation
    test1 = 2
    do while (test1 .gt. 1) ! effective search for sin/cos(phi)...
      call random_number(rand)
      w1 = 1. - 2.*rand
      call random_number(rand)
      w2 = 1. - 2.*rand ! ... if phi is random
      test1 = w1*w1 + w2*w2 
    end do
    help1 = sqrt(test1)
    cos_rot = w1/help1
    sin_rot = w2/help1
  end subroutine random_rotation

  subroutine PRTA ! Ray Tracing for external rays 
  integer(8) i, j, k

    call tilt ! oder w = ev31*fn1(utr) + ev32*fn2(utr) + ev33*fn3(utr) < 0!!
    if (w .lt. -1.) w = -1. ! numerical reasons
    w = -abs(w) ! numerical reasons
    angle_i = int(acos((abs(w)))*rayconst1) ! angel bin, abs(w) instead -w for numerical reasons.  
    call reflection_vectors ! prop and e_perp for reflection: vk9, vk7
    help2 = cos_theta_t_a(angle_i) ! -> refraction_vector
    wf = wfa(angle_i)
    call refraction_vectors ! prop for refraction: vk12, e_perp = vk7
    call rotation_matrix(vk71, vk72, vk73, ev11, ev12, ev13, ev31, ev32, ev33)
    new_in_stokes = matmul(rot_stokes,in_stokes)

    do i = 1, 4
      do j = 1, 4
        trans_r(i,j) = ra(i,j,angle_i)
        trans_t(i,j) = ta(i,j,angle_i)
      end do
    end do

    refl_stokes = matmul(trans_r, new_in_stokes)
    trans_stokes = matmul(trans_t, new_in_stokes)

    test1 = in_stokes(1,1) - refl_stokes(1,1) - trans_stokes(1,1) ! ### 
    energy_violation = energy_violation + abs(test1) ! ### 

    rayhaben = rayhaben + 1 ! add reflected ray to ray tree 
    raytree(rayhaben)%in_stokes = refl_stokes
    raytree(rayhaben)%es1 = vk71
    raytree(rayhaben)%es2 = vk72
    raytree(rayhaben)%es3 = vk73
call cross_product(vk91, vk92, vk93, vk71, vk72, vk73)
    raytree(rayhaben)%ep1 = cross_1
    raytree(rayhaben)%ep2 = cross_2
    raytree(rayhaben)%ep3 = cross_3
    raytree(rayhaben)%ek1 = vk91
    raytree(rayhaben)%ek2 = vk92
    raytree(rayhaben)%ek3 = vk93
    raytree(rayhaben)%utrold = utr
    raytree(rayhaben)%xold = xs1(utr)
    raytree(rayhaben)%yold = xs2(utr)
    raytree(rayhaben)%zold = xs3(utr)
    raytree(rayhaben)%aorb = 'A'
    raytree(rayhaben)%recdeep = raytree(raysoll)%recdeep + 1
    raytree(rayhaben)%totrefl = raytree(raysoll)%totrefl

    rayhaben = rayhaben + 1 ! add refracted ray to ray-tree 
    raytree(rayhaben)%in_stokes = trans_stokes
    raytree(rayhaben)%es1 = vk71
    raytree(rayhaben)%es2 = vk72
    raytree(rayhaben)%es3 = vk73
call cross_product(vk121, vk122, vk123, vk71, vk72, vk73)
    raytree(rayhaben)%ep1 = cross_1
    raytree(rayhaben)%ep2 = cross_2
    raytree(rayhaben)%ep3 = cross_3
    raytree(rayhaben)%ek1 = vk121
    raytree(rayhaben)%ek2 = vk122
    raytree(rayhaben)%ek3 = vk123
    raytree(rayhaben)%utrold = utr
    raytree(rayhaben)%xold = xs1(utr)
    raytree(rayhaben)%yold = xs2(utr)
    raytree(rayhaben)%zold = xs3(utr)
    raytree(rayhaben)%aorb = 'B'
    raytree(rayhaben)%recdeep = raytree(raysoll)%recdeep + 1
    raytree(rayhaben)%totrefl = raytree(raysoll)%totrefl
  end subroutine PRTA
  
  subroutine PRTB ! Ray Tracing for internal rays 
  integer(8) i, j, k

    call tilt ! oder w = ev31*fn1(utr) + ev32*fn2(utr) + ev33*fn3(utr) > 0!!
    if (w .gt. 1.) w = 1. ! numerical reasons
    w = abs(w)    
    angle_i = int(dacos(abs(w))*rayconst1)!  angel bin
    call reflection_vectors ! prop and e_perp for reflection: vk9, vk7
    if (angle_i .lt. angle_g) then ! if no total reflection then 
      help2 = cos_theta_t_b(angle_i) ! -> refraction_vector
      wf = wfb(angle_i) ! -> refraction_vector
      call refraction_vectors ! prop for refraction: vk12, e_perp = vk7
      call rotation_matrix(vk71, vk72, vk73, ev11, ev12, ev13, ev31, ev32, ev33)
      new_in_stokes = matmul(rot_stokes, in_stokes)
      do i = 1, 4
        do j = 1, 4
          trans_r(i,j) = rb(i,j,angle_i)
          trans_t(i,j) = tb(i,j,angle_i)
        end do
      end do
      refl_stokes = matmul(trans_r, new_in_stokes)
      trans_stokes = matmul(trans_t, new_in_stokes)

      test1 = in_stokes(1,1)  - refl_stokes(1,1) - trans_stokes(1,1)! ### 
      energy_violation = energy_violation + abs(test1) ! ### 
      energy_abs = energy_abs + (1 - absfak)*in_stokes(1,1)

      rayhaben = rayhaben + 1 ! add refracted ray to raytree 
      raytree(rayhaben)%in_stokes = absfak*trans_stokes
      raytree(rayhaben)%es1 = vk71
      raytree(rayhaben)%es2 = vk72
      raytree(rayhaben)%es3 = vk73
call cross_product(vk121, vk122, vk123, vk71, vk72, vk73)
      raytree(rayhaben)%ep1 = cross_1
      raytree(rayhaben)%ep2 = cross_2
      raytree(rayhaben)%ep3 = cross_3
      raytree(rayhaben)%ek1 = vk121
      raytree(rayhaben)%ek2 = vk122
      raytree(rayhaben)%ek3 = vk123
      raytree(rayhaben)%utrold = utr
      raytree(rayhaben)%xold = xs1(utr)
      raytree(rayhaben)%yold = xs2(utr)
      raytree(rayhaben)%zold = xs3(utr)
      raytree(rayhaben)%aorb = 'A'
      raytree(rayhaben)%recdeep = raytree(raysoll)%recdeep + 1
      raytree(rayhaben)%totrefl = raytree(raysoll)%totrefl

      rayhaben = rayhaben + 1 ! add reflected ray to raytree 
      raytree(rayhaben)%in_stokes = absfak*refl_stokes
      raytree(rayhaben)%es1 = vk71
      raytree(rayhaben)%es2 = vk72
      raytree(rayhaben)%es3 = vk73
call cross_product(vk91, vk92, vk93, vk71, vk72, vk73)
      raytree(rayhaben)%ep1 = cross_1
      raytree(rayhaben)%ep2 = cross_2
      raytree(rayhaben)%ep3 = cross_3
      raytree(rayhaben)%ek1 = vk91
      raytree(rayhaben)%ek2 = vk92
      raytree(rayhaben)%ek3 = vk93
      raytree(rayhaben)%utrold = utr
      raytree(rayhaben)%xold = xs1(utr)
      raytree(rayhaben)%yold = xs2(utr)
      raytree(rayhaben)%zold = xs3(utr)
      raytree(rayhaben)%aorb = 'B'
      raytree(rayhaben)%recdeep = raytree(raysoll)%recdeep + 1
      raytree(rayhaben)%totrefl = raytree(raysoll)%totrefl  

    else ! Totalreflection 

      call rotation_matrix(vk71, vk72, vk73, ev11, ev12, ev13, ev31, ev32, ev33)
      new_in_stokes = matmul(rot_stokes, in_stokes)
      do i = 1, 4
        do j = 1, 4
          trans_r(i,j) = rb(i,j,angle_i)
        end do
      end do
      refl_stokes = matmul(trans_r, new_in_stokes)

      test1 = in_stokes(1,1)  - refl_stokes(1,1) ! ### 
      energy_violation = energy_violation + abs(test1) ! ### 
      energy_abs = energy_abs + (1 - absfak)*in_stokes(1,1)

      rayhaben = rayhaben + 1
      raytree(rayhaben)%in_stokes = absfak*refl_stokes
      raytree(rayhaben)%es1 = vk71 ! add 'total' reflected ray to raytree 
      raytree(rayhaben)%es2 = vk72
      raytree(rayhaben)%es3 = vk73
call cross_product(vk91, vk92, vk93, vk71, vk72, vk73)
      raytree(rayhaben)%ep1 = cross_1
      raytree(rayhaben)%ep2 = cross_2
      raytree(rayhaben)%ep3 = cross_3
      raytree(rayhaben)%ek1 = vk91
      raytree(rayhaben)%ek2 = vk92
      raytree(rayhaben)%ek3 = vk93
      raytree(rayhaben)%utrold = utr
      raytree(rayhaben)%xold = xs1(utr)
      raytree(rayhaben)%yold = xs2(utr)
      raytree(rayhaben)%zold = xs3(utr)
      raytree(rayhaben)%aorb = 'B'
      raytree(rayhaben)%recdeep = raytree(raysoll)%recdeep
! no count-up of ray-depth for total internal reflection !!! 
      raytree(rayhaben)%totrefl = raytree(raysoll)%totrefl + 1

    end if ! total reflection

  end subroutine PRTB

  subroutine PZS            ! store intermediate status 
  real(8) ext_cs, ext_eff, abs_cs, abs_eff, scat_cs, scat_eff, r_eq
  integer(8) i, j, k
  character(8) date
  character(10) time
  character(5) zone
  integer values(8)

    call date_and_time(date, time, zone, values)
    call system_clock(end_time)

    open(unit = 12, file = ifn(1:pure_str)//'.save', form = 'unformatted')
    do i = 0, 180
      do j  = 1, nstokes 
        do k  = 1, 4
           write(12) rtd(i,j,k)
        end do
      end do
    end do
    close(12)

    ext_eff = 2. ! extinction efficiency = 2, what else ...

    if (no_diff .eq. 0) then 
      gcs = gcs_mean_poly
    else 
      gcs = gcs_mean 
    end if
    r_eq = sqrt(gcs/pi) ! radius of surface eq. sphere 
    ext_cs = ext_eff*gcs ! extinction cross section 
    scat_eff = 1 + (energy_in - energy_abs)/energy_in ! ### -> S. 226
    scat_cs  = gcs*scat_eff ! scattering cross section
    abs_eff = ext_eff - scat_eff ! absorption efficiency 
    abs_cs = abs_eff*gcs ! absorption cross section 
    
    w_0 = scat_eff/ext_eff ! single scattering albedo 
    nt = 0
    open(12, file = ifn(1:pure_str)//'.status', form = 'formatted')
    write(12,*) ' Inputfile   : ', ifn
    write(12,*) ' Date & time : ', date
    write(12,*) ' Elapsed time: ', int((end_time - start_time)/1.e6), ' sec.' ! time
    write(12,*) ' Random seed : ', random_arg
    write(12,*) ' '
    write(12,*) orient-1, ' out of ', orientmax, ' orientations'
    write(12,*) poly_orient,' good orientations'
    write(12,*) ntotal, ' hits out of ', ntm_sum, ' incoming rays'
    write(12,*) ' '
    write(12,*) 'Incoming energy  : ', energy_in
    write(12,*) 'Scattered energy : ', energy_out
    write(12,*) 'Real abs. energy : ', energy_abs
help1 = energy_in - energy_out - energy_abs
    write(12,*) 'Missing energy   : ',help1,' = ',help1/energy_in*100, '%'
    write(12,*) 'FreSnell-Error   : ', energy_violation,' = ',energy_violation/energy_in*100, '%'
    write(12,*) 'Finite tot. refl.: ', lost_tr,' = ',lost_tr/energy_in*100, '%'
    write(12,*) 'Finite ray recur.: ', lost_rm,' = ',lost_rm/energy_in*100, '%'
help1 = help1 - lost_tr - lost_rm
    write(12,*) 'Still missing    : ', help1,' = ',help1/energy_in*100, '%'
    write(12,*) ' '
    write(12,*) '''scanning'' cross section: ', gcs_mean
    write(12,*) '''polygon '' cross section: ', gcs_mean_poly
    write(12,*) '              extinction, scattering, absorption : '
    write(12,*) 'cross sec.: ',ext_cs,' ',scat_cs,' ',abs_cs
    write(12,*) 'effic.''s : ',ext_eff,' ',scat_eff,' ',abs_eff
    write(12,*) ' '
    write(12,*) 'w_0    : ', w_0
    write(12,*) 'g_rt   : ', g_rt
    write(12,*) 'g_diff (g_diff_discrete) : ', g_diff, g_diff_discrete
    write(12,*) 'g_total: ', g
    write(12,*) ' '
    write(12,*) 'delta forward peak : ', e_0
    write(12,*) '      backward "   : ', e_pi
    close(12)

! save it again to small binary file for easy analysis
    open(12, file = ifn(1:pure_str)//'.status.bin', form = 'unformatted')
    write(12) energy_in
    write(12) energy_out
    write(12) energy_abs
    write(12) energy_violation
    write(12) lost_tr
    write(12) lost_rm
    write(12) gcs_mean
    write(12) gcs_mean_poly
    write(12) ext_eff, scat_eff, abs_eff
    write(12) w_0
    write(12) g_rt
    write(12) g_diff
    write(12) g
    close(12)

  end subroutine PZS

  subroutine PRTNORM ! final normalization and output 
  
  real(8) ptotal, rt11, diff, inc_fac, x1, x2, y1, y2, x, y
  integer(8) i, j, k

    help1 = rad/resol_fak ! delta theta 
    help2 = rad*0.25d0/resol_fak
    help3 = rad*(180.d0 - 0.25d0/resol_fak) 
    help5 = 2.d0*w_0 - 1.d0
    help6 = 2.d0*w_0
    help7 = help5/help6

    if (no_diff .eq. 0) then 

      do i = 1, theta_diff
        diff_total(theta_diff + i) = diff_total(theta_diff - i)
      end do
! The following do-loop is based on guesses and requires further proof
      do i = 0, rtddim ! multiply "inclination factor" 
        angle = dble(i)*rad2resol_fak
        inc_fac = (0.5*(1. + cos(angle)))**2
        diff_total(i) = diff_total(i)*inc_fac
      end do
      diff_total(rtddim) = diff_total(rtddim - 1) ! avoid "0" at 180 degree

      call diff_discrete ! normalize and write "discrete" diffraction phase function

! Interpolate diffraction at center angle of bin 1 ...
      y1 = diff_total(0)
      y2 = diff_total(1)
      x1 = 0.
      x2 = 1./resol_fak
      x = help2*deg
      y = (y1 - y2)/(x1 - x2)*(x - x1) + y1
      diff_total(0) = y ! ... and replace exact forward diff. with this
! Interpolate diffraction at center angle of bin rtddim ...
      y1 = diff_total(rtddim - 1)
      y2 = diff_total(rtddim)
      x1 = (rtddim - 1)*help1*deg
      x2 = 180
      x = help3*deg
      y = (y1 - y2)/(x1 - x2)*(x - x1) + y1
      diff_total(rtddim) = y ! ... and replace exact backward diff. with this
! all additional diffraction values can be regarded as bin means. 

      sum1 = 0 ! normalize diffraction ...
      do i = 0, rtddim 
        d_o1 = 2*sin(i*help1)*sin(help1/2) ! solid angle interval 
        if (i .eq. 0) d_o1 = 2*sin(help2)*sin(help2)
        if (i .eq. rtddim) d_o1 = 2*sin(help3)*sin(help2)
        sum1 = sum1 + diff_total(i)*d_o1 ! dE/dOmega * dOmega = dE
      end do
      do i = 0, rtddim
        diff_total(i) = diff_total(i)/sum1*2.
      end do

      open(63, file = ifn(1:pure_str)//'.diff')  ! Diffraction pattern 
      g_diff = 0.
      do i = 0, rtddim 
        angle = i*help1*deg
        d_o1 = 2*sin(i*help1)*sin(help1/2) ! solid angle interval 
        if (i .eq. 0) then 
          angle = help2*deg
          d_o1 = 2*sin(help2)*sin(help2)
        end if
        if (i .eq. rtddim) then
          angle = help3*deg
          d_o1 = 2*sin(help3)*sin(help2)
        end if
        write(63, 1004) angle, diff_total(i)
        help4 = cos(angle*rad)*d_o1
        g_diff = g_diff + diff_total(i)*help4
      end do
      g_diff = g_diff/2.
      close(63)
    end if ! (no_diff = 0) if

    sum2 = 0. ! ... and ray tracing
    do i = 0, rtddim 
      sum2 = sum2 + rtd(i,1,1) ! == dE
    end do

    open(61, file = ifn(1:pure_str)//'.1dr') ! ray tracing phase function: rt11 
    do i = 0, rtddim ! divide dE by solid angle interval -> phase function
      d_o1 = 2*sin(i*help1)*sin(help1/2) ! solid angle interval 
      angle = i*help1*deg
      if (i .eq. 0) then 
        d_o1 = 2*sin(help2)*sin(help2)
        angle = help2*deg
      end if
      if (i .eq. rtddim) then
        d_o1 = 2*sin(help3)*sin(help2)
        angle = help3*deg
      end if
      help4 = cos(angle*rad)*d_o1

      do j = 1, 4
        do k = 1, 4
          rtd(i,j,k) = rtd(i,j,k)/d_o1/sum2*2.
        end do
      end do
      write(61, 1004) angle, rtd(i,1,1)
      g_rt = g_rt + rtd(i,1,1)*help4
    end do
    close(61)
    g_rt = g_rt/2. 
    e_0 = e_0/sum2 ! RT energy fraction for delta forward ...
    e_pi = e_pi/sum2 ! ... and delta backward peak


 1003 format(f8.4,1x,f7.4)
 1004 format(f8.4,1x,e9.3)

    if (diff_seperate .eq. 0 .and. no_diff .eq. 0) then 

      open(12, file = ifn(1:pure_str)//'.p12')  
!     open(21, file = ifn(1:pure_str)//'.p21') ! only for testing purposes
      open(22, file = ifn(1:pure_str)//'.p22')
      open(33, file = ifn(1:pure_str)//'.p33')
      open(34, file = ifn(1:pure_str)//'.p34')
!     open(43, file = ifn(1:pure_str)//'.p43') ! only for testing purposes
      open(44, file = ifn(1:pure_str)//'.p44')
!     open(25, file = ifn(1:pure_str)//'.pij.bin') ! only for testing purposes
      open(64, file = ifn(1:pure_str)//'.p11') ! total phase function 

      g = 0. ! asymmetry parameter 
      do i = 0, rtddim 
        angle = i*help1*deg
        d_o1 = 2*sin(i*help1)*sin(help1/2) ! solid angle interval 
        if (i .eq. 0) then 
          angle = help2*deg
          d_o1 = 2*sin(help2)*sin(help2)
        end if
        if (i .eq. rtddim) then
          angle = help3*deg
          d_o1 = 2*sin(help3)*sin(help2)
        end if

        rt11 = rtd(i,1,1)
        diff = diff_total(i)
        ptotal = (diff + help5*rt11)/help6

        do j = 1, 4 ! add diffraction to diagonal elements
          do k = 1, 4
            if (j .eq. k) then 
              pij(i,j,k) = (diff + help5*rtd(i,j,k))/help6
            else 
              pij(i,j,k) = rtd(i,j,k)*help7
            end if
          end do
        end do
        write(64, 1004) angle, ptotal


!       write(25) angle, ((pij(i,j,k)/pij(i,1,1), j = 1,4), k = 1,4)
        write(12, 1003) angle, pij(i,1,2)/pij(i,1,1)
!       write(21, 1003) angle, pij(i,2,1)/pij(i,1,1)
        write(22, 1003) angle, pij(i,2,2)/pij(i,1,1)
        write(33, 1003) angle, pij(i,3,3)/pij(i,1,1)
        write(34, 1003) angle, pij(i,3,4)/pij(i,1,1)
!       write(43, 1003) angle, pij(i,4,3)/pij(i,1,1)
        write(44, 1003) angle, pij(i,4,4)/pij(i,1,1)

        help4 = cos(angle*rad)*d_o1
        g = g + ptotal*help4
      end do
!     close(25)
      close(12)
!     close(21)
      close(22)
      close(33)
      close(34)
!     close(43)
      close(44)
      close(64)
      g = g/2.0

    else ! pure ray tracing

      open(12, file = ifn(1:pure_str)//'.r12')  
!     open(21, file = ifn(1:pure_str)//'.r21') ! only for testing purposes
      open(22, file = ifn(1:pure_str)//'.r22')
      open(33, file = ifn(1:pure_str)//'.r33')
      open(34, file = ifn(1:pure_str)//'.r34')
!     open(43, file = ifn(1:pure_str)//'.r43') ! only for testing purposes
      open(44, file = ifn(1:pure_str)//'.r44')
!     open(25, file = ifn(1:pure_str)//'.rij.bin') ! only for testing purposes
      do i = 0, rtddim
        angle = i*help1*deg
        if (i .eq. 0) angle = help2*deg
        if (i .eq. rtddim) angle = help3*deg
!       write(25) angle, ((rtd(i,j,k)/rtd(i,1,1), j = 1,4), k = 1,4)
        write(12, 1003) angle, rtd(i,1,2)/rtd(i,1,1)
!       write(21, 1003) angle, rtd(i,2,1)/rtd(i,1,1)
        write(22, 1003) angle, rtd(i,2,2)/rtd(i,1,1)
        write(33, 1003) angle, rtd(i,3,3)/rtd(i,1,1)
        write(34, 1003) angle, rtd(i,3,4)/rtd(i,1,1)
!       write(43, 1003) angle, rtd(i,4,3)/rtd(i,1,1)
        write(44, 1003) angle, rtd(i,4,4)/rtd(i,1,1)      
      end do

!     close(25)
      close(12)
!     close(21)
      close(22)
      close(33)
      close(34)
!     close(43)
      close(44)

    end if

  end subroutine PRTNORM  

  subroutine diff_discrete
  integer(8) i, j
  real(8) diff_total_discrete(0:rtddim_max)    ! discrete diffraction function 

! The following would be a good integration of the diffraction function. 
! Unfortunately, it does not fit to the 'binned' ray tracing results.
! We do it anyway and save it under inf + '.diff.discrete' 
! The corresponding asymmetry parameter -> g_diff_discrete

    help1 = rad/resol_fak ! delta theta 
    help2 = rad*0.25/resol_fak
    help3 = rad*(180 - 0.25/resol_fak) 

    simpson = 2.*pi/dble(rtddim)/6.
    simpson_fak(0) = 1.
    simpson_fak(1) = 4.
    simpson_fak(2) = 1.

    sum = 0.
    do i = 0, rtddim - 2, 2! test diffraction phase function normalization
      do j = 0, 2
        iangle = i + j
        angle = dble(iangle)*rad2resol_fak
        sintheta = sin(angle)
        sum = sum + diff_total(iangle)*sintheta*simpson_fak(j)
      end do
    end do
    sum = sum*simpson ! seems to be = 2pi !!!

    do i = 0, rtddim ! renormalize diffraction phase function
      diff_total_discrete(i) = diff_total(i)/sum*2.
    end do

    g_diff_discrete = 0
    do i = 0, rtddim - 2, 2 ! calculate asymmetry parameter
      do j = 0, 2
        iangle = i + j
        angle = dble(iangle)*rad2resol_fak
        sintheta = sin(angle)
        costheta = cos(angle)
        g_diff_discrete = g_diff_discrete + diff_total_discrete(iangle)*costheta*sintheta*simpson_fak(j)
      end do
    end do
    g_diff_discrete = g_diff_discrete*simpson/2.

    open(10, file = ifn(1:pure_str)//'.diff.discrete')
    do i = 0, rtddim 
      angle = i*help1*deg
      write(10,1004) angle, diff_total_discrete(i)
    end do
    close(10)
 1004 format(f8.4,1x,e9.3)

  end subroutine diff_discrete
  

  subroutine PDAL
! reads crystal coordinates from file 'afn'. The structure of afn is
! as follows:
! um               : number of plane surfaces
! vm(1)            : number of edge points in surface 1
!  .
!  .
! vm(um)           :               "                  um
! p1(1,1),p2(1,1),p3(1,1)      : x,y,z coordinates of point 1 in plane 1
!          .
!          .
! p1(1,vm(1)),p2(1,vm(1)),p3(1,vm(1)) :        "            vm(1)  "   1
! p1(2,1),p2(2,1),p3(2,1)      :               "            1          2
!          .
!          .
! p1(2,vm(1)),p2(2,vm(1)),p3(2,vm(1)) :        "            vm(1)  "   2
!          
! and so on :.
!
! NOTE: The crystal edge points must be given in clockwise order
!       when looking from outside !!!

  
  integer(8) i, j
    read(10, fmt = '(A100)') afn
    open(unit = 11, file = afn, status = 'old', form = 'formatted')
    read(11, *) um ! number of plane surfaces 
    do j  = 1, um
      read(11, *) vm(j) ! number of points per surface 
    end do
    do i = 1, um
      do j  = 1, vm(i)
        read(11, *) p1(i, j), p2(i, j), p3(i, j) ! x,y,z coordinates
       end do
    end do 
    do i = 1, um
      do j = vm(i) + 1, vn
        p1(i, j) = 0
        p2(i, j) = 0
        p3(i, j) = 0
      end do
    end do
    do i = 1, um
      do j  = 1, vm(i)
        p4(i, j) = p1(i, j)
        p5(i, j) = p2(i, j)
        p6(i, j) = p3(i, j)
        p7(i, j) = p1(i, j)
        p8(i, j) = p2(i, j)
        p9(i, j) = p3(i, j)
      end do
    end do
    close(11)
  end subroutine PDAL

  subroutine PDAS
! writes crystal coordinates to file 'rot.crystal'. 
!  For testing purposes only, usually disabled
  
!  integer(8) i, j
    open(11, file = 'test.crystal')
    write(11,*) um
    do j_sub  = 1, um
      write(11, *) vm(j_sub)
    end do
    do i_sub = 1, um
      do j_sub = 1, vm(i_sub)
        write(11,*) p1(i_sub, j_sub), p2(i_sub, j_sub), p3(i_sub, j_sub)
      end do
    end do
    close(11)
  end subroutine PDAS
  
  subroutine SDATIN ! Input of ray tracing parameter 
  
    read(10,*) la, rbi, ibi ! wavelength, real- and im. part of refra. 
    read(10,*) rays ! Number of incident rays per crystal orientation 
    read(10,*) orientmax ! Number of crystal orientations 
    read(10,*) recdeepMax ! Number of ray recursions
    read(10,*) totreflmax ! # of total internal reflections 
    read(10,*) distortion ! degree of crystal distortion, 0 < distortion < 1 
    read(10,*) rtddim
    read(10,*) diff_seperate
    read(10,*) no_diff
    close(10) 
  end subroutine SDATIN

  subroutine PEX ! determines extreme coordinates 
  
  integer(8) i, j
    xmax = p1(1, 1)
    ymax = p2(1, 1)
    zmax = p3(1, 1)
    xmin = p1(1, 1)
    ymin = p2(1, 1)
    zmin = p3(1, 1)
    do i = 1, um
      do j  = 1, vm(i)
        if (xmax .lt. p1(i, j)) xmax = p1(i, j) ! use maxval, minval 
        if (ymax .lt. p2(i, j)) ymax = p2(i, j)
        if (zmax .lt. p3(i, j)) zmax = p3(i, j)
        if (xmin .gt. p1(i, j)) xmin = p1(i, j)
        if (ymin .gt. p2(i, j)) ymin = p2(i, j)
        if (zmin .gt. p3(i, j)) zmin = p3(i, j)
      end do
    end do
  end subroutine PEX

  subroutine PSP ! center of gravity
   
  integer(8) n, su, sv
    sp1 = 0
    sp2 = 0
    sp3 = 0
    n = 0
    do su  = 1, um 
      do sv  = 1, vm(su) 
        n = n + 1
        sp1 = p1(su, sv) + sp1
        sp2 = p2(su, sv) + sp2
        sp3 = p3(su, sv) + sp3
      end do
    end do
    sp1 = sp1 / n
    sp2 = sp2 / n
    sp3 = sp3 / n
  end subroutine PSP
 
  subroutine PNK ! coefficients of plane equation   
!  nk1(i)*x + nk2(i)*y + nk3(i)*z = 1, i = 1,..,um
  
    co11 = p1(mainu, 1)
    co12 = p2(mainu, 1)
    co13 = p3(mainu, 1)
    co14 = 1
    co21 = p1(mainu, 2)
    co22 = p2(mainu, 2)
    co23 = p3(mainu, 2)
    co24 = 1
    co31 = p1(mainu, 3)
    co32 = p2(mainu, 3)
    co33 = p3(mainu, 3)
    co34 = 1
    call GLS
    nk1(mainu) = x1
    nk2(mainu) = x2
    nk3(mainu) = x3
  end subroutine PNK

  subroutine PPD ! normal vectors with respect to surface segments 
  
  integer(8) v, forlim, um1, endp
  real(8) diff_x, diff_y, diff_z
  logical inside
    um1 = mainu
    inside = .true.
    vmu = vm(um1)
    forlim = vmu
    do v  = 1, forlim 
      endp = v + 1
      if (endp .eq. vmu + 1) endp = 1
      diff_x = p1(um1,endp)-p1(um1,v) ! vector in the polygon plane 
      diff_y = p2(um1,endp)-p2(um1,v) ! joining p_n -> p_(n+1) 
      diff_z = p3(um1,endp)-p3(um1,v)
      kano1(um1,v) = (fn2(um1)*diff_z-fn3(um1)*diff_y) ! vector in the plane 
      kano2(um1,v) = (fn3(um1)*diff_x-fn1(um1)*diff_z) ! perp. to diff 
      kano3(um1,v) = (fn1(um1)*diff_y-fn2(um1)*diff_x) ! pointing outward 
    end do
  end subroutine PPD

  subroutine PFN ! outward directed normal vectors of crystal surfaces 
  
  real(8) avek1, avek2, avek3, bvek1, bvek2, bvek3, cvek1, cvek2, cvek3, laenge
    avek1 = p1(mainu, 2) - p1(mainu, 1)
    avek2 = p2(mainu, 2) - p2(mainu, 1)
    avek3 = p3(mainu, 2) - p3(mainu, 1)
    bvek1 = p1(mainu, 3) - p1(mainu, 2)
    bvek2 = p2(mainu, 3) - p2(mainu, 2)
    bvek3 = p3(mainu, 3) - p3(mainu, 2)
    cvek1 = avek2*bvek3 - avek3*bvek2 ! c = a x b 
    cvek2 = avek3*bvek1 - avek1*bvek3 ! points inward 
    cvek3 = avek1*bvek2 - avek2*bvek1
    laenge = Sqrt(cvek1*cvek1 + cvek2*cvek2 + cvek3*cvek3)
    fn1(mainu) = -cvek1/laenge ! normalized & points outward 
    fn2(mainu) = -cvek2/laenge
    fn3(mainu) = -cvek3/laenge
  end subroutine PFN

  subroutine PROT           ! crystal rotation 
  integer(8) i, j
  real(8) r11, r12, r13, &
       r21, r22, r23, &
       r31, r32, r33, & ! Euler matrix do rotation 
       s1, s2, s3, c1, c2, c3 ! Sin, Cos of Euler angles 

    call random_number(rand)
    alpha_euler = twopi*(rand)
    call random_number(rand)
    beta_euler  = dacos(1.0 - 2.0*dble(rand))
    call random_number(rand)
    gamma_euler = twopi*rand

    s1 = dsin(alpha_euler) 
    s2 = dsin(beta_euler) 
    s3 = dsin(gamma_euler)  
    c1 = dcos(alpha_euler) 
    c2 = dcos(beta_euler) 
    c3 = dcos(gamma_euler)
 
    r11 = -c2*s1*s3 + c1*c3
    r12 = -c2*s1*c3 - c1*s3
    r13 =  s2*s1
    r21 =  c2*c1*s3 + s1*c3
    r22 =  c2*c1*c3 - s1*s3
    r23 = -s2*c1
    r31 =  s2*s3
    r32 =  s2*c3
    r33 =  c2

    do i = 1, um
      do j = 1, vm(i)            ! translation to the center of grav. 
        p4(i, j) = p1(i, j) - sp1
        p5(i, j) = p2(i, j) - sp2
        p6(i, j) = p3(i, j) - sp3
      end do
    end do

    do i = 1, um
      do j = 1, vm(i) 
        p1(i, j) =  p4(i, j)*r11 + p5(i, j)*r12 + p6(i, j)*r13
        p2(i, j) =  p4(i, j)*r21 + p5(i, j)*r22 + p6(i, j)*r23
        p3(i, j) =  p4(i, j)*r31 + p5(i, j)*r32 + p6(i, j)*r33 
      end do
    end do
    do i = 1, um
      do j = 1, vm(i)
        p4(i, j) = p1(i, j)
        p5(i, j) = p2(i, j)
        p6(i, j) = p3(i, j)
      end do
    end do
    call PEX

    do i = 1, um
      do j = 1, vm(i) 
        p6(i, j) = p6(i, j) - zmax - 1
        p3(i, j) = p6(i, j)
      end do
    end do
  end subroutine PROT

  subroutine PINIT ! initializes this and that 
  
  integer(8) i, j, k
  real(4) seed
  real(8) nr, ni, nr1, ni1, nr2, ni2, k1, k2, krel, ref1, ref2, ref3
  real(8) arg, fac1, theta_i, theta_t, ref6, nrel
  real(8) r11, r12, r33, r34, t11, t12, t33, t34
  character(8) date
  character(10) time
  character(5) zone
  integer values(8)

    call system_clock(start_time)
    nerror = 0
    pi = 4*datan(1.d0) 
    twopi = 2*pi 
    pitwo = pi/2  
    deg = 180/pi 
    rad = pi/180
    gwu = 1.0E-11 ! math. epsilon to avoid numerical errors 
    epsilon = 1.0d-7 ! ### for delta peaks    
!    process_id = getpid_() ! IBM ###
    call random_seed()    
! ### maschine dependent, seeds random number generator 

    resol_fak = dble(rtddim)/180. ! 1/resol_fak = angular resolution in degree 
    theta_diff = rtddim/2 ! -> polydiff
    rad2resol_fak = rad/resol_fak

    rtd = 0 ! init. array rtd 
    utr = 0
    energy_in = 0. 
    energy_out = 0.
    energy_abs = 0.
    nt = 0
    ntotal = 0
    ntm_sum = 0
    gcs_sum = 0.
    poly_orient = 0

    mabs = -4.*Pi/la*ibi ! absorption coefficient 
    myfunc_arg = 1/rbi*1/Sqrt( 1 - (ibi/rbi)*(ibi/rbi))
    if (myfunc_arg .gt. 1.0) myfunc_arg = 1.0 ! ###  
    grenzw = Dasin(myfunc_arg) ! angle of total reflection and 
    angle_g = dint(grenzw*2/pi*angle_max) ! corresp. bin 
    cosgrenzw = dcos(grenzw)
    waveno = twopi/la
    diffconst1 = 1./(2.*waveno*waveno )
    diffconst2 = waveno*waveno/4./pi/pi 
    rayconst1 = 2./pi*angle_max
    la2 = la*la

    help1 = 0.5*rad ! avoid pi/2, pi, 3/2*pi, 2pi in azimuth integration
    do i = 1, n_phi_diff ! azimuth integration in polydiff
      phi_diff(i) = twopi/(n_phi_diff - 1)*(i - 1) + help1
    end do
    simpson = (phi_diff(n_phi_diff) - phi_diff(1))/n_phi_diff/3. ! Simpson integ.
    simpson_fak(0) = 1.
    simpson_fak(1) = 4.
    simpson_fak(2) = 1.

    rot_stokes = 0.
    ra = 0.
    ta = 0.
    rb = 0.
    tb = 0.
    rot_stokes(1,1) = 1. 
    rot_stokes(4,4) = 1.
    id_mat = 0.
    do i = 1, 4
      id_mat(i,i) = 1.
    end do

    nr1 = 1.0
    ni1 = 0.0
    nr2 = rbi
    ni2 = ibi

    k1 = ni1/nr1
    k2 = ni2/nr2
    krel = (k2 - k1)/(1 + k1*k2)
    nrel = nr2/nr1*(1 + k1*k2)/(1 + k1*k1)
    ref1 = nrel*nrel
    ref2 = krel*krel
    ref3 = (1 + ref2)*(1 + ref2)
    ref6 = ref1*ref3/(1 + krel*k2)/(1 + krel*k2)
    fac1 = pi/2./angle_max

    do i = 0, angle_max - 1
      theta_i = (i + 0.5)*fac1
      if (theta_i .le. pi/2) then 
        call refraction(ref1, ref2, ref3, ref6, theta_i, theta_t, krel, k2)
call fresnel(nr1, ni1, nr2, ni2, theta_i, r11, r12, r33, r34, t11, t12, t33, t34)
        refrwa(i) = theta_t
        help5 = 0.5 
        help6 = 0.5*cos(theta_t)*nr2/nr1/cos(theta_i)

        ra(1,1,i) = help5*r11 
        ra(1,2,i) = help5*r12
        ra(3,3,i) = help5*r33
        ra(3,4,i) = help5*r34 
        ta(1,1,i) = help6*t11
        ta(1,2,i) = help6*t12
        ta(3,3,i) = help6*t33
        ta(3,4,i) = help6*t34
        ra(2,2,i) = ra(1,1,i)
        ra(2,1,i) = ra(1,2,i)
        ra(4,4,i) = ra(3,3,i)
        ra(4,3,i) =-ra(3,4,i)
        ta(2,2,i) = ta(1,1,i)
        ta(2,1,i) = ta(1,2,i)
        ta(4,4,i) = ta(3,3,i)
        ta(4,3,i) =-ta(3,4,i)

        w = dcos(theta_i)
        help1 = sqrt(1 - w*w)
        if (help1 .eq. 0) then 
          wfa(i) = 0 ! ###
        else
          wfa(i) = dsin(refrwa(i))/help1
        endif
        cos_theta_t_a(i) = dcos(refrwa(i))

      endif
    enddo

    nr1 = rbi
    ni1 = ibi
    nr2 = 1.0
    ni2 = 0.0
    k1 = ni1/nr1
    k2 = ni2/nr2
    krel = (k2 - k1)/(1 + k1*k2)
    nrel = nr2/nr1*(1 + k1*k2)/(1 + k1*k1)
    ref1 = nrel*nrel
    ref2 = krel*krel
    ref3 = (1 + ref2)*(1 + ref2)
    ref6 = ref1*ref3/(1 + krel*k2)/(1 + krel*k2)
    do i = 0, angle_max - 1
      theta_i = (i + 0.5)*fac1
      if (theta_i .le. grenzw) then ! Refrac. angle & Fresnel-mat.: Ref. & Trans  
        call refraction(ref1, ref2, ref3, ref6, theta_i, theta_t, krel, k2)

call fresnel(nr1, ni1, nr2, ni2, theta_i, r11, r12, r33, r34, t11, t12, t33, t34)

        refrwb(i) = theta_t  
        help5 = 0.5 
        help6 = 0.5*cos(theta_t)*nr2/nr1/cos(theta_i)
 
        rb(1,1,i) = help5*r11
        rb(1,2,i) = help5*r12
        rb(3,3,i) = help5*r33
        rb(3,4,i) = help5*r34
        tb(1,1,i) = help6*t11
        tb(1,2,i) = help6*t12
        tb(3,3,i) = help6*t33
        tb(3,4,i) = help6*t34

        rb(2,2,i) = rb(1,1,i)
        rb(2,1,i) = rb(1,2,i)
        rb(4,4,i) = rb(3,3,i)
        rb(4,3,i) =-rb(3,4,i)
        tb(2,2,i) = tb(1,1,i)
        tb(2,1,i) = tb(1,2,i)
        tb(4,4,i) = tb(3,3,i)
        tb(4,3,i) =-tb(3,4,i)

        w = dcos(theta_i)
        help1 = sqrt(1 - w*w)
        if (help1 .eq. 0) then 
          wfb(i) = 0
        else
          wfb(i) = dsin(refrwb(i))/help1
        endif
        cos_theta_t_b(i) = dcos(refrwb(i))

      else ! Fresnel-matr. for total internal reflection

call fresnel(nr1, ni1, nr2, ni2, theta_i, r11, r12, r33, r34, t11, t12, t33, t34)
        help5 = 0.5
        rb(1,1,i) = help5*r11
        rb(1,2,i) = help5*r12
        rb(3,3,i) = help5*r33
        rb(3,4,i) = help5*r34
        rb(2,2,i) = rb(1,1,i)
        rb(2,1,i) = rb(1,2,i)
        rb(4,4,i) = rb(3,3,i)
        rb(4,3,i) =-rb(3,4,i)
      end if
    end do


    call PEX ! -> lq_max 
    lq_max = (xmax - xmin)*(xmax - xmin) + & 
             (ymax - ymin)*(ymax - ymin) + &
             (zmax - zmin)*(zmax - zmin)

  end subroutine PINIT

  subroutine refraction(ref1, ref2, ref3, ref6, theta_i, theta_t, krel, k2)

  real(8) ref1, ref2, ref3, ref6, theta_i, theta_t, krel, k2
  real(8) sintiq, ref4, ref5, q4, q2, g, test1, test2
  real(8) ref7, rnstar

!  Calculates the angle of refraction for transmission from 
!  (nr1, ni1) -> (nr2,ni2)

    sintiq = sin(theta_i)*sin(theta_i)
    ref4 = 1 - (1 - ref2)/ref1/ref3*sintiq
    ref5 = 2*krel/ref1/ref3*sintiq
    q4 = ref4*ref4 + ref5*ref5
    q2 = sqrt(q4)
    g = asin(ref5/q2)/2
    test1 = acos(ref4/q2)/2
    test2 = atan(ref5/ref4)/2
    g = test1
    ref7 = (cos(g) - k2*sin(g))*(cos(g) - k2*sin(g))
    rnstar = sqrt(sintiq + ref6*q2*ref7)        
    theta_t = asin(sin(theta_i)/rnstar)

  end subroutine refraction

  subroutine fresnel(nr1, ni1, nr2, ni2, theta_i,r11, r12, r33, r34, t11, t12, t33, t34)

!   input: complex refractive index of medium 1: cn1 = nr1 + i ni1
!                         "                   2: cn2 = nr2 + i ni2
!          angle of incidence in rad. : theta_i
!
!   output: fresnel coefficients for transforming field amplitudes
!           ats,..,brp , where
!   a.. == real part, b.. == imaginary part   
!   .r. == reflection, .t. == transmission
!   ..p == parallel, ..s = perpendicular

  real(8) nr1, ni1, nr2, ni2, theta_i, sinti
  real(8) r11, r12, r33, r34, t11, t12, t33, t34
  complex(8) ci, cn1, cn2, csintiq, costi, chelp1, cts, crs, ctp, crp

     ci = dcmplx(0,-1) 
     cn1 = dcmplx(nr1, ni1)
     cn2 = dcmplx(nr2, ni2)
     sinti = dsin(theta_i)
     csintiq = dcmplx(sinti*sinti)
     costi = dcmplx(cos(theta_i))
     chelp1 = cdsqrt(cn2*cn2 - cn1*cn1*csintiq)
       
     cts = 2.*cn1*costi/(cn1*costi + chelp1)
     crs = (cn1*costi - chelp1)/(cn1*costi + chelp1)
     ctp = 2.*cn1*cn2*costi/(cn2*cn2*costi + cn1*chelp1)
     crp = (cn2*cn2*costi - cn1*chelp1)/(cn2*cn2*costi + cn1*chelp1)

     r11 = Real(crp*Conjg(crp) + crs*Conjg(crs))
     r12 = Real(crp*Conjg(crp) - crs*Conjg(crs))
     r33 = 2.*Real(crp*Conjg(crs))
     r34 = 2.*Aimag(crp*Conjg(crs))
     t11 = Real(ctp*Conjg(ctp) + cts*Conjg(cts))
     t12 = Real(ctp*Conjg(ctp) - cts*Conjg(cts))
     t33 = 2.*Real(ctp*Conjg(cts))
     t34 = 2.*Aimag(ctp*Conjg(cts))

  end subroutine fresnel


  subroutine PGCS
! det. exact geoemtrical cross section, given by a closed polygon.
!  Note: The subroutine works, but should be rewritten in a more 
!        transparent way. Maybe later?! 
  
  integer(8) kand(1:1000), i, j, kand_max, s_start, p_start, endp
  real(8) xmin_s, t1, t2, x1_s, y1, x2_s, y2, z, next_kand_x, next_kand_y, &
  len_min, n, p1x, p2x, p3x, p4x, p1y, p2y, p3y, p4y, w_max, eps, angle
  logical treffer

    eps = 1.0e-5 ! ###

    Kand_max = 0 ! select those planes which 'face' upward 
    do i = 1, um 
      if (fn3(i) .gt. 0) then
        kand_max = kand_max + 1
        kand(kand_max) = i ! kand cont. surface numbers 
      end if
    end do
  
    p1x = 0 
    p1y = 0 
    p2x = 0 
    p2y = -1
    xmin_s = 0 ! search do a save candidate to start with 
    do i = 1, kand_max 
      do j  = 1, vm(kand(i))       
        if (p1(kand(i),j) .lt. xmin_s) then
          xmin_s = p1(kand(i),j) 
          s_start = kand(i) ! surface with miminum x 
          p_start = j ! corresponding point 
          if (j .eq. vm(kand(i))) then 
            endp = 1 
          else 
            endp = j + 1
          end if
          p3x = p1(kand(i),j)
          p3y = p2(kand(i),j)
          p4x = p1(kand(i),endp)
          p4y = p2(kand(i),endp)  
          x1_s = p2x - p1x
          y1 = p2y - p1y
          nf = sqrt(x1_s*x1_s + y1*y1)
          x1_s = x1_s/nf 
          y1 = y1/nf ! direction of X_1 
          x2_s = p4x - p3x 
          y2 = p4y - p3y
          nf = sqrt(x2_s*x2_s + y2*y2)
          x2_s = x2_s/nf 
          y2= y2/nf ! direction of X_2 
          test1 = x2_s*y1 - y2*x1_s ! -> S. 123 
          w = -x1_s*x2_s - y1*y2 ! dcos(Winkel zw. Richtungsvektoren) 
          if (test1 .ge. 0) angle = dacos(w)
          if (test1 .lt. 0) angle = twopi - dacos(w)
          w_max = angle      
        end if
        if (p1(kand(i),j) .eq. xmin_s) then     
          if (j .eq. vm(kand(i))) then 
            endp = 1 
          else 
            endp = j + 1
          end if
          p3x = p1(kand(i),j)
          p3y = p2(kand(i),j)
          p4x = p1(kand(i),endp)
          p4y = p2(kand(i),endp)  
          x1_s = (p2x - p1x) 
          y1 = (p2y - p1y)
          nf = sqrt(x1_s*x1_s + y1*y1)
          x1_s = x1_s/nf 
          y1 = y1/nf ! direction of X_1 
          x2_s = (p4x - p3x) 
          y2 = (p4y - p3y) 
          nf = sqrt(x2_s*x2_s + y2*y2)
          x2_s = x2_s/nf 
          y2 = y2/nf ! direction of X_2 
          test1 = x2_s*y1 - y2*x1_s ! -> S. 123 
          w = -x1_s*x2_s - y1*y2 ! dcos(Winkel zw. Richtungsvektoren) 
          if (test1 .ge. 0) then 
            angle = dacos(w)
          else
            angle = twopi - dacos(w)
          end if
          if (angle .gt. w_max) then 
            w_max = angle      
            s_start = kand(i) ! surface with miminum x 
            p_start = j ! corresponding point 
          end if
        end if
      end do ! j
    end do ! i

    poly_x(1) = p1(s_start,p_start) ! first polygon point 
    poly_y(1) = p2(s_start,p_start)
    if (p_start .eq. vm(s_start)) then
      endp = 1 
    else 
      endp = p_start + 1
    end if
    next_kand_x = p1(s_start, endp)           ! next canditate right 
    next_kand_y = p2(s_start, endp)        ! handed to the prev. one 

    n_poly = 0
    do while (n_poly .eq. 0 .or. & ! Herrjeh!!
              (poly_x(n_poly + 1) .ne. poly_x(1) .or. & 
              poly_y(n_poly + 1) .ne. poly_y(1)) .and. &
              n_poly .lt. 900) ! ### numerical reasons

      n_poly = n_poly + 1

      p1x = poly_x(n_poly) 
      p1y = poly_y(n_poly)
      p2x = next_kand_x
      p2y = next_kand_y
      treffer = .false.
      len_min = 1  
      w_max = 0

      do i = 1, kand_max 
        do j  = 1, vm(kand(i)) 
          if (j .eq. vm(kand(i))) then 
            endp = 1 
          else 
            endp = j + 1
          end if
          p3x = p1(kand(i),j)
          p3y = p2(kand(i),j)
          p4x = p1(kand(i),endp)
          p4y = p2(kand(i),endp)

          n = (p2x - p1x)*(p4y - p3y) - (p2y - p1y)*(p4x - p3x) ! S.123 
          if (dabs(n) .gt. eps) then  ! avoid parallel lines 
            z = (p1y - p3y)*(p4x - p3x) - (p1x - p3x)*(p4y - p3y) ! " 
            t1 = z/n ! -> intersection at test line 
            z = (p1y - p3y)*(p2x - p1x) - (p1x - p3x)*(p2y - p1y)
            t2 = z/n
            if (eps .lt. t1 .and. t1 .lt. 1 - eps) then ! t1 aus )0,1( 
              if (eps .lt. t2 .and. t2 .lt. 1 - eps) then ! t2 aus )0,1( 
                if (t1 .lt. len_min) then ! search do smallest t1 
                  treffer = .true.
                  len_min = t1 ! new mimimum length  
                  poly_x(n_poly + 1) = t1*p2x + (1 - t1)*p1x
                  poly_y(n_poly + 1) = t1*p2y + (1 - t1)*p1y
                  next_kand_x = p4x
                  next_kand_y = p4y
                end if 
              end if
            end if ! case (iv)a on page 122 
            if (dabs(t1 - 1) .lt. eps .and. dabs(t2) .lt. eps) then
!   numerically for t1 = 1 and t2 = 0 
              if (treffer .eqv. .false.) then ! treffer hat Vorrang! 
                x1_s = (p2x - p1x) 
                y1 = (p2y - p1y)
                nf = sqrt(x1_s*x1_s + y1*y1)
                x1_s = x1_s/nf 
                y1 = y1/nf ! direction of X_1 
                x2_s = (p4x - p3x) 
                y2 = (p4y - p3y) 
                nf = sqrt(x2_s*x2_s + y2*y2)
                x2_s = x2_s/nf 
                y2= y2/nf ! direction of X_2 
                test1 = x2_s*y1 - y2*x1_s ! -.gt. S. 123 
                w = -x1_s*x2_s - y1*y2 ! dcos(Winkel zw. Richtungsvektoren) 

                if (dabs(w) .gt. 1) w = w/dabs(w) ! numerically reasons 
                if (test1 .ge. 0) then 
                  angle = dacos(w)
                else 
                  angle = twopi - dacos(w)
                endif
                if (angle .gt. w_max) then ! search do largest angle 
                  w_max = angle
                  poly_x(n_poly + 1) = p3x ! oder p2x, egal! 
                  poly_y(n_poly + 1) = p3y
                  next_kand_x = p4x
                  next_kand_y = p4y                
                end if ! search do largest angle 
              end if ! treffer hat Vorrang! 
            end if ! t2 = 1 and t2 = 0 
          end if ! avoid parallel lines 
        end do ! point loop 
      end do ! surface loop 
    
    end do ! while

    if (n_poly .lt. 900) then ! ###
      area = 0
      do i = 1, n_poly 
        area = area + poly_x(i)*poly_y(i+1) - poly_x(i+1)*poly_y(i)
      end do
      area = area/2.
      poly_sum = poly_sum - area
      poly_orient = poly_orient + 1
      gcs_mean_poly = poly_sum/dble(poly_orient) 
    end if

  end subroutine PGCS

  subroutine polydiff ! diffraction at crystal projection <- PGCS 
  
  real(8) dx(1:100), dy(1:100), theta, sint, sint2, sinp, cosp, fac1, l, &
       sumre, sumim, sum, qj, qj1, uj, rj, re, im, phi_av
  integer(8) i, j, it, ip, jp
  if (n_poly .lt. 900) then ! ###
    do i = 1, n_poly 
      dx(i) = poly_x(i+1) - poly_x(i)
      dy(i) = poly_y(i+1) - poly_y(i)
    end do
    fac1 = -diffconst2/area
    do it = 1, theta_diff ! no forward diffraction here 
      theta = dble(it)*rad2resol_fak ! diffraction exactly at angel, NO BIN AVER.
      sint = dsin(theta) 
      sint2 = sint*sint
      l = diffconst1/sint2
      phi_av = 0
      do ip = 1, n_phi_diff 
        sinp = dsin(phi_diff(ip)) 
        cosp = dcos(phi_diff(ip))
        sumre = 0 
        sumim = 0
        qj = waveno*(poly_x(1)*cosp + poly_y(1)*sinp)*sint
        do j = 1, n_poly 
          rj = 1/(dx(j)*cosp + dy(j)*sinp)
          uj = dx(j)/sinp - dy(j)/cosp
          qj1=waveno*(poly_x(j+1)*cosp + poly_y(j+1)*sinp)*sint
          sumre = sumre + rj*(dcos(qj1) - dcos(qj))*uj
          sumim = sumim + rj*(dsin(qj1) - dsin(qj))*uj
          qj = qj1
        end do
        re = l*sumre
        im = l*sumim
        sqrint(ip) = re*re + im*im
      end do
      sum = 0
      do ip = 1, n_phi_diff - 2, 2 ! Integrate along phi
        do jp = 0, 2 ! Simpson
          sum = sum + sqrint(ip + jp)*simpson_fak(jp)
        end do
      end do
      phi_av = sum*simpson*fac1/twopi
      diff_total(it) = diff_total(it) + phi_av
    end do ! theta-loop 
    diff_total(0) = diff_total(0) - area/la2 ! add direct forward peak
  end if
  end subroutine polydiff

end program rt! program rt containing all subroutines



