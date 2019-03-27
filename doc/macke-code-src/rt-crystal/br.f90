program bullet_rosettes
implicit none

  integer, parameter ::  MaxFlaechenzahl = 100
  integer, parameter ::  MaxEckenzahl = 10   ! Max. Zahl der Ecken pro Flaeche 
  integer, parameter ::  MaxPunktezahl = 100   ! Max. Zahl der Eckp. im Kristall 
  real p1(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real p2(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real p3(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real p4(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real p5(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real p6(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real pf1(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real pf2(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real pf3(1:MaxFlaechenzahl, 1:MaxEckenzahl)
  real fnsign(1:MaxFlaechenzahl)
  real ek1(1:MaxPunktezahl)
  real ek2(1:MaxPunktezahl)
  real ek3(1:MaxPunktezahl)

  integer vm(1:MaxFlaechenzahl)
  integer vmf(1:MaxFlaechenzahl)

  real xmax, ymax, zmax, xmin, ymin, zmin, x_rot, y_rot, z_rot, sp1, sp2, sp3, r, L, h, pi

  integer pm, pn, um, u, v, i, j, bullet_max, um_total, n
  character(1)   cha, chb, chc
  character(100) afn, ifn, vfn


  print *, ' Builds coordinate-file for a bullet rosette. '
  print *, ' You have to have an input file of the form:'
  print *, ' bullet_max ! number of branches, and for each branch'
  print *, ' L, r, h, x_rot, y_rot, z_rot  which means '
  print *, ' length, radius, pyramide-height, Rotation angles [degree].'
  print *, ' Original direction of L-axis is (0,0,1) '
  print *, ' '
  print *, ' You haven''t? Then type nofile'
  print *, 'bullet_rosette> Inputfile : '
  read(5, fmt = '(A100)') ifn
  if (ifn(1:6) .eq. 'nofile') then 
    print *, 'Number of branches: '
    read *, bullet_max
    do i = 1, bullet_max
      print *, 'bullet ', i,': Give L, r, h, x_rot, y_rot, z_rot'
      read *, L, r, h, x_rot, y_rot, z_rot
      call Make_Bullet
      call Rotate_Bullet
      call Translate_Bullet
      call Add_Bullet_to_Rosette
    end do
    call Write_Bullet_Rosett
  else 
    open(unit = 10, file = ifn, status = 'old', form = 'formatted')
    read(10,*) bullet_max ! Anzahl der Bullets 
    call Init_All
    do i = 1, bullet_max
      read(10,*) L, r, h, x_rot, y_rot, z_rot
      call Make_Bullet
      call Rotate_Bullet
      call Translate_Bullet
      call Add_Bullet_to_Rosette
    end do
    call Write_Bullet_Rosett
  end if

  contains 

  subroutine PEX            ! Bestimmung der Extrema 
  integer i,j

    xmax = p1(1, 1)
    ymax = p2(1, 1)
    zmax = p3(1, 1)
    xmin = p1(1, 1)
    ymin = p2(1, 1)
    zmin = p3(1, 1)
    do i = 1, um
      do j = 1, vm(i)    
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
 
  subroutine Make_Bullet
  real f1, z1, z2
  integer u, v, i, j
 
    f1 = Sqrt(3.)/2.
    z1 = -L - h
    z2 = -h
  
    ek1( 1) = 0     
    ek2( 1) = -r   
    ek3( 1) = z1
    ek1( 2) = r*f1
    ek2( 2) = -r/2
    ek3( 2) = z1
    ek1( 3) = r*f1
    ek2( 3) = r/2
    ek3( 3) = z1
    ek1( 4) = 0
    ek2( 4) = r
    ek3( 4) = z1
    ek1( 5) = -r*f1
    ek2( 5) = r/2
    ek3( 5) = z1
    ek1( 6) = -r*f1
    ek2( 6) = -r/2
    ek3( 6) = z1
  
    ek1( 7) = 0
    ek2( 7) = -r
    ek3( 7) = z2
    ek1( 8) = r*f1
    ek2( 8) = -r/2
    ek3( 8) = z2
    ek1( 9) = r*f1
    ek2( 9) = r/2
    ek3( 9) = z2
    ek1(10) = 0
    ek2(10) = r 
    ek3(10) = z2
    ek1(11) = -r*f1
    ek2(11) = r/2
    ek3(11) = z2
    ek1(12) = -r*f1
    ek2(12) = -r/2 
    ek3(12) = z2
  
    ek1(13) = 0     
    ek2(13) = 0    
    ek3(13) = 0 
  
    open(30, file = 'ssb.Points', status = 'old')
  
    u = 0
    chb = 'y'
    do while (chb .ne. 'n')
      u = u + 1
      v = 0
      pn = 1
      do while (pn .ne. 0)
        v = v + 1
        read(30,*) pn
        if (pn .gt. 0) then
          p1(u, v) = ek1(pn)
          p2(u, v) = ek2(pn)
          p3(u, v) = ek3(pn)
        end if
      end do
      vm(u) = v - 1
      READ(30,fmt = '(A1)') chb
    end do
    um = u
    close(30)

    do i = 1, um
      do j = 1, vm(i)
        p4(i, j) = p1(i, j)
        p5(i, j) = p2(i, j)
        p6(i, j) = p3(i, j)
      end do
    end do 

  end subroutine Make_Bullet

  subroutine Rotate_Bullet
  integer i,j
  real zmaxa, siga, coga, sild, cold, sips, cops

    call PEX
    zmaxa = zmax
    do i = 1, um
      do j = 1, vm(i) ! Translation zum Schwerpunkt 
        p4(i, j) = p4(i, j) - sp1
        p5(i, j) = p5(i, j) - sp2
        p6(i, j) = p6(i, j) - sp3
      end do
    end do

    y_rot = y_rot * pi / 180.0 ! Drehung um y - Achse 
    siga = SIN(y_rot)
    coga = COS(y_rot)
    do i = 1, um
      do j = 1, vm(i) 
        p1(i, j) = p4(i, j) * coga + p6(i, j) * siga
        p2(i, j) = p5(i, j)
        p3(i, j) = -p4(i, j) * siga + p6(i, j) * coga
     end do
    end do
    do i = 1, um
      do j = 1, vm(i)
        p4(i, j) = p1(i, j)
        p5(i, j) = p2(i, j)
        p6(i, j) = p3(i, j)
      end do
    end do

    x_rot = x_rot * pi / 180.0  ! Drehung um x - Achse 
    sild = SIN(x_rot)
    cold = COS(x_rot)
    do i = 1, um
      do j = 1, vm(i)
        p1(i, j) =  p4(i, j)
        p2(i, j) =  p5(i, j) * cold + p6(i, j) * sild
        p3(i, j) = -p5(i, j) * sild + p6(i, j) * cold
      end do
    end do
    do i = 1, um
      do j = 1, vm(i)
        p4(i, j) = p1(i, j)
        p5(i, j) = p2(i, j)
        p6(i, j) = p3(i, j)
      end do
    end do

    z_rot = z_rot * pi / 180.0
    sips = SIN(z_rot) ! Drehung um z-Achse 
    cops = COS(z_rot)
    do i = 1, um
      do j = 1, vm(i)
        p1(i, j) =  p4(i, j) * cops + p5(i, j) * sips
        p2(i, j) = -p4(i, j) * sips + p5(i, j) * cops
        p3(i, j) =  p6(i, j)
      end do
    end do
    do i = 1, um
      do j = 1, vm(i)
        p4(i, j) = p1(i, j)
        p5(i, j) = p2(i, j)
        p6(i, j) = p3(i, j)
      end do
    end do

  end subroutine Rotate_Bullet

  subroutine Translate_Bullet
  real T1, T2, T3
  integer i, j

    T1 = p1(13,2) 
    T2 = p2(13,2) 
    T3 = p3(13,2)

    do i = 1, um  ! Translation zurueck 
      do j = 1, vm(i)
        p4(i, j) = p4(i, j) - T1
        p5(i, j) = p5(i, j) - T2
        p6(i, j) = p6(i, j) - T3
        p1(i, j) = p4(i, j)
        p2(i, j) = p5(i, j)
        p3(i, j) = p6(i, j)
      end do
    end do

  end subroutine Translate_bullet

  subroutine Add_Bullet_to_Rosette
  integer i, j

    um_total = um_total + um
    do i = 1, um
      vmf(um_total - um + i) = vm(i)
      do j = 1, vm(i)
        pf1(um_total - um + i, j) = p1(i,j)
        pf2(um_total - um + i, j) = p2(i,j)
        pf3(um_total - um + i, j) = p3(i,j)
      end do
    end do

  end subroutine Add_Bullet_to_Rosette 

  subroutine Write_Bullet_Rosett
  integer i, j, k

    print *, 'bullet> Give output file:  '
    READ(6, fmt = '(A100)') afn
    open(20, file = afn, status = 'unknown')
    WRITE(20,*)  um_total
    do i = 1, um_total
      WRITE(20, *) vmf(i)
    end do
    do i = 1, um_total
      do j = 1, vmf(i)
        WRITE(20, *) pf1(i, j), pf2(i, j), pf3(i, j)
      end do
    end do
    CLOSE(20)

  end subroutine Write_Bullet_Rosett

  subroutine Init_All
    um_total = 0
    pi = 4.*atan(1.)
  end  subroutine Init_All


end program bullet_rosettes
