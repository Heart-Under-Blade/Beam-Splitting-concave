program column_plates
implicit none

  real, parameter ::    pi = 3.141592653589
  integer, parameter ::  MaxFlaechenzahl = 8
  integer, parameter ::  MaxEckenzahl = 6   ! Max. Zahl der Ecken pro Flaeche 
  integer, parameter ::  MaxPunktezahl = 24   ! Max. Zahl der Eckp. im Kristall 

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
		
  real xmax, ymax, zmax, xmin, ymin, zmin, f1, z1, z2
  real x_rot, y_rot, z_rot, sp1, sp2, sp3, r, L, h

  integer vm(1:MaxFlaechenzahl)
  integer vmf(1:MaxFlaechenzahl)

  integer pm, pn, um, u, v, i, j, bullet_max, um_total, n
  character(1) cha, chb, chc
  character(100)  afn, ifn, vfn


  print *, 'colpla> Length, radius: '
  read(5,*) l, R
  call Make_PLate
  call Write_Plate

  contains

  subroutine Make_Plate 

    f1 = sqrt(3.)/2
    z1 = 0
    z2 = -L
  
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
  
    open(10,file = 'Hex.Points', form = 'formatted')
  
    u = 0
    chb = 'y'
    do while (chb .ne. 'n')
      u = u + 1
      v = 0
      pn = 1
      do while (pn .ne. 0)
        v = v + 1
        read(10,*) pn
        if (pn .gt. 0) then
          p1(u, v) = ek1(pn)
          p2(u, v) = ek2(pn)
          p3(u, v) = ek3(pn)
        end if
      end do
      vm(u) = v - 1
      READ(10,fmt = '(A1)') chb
    end do
    um = u
    close(10)

  end subroutine Make_Plate 


  subroutine Write_Plate
    print *, 'colpla> Crystal filename (extension is .crystal): '
    read(5,*) afn
    open(10,file = afn)
    write(10,*) um
    do j = 1, um
      write(10,*)  vm(j)
    end do
    do i = 1, um
      do j = 1, vm(i)
        write(10,*) p1(i, j), p2(i, j), p3(i, j)
      end do
    end do
    CLOSE(10)

  end subroutine Write_Plate

end program column_plates

