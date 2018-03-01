program fractal3D
implicit none

integer, parameter ::    fl_max = 10000 ! max # of surfaces 
integer, parameter ::    pt_max = 3 ! max # of points per surface 
    
real main_p1(1:fl_max, 1:pt_max), main_p2(1:fl_max, 1:pt_max)
real main_p3(1:fl_max, 1:pt_max), sub_p1(1:fl_max, 1:pt_max)
real sub_p2(1:fl_max, 1:pt_max), sub_p3(1:fl_max, 1:pt_max)
real fn1(1:fl_max), fn2(1:fl_max), fn3(1:fl_max)
   
real sp1, sp2, sp3  ! surface center of gravity 
real sph1, sph2, sph3  ! local pyramid height 
real len ! length of triangle side 
real help1, help2, help3, rdummy
real height ! pyramid height  
real random_mag, random_arg ! distortion magnitude, seed of RNG 
real true_length !
real rand ! arguemt of random_number
    
integer n_gen ! # of iterations 
integer i, j, n_sub, n_max
integer(4) process_id, getpid_ ! seed of random number generator
  
character(100) afn, ofn, opop
character(1) random_char

  random_char = 'n'
  print *, ' Generation: '; read *, n_gen ! # of generations 
  print *, ' True size: '; read *, true_length 
  print *, ' Filename: '; read *, ofn ! outputfilename 
  print *, ' Mit Random-Verzerrungen (y/n) ? : '; read *, random_char
  call PDAL ! read Initiator 
   
  help1 = main_p1(1,1)-main_p1(1,2)
  help2 = main_p2(1,1)-main_p2(1,2)
  help3 = main_p3(1,1)-main_p3(1,2)
  len = Sqrt(help1*help1 + help2*help2 + help3*help3)            
             
  if (random_char .eq. 'y') then 
    print *,' Verzerrungsmass (0..1): '; read *, random_mag
!    process_id = getpid_() ! IBM ###
!    random_arg = abs(sin(log(real(process_id))))
    call random_seed
  else 
    random_mag = 0
  end if

  do i = 1, n_gen 
    len = len/2
    n_sub = 0
    do j = 1, n_max 
      call PFN ! surface normal
    end do 
    do j = 1, n_max
      call MakeNextGeneration
      n_sub = n_sub + 6
    end do
    n_max = n_sub
    call ReplaceOldGeneration
  end do
   
  call WriteFractal

  contains 


subroutine PDAL           ! Kristalldaten lesen 
integer i,j, dummy

  open(10, file = 'tetraeder', status = 'old')
  read(10,*) n_max
  do i = 1, n_max
    read(10,*) dummy
  end do
  do i = 1, n_max
    do j = 1, 3
      read(10,*) main_p1(i, j), main_p2(i, j), main_p3(i, j)
    end do
  end do
  close(10)

  do i = 1, n_max
    do j = 1, 3
      main_p1(i,j) = main_p1(i,j)*true_length
      main_p2(i,j) = main_p2(i,j)*true_length
      main_p3(i,j) = main_p3(i,j)*true_length
    end do
  end do

end subroutine PDAL
 
subroutine MakeNextGeneration

real  p11, p12, p13, p21, p22, p23, p31, p32, p33
real  h11, h12, h13, h21, h22, h23, h31, h32, h33
real  x11, x12, x13, x21, x22, x23, h1, h2, h3
real  a11, a12, a13, a21, a22, a23, a31, a32, a33
real  norm, fak1, fak2, zerr

  p11 = main_p1(j,1) ! Eckpunkte des 'alten' Dreiecks 
  p12 = main_p2(j,1) ! keine Verzerrung! 
  p13 = main_p3(j,1)
  p21 = main_p1(j,2) 
  p22 = main_p2(j,2)
  p23 = main_p3(j,2)
  p31 = main_p1(j,3) 
  p32 = main_p2(j,3)
  p33 = main_p3(j,3)

  call random_number(rand)
  zerr = 2*random_mag*rand + 0.5 - random_mag 
  h11 = p11 + zerr*(p21 - p11) 
  h12 = p12 + zerr*(p22 - p12)
  h13 = p13 + zerr*(p23 - p13)

  call random_number(rand)
  zerr = 2*random_mag*rand + 0.5 - random_mag 
  h21 = p21 + zerr*(p31 - p21) 
  h22 = p22 + zerr*(p32 - p22)
  h23 = p23 + zerr*(p33 - p23)

  call random_number(rand)
  zerr = 2*random_mag*rand + 0.5 - random_mag 
  h31 = p31 + zerr*(p11 - p31) 
  h32 = p32 + zerr*(p12 - p32)
  h33 = p33 + zerr*(p13 - p33)

  a11 = (p11 + p21)/2 
  a12 = (p12 + p22)/2
  a13 = (p13 + p23)/2
  a21 = (p21 + p31)/2
  a22 = (p22 + p32)/2
  a23 = (p23 + p33)/2
  a31 = (p11 + p31)/2
  a32 = (p12 + p32)/2 
  a33 = (p13 + p33)/2  

  x11 = (a11 + a21)/2; x12 = (a12 + a22)/2; x13 = (a13 + a23)/2
  
  x21 = a31 - x11; x22 = a32 - x12; x23 = a33 - x13
  norm = sqrt(x21*x21 + x22*x22 + x23*x23)
  x21 = x21/norm; x22 = x22/norm; x23 = x23/norm

  fak1 = 1/2/sqrt(3.)*len; fak2 = sqrt(2./3.)*len

  call random_number(rand)
  zerr = (rand - 0.5)*random_mag*len/4
  h1 = x11 + fak1*x21 + fak2*fn1(j) + zerr
  call random_number(rand)
  zerr = (rand - 0.5)*random_mag*len/4
  h2 = x12 + fak1*x22 + fak2*fn2(j) + zerr
  call random_number(rand)
  zerr = (rand - 0.5)*random_mag*len/4
  h3 = x13 + fak1*x23 + fak2*fn3(j) + zerr
  
  sub_p1(n_sub + 1, 1) = h11
  sub_p2(n_sub + 1, 1) = h12
  sub_p3(n_sub + 1, 1) = h13
  
  sub_p1(n_sub + 1, 2) = p21
  sub_p2(n_sub + 1, 2) = p22
  sub_p3(n_sub + 1, 2) = p23
  
  sub_p1(n_sub + 1, 3) = h21
  sub_p2(n_sub + 1, 3) = h22
  sub_p3(n_sub + 1, 3) = h23

  
  sub_p1(n_sub + 2, 1) = h21
  sub_p2(n_sub + 2, 1) = h22
  sub_p3(n_sub + 2, 1) = h23
  
  sub_p1(n_sub + 2, 2) = p31
  sub_p2(n_sub + 2, 2) = p32
  sub_p3(n_sub + 2, 2) = p33
  
  sub_p1(n_sub + 2, 3) = h31
  sub_p2(n_sub + 2, 3) = h32
  sub_p3(n_sub + 2, 3) = h33
  
  
  sub_p1(n_sub + 3, 1) = h31
  sub_p2(n_sub + 3, 1) = h32
  sub_p3(n_sub + 3, 1) = h33
  
  sub_p1(n_sub + 3, 2) = p11
  sub_p2(n_sub + 3, 2) = p12
  sub_p3(n_sub + 3, 2) = p13
  
  sub_p1(n_sub + 3, 3) = h11
  sub_p2(n_sub + 3, 3) = h12
  sub_p3(n_sub + 3, 3) = h13
  
  
  sub_p1(n_sub + 4, 1) = h11
  sub_p2(n_sub + 4, 1) = h12
  sub_p3(n_sub + 4, 1) = h13
  
  sub_p1(n_sub + 4, 2) = h21
  sub_p2(n_sub + 4, 2) = h22
  sub_p3(n_sub + 4, 2) = h23
  
  sub_p1(n_sub + 4, 3) = h1
  sub_p2(n_sub + 4, 3) = h2
  sub_p3(n_sub + 4, 3) = h3
 
  
  sub_p1(n_sub + 5, 1) = h21
  sub_p2(n_sub + 5, 1) = h22
  sub_p3(n_sub + 5, 1) = h23
  
  sub_p1(n_sub + 5, 2) = h31
  sub_p2(n_sub + 5, 2) = h32
  sub_p3(n_sub + 5, 2) = h33
  
  sub_p1(n_sub + 5, 3) = h1
  sub_p2(n_sub + 5, 3) = h2
  sub_p3(n_sub + 5, 3) = h3
  
  
  sub_p1(n_sub + 6, 1) = h31
  sub_p2(n_sub + 6, 1) = h32
  sub_p3(n_sub + 6, 1) = h33
  
  sub_p1(n_sub + 6, 2) = h11
  sub_p2(n_sub + 6, 2) = h12
  sub_p3(n_sub + 6, 2) = h13
  
  sub_p1(n_sub + 6, 3) = h1
  sub_p2(n_sub + 6, 3) = h2
  sub_p3(n_sub + 6, 3) = h3
  
end subroutine MakeNextGeneration
  
  
subroutine ReplaceOldGeneration
integer k, l

  do k = 1, n_max
    do l = 1, 3
      main_p1(k,l) = sub_p1(k, l)
      main_p2(k,l) = sub_p2(k, l)
      main_p3(k,l) = sub_p3(k, l)
    end do
  end do
      
end subroutine ReplaceOldGeneration
 
subroutine PFN            ! Bestimmung der Flaechennormalen  
! Vorausgesetzt, die Eckpunke sind rechtssinnig angeordnet !!! 

real  avek1, avek2, avek3, bvek1, bvek2, bvek3, cvek1, cvek2, cvek3, laenge
integer  h, i

  avek1 = main_p1(j, 2) - main_p1(j, 1)
  avek2 = main_p2(j, 2) - main_p2(j, 1)
  avek3 = main_p3(j, 2) - main_p3(j, 1)
  bvek1 = main_p1(j, 3) - main_p1(j, 2)
  bvek2 = main_p2(j, 3) - main_p2(j, 2)
  bvek3 = main_p3(j, 3) - main_p3(j, 2)

  cvek1 = avek2*bvek3 - avek3*bvek2 ! c = a x b 
  cvek2 = avek3*bvek1 - avek1*bvek3
  cvek3 = avek1*bvek2 - avek2*bvek1

  laenge = Sqrt(cvek1*cvek1 + cvek2*cvek2 + cvek3*cvek3)

  fn1(j) = -cvek1/laenge
  fn2(j) = -cvek2/laenge
  fn3(j) = -cvek3/laenge

end subroutine PFN 
   
subroutine WriteFractal
integer k, l

  open(10, file = ofn)
  write(10,*) n_max
  do k = 1, n_max 
    write(10,*) 3
  end do
  do k = 1, n_max
    do l = 1, 3
      write(10,*) main_p1(k, l),' ', main_p2(k,l),' ', main_p3(k,l)
    end do
  end do
  write(10,*) random_arg
  close(10)
end subroutine WriteFractal

 
end program fractal3D
       
 


