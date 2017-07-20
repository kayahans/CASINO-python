!-------------------------------------------------------------------------!
! MAKE_P_STARS,  NDD,  October 2006.                                      !
! Utility for generating Jastrow p terms in which G vectors with the same !
! |G| have the same label.                                                !
!-------------------------------------------------------------------------!

MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)


CONTAINS


 SUBROUTINE invert(A,Ainv)
!-----------------------!
! Inverts 3x3 matrices. !
!-----------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: A(3,3)
  REAL(dp),INTENT(out) :: Ainv(3,3)
  REAL(dp) d
  d=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3)) &
   &+A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  if(d==0.d0)then
   write(6,*)'Trying to invert a singular determinant.' ; stop
  endif ! d=0
  d=1.d0/d
  Ainv(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))*d
  Ainv(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))*d
  Ainv(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))*d
  Ainv(2,1)=(A(3,1)*A(2,3)-A(2,1)*A(3,3))*d
  Ainv(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))*d
  Ainv(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))*d
  Ainv(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))*d
  Ainv(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))*d
  Ainv(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))*d
 END SUBROUTINE invert


 SUBROUTINE create_index(y,x_index)
!-----------------------------------------------------------------------------!
! This subroutine creates an index array x_index for the n items of data in   !
! the array y.  Adapted from Numerical Recipes.                               !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: y(:)
  INTEGER,INTENT(out) :: x_index(:)
  INTEGER,PARAMETER :: ins_sort_thresh=7,stacksize=80
  INTEGER n,i,x_indexj,ir,itemp,j,jstack,k,l,lp1,istack(stacksize)
  REAL(dp) yj
  n=size(x_index)
  do j=1,n
   x_index(j)=j
  enddo ! j
  if(n<=1)return
  jstack=0
  l=1
  ir=n
  do
   if(ir-l<ins_sort_thresh)then
    jloop: do j=l+1,ir
     x_indexj=x_index(j) ; yj=y(x_indexj)
     do i=j-1,l,-1
      if(y(x_index(i))<=yj)then
       x_index(i+1)=x_indexj
       cycle jloop
      endif ! y(x_index(i))<=yj
      x_index(i+1)=x_index(i)
     enddo ! i
     x_index(l)=x_indexj
    enddo jloop ! j
    if(jstack==0)return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
   else
    k=(l+ir)/2
    lp1=l+1
    itemp=x_index(k)    ; x_index(k)=x_index(lp1)  ; x_index(lp1)=itemp
    if(y(x_index(l))>y(x_index(ir)))then
     itemp=x_index(l)   ; x_index(l)=x_index(ir)   ; x_index(ir)=itemp
    endif
    if(y(x_index(lp1))>y(x_index(ir)))then
     itemp=x_index(lp1) ; x_index(lp1)=x_index(ir) ; x_index(ir)=itemp
    endif
    if(y(x_index(l))>y(x_index(lp1)))then
     itemp=x_index(l)   ; x_index(l)=x_index(lp1)  ; x_index(lp1)=itemp
    endif
    i=lp1
    j=ir
    x_indexj=x_index(lp1)
    yj=y(x_indexj)
    do
     do
      i=i+1
      if(y(x_index(i))>=yj)exit
     enddo ! i
     do
      j=j-1
      if(y(x_index(j))<=yj)exit
     enddo ! j
     if(j<i)exit
     itemp=x_index(i) ; x_index(i)=x_index(j) ; x_index(j)=itemp
    enddo
    x_index(lp1)=x_index(j)
    x_index(j)=x_indexj
    jstack=jstack+2
    if(jstack>stacksize)then
     write(6,*)'stacksize is too small.'
     stop
    endif ! jstack>stacksize
    if(ir-i+1>=j-l)then
     istack(jstack)=ir
     istack(jstack-1)=i
     ir=j-1
    else
     istack(jstack)=j-1
     istack(jstack-1)=l
     l=i
    endif ! ir-i+1>=j-l
   endif ! ir-l<ins_sort_thresh
  enddo
 END SUBROUTINE create_index


 SUBROUTINE wordwrap(text,unit_in,linelength_in)
!-------------------------------------------------------------------------!
! This subroutine prints out the contents of the character string 'text', !
! ensuring that line breaks only occur at space characters.  The output   !
! is written to unit unit_in if this parameter is supplied; otherwise the !
! output is written to unit o.  The maximum length of each line is given  !
! by linelength_in if this is supplied; otherwise the default line length !
! is 79 characters.                                                       !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in),OPTIONAL :: unit_in,linelength_in
  CHARACTER(*),INTENT(in) :: text
  CHARACTER(260) :: temp
  INTEGER i,unit,lentext,startpoint,stoppoint,lastpos,linelength
  if(present(unit_in))then
   unit=unit_in
  else
   unit=6
  endif ! unit supplied.
  lentext=len(trim(text))
  if(lentext<1)then
   write(unit,*)
   return
  endif ! No text
  if(present(linelength_in))then
   if(linelength_in>=2)then
    linelength=linelength_in
   else
    linelength=2
   endif ! sensible line-length supplied.
  else
   linelength=79
  endif ! linelength present.
  startpoint=1
  i=0
  do
   i=i+1
   stoppoint=startpoint+linelength-1
   if(stoppoint<=lentext)then
    lastpos=index(trim(text(startpoint:stoppoint))," ",.true.)
    if(lastpos>0)stoppoint=startpoint+lastpos-1
   else
    stoppoint=lentext
   endif ! stoppoint <= length of text
   if(i==1)then
! Allow the user to indent the first line, if (s)he wishes.
    temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
    write(unit,*)trim(temp)
   else
    temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
    write(unit,*)trim(adjustl(temp))
   endif ! i=1
   if(stoppoint==lentext)then
    exit
   else
    startpoint=stoppoint+1
   endif ! Finished text?
  enddo ! Loop over lines.
 END SUBROUTINE wordwrap


END MODULE utils


PROGRAM make_p_stars
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE utils
 IMPLICIT NONE
 INTEGER ierr,no_stars,ig,igp,i,j,k,istar,no_shells(3),no_G_tot,jmin,kmin, &
  &ngvec,dimensionality,ilayer,nlayer,label
 INTEGER,ALLOCATABLE :: ira(:),gvec_int(:,:)
 REAL(dp) latt_vect(3,3),rec_vect(3,3),gvec(3),old_G2,tol_G2
 REAL(dp),ALLOCATABLE :: gvec_sq(:)
 REAL(dp),PARAMETER :: tol=1.d-8
 CHARACTER(1) latt_type

 write(6,*)
 write(6,*)'MAKE_P_STARS'
 write(6,*)'============'
 write(6,*)

! Get user to enter the dimensionality.
 do
  write(6,*)'Please enter the dimensionality (3, 2 or 1).'
  read(5,*,iostat=ierr)dimensionality
  if(ierr/=0)dimensionality=-1
  if(dimensionality<1.or.dimensionality>3)then
   write(6,*)'Please try again.'
  else
   exit
  endif ! dimensionality other than 1-3.
 enddo

! Obtain lattice vectors (and number of layers for 2D systems).
 nlayer=1
 if(dimensionality==3)then
  do
   write(6,*)'Please enter the lattice type (F=FCC; S=SC; B=BCC; O=Other).'
   read(5,*,iostat=ierr)latt_type
   if(ierr/=0)latt_type=''
   if(latt_type=='f'.or.latt_type=='F')then
! FCC lattice
    latt_vect(1:3,1)=(/0.d0,1.d0,1.d0/)
    latt_vect(1:3,2)=(/1.d0,0.d0,1.d0/)
    latt_vect(1:3,3)=(/1.d0,1.d0,0.d0/)
    exit
   elseif(latt_type=='s'.or.latt_type=='S')then
! SC lattice
    latt_vect(1:3,1)=(/1.d0,0.d0,0.d0/)
    latt_vect(1:3,2)=(/0.d0,1.d0,0.d0/)
    latt_vect(1:3,3)=(/0.d0,0.d0,1.d0/)
    exit
   elseif(latt_type=='b'.or.latt_type=='B')then
! BCC lattice
    latt_vect(1:3,1)=(/-1.d0,1.d0,1.d0/)
    latt_vect(1:3,2)=(/1.d0,-1.d0,1.d0/)
    latt_vect(1:3,3)=(/1.d0,1.d0,-1.d0/)
    exit
   elseif(latt_type=='o'.or.latt_type=='O')then
! General lattice.  Get user to enter lattice vectors.
    do
     write(6,*)'Please enter the first lattice vector (a_x,a_y,a_z):'
     read(5,*,iostat=ierr)latt_vect(1:3,1)
     if(ierr/=0)then
      write(6,*)'Please try again (and get it right this time).'
     else
      exit
     endif ! ierr/=0
    enddo ! Get first lattice vector
    do
     write(6,*)'Please enter the second lattice vector (a_x,a_y,a_z):'
     read(5,*,iostat=ierr)latt_vect(1:3,2)
     if(ierr/=0)then
      write(6,*)'Please try again (and get it right this time).'
     else
      exit
     endif ! ierr/=0
    enddo ! Get second lattice vector
    do
     write(6,*)'Please enter the third lattice vector (a_x,a_y,a_z):'
     read(5,*,iostat=ierr)latt_vect(1:3,3)
     if(ierr/=0)then
      write(6,*)'Please try again (and get it right this time).'
     else
      exit
     endif ! ierr/=0
    enddo ! Get third lattice vector
    exit
   else
    write(6,*)'Please try again'
   endif ! latt_type
  enddo ! Get lattice type
 elseif(dimensionality==2)then
  do
   write(6,*)'Please enter the lattice type (S=Square; H=Hexagonal; O=Other).'
   read(5,*,iostat=ierr)latt_type
   if(ierr/=0)latt_type=''
   if(latt_type=='s'.or.latt_type=='S')then
! Square lattice
    latt_vect(1:2,1)=(/1.d0,0.d0/)
    latt_vect(1:2,2)=(/0.d0,1.d0/)
    exit
   elseif(latt_type=='h'.or.latt_type=='H')then
! Hexagonal lattice
    latt_vect(1:2,1)=(/1.d0,0.d0/)
    latt_vect(1:2,2)=(/0.5d0,0.5d0*sqrt(3.d0)/)
    exit
   elseif(latt_type=='o'.or.latt_type=='O')then
! General lattice.  Get user to enter lattice vectors.
    do
     write(6,*)'Please enter the first lattice vector (a_x,a_y):'
     read(5,*,iostat=ierr)latt_vect(1:2,1)
     if(ierr/=0)then
      write(6,*)'Please try again (and get it right this time).'
     else
      exit
     endif ! ierr/=0
    enddo ! Get first lattice vector
    do
     write(6,*)'Please enter the second lattice vector (a_x,a_y):'
     read(5,*,iostat=ierr)latt_vect(1:2,2)
     if(ierr/=0)then
      write(6,*)'Please try again (and get it right this time).'
     else
      exit
     endif ! ierr/=0
    enddo ! Get second lattice vector
    exit
   else
    write(6,*)'Please try again'
   endif ! latt_type
  enddo ! Get lattice type
  latt_vect(3,1)=0.d0 ; latt_vect(3,2)=0.d0
  latt_vect(1:2,3)=0.d0 ; latt_vect(3,3)=1.d0
  do
   call wordwrap('Please enter the number of layers if you are looking at a &
    &real system with atoms, e.g., enter "2" for bilayer graphene.  Otherwise, &
    &if you have a strictly 2D model system (even one with multiple layers), &
    &please choose "1".')
   read(5,*,iostat=ierr)nlayer
   if(ierr/=0)nlayer=-1
   if(nlayer<1)then
    write(6,*)'Please enter a positive integer.'
   else
    exit
   endif ! nlayer<1
  enddo
 else ! 1D system
  latt_vect=0.d0
  latt_vect(1,1)=1.d0 ; latt_vect(2,2)=1.d0 ; latt_vect(3,3)=1.d0
 endif

! Get user to enter number of stars of G vectors to include in p.
 do
  write(6,*)'Please enter the number of stars of G vectors to include in p.'
  read(5,*,iostat=ierr)no_stars
  if(ierr/=0)no_stars=-1
  if(no_stars<=1)then
   write(6,*)'Need two or more stars of G vectors.'
  else
   exit
  endif ! no_stars<=0
 enddo
 write(6,*)

! Get reciprocal lattice vectors.
 call invert(latt_vect,rec_vect)
 rec_vect=transpose(rec_vect)

! Size of G vector grid to study.  Grid is -no_shells(i) .. no_shells(i) in
! each direction i.
 no_shells=0
 no_shells(1:dimensionality)=no_stars-1

! Number of G vectors.  NB, only one out of each (+/-)G is used.
 no_G_tot=(2*no_shells(1)+1)*(2*no_shells(2)+1)*(2*no_shells(3)+1)/2+1
 allocate(gvec_int(3,no_G_tot),gvec_sq(no_G_tot),ira(no_G_tot),stat=ierr)
 if(ierr/=0)then
  write(6,*)'Allocation error.' ; stop
 endif ! ierr/=0

! Loop over G vectors in grid, storing integer coordinates & |G|^2.
 ig=0
 do i=0,no_shells(1)
  if(i==0)then
   jmin=0
  else
   jmin=-no_shells(2)
  endif ! i=0
  do j=jmin,no_shells(2)
   if(i==0.and.j==0)then
    kmin=0
   else
    kmin=-no_shells(3)
   endif ! i=j=0
   do k=kmin,no_shells(3)
    ig=ig+1
    gvec(1:3)=dble(i)*rec_vect(1:3,1)+dble(j)*rec_vect(1:3,2) &
     &+dble(k)*rec_vect(1:3,3)
    gvec_sq(ig)=dot_product(gvec,gvec)
    gvec_int(1,ig)=i ; gvec_int(2,ig)=j ; gvec_int(3,ig)=k
   enddo ! k
  enddo ! j
 enddo ! i

! Create array for indexing G vectors in order of ascending |G|^2.
 call create_index(gvec_sq,ira)

! Tolerance for comparing G^2.
 if(dimensionality==3)then
  tol_G2=tol*max(dot_product(rec_vect(:,1),rec_vect(:,1)), &
   &dot_product(rec_vect(:,2),rec_vect(:,2)), &
   &dot_product(rec_vect(:,3),rec_vect(:,3)))
 elseif(dimensionality==2)then
  tol_G2=tol*max(dot_product(rec_vect(:,1),rec_vect(:,1)), &
   &dot_product(rec_vect(:,2),rec_vect(:,2)))
 else ! 1D system
  tol_G2=tol*abs(rec_vect(1,1))
 endif ! dimensionality

! Count number of G vectors to be written out.
 istar=1
 old_G2=-100.d0
 ngvec=0
 do ig=2,no_G_tot
  igp=ira(ig)
  if(gvec_sq(igp)-old_G2>tol_G2)then
   istar=istar+1
   if(istar>no_stars)then
    ngvec=ig-2 ! Number of G vectors to write out.
    exit
   endif ! Time to exit
   old_G2=gvec_sq(igp)
  endif ! new star
 enddo ! ig
 if(ngvec==0)ngvec=no_G_tot-1
! Increase the number of G vectors for subsequent layers.
 ngvec=ngvec+(nlayer-1)*(2*ngvec+1)

 write(6,*)'********** START OF OUTPUT **********'
 write(6,*)'START P TERM'
 write(6,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud)'
 write(6,*)0
 write(6,*)'Number of simulation-cell G-vectors (NB, cannot have both G & -G)'
 write(6,*)ngvec
 write(6,*)'G-vector (in terms of rec latt vects) ; label'
 label=0
 do ilayer=1,nlayer
  if(ilayer>1)then
   label=label+1
   write(6,'("   ",3(i5," "),"   ",i4)')0,0,ilayer-1,label
  endif ! ilayer>1
  istar=1
  old_G2=-100.d0
  do ig=2,no_G_tot
   igp=ira(ig)
   if(gvec_sq(igp)-old_G2>tol_G2)then
    istar=istar+1
    if(istar>no_stars)exit
    label=label+1
    old_G2=gvec_sq(igp)
   endif ! new star
   if(ilayer==1)then
    write(6,'("   ",3(i5," "),"   ",i4)')gvec_int(1:3,igp),label
   else
    write(6,'("   ",3(i5," "),"   ",i4)')gvec_int(1:2,igp),ilayer-1,label
    write(6,'("   ",3(i5," "),"   ",i4)')-gvec_int(1:2,igp),ilayer-1,label
   endif ! ilayer=1
  enddo ! ig
 enddo ! ilayer
 write(6,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
 write(6,*)'END P TERM'
 write(6,*)'*********** END OF OUTPUT ***********'

 deallocate(gvec_int,gvec_sq,ira)

 write(6,*)
 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM make_p_stars
