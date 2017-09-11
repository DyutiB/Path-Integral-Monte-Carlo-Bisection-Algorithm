program pimc_bisection
include 'param.dat'
!! this is the main program for bisection algorithm. All the other subroutines are included with this file. This file will 
!! generate the stable/equilibrium positions and positional and orientational correlations along with lindemann ratio. 
!! No of interacting particles and quantum fluctions are already included in the param.dat file. The quantum Boltzmann 
!! particles are interacting via Coulomb potential, changing it to other potential is easily doable.   
!! the interacting particles are confined by irregular potential (the parameters are in Votparam.dat)
!! Changing the potential is easily doable. 

integer :: i,j,ij,ii,t,ipt,jj,indx,n
integer :: bulk_counter, imtcle
integer :: iptcle,kbead
double precision :: accept
double precision :: Centbin(Ndim,Nbin2),Gvalbin(Nbin2), xmax,xmin, &
& ymax,ymin,incx,incy,poscm(Ndim,Npt), cmaccept 
double precision :: dist1p(Npt,Npt), cbondbin(Nbin), Ampbin(Nbin), & 
& AmpGbr(Nbin),Amp(Npt),Ph(Npt),inc,div,dia,grbin(Nbin)
double precision :: avgr(Npt,tbead), avgr2(Npt,tbead),flucr,numb
double complex :: psi,mpsi(Npt)
 character(len=20) :: outname,denstrng,boostrng,bocfstrng 

dia=0.d0; inc=0.d0; div=0.d0
!	 dia = diameter ; inc = binning increment for gofr,bocf,
!	 div= binning increment for boo

 write(outname,30) nrsb ! nrsb = quantum parameter (n in paper)
 30 format(f6.4)
 denstrng = 'densityN57n'//trim(outname)
 boostrng = 'boon'//trim(outname)
 bocfstrng = 'ampgbrn'//trim(outname)

open(9,file=denstrng,status='unknown')                  !! file to write stable/equilibrium positions
open(10,file=boostrng,status='unknown')                 !! file to write the Bond-Orientational order after binning
open(11,file=bocfstrng,status='unknown')                !! file to write the orientational correlation after binning. 
open(12,file='nrsb_lindamann',status='unknown',position='append')
                                                        !! file to write Lindemann ratio with changing 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initial position of Npt particles and their center of mass(cm) and the distances. In this initial position, position
!!! of each beads are given the same position as the position of center of mass.  
call Initpos(pos_old,cmpos_old)
call distances()
!! calculates Npair distance between ptcls both bead-wise and cm r_pp(tbead,npair), cmr_pp(Npair)
!! pos_old(Ndim,tbead,Npt) = updated positions on total n of beads	
!! po(Ndim,nb,Npt) = saved positions of `nb' bead before multi-slice (nb-(time)slice) beads movement 
!! pn(Ndim,nb,Npt) = updated positions of nb bead after bead movement	 
!!         nb = total no of beads reorganised by multi-slice move	
 	
!-------------density calculation------------------------------
!---bining for position density calculation below--------------
xmax = 2.50d0; xmin = -2.5d0
ymax = 2.5d0; ymin = -2.5d0
incx = (xmax-xmin)/dfloat(Nbin)
incy = (ymax-ymin)/dfloat(Nbin)

ij = 0
do i=1,Nbin
do j=1,Nbin
  ij = ij+1

  centbin(1,ij) = xmin + (dfloat(i)-0.5d0)*incx
  centbin(2,ij) = ymin + (dfloat(j)-0.5d0)*incy
  Gvalbin(ij) = 0.d0

enddo
enddo 
!-------------------------------------------------------------
 accept = 0.d0  ! acceptance ratio of bead movement (bisection/staging)
 cmaccept = 0.d0  ! acceptance ratio of center of mass movement
 indx = 0

!!!!!! In the final iteration, where both the beads and center of mass of the beads are to be moved together. 
!!!!!! Before the final iteration, we do some iteration over the beads positions as initially every beads are allowed 
!!!!!! to have the same initial positions. 
!---------------------------------------------------------------
!-------before equilibrium run----------------------------------
do ii=1,nblkeq*nstep
  call randombead(kbead)

   do iptcle = 1,Npt
    call pimc_metro(iptcle,kbead,accept)  !bisection MC move
   enddo

enddo 
!----------------------equilibrium run------------------------
!-------------------------------------------------------------
accept = 0.d0 

do ii=1,nstep*const
!  print*,ii
 call cmmetro(cmaccept) ! cm movement
!	(displacement move -- without altering
!	"trajectories" or "polymers")

  do jj=1,nblkeq

   call randombead(kbead)
!	selects kbead to (kbead+nb) timeslices to reorganise by bisection

   do iptcle = 1,Npt
    call pimc_metro(iptcle,kbead,accept)  !bisection MC move
   enddo

 enddo
enddo
!-------------------------------------------------------------
print*,'equlibrium done'

 accept = 0.d0
 cmaccept = 0.d0
 indx = 0

 avgr = 0.d0; avgr = 0.d0  ! re-initialize after equilibriation

 call diameter(dia)
 
 dia = dia
 inc = (dia*1.05d0)/dfloat(Nbin)
 div = 1.5d0/dfloat(Nbin)

 do n=1,Nbin
  grbin(n)=(dfloat(n)-0.5d0)*inc 
  cbondbin(n)=(dfloat(n)-0.5d0)*div
  AmpGbr(n)=0.d0
  Ampbin(n)=0.d0
 enddo 
 
do ii=1,const*nstep  ! const is an integer factor larger than 1
!	ensures that that No. of MC passes are cont times eqlb passes
 
 call cmmetro(cmaccept)  ! cm move of the "ring polymers"

 do jj=1,nblk
 
  indx = indx + 1 

   call randombead(kbead) ! usage similar to what was done during eqlb run

   do iptcle = 1,Npt
    call pimc_metro(iptcle,kbead,accept)
   enddo
!	Note that the bisection move in one MC pass must occur over all the
!	Nptcls. Here theloop over Npt is kept outside, could be made a
!	part of bisection.

enddo ! do jj=1,nblk ends here

!! if(mod(indx,40)==0) then
  call dsG(incx,incy,tbead,pos_old,centbin,Gvalbin)
!! endif

  call dist1p_wothers(dist1p) ! dist1p(Npt,Npt) array calcualtion
!	calculates the dist of iptcle with all the Nptcles (dist with
!	itself is trivially zero), and then take a loop over iptcle from
!	1 to Npt (ideally there will be Npair entries,
!	but here there will be additional Npt number of zeros)

   bulk_counter = 0 ! Checks the no. of ptcles in bulk
!	by removing those on the boundary

 if(wrtbnd) then ! flagtocalculate bond-orintational corr fn

  do iptcle=1,Npt ! Similar to classical MC
     psi=(0.d0,0.d0)
     call nbr_angle(iptcle,dist1p(iptcle,:),bulk_counter,psi)
     mpsi(iptcle)=psi
     Amp(iptcle)=cdabs(psi)
     Ph(iptcle)=datan2(dimag(psi),dble(psi))
  enddo     
 
  call pairG(div,Npt,Amp,Nbin,cbondbin,Ampbin,bulk_counter)
    
  do imtcle=1,Npt
    call pairGBr(imtcle,inc,grbin,AmpGbr,mpsi,dist1p) 
  enddo

 endif ! if(wrtbnd) ends here

   call lindemann(avgr,avgr2)


enddo ! do ii=1,nstep ends here

print*, accept/dfloat(const*nstep*nblk)
print*, cmaccept/dfloat(const*nstep)

 ij = 0
 do i=1,Nbin
  write (9,*)
  do j=1,Nbin
   ij = ij + 1
   write(9,*) centbin(1,ij), centbin(2,ij), Gvalbin(ij)/dfloat(const*nstep)
  enddo 
 enddo

do n=1,Nbin
 write(10,*) cbondbin(n), Ampbin(n)/dfloat(const*nstep)
enddo 

do n=1,Nbin
 write(11,*) grbin(n), Ampgbr(n)/dfloat(const*nstep)
enddo 

flucr = 0.d0
do n=1,Npt
 do t=1,tbead
    numb = (avgr2(n,t)/dfloat(const*nstep)) - &
 &   (avgr(n,t)/dfloat(const*nstep))**2

   flucr = flucr + numb 
 enddo 
enddo 
 write(12,*) nrsb,flucr/dfloat(Npt)

stop
end program pimc_bisection
!--------------------------------------------------------------------
subroutine pimc_metro(iptcle,kbead,accept)
include 'param.dat'
integer :: i,j,k,kbead,t,kk
integer :: iptcle,lev,ilev,irg,irg2,nin,lac,lin
double precision :: dirg2,alpha,potold,accept,eta1,eta2
double precision :: Prob_old,Prob2,AV(nlev,Ndim),Uprt(nlev), &
& du,potnew,T_acc,rand
double precision, parameter :: two = 2
logical :: lev_run !flag to know if bisection run in progress

lev_run = .true.
Prob_old = 1.d0  !initial prob. set to 1, it will be updated as and when the move at each level is accepted

do lev = nlev,2,-1
if(lev_run.eqv.(.true.)) then
     ilev = lev - 1   !level number
     irg = 2**(ilev)  !distance between the beads being moved in level ilev
     irg2 = 2**(ilev-1)  !the starting bead to move counting from r_0
     nin = 2**(nlev-ilev-1)  !no of bead moving in level ilev
     dirg2 = (2.d0)**(ilev-2)  !related to irg2 but made double precision,
!	as it enters in the factor multiplying gaussian random number
     alpha = dsqrt(dble(dirg2)*tau*lambda) !factor multiplying eta

 call poten1b(iptcle,lev,potold)

       lac = irg2 + 1  !setting r_0 = r_1 (convention)
!print*,lev,irg,irg2,nin,lac
 do lin=1,nin  !index to move the appropriate beads at a given level
  call boxmuller(eta1,eta2)
   pn(1,lac,iptcle) = 0.5d0*(pn(1,lac+irg2,iptcle)+pn(1,lac-irg2,iptcle)) &
         &  + eta1*alpha
   pn(2,lac,iptcle) = 0.5d0*(pn(2,lac+irg2,iptcle)+pn(2,lac-irg2,iptcle)) &
         &  + eta2*alpha
!	implementing Eq (5.51)

	   lac = lac + irg
!	 finding next bead to move in the same level of bisection
 enddo ! do lin=1,nin ends here

 call poten1b(iptcle,lev,potnew)
   du = potnew - potold

 call return_exp(du,Prob2)
    T_acc = Prob2/Prob_old

 call random_number(rand)

!    print*,T_acc,rand

 if( T_acc.gt.rand) then
    Prob_old = Prob2
    lev_run = .true.

 !   print*,lev
 else
    lev_run = .false.

    do t=1,nb
     pn(1,t,iptcle) = po(1,t,iptcle)
     pn(2,t,iptcle) = po(2,t,iptcle)
    enddo 

 endif ! if(T_acc...) ends here

if((lev.eq.two).and.(lev_run.eqv.(.true.))) then

   call savepos_update(iptcle,kbead)
   call cmdist_update(iptcle)
   call beaddist_update(iptcle,kbead)

    accept = accept + 1.d0/dfloat(Npt)
endif ! if((lev.eq.two)...) ends here


endif ! if(lev_run...ends here)
enddo ! do lev=nlev,2,-1 ends here

return
end subroutine pimc_metro
!--------------------------------------------------------------------
subroutine beaddist_update(iptcle,kbead)
include 'param.dat'
integer :: iptcle,kbead,i,jj,ij,j,k,k0,l
double precision :: dx,dy,dist_ij

do jj=1,Npt
 if( jj /= iptcle) then

  if(jj < iptcle) then
    i = iptcle
    j = jj
  else
    i = jj
    j = iptcle
  endif 

   ij = j+((i-1)*(i-2))/2

   k0 = kbead
    l = 1
  do k=1,nb 
    if((k0+k-1).gt.tbead) then
       k0 = 1
        l = k
    endif 

   dx = pn(1,k,i) - pn(1,k,j)
   dy = pn(2,k,i) - pn(2,k,j)
   dist_ij = dx*dx + dy*dy
   dist_ij = sqrt(dist_ij) 

   r_pp(k0+k-l,ij) = dist_ij

  enddo ! do k=1,nb ends here

 endif 
enddo 

return
end subroutine beaddist_update
!--------------------------------------------------------------------
subroutine cmdist_update(iptcle)
include 'param.dat'
integer :: iptcle,i,j,ij,jj,k

do jj=1,Npt
if(jj /= iptcle) then
 
 if( jj < iptcle) then
   i=iptcle
   j=jj
 else
   i=jj
   j=iptcle
 endif 
   ij=j+((i-1)*(i-2))/2

 do k=1,Ndim
   cmrdist(k,ij)=cmpos_new(k,i)-cmpos_new(k,j)
   cmrdist_save(k,ij) = cmrdist(k,ij)
 enddo 

 cmr_pp(ij)=0.d0
 do k=1,Ndim
 cmr_pp(ij) = cmr_pp(ij)+cmrdist(k,ij)*cmrdist(k,ij)
 enddo
 cmr_pp(ij) = dsqrt(cmr_pp(ij))
 cmr_ppsave(ij) = cmr_pp(ij)

endif
enddo

return
end subroutine cmdist_update
!--------------------------------------------------------------------
subroutine savepos_update(iptcle,kbead)
include 'param.dat'
integer :: i,j,k,k0,kbead,iptcle,l,t,kk

 k0 = kbead
 l = 1
do k=1,nb
 if ((k0+k-1).gt.tbead) then
    k0 = 1
     l = k
 endif

  pos_old(1,k0+k-l,iptcle) = pn(1,k,iptcle)
  pos_old(2,k0+k-l,iptcle) = pn(2,k,iptcle)

enddo 

do k=1,Ndim
  cmpos_old(k,iptcle) = 0.d0
  cmpos_new(k,iptcle) = 0.d0
enddo

do t=1,tbead
 do k=1,Ndim
   cmpos_old(k,iptcle) = cmpos_old(k,iptcle) + pos_old(k,t,iptcle)/dfloat(tbead)
   cmpos_new(k,iptcle) = cmpos_new(k,iptcle) + pos_old(k,t,iptcle)/dfloat(tbead)
 enddo 
enddo

return
end subroutine savepos_update
!--------------------------------------------------------------------
subroutine return_exp(du,Prob2)
include 'param.dat'
double precision :: du,Prob2

if ( du.lt.(-10.0)) then
   Prob2 = 22026.460d0 
!	 if du is very small, must move, hence choose large prob
else
   Prob2 = Exp(-du)
endif

return
end subroutine return_exp
!---------------------------------------------------------------------
subroutine randombead(kbead)
include 'param.dat'

integer :: i,j,k,k0,l,kbead,iptcle,kk,ij
double precision :: rann

kbead = 0
 call Prepseed()
 call random_number(rann)
 kbead = int(rann*tbead+1.d0)

! write(11,*) 
! write(11,*)'#', kbead
! write(12,*) 
! write(12,*) kbead

do iptcle = 1,Npt
  k0 = kbead
   l = 1
!  write(11,*) '#particle no=',iptcle
!  write(12,*) '#particle no=',iptcle

 do k=1,nb
    if((k0+k-1).gt.tbead) then
       k0 = 1
        l = k
    endif

   do kk=1,Ndim
     pn(kk,k,iptcle) = pos_old(kk,k0+k-l,iptcle)
     po(kk,k,iptcle) = pos_old(kk,k0+k-l,iptcle)   
   enddo
!  write(2,*) pn(1,k,iptcle), pn(2,k,iptcle)
!   write(11,*) pn(1,k,iptcle),pn(2,k,iptcle)
!   write(12,*) po(1,k,iptcle),po(2,k,iptcle)

  enddo  ! do k=1,nb ends here
enddo ! do iptcle ends here


return
end subroutine randombead
!---------------------------------------------------------------------
subroutine boxmuller(eta1,eta2)
include 'param.dat'
! box-muller random number generator..used to move the bead by a random no
! whose mean is zero and standard deviation 1
! always gives two random no in single calling

double precision :: r1,r2,term,theta,r,eta1,eta2
double precision,parameter :: pi = 4.d0*datan(1.d0)

 call Prepseed()
 call random_number(r1)
 call random_number(r2)

 term = -2.d0*dlog(1.d0-r1)
 r = dsqrt(term)
 theta = 2.d0*pi*r2

 eta1 = r*dcos(theta)
 eta2 = r*dsin(theta)

return
end subroutine boxmuller
!---------------------------------------------------------------------
subroutine poten1b(iptcle,lev,poten)
include 'param.dat'
integer :: i,jj
integer :: iptcle,lev,ilev,irg,irg2,lin,nin,lac
double precision :: poten, irregpot, alpha, dirg2, &
   & posn(Ndim), dx, dy, dist_ij

   ilev = lev - 1
   irg = 2**(ilev)
   irg2 = 2**(ilev-1)  
   nin = 2**(nlev-ilev-1)
   dirg2 = (2.d0)**(ilev-1)
   alpha = dirg2*tau

   poten = 0.d0
   do jj=1,Npt
   if( jj /= iptcle) then
      lac = irg2 + 1

    do lin=1,nin

	dx = pn(1,lac,iptcle) - pn(1,lac,jj)
        dy = pn(2,lac,iptcle) - pn(2,lac,jj)
        dist_ij = dx*dx + dy*dy
        dist_ij = dsqrt(dist_ij)

        poten = poten + 1.d0/dist_ij
 
          lac = lac + irg
 
    enddo  ! do lin=1,nin ends here

   endif  ! if(..) ends here
   enddo  ! do jj...ends here

    lac = irg2 + 1
    do i=1,nin
    posn(1) = pn(1,lac,iptcle); posn(2) = pn(2,lac,iptcle)
    poten = poten + irregpot(posn)

          lac = lac + irg 
    enddo 

        poten = alpha*poten 

return
end subroutine poten1b
!--------------------------------------------------------------------
subroutine Initpos(pos0,pos)
include 'param.dat'
integer :: i,t,n,j,indx
double precision :: pos(Ndim,Npt),pos0(Ndim,tbead,Npt), &
&  posn0(Ndim,Npt)
double precision :: site_r,site_theta,randR(Npt),randTh(Npt),pi
 character(len=20) :: outnrsb, outname
 character(len=50) :: outputstrng

!write(outnrsb,10) nrsb
!10 format(f5.1)
pi=4.d0*datan(1.d0)

! call PrepSeed()            !prepares seed for random number generator
 call random_number(randR)  !for random radial position of particles
 call random_number(randTh) !for random angular position of particles

indx=0           !keeps track of how many particles already placed
do i=1,Npt
     indx=indx+1
                           !Rad_ini is the scale of max radius of particles
101  site_r=Rad_ini*dsqrt(randR(i))  !randomly choosen radial position
     site_theta=2.d0*pi*randTh(i)    !randomly chosen angular position

     posn0(1,i)=site_r*dcos(site_theta)   !x-position
     posn0(2,i)=site_r*dsin(site_theta)   !y-position

	 do j=i-1,1, -1  !if posn of any 2 ptcles are same, within 10^{-10} 
	   if ( (dabs(posn0(1,i)-posn0(1,j)).lt.(1.0e-10)).and.&
     &      (dabs(posn0(2,i)-posn0(2,j)).lt.(1.0e-10)) ) then
         goto 101
       endif
     enddo
enddo  ! do i=1,Npt ends here

do i=1,Npt
  pos(1,i) = posn0(1,i)     ! initiating the position of cm 
  pos(2,i) = posn0(2,i) 

  do n = 1,tbead    ! all beads are given same initial position as cm
   pos0(1,n,i) = posn0(1,i)   
   pos0(2,n,i) = posn0(2,i)
  enddo ! do n=1,tbead ends here
enddo ! do i=1,Npt ends here

!endif 

return
end subroutine Initpos
!-----------------------------------------------------------------
subroutine distances()
include 'param.dat'

integer :: ij,i,j,n
double precision :: dx,dy,dist

ij = 0
do i=2,Npt
 do j=1,i-1

  ij = ij + 1
  dx = cmpos_old(1,i) - cmpos_old(1,j)
  dy = cmpos_old(2,i) - cmpos_old(2,j)
   dist = dx*dx + dy*dy
   dist = dsqrt(dist)
 
    cmr_pp(ij) = dist
    
   do n=1,tbead
    r_pp(n,ij) = dist
   enddo

 enddo  ! do j=1,i-1 ends here
enddo  ! do i=2,Npt ends here

return
end subroutine distances
!-----------------------------------------------------------------
subroutine PrepSeed()
implicit none

integer :: i,sdsize
integer :: clockno
integer, allocatable :: seed(:)
real(8), allocatable :: vals(:)

 call random_seed(sdsize)
allocate (seed(sdsize))
allocate (vals(sdsize))

 call random_number(vals)

 call system_clock(count=clockno)

do i=1,sdsize
     seed(i)=int(5000.d0*vals(i))+clockno
enddo
 call random_seed(put=seed)

end subroutine PrepSeed
!------------------------------------------------------------------
function irregpot(x)
include 'votparam.dat'
!implicit none

double precision :: x(2), irregpot

irregpot = x(1)*x(1)*x(1)*x(1)/b + b*x(2)*x(2)*x(2)*x(2)  &
  &  -2.d0*lam*x(1)*x(1)*x(2)*x(2) + gam*( x(1)-x(2) )  &
  &  *x(1)*x(2)*dsqrt( x(1)*x(1) + x(2)*x(2) )
irregpot = a*irregpot

!irregpot = x(1)*x(1) + x(2)*x(2)

return
end function irregpot
!-------------------------------------------------------------------
subroutine cmmetro(cmaccept)
include 'param.dat'

integer :: ipt,k
double precision :: poten1po, poten1pn, dx(2), rann, &
 & rand, cmaccept, p
double precision :: half=0.5d0, one=1.d0

do ipt = 1,Npt

 do k=1,Ndim
   dx(k)=0.d0
 enddo 

 call poten_1cminit(ipt,poten1po)

 do k=1,Ndim
 call random_number(rand)
 dx(k)=2.d0*(rand-half)*delta
 cmpos_new(k,ipt) = cmpos_old(k,ipt) + dx(k)
 enddo 

 call poten_1cmmove(ipt,dx,poten1pn)

 if( poten1pn < poten1po ) then
      p = 1.d0
  else
      p = dexp(-(poten1pn-poten1po)/(redt*tbead))
 endif  

 call random_number(rann) 

 if ( p > rann ) then

  call cmupdate(ipt)
  call beadupdate(ipt)

   cmaccept = cmaccept + one/dfloat(Npt)
   cmacceptpar(ipt) = cmacceptpar(ipt) + one

 else

   call posrestore(ipt)

 endif ! if( p > rann) ends here
enddo ! do ipt=1,Npt ends here

return
end subroutine cmmetro
!-------------------------------------------------------------------
subroutine beadupdate(iptcle)
include 'param.dat'
integer :: iptcle,t,jj,i,j,ij

do jj=1,Npt
 if ( jj /= iptcle) then

  if(jj < iptcle) then
   i=iptcle
   j=jj
  else
   i=jj
   j=iptcle
 endif
 ij=j+((i-1)*(i-2))/2

 do t=1,tbead
  r_pp(t,ij) = r_pptrial(t,jj)
 enddo

 endif 
enddo 

return
end subroutine beadupdate
!-------------------------------------------------------------------
subroutine cmupdate(iptcle)
include 'param.dat'

integer :: iptcle,k,jj,i,j,ij

do k=1,Ndim
  cmpos_old(k,iptcle) = cmpos_new(k,iptcle)
enddo

do jj=1,Npt
 if( jj /= iptcle) then

   if( jj < iptcle) then
       i = iptcle
       j = jj
   else
       i = jj
       j = iptcle
   endif
 
   ij = j + ((i-1)*(i-2))/2

   do k=1,Ndim
    cmrdist(k,ij) = cmpos_new(k,i)-cmpos_new(k,j)
   enddo 

   cmr_pp(ij) = 0.d0
   do k=1,Ndim
   cmr_pp(ij) = cmr_pp(ij) + cmrdist(k,ij)*cmrdist(k,ij)
   enddo

   cmr_pp(ij) = dsqrt(cmr_pp(ij))

 endif 
enddo ! do jj=1,Npt ends here

return
end subroutine cmupdate
!-------------------------------------------------------------------
subroutine posrestore(iptcle)
include 'param.dat'

integer :: t,iptcle

 cmpos_new(1,iptcle) = cmpos_old(1,iptcle)
 cmpos_new(2,iptcle) = cmpos_old(2,iptcle)

do t=1,tbead
  pos_old(1,t,iptcle) = pos_save(1,t,iptcle)
  pos_old(2,t,iptcle) = pos_save(2,t,iptcle)
enddo 

return
end subroutine posrestore
!-------------------------------------------------------------------
subroutine poten_1cminit(iptcle,poten)
include 'param.dat'

integer :: i,j,ij,jj,t,iptcle
double precision :: poten, posn(2), irregpot
double precision, parameter :: one=1.d0

poten = 0.d0
do jj=1,Npt
 if( jj /= iptcle) then

  if( jj < iptcle) then
     i = iptcle
     j = jj
  else
     i = jj
     j = iptcle
  endif 
     ij = j+((i-1)*(i-2))/2

  do t=1,tbead
   poten = poten + one/r_pp(t,ij)
  enddo 

 endif 
enddo ! do jj=1,Npt ends here

do t=1,tbead
 posn(1) = pos_old(1,t,iptcle); posn(2) = pos_old(2,t,iptcle)
 poten = poten + irregpot(posn)
enddo

return
end subroutine poten_1cminit
!------------------------------------------------------------------
subroutine poten_1cmmove(iptcle,dx,poten)
include 'param.dat'

integer :: n,iptcle,jj,t
double precision :: dx(2),rdx,rdy,poten,posn(2),irregpot
double precision, parameter :: one = 1.d0

do n=1,tbead
  pos_save(1,n,iptcle) = pos_old(1,n,iptcle)
  pos_save(2,n,iptcle) = pos_old(2,n,iptcle)

  pos_old(1,n,iptcle) = pos_old(1,n,iptcle) + dx(1)
  pos_old(2,n,iptcle) = pos_old(2,n,iptcle) + dx(2)
enddo 

poten=0.d0
do jj=1,Npt
if ( jj /= iptcle) then

  do t=1,tbead
   rdx = pos_old(1,t,iptcle)-pos_old(1,t,jj)
   rdy = pos_old(2,t,iptcle)-pos_old(2,t,jj)
     r_pptrial(t,jj) = rdx*rdx + rdy*rdy
     r_pptrial(t,jj) = dsqrt(r_pptrial(t,jj))
     poten = poten + one/r_pptrial(t,jj)
  enddo 

endif 
enddo

do t=1,tbead
  posn(1) = pos_old(1,t,iptcle); posn(2) = pos_old(2,t,iptcle)
  poten = poten + irregpot(posn)
enddo

return
end subroutine poten_1cmmove
!------------------------------------------------------------------
subroutine lindemann(avgr,avgr2)
include 'param.dat'
integer :: n,t
double precision :: avgr(Npt,tbead),avgr2(Npt,tbead)

do n=1,Npt
 do t=1,tbead
     avgr(n,t) = avgr(n,t) + dsqrt(pos_old(1,t,n)*pos_old(1,t,n)+&
 &   pos_old(2,t,n)*pos_old(2,t,n))
     avgr2(n,t) = avgr2(n,t) + (pos_old(1,t,n)*pos_old(1,t,n)+&
 &   pos_old(2,t,n)*pos_old(2,t,n))
 enddo 
enddo ! do n=1,Npt ends here

return
end subroutine lindemann
!---------------------------------------------------------------------------
subroutine dsG(incrx,incry,en,dat,cent,Gval)
include 'param.dat'
integer :: ij,kbin,i,t,en
double precision :: incrx,incry,halfbinx,halfbiny, &
 & cent(Ndim,Nbin2),Gval(Nbin2),dat(Ndim,en,Npt)
double precision,parameter :: one=1.d0 
logical :: run_bin

halfbinx = 0.5d0*incrx
halfbiny = 0.5d0*incry

do i=1,Npt
do t = 1,en
 if((dabs(dat(1,t,i)).gt.(1.0e-10)).and.(dabs(dat(2,t,i)).gt.(1.0e-10))) then
      run_bin=.true.  

   do kbin = 1,Nbin2
   if (run_bin.eqv.(.true.)) then
      if ((dat(1,t,i) > (cent(1,kbin)-halfbinx)).and. &
      & (dat(1,t,i) <= (cent(1,kbin)+halfbinx)).and. &
      & (dat(2,t,i) > (cent(2,kbin)-halfbiny)).and. &
      & (dat(2,t,i) <= (cent(2,kbin)+halfbiny)) ) then
 
      Gval(kbin) = Gval(kbin) + one
           run_bin = .false.     
      endif
   endif
   enddo ! do kbin =1,Nbin2 ends here

 endif 
enddo ! do t=1,tbead ends here
enddo ! do i=1,Npt ends here

return
end subroutine dsG
!---------------------------------------------------------------------------
SUBROUTINE qsortd(x,ind,n)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-18  Time: 11:55:47

IMPLICIT NONE
INTEGER :: n
double precision :: x(n)
INTEGER :: ind(n)

!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

INTEGER   :: iu(21), il(21)
INTEGER   :: m, i, j, k, l, ij, it, itt, indx
REAL      :: r
double precision :: t

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO  i = 1, n
  ind(i) = i
END DO
m = 1
i = 1
j = n
r = .375

! TOP OF LOOP

20 IF (i >= j) GO TO 70
IF (r <= .5898437) THEN
  r = r + .0390625
ELSE
  r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = i + r*(j-i)
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
  ind(ij) = indx
  ind(i) = it
  it = indx
  t = x(it)
END IF

! INITIALIZE L

l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

indx = ind(j)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 80
END IF

il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) RETURN
i = il(m)
j = iu(m)

80 IF (j-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90
END SUBROUTINE qsortd
!-----------------------------------------------------------
subroutine dist1p_wothers(r1p_ij)
include 'param.dat'
integer :: iptcle,ij,jj,k,i,j,n
double precision :: r1p_ij(Npt,Npt)

do iptcle=1,Npt
  do n=1,Npt
    r1p_ij(iptcle,n)=0.d0
  enddo
enddo

do iptcle=1,Npt

  do jj=1,Npt

    if (jj /= iptcle) then

      if (jj < iptcle) then  !this if-construct is to avoid
        i=iptcle             !double counting of pairs
        j=jj
      else
        i=jj
        j=iptcle
      endif

      ij=j+((i-1)*(i-2))/2  !index ij of npair corr to i & j of N_plcle

      r1p_ij(iptcle,jj)=cmr_pp(ij)

    endif   !(jj /= iptcle)
  enddo   !jj=1,N_ptcle
enddo   !iptcle=1,N_ptcle

return
end subroutine dist1p_wothers
!-----------------------------------------------------------
subroutine nbr_angle(iptcle,r1p_dist,blkounter,psi)
include'param.dat'

integer :: i,k,iptcle,srt_ind(Npt),srt_ind2(nghbr),en, blkounter
double precision :: r1p_dist(Npt), pos_old_nbr(Ndim,nghbr),nbr_vec(Ndim,nghbr),angle(nghbr),  & 
  &  ang_vec,minang,angdiff(nghbr)
double complex :: psi
double precision, parameter :: pi=4.d0*datan(1.d0)

en = Npt
 call qsortd(r1p_dist,srt_ind,en)

  do k=1,nghbr
    pos_old_nbr(1,k)=cmpos_old(1,srt_ind(k+1)); pos_old_nbr(2,k)=cmpos_old(2,srt_ind(k+1))
 nbr_vec(1,k)=pos_old_nbr(1,k)-cmpos_old(1,iptcle); nbr_vec(2,k)=pos_old_nbr(2,k)-cmpos_old(2,iptcle)
  enddo

angle=0.d0

do k=1,nghbr
  angle(k)=datan2(nbr_vec(2,k), nbr_vec(1,k))
  if (angle(k) < 0.d0) angle(k) = 2*pi + angle(k)
enddo

minang=minval(angle)
do k=1,nghbr
  angle(k) = ( angle(k) - minang )
enddo

 call qsortd(angle,srt_ind2,nghbr)

do k=2,nghbr
  angdiff(k-1)=angle(srt_ind2(k))-angle(srt_ind2(k-1))
enddo
angdiff(6)=2*pi-angle(srt_ind2(nghbr))


if (maxval(angdiff) < 1.2d0*(2.d0*pi/3.d0)) then

  do k=1,nghbr
     psi=psi+exp(6.d0*(0,1)*angle(k))/dfloat(nghbr) 
   !    psi=psi+exp(6.d0*(0,1)*angdiff(k))/dfloat(nghbr) 
  enddo
 blkounter = blkounter + 1

endif

return
end subroutine nbr_angle
!-------------------------------------------------------------------------------
subroutine pairGBr(iptcle,incr,cent,Gval,tpsi,r1p_ij)
include 'param.dat'

integer :: ij,kbin,run_bin,i
integer :: srt_ind(Npt),imtcle,iptcle
double precision :: incr, halfbin,r1p_ij(Npt,Npt)
double precision :: cent(Nbin), Gval(Nbin),Pval(Nbin),dat(Npt)
 double complex :: psii,tpsi(Npt),value

halfbin=0.5d0*incr
 do i=1,Npt
  dat(i)=r1p_ij(iptcle,i)
 enddo

 call qsortd(dat,srt_ind,Npt)

run_bin = 1

do ij=1,Npt
    imtcle=0
    value=0.d0

  do kbin=run_bin,Nbin

    if ( (dat(srt_ind(ij)) > 1.0e-10) .and.  &
  & (dat(srt_ind(ij)) > (cent(kbin)-halfbin)).and. &
  & (dat(srt_ind(ij)) <= (cent(kbin)+halfbin)) ) then
     
    imtcle = srt_ind(ij)
    value = tpsi(imtcle)*conjg(tpsi(iptcle))     

    Gval(kbin) = Gval(kbin) + cdabs(value)
 
!   if(wrtph)then
!    Pval(kbin) = Pval(kbin) + datan2(dimag(value),dble(value))
!   endif   
      run_bin=kbin
    endif
  enddo

enddo

return
end subroutine pairGBr
!----------------------------------------------------------------------------
subroutine pairG(incr,en,dat,Nbin,cent,Gval,totpt)
implicit none
integer :: ij,kbin,run_bin,en,Nbin,totpt
integer :: srt_ind(en)
double precision :: incr, halfbin,dat(en)
double precision :: cent(Nbin), Gval(Nbin)
double precision, parameter :: one=1.d0

halfbin=0.5d0*incr 
 call qsortd(dat,srt_ind,en)

run_bin = 1

do ij=1,en

  do kbin=run_bin,Nbin
    if ( (dat(srt_ind(ij)) > 1.0e-10) .and.  &
  & (dat(srt_ind(ij)) > (cent(kbin)-halfbin)).and. &
  & (dat(srt_ind(ij)) <= (cent(kbin)+halfbin)) ) then
      Gval(kbin) = Gval(kbin) + one/dfloat(totpt)
      run_bin=kbin
    endif
  enddo

enddo

return
end subroutine pairG
!----------------------------------------------------------------------------
subroutine diameter(diam)
include 'param.dat'

integer :: i
double precision :: dat(Npt),diam

diam = 0.d0

do i=1,Npt
 dat(i) = dsqrt(cmpos_old(1,i)*cmpos_old(1,i)+&
 & cmpos_old(2,i)*cmpos_old(2,i))
enddo 

diam = maxval(dat)

return
end subroutine diameter
!-----------------------------------------------------------------------------





