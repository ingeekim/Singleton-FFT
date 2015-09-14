#line 1 "primfact.F"
      subroutine primfact(n,nfact,maxf,multp) 
!!  Determines the prime factors of the input number n.
!!  Returns its numer of prime factors nfac,
!!    the maximum prime factor maxf,
!!    and the multipls of its square-free prime factors multp.
!!  I.G. Kim, Oct. 2003
!!  Eliminated 'go to' statemens    I.G. Kim, May. 2006
! internal variables
      use kinds, only: ikind
      use mpicfg, only: inode
      implicit none
      logical, dimension(:), allocatable :: prime
      integer(kind=ikind), intent(in) :: n
      integer(kind=ikind), intent(out) :: nfact,maxf,multp
      !! local variables
      integer(kind=ikind), dimension(:), allocatable :: prbelow,cfact
      integer(kind=ikind), dimension(:,:), allocatable :: prnfact
      integer(kind=ikind) :: i,j,m,nfac,p,pcount
!        10        20        30        40        50        60        70
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      if(n.le.2) then
        if(inode.eq.0) then
          write(*,*) 'SUBROUTINE primfact:'
          write(*,*) 'n should be greater than 2!'
        endif  
        stop 'primfact1'
      end if
      allocate(prime(n))

      ! set all numbers are prime for the beginning
      prime(:)=.true.

      ! find prime numbers
      pcount=0
      do i=2,n                 ! 2 is a prime number
        if(prime(i)) then      ! if i is a prime number,
            do j=i+i,n,i       ! eliminate its multiples
              prime(j)=.false.
            end do
            pcount=pcount+1
        end if
      end do

      ! put the found prime numbers into prbelow
      allocate(prbelow(pcount))
      prbelow(:)=0
      j=1
      do i=2,n
        if(prime(i)) then
          prbelow(j)=i
          j=j+1
        end if
      end do
      deallocate(prime)  !! clear, no more used

      !
      ! begin factorizing the number n into prime numbers
      !
      allocate(cfact(pcount))
      cfact(:)=0
      m=n
      factoring: do
        i=1
        factor: do
          p=prbelow(i)
          if(mod(m,p).ne.0) then
            i=i+1
            cycle factor
          else
            cfact(i)=cfact(i)+1
            m=m/p
            if(m.gt.1) cycle factoring
          end if
          exit factor
        end do factor
        exit factoring
      end do factoring

      !
      ! count the number of prime factors
      !
      nfac=0  !! the number of kinds of prime factors
      nfact=0  !! the number of prime factors
      do i=1,pcount
        if(cfact(i).gt.0) then
          nfac=nfac+1
          nfact=nfact+cfact(i)
        end if
      end do

      !
      ! pack the prime factors of n and
      ! its number occurred into prnfact(nfac,2)
      !
      allocate(prnfact(2,nfac))
      prnfact=0
      j=1
      do i=1,pcount
        if(cfact(i).ne.0) then
          prnfact(1,j)=prbelow(i)
          prnfact(2,j)=cfact(i)
          j=j+1
        end if
      end do
      deallocate(prbelow)  ! no more used
      deallocate(cfact)    ! no more used

      ! the largest prime factor of n
      maxf=prnfact(1,nfac)

      ! the multiples of square-free factors of the prime factors of n
      multp=1
      do i=1,nfac
        j=mod(prnfact(2,i),2)
        if (j.ne.0) multp=multp*prnfact(1,i)
      end do 
      deallocate(prnfact)  !! clear
      return
      end subroutine primfact
!        10        20        30        40        50        60        70
!23456789012345678901234567890123456789012345678901234567890123456789012
