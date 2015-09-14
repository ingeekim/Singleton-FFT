#line 1 "cfft.F"
      subroutine cfft(a,b,ntot,n,nspan,isgn)
!        10        20        30        40        50        60        70
!23456789012345678901234567890123456789012345678901234567890123456789012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Multivariate complex fourier transform, computed in place
!     using mixed-radix fast fourier transform algorithm.
!  By r. c. singleton, stanford research institute, oct. 1968
!  Arrays a and b originally hold the real and imaginary
!     components of the data, and return the real and
!     imaginary components of the resulting fourier coefficients.
!  Multivariate data is indexed according to the fortran
!     array element successor function, without limit
!     on the number of implied multiple subscripts.
!  The subroutine is called once for each variate.
!     the calls for a multivariate transform may be in any order.
!
!     ntot is the total number of complex data values.
!     n is the dimension of the current variable.
!     nspan/n is the spacing of consucutive data values
!
!     while indexing the current variable.
!  The sign of isn determines the sign of the complex
!     exponential, and the magnitude of isn is normally one.
!  For a single-variate transform,
!     ntot = n = nspan = (number of complex data values), f.g.
!     call cft(a,b,n,n,n,1)
!
!  A tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
!     is computed by
!     call cft(a,b,n1*n2*n3,n1,n1,1)
!     call cft(a,b,n1*n2*n3,n2,n1*n2,1)
!     call cft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
!
!  The data may alternatively be stored in a single complex
!     array a, then the magnitude of isn changed to two to
!     give the correct indexing increment and the second parameter
!     used to pass the initial address for the sequence of
!     imaginary values, e.g.
!        real s(2)
!        equivalence (a,s)
!        ....
!        ....
!        call cft(a,s(2),ntot,n,nspan,2)
!  Arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
!     are used for temporary storage. if the available storage
!     is insufficient, the program is terminated by a stop.
!    
!     maxn must be .ge. the number of factors of n.
!     maxf must be .ge. the maximum prime factor of n.
!     maxp must be .gt. the number of prime factors of n.
!  In addition, maxn is assumed to be odd.  
!  If the square-free portion k of n has two or more prime factors, 
!     then maxp must be .ge. k-1.
!  Array storage in nfac for a maximum of 11 factors of n. !! eliminated
!  If n has more than one square-free factor, the product of the
!     square-free factors must be .le. 210  !! eliminated
!
!! The temporary storages are allocated dynamically, so that there is
!!    no limit on the array size of n, except the integer limit.
!! I.G.Kim 2003
!! 
!! A function float(i) is replaced by real(i,fp) by i.g.kim Apr.2005
!!
!! Completely redesigned to eliminate 'go to' statements
!!                                             by i.g.kim May. 2006
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use kinds, only: ikind,fp,ofh
      use constants, only: pi,tpi
      use mpicfg, only: inode
      implicit none
      !! Define the interface to parent procedures
      integer(kind=ikind), intent(in) :: ntot,n,nspan
      integer(kind=ikind), intent(inout) :: isgn
      real(kind=fp), dimension(ntot), intent(inout) :: a,b
      !! local variables
      !! Define the interface to the primfact subroutine
      interface
         subroutine primfact(n,nfact,maxf,multp) 
         use kinds, only: ikind
         implicit none
         integer(kind=ikind), intent(in) :: n
         integer(kind=ikind), intent(out) :: nfact,maxf,multp
         end subroutine primfact
      end interface
!! Define allocation related variables
      integer(kind=ikind), dimension(:), allocatable :: nfac,np
      real(kind=fp), dimension(:), allocatable :: at,ck,bt,sk
      integer(kind=ikind) :: maxn,maxf,maxp
!! Define internal variables
      integer(kind=ikind) :: i,ii,inc,j,jc,jf,jj,jx
      integer(kind=ikind) :: k,k1,k2,k3,k4,kk,kkx,ks,kspan,kspnn,kt
      integer(kind=ikind) :: m,nn,nt
      real(kind=fp) :: rad,radf,s72,c72,s120,sd
      real(kind=fp) :: aa,aj,ajm,ajp,ak,akp,akm,bb,bj,bjm,bjp,bk,bkm,bkp
      real(kind=fp) :: c1,c2,c3,cd,s1,s2,s3
!        10        20        30        40        50        60        70
!23456789012345678901234567890123456789012345678901234567890123456789012
      equivalence(i,ii)
      !!
      !! Check the size of array, and determine the size of
      !!     the temporary storages
      !!

      !! The origianl memory sizes
      call primfact(n,maxn,maxf,maxp)
      if(maxp.le.maxn) maxp=maxn+1
      if(mod(maxn,2).eq.0) maxn=maxn+1  ! make maxn odd

      !! Now allocate the temporary storages
      allocate(nfac(maxn+1))
      allocate(at(maxf))     
      allocate(ck(maxf))     
      allocate(bt(maxf))     
      allocate(sk(maxf))     
      allocate(np(maxp))
      if((n-1).lt.0) then
        if(inode.eq.0) then
          write(ofh,*) 'cfft array size n is less than 1'
        endif  
        stop 'cfft1'
      end if

      inc=isgn
      !! the following constants are rad = 2.*pi, s72 = sin(0.4*pi),
      !! c72 = cos(0.4*pi) and s120 = sqrt(0.75)
      rad=tpi
      s72=sin(0.4_fp*pi)
      c72=cos(0.4_fp*pi)
      s120=sqrt(0.75_fp)
      if(isgn.lt.0) then
        s72=-s72
        s120=-s120
        rad=-rad
        inc=-inc
      end if

      !! size parameters of the arrays, a and b
      nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*real(jc,fp)*0.5_fp
      i=0
      jf=0

      !!
      !! determine the factors of n
      !!
      m=0
      k=n

      do while(k-(k/16)*16.eq.0)
        m=m+1
        nfac(m)=4
        k=k/16
      end do

      j=3
      jj=9
      det_fac3: do
        do while(mod(k,jj).eq.0)
          m=m+1
          nfac(m)=j
          k=k/jj
        end do
        j=j+2
        jj=j*j
        if(jj.gt.k) exit det_fac3
      end do det_fac3

      det_fac_gt4: if(k.gt.4) then
        if(k-(k/4)*4.eq.0) then
          m=m+1
          nfac(m)=2
          k=k/4
        end if
        kt=m
        j=2
        det_fac_odd: do
          if(mod(k,j).eq.0) then
            m=m+1
            nfac(m)=j
            k=k/j
          end if
          j=((j+1)/2)*2+1
          if(j.gt.k) exit det_fac_odd
        end do det_fac_odd
      else det_fac_gt4
        kt=m
        nfac(m+1)=k
        if(k.ne.1) m=m+1
      end if det_fac_gt4

      if(kt.ne.0) then
        j=kt
        do jx=j,1,-1
          m=m+1
          nfac(m)=nfac(jx)
        end do ! jx_loop
      end if

      !!
      !!===> now compute fourier transform
      !!
      fourier_transform: do
        sd=radf/real(kspan,fp)
        cd=2.0_fp*sin(sd)**2
        sd=sin(sd+sd)
        kk=1
        i=i+1
        !! the others rather than factor of 2 are managed 
        !! in factor of 4 parts
        is_factor4: if(nfac(i).ne.2) then
          !---> transform for factor of 4
          is_odd_factor: if(nfac(i).ne.4) then
            !---> transform for odd factors
            k=nfac(i)
            kspnn=kspan
            kspan=kspan/k
            which_odds: select case (k)
              case (3) which_odds
                !---> transform for factor of 3 (optional code)
                factor3: do
                  k1=kk+kspan
                  k2=k1+kspan
                  ak=a(kk)
                  bk=b(kk)
                  aj=a(k1)+a(k2)
                  bj=b(k1)+b(k2)
                  a(kk)=ak+aj
                  b(kk)=bk+bj
                  ak=-0.5_fp*aj+ak
                  bk=-0.5_fp*bj+bk
                  aj=(a(k1)-a(k2))*s120
                  bj=(b(k1)-b(k2))*s120
                  a(k1)=ak-bj
                  b(k1)=bk+aj
                  a(k2)=ak+bj
                  b(k2)=bk-aj
                  kk=k2+kspan
                  if(kk.lt.nn) cycle factor3
                  kk=kk-nn
                  if(kk.gt.kspan) exit factor3
                end do factor3
                !! the factor-of-3 finished, do multiply rotation factors
              case (5) which_odds
                !---> transform for factor of 5 (optional code)
                c2=c72**2-s72**2
                s2=2.0_fp*c72*s72
                factor5: do
                  k1=kk+kspan
                  k2=k1+kspan
                  k3=k2+kspan
                  k4=k3+kspan
                  akp=a(k1)+a(k4)
                  akm=a(k1)-a(k4)
                  bkp=b(k1)+b(k4)
                  bkm=b(k1)-b(k4)
                  ajp=a(k2)+a(k3)
                  ajm=a(k2)-a(k3)
                  bjp=b(k2)+b(k3)
                  bjm=b(k2)-b(k3)
                  aa=a(kk)
                  bb=b(kk)
                  a(kk)=aa+akp+ajp
                  b(kk)=bb+bkp+bjp
                  ak=akp*c72+ajp*c2+aa
                  bk=bkp*c72+bjp*c2+bb
                  aj=akm*s72+ajm*s2
                  bj=bkm*s72+bjm*s2
                  a(k1)=ak-bj
                  a(k4)=ak+bj
                  b(k1)=bk+aj
                  b(k4)=bk-aj
                  ak=akp*c2+ajp*c72+aa
                  bk=bkp*c2+bjp*c72+bb
                  aj=akm*s2-ajm*s72
                  bj=bkm*s2-bjm*s72
                  a(k2)=ak-bj
                  a(k3)=ak+bj
                  b(k2)=bk+aj
                  b(k3)=bk-aj
                  kk=k4+kspan
                  if(kk.lt.nn)  cycle factor5
                  kk=kk-nn
                  if(kk.gt.kspan) exit factor5
                end do factor5      
                !! the factor-of-5 finished, do multiply rotation factors
              case default which_odds
                if(k.ne.jf) then
                  jf=k
                  s1=rad/real(k,fp)
                  c1=cos(s1)
                  s1=sin(s1)
                  if(jf.gt.maxf) then
                    if(inode.eq.0) then
                      write(ofh,*) 'cfft array size n has a prime'//
     &                             ' factor larger than maxf'
                     write(ofh,*) '  or the square-free portion of n is'
                      write(ofh,*) '  greater than maxp+1'
                      write(ofh,*) 'The array size n is ',n
                    endif  
                    stop 'cfft2'
                  end if
                  ck(jf)=1.0_fp
                  sk(jf)=0.0_fp
                  j=1
                  !!!!! THIS PART IS SLIGHTLY CHANGED COMPARED TO 'cfft.f' !!!
                  init_odd_factor: do
                    ck(j)=ck(k)*c1+sk(k)*s1
                    sk(j)=ck(k)*s1-sk(k)*c1
                    ck(k-1)=ck(j)
                    sk(k-1)=-sk(j)
                    k=k-1
                    j=j+1
                    if(j.ge.k) exit init_odd_factor
                  end do init_odd_factor
                end if
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                other_odd_factors: do
                  k1=kk
                  k2=kk+kspnn
                  aa=a(kk)
                  bb=b(kk)
                  ak=aa
                  bk=bb
                  j=1
                  k1=k1+kspan
                  other_part1: do
                    k2=k2-kspan
                    j=j+1
                    at(j)=a(k1)+a(k2)
                    ak=at(j)+ak
                    bt(j)=b(k1)+b(k2)
                    bk=bt(j)+bk
                    j=j+1
                    at(j)=a(k1)-a(k2)
                    bt(j)=b(k1)-b(k2)
                    k1=k1+kspan
                    if(k1.ge.k2) exit other_part1
                  end do other_part1
                  a(kk)=ak
                  b(kk)=bk
                  k1=kk
                  k2=kk+kspnn
                  j=1
                  other_part2: do
                    k1=k1+kspan
                    k2=k2-kspan
                    jj=j
                    ak=aa
                    bk=bb
                    aj=0.0_fp
                    bj=0.0_fp
                    k=1
                    sub_other_part2: do
                      k=k+1
                      ak=at(k)*ck(jj)+ak
                      bk=bt(k)*ck(jj)+bk
                      k=k+1
                      aj=at(k)*sk(jj)+aj
                      bj=bt(k)*sk(jj)+bj
                      jj=jj+j
                      if(jj.gt.jf) jj=jj-jf
                      if(k.ge.jf) exit sub_other_part2
                    end do sub_other_part2
                    k=jf-j
                    a(k1)=ak-bj
                    b(k1)=bk+aj
                    a(k2)=ak+bj
                    b(k2)=bk-aj
                    j=j+1
                    if(j.ge.k) exit other_part2
                  end do other_part2
                  kk=kk+kspnn
                  if(kk.le.nn) cycle other_odd_factors
                  kk=kk-nn
                  if(kk.gt.kspan) exit other_odd_factors
                end do other_odd_factors 
              !! end of the select
            end select which_odds
            !
            !---> multiply by rotation factor (except for factors of 2 and 4)
            !
            !! factor-of-2 and 4 are finished, do permutation
            if(i.eq.m) exit fourier_transform
            kk=jc+1
            multiply_rot_factor: do
              c2=1.0_fp-cd
              s1=sd
              sub_mult_rot_factor: do
                c1=c2
                s2=s1
                kk=kk+kspan
                multiply: do
                  ak=a(kk)
                  a(kk)=c2*ak-s2*b(kk)
                  b(kk)=s2*ak+c2*b(kk)
                  kk=kk+kspnn
                  if(kk.le.nt) cycle multiply
                  ak=s1*s2
                  s2=s1*c2+c1*s2
                  c2=c1*c2-ak
                  kk=kk-nt+kspan
                  if(kk.gt.kspnn) exit multiply
                end do multiply 
                c2=c1-(cd*c1+sd*s1)
                s1=s1+(sd*c1-cd*s1)
                ! the following three statements compensate
                ! for truncation error. 
                ! if rounded arithmetic is used, they may be deleted.
                ! c1=0.5/(c2**2+s1**2)+0.5
                ! s1=c1*s1
                ! c2=c1*c2
                kk=kk-kspnn+jc
                if(kk.gt.kspan) exit sub_mult_rot_factor
              end do sub_mult_rot_factor
              kk=kk-kspan+jc+inc
              if(kk.gt.(jc+jc)) exit multiply_rot_factor
            end do multiply_rot_factor
            cycle fourier_transform 
            !! end of the odd factor fourier transforms
          end if is_odd_factor
          kspnn=kspan
          kspan=kspan/4
          ! factor4 loops
          factor4: do
            c1=1.0_fp
            s1=0.0_fp
            sub_factor4: do
              k1=kk+kspan
              k2=k1+kspan
              k3=k2+kspan
              akp=a(kk)+a(k2)
              akm=a(kk)-a(k2)
              ajp=a(k1)+a(k3)
              ajm=a(k1)-a(k3)
              a(kk)=akp+ajp
              ajp=akp-ajp
              bkp=b(kk)+b(k2)
              bkm=b(kk)-b(k2)
              bjp=b(k1)+b(k3)
              bjm=b(k1)-b(k3)
              b(kk)=bkp+bjp
              bjp=bkp-bjp
              direction4: if(isgn.lt.0) then
                akp=akm+bjm
                akm=akm-bjm
                bkp=bkm-ajm
                bkm=bkm+ajm
              else direction4
                akp=akm-bjm
                akm=akm+bjm
                bkp=bkm+ajm
                bkm=bkm-ajm
              end if direction4
              size4: if(s1.eq.0.0_fp) then
                a(k1)=akp
                b(k1)=bkp
                a(k2)=ajp
                b(k2)=bjp
                a(k3)=akm
                b(k3)=bkm
              else size4
                a(k1)=akp*c1-bkp*s1
                b(k1)=akp*s1+bkp*c1
                a(k2)=ajp*c2-bjp*s2
                b(k2)=ajp*s2+bjp*c2
                a(k3)=akm*c3-bkm*s3
                b(k3)=akm*s3+bkm*c3
              end if size4
              kk=k3+kspan
              if(kk.le.nt) cycle sub_factor4
              c2=c1-(cd*c1+sd*s1)
              s1=(sd*c1-cd*s1)+s1
              ! the following three statements compensate for truncation
              ! error. if rounded arithmetic is used, they may be deleted.
              ! c1=0.5/(c2**2+s1**2)+0.5
              ! s1=c1*s1
              ! c1=c1*c2
              ! next statement should be deleted
              ! if non-rounded arithmetic is used
              c1=c2
              c2=c1**2-s1**2
              s2=2.0_fp*c1*s1
              c3=c2*c1-s2*s1
              s3=c2*s1+s2*c1
              kk=kk-nt+jc
              if(kk.gt.kspan) exit sub_factor4
            end do sub_factor4
            kk=kk-kspan+inc
            if(kk.gt.jc) exit factor4
          end do factor4
          if(kspan.eq.jc) exit fourier_transform
          cycle fourier_transform
        end if is_factor4

        !!---> transform for factor of 2 (including rotation factor)
        !---> factor2 part1
        kspan=kspan/2
        k1=kspan+2
        init_factor2: do
          k2=kk+kspan
          ak=a(k2)
          bk=b(k2)
          a(k2)=a(kk)-ak
          b(k2)=b(kk)-bk
          a(kk)=a(kk)+ak
          b(kk)=b(kk)+bk
          kk=k2+kspan
          if(kk.le.nn) cycle init_factor2
          kk=kk-nn
          if(kk.gt.jc) exit init_factor2
        end do init_factor2
        ! when radix-2 FFT finished, start the permutation
        if(kk.gt.kspan) exit fourier_transform
        !-->  factor2 part2
        factor2: do
          c1=1.0_fp-cd
          s1=sd
          sub_factor2: do
            k2=kk+kspan
            ak=a(kk)-a(k2)
            bk=b(kk)-b(k2)
            a(kk)=a(kk)+a(k2)
            b(kk)=b(kk)+b(k2)
            a(k2)=c1*ak-s1*bk
            b(k2)=s1*ak+c1*bk
            kk=k2+kspan
            if(kk.lt.nt) cycle sub_factor2
            k2=kk-nt
            c1=-c1
            kk=k1-k2
            if(kk.gt.k2) cycle sub_factor2
            ak=c1-(cd*c1+sd*s1)
            s1=(sd*c1-cd*s1)+s1
            ! the following three statements compensate for truncation
            ! error. if rounded arithmetic is used, they may be deleted.
            ! c1=0.5/(ak**2+s1**2)+0.5
            ! s1=c1*s1
            ! c1=c1*ak
            ! next statement should be deleted
            ! if non-rounded arithmetic is used
            c1=ak
            kk=kk+jc
            if(kk.ge.k2) exit sub_factor2
          end do sub_factor2
          k1=k1+inc+inc
          kk=(k1-kspan)/2+jc
          if(kk.gt.(jc+jc)) exit factor2
        end do factor2
      end do fourier_transform

      !!
      !!===> permute the results to normal order---done in two stages
      !!
      !---> permutation for square factors of n
      np(1)=ks
      is_permute_square: if(kt.ne.0) then
        k=kt+kt+1
        if(m.lt.k) k=k-1
        j=1
        np(k+1)=jc
        do 
          np(j+1)=np(j)/nfac(j)
          np(k)=np(k+1)*nfac(j)
          j=j+1
          k=k-1
          if(j.ge.k) exit
        end do
        k3=np(k+1)
        kspan=np(2)
        kk=jc+1
        k2=kspan+1
        j=1
        is_single: if(n.eq.ntot) then
          ! permutation for single-variate transform (optional code)
          permute_single: do
            do
              ak=a(kk)
              a(kk)=a(k2)
              a(k2)=ak
              bk=b(kk)
              b(kk)=b(k2)
              b(k2)=bk
              kk=kk+inc
              k2=kspan+k2
              if(k2.ge.ks) exit
            end do
            sub_permute_single: do
              k2=k2-np(j)
              j=j+1
              k2=np(j+1)+k2
              if(k2.gt.np(j)) cycle sub_permute_single
              j=1
              do
                if(kk.lt.k2) cycle permute_single
                kk=kk+inc
                k2=kspan+k2
                if(k2.ge.ks) exit
              end do
              if(kk.ge.ks) exit sub_permute_single
            end do sub_permute_single
            exit permute_single
          end do permute_single
          jc=k3
        else is_single
          ! permutation for multivariate transform
          permute_multi: do
            k=kk+jc
            do kkx=kk,k-1,inc
              ak=a(kkx)
              a(kkx)=a(k2)
              a(k2)=ak
              bk=b(kkx)
              b(kkx)=b(k2)
              b(k2)=bk
              k2=k2+inc
            end do
            kk=kkx+ks-jc
            k2=k2+ks-jc
            if(kk.lt.nt) cycle permute_multi
            k2=k2-nt+kspan
            kk=kk-nt+jc
            if(k2.lt.ks) cycle permute_multi
            sub_permute_multi: do
              k2=k2-np(j)
              j=j+1
              k2=np(j+1)+k2
              if(k2.gt.np(j)) cycle sub_permute_multi
              j=1
              do
                if(kk.lt.k2) cycle permute_multi
                kk=kk+jc
                k2=k2+kspan
                if(k2.ge.ks) exit
              end do
              if(kk.ge.ks) exit sub_permute_multi
            end do sub_permute_multi
            exit permute_multi
          end do permute_multi
          jc=k3
        end if is_single
      end if is_permute_square

      if((2*kt+1).ge.m) then
         deallocate(nfac)
         deallocate(at)     
         deallocate(ck)     
         deallocate(bt)     
         deallocate(sk)     
         deallocate(np)     
         return
      end if
      kspnn=np(kt+1)

      !---> permutation for square-free factors of n
      j=m-kt
      nfac(j+1)=1
      do
        nfac(j)=nfac(j)*nfac(j+1)
        j=j-1
        if(j.eq.kt) exit
      end do

      kt=kt+1
      nn=nfac(kt)-1
      if(nn.gt.maxp) then
        ! error control
        if(inode.eq.0) then
          write(ofh,*) 'cfft array size n has a prime factor'//
     &               'larger than maxf'
          write(ofh,*) '  or the square-free portion of n is'
          write(ofh,*) '  greater than maxp+1'
          write(ofh,*) 'The array size n is ',n
        endif  
        stop 'cfft2'
      end if
      jj=0
      j=0
      prep_square_free: do
        k2=nfac(kt)
        k=kt+1
        kk=nfac(k)
        j=j+1
        if(j.gt.nn) exit prep_square_free
        sub_prep_squ_free: do
          jj=kk+jj
          if(jj.lt.k2) then
            np(j)=jj
            cycle prep_square_free
          else
            jj=jj-k2
            k2=kk
            k=k+1
            kk=nfac(k)
            cycle sub_prep_squ_free
          end if
        end do sub_prep_squ_free
      end do prep_square_free
      j=0
      !---> determine the permutation cycles of length greater than 1
      det_permute_cycles: do
        j=j+1
        kk=np(j)
        if(kk.lt.0) cycle det_permute_cycles
        if(kk.ne.j) then
          sub_det_perm_cycle: do
            k=kk
            kk=np(k)
            np(k)=-kk
            if(kk.eq.j) exit sub_det_perm_cycle
          end do sub_det_perm_cycle
          k3=kk
          cycle det_permute_cycles
        end if
        np(j)=-j
        if(j.eq.nn) exit det_permute_cycles
      end do det_permute_cycles
      maxf=inc*maxf
      !---> reorder a and b, following the permutation cycles
      permute_square_free: do
        j=k3+1
        nt=nt-kspnn
        ii=nt-inc+1
        if(nt.lt.0) exit permute_square_free
        permute_squ_free: do
          j=j-1
          if(np(j).lt.0) cycle permute_squ_free
          jj=jc
          sub_permute_squ_free: do
            kspan=jj
            if(jj.gt.maxf) kspan=maxf
            jj=jj-kspan
            k=np(j)
            kk=jc*k+ii+jj
            k1=kk+kspan
            k2=0
            backup: do
              k2=k2+1
              at(k2)=a(k1)
              bt(k2)=b(k1)
              k1=k1-inc
              if(k1.eq.kk) exit backup
            end do backup 
            actual_permute_squ_free: do
              k1=kk+kspan
              k2=k1-jc*(k+np(k))
              k=-np(k)
              sub_actual_permute: do
                a(k1)=a(k2)
                b(k1)=b(k2)
                k1=k1-inc
                k2=k2-inc
                if(k1.eq.kk) exit sub_actual_permute 
              end do sub_actual_permute
              kk=k2
              if(k.eq.j) exit actual_permute_squ_free
            end do actual_permute_squ_free
            k1=kk+kspan
            k2=0
            restore: do
              k2=k2+1
              a(k1)=at(k2)
              b(k1)=bt(k2)
              k1=k1-inc
              if(k1.eq.kk) exit restore
            end do restore
          if(jj.eq.0) exit sub_permute_squ_free
          end do sub_permute_squ_free
          if(j.eq.1) exit permute_squ_free
        end do permute_squ_free
      end do permute_square_free
      !---> all done, clear variables
      deallocate(nfac)
      deallocate(at)     
      deallocate(ck)     
      deallocate(bt)     
      deallocate(sk)     
      deallocate(np)     
      return
      end subroutine cfft
