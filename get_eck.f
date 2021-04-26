      subroutine get_eck(x)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      parameter (natoms = natom)
      common/masses/um(natoms),a0(3,natoms),utot
      dimension f(3,3),finv(3,3),transf(3,3),eig(3),tmp(3),f0(3,3)
      dimension eck(3,3),vec(3,3),geom(3,natoms),com(3)
      dimension x(3,natoms),dip(3)
C	  Rotate molecule into an eckart frame
C	  Inputs:
C	  x = geometry of molecule (will lose this information)
C	  Outputs:
C	  x = geometry rotated to eckart frame (overwritten from input)
      com = 0.d0
      do i = 1,natom
         do j = 1,3
        geom(j,i) = x(j,i)
            com(j) = com(j) + um(i)*geom(j,i)
         enddo
      enddo
      vec = 0.d0
      do j = 1,3
         com(j) = com(j)/utot
         do i = 1,natom
            g0 = geom(j,i)
            geom(j,i) = geom(j,i)-com(j)
            do k = 1,3
               vec(j,k) = vec(j,k) + um(i)*geom(j,i)*a0(k,i)
            enddo
         enddo
      enddo
c
      f = 0.d0
      do i = 1,3
         do j = 1,3
            do k = 1,3
               f(j,i) = f(j,i) + vec(k,j)*vec(k,i)
            enddo
         enddo
      enddo
      f0 = f
      call house(f,3,3,eig,tmp)
      do j = 1,3
         eig(j) = 1./sqrt(eig(j))
      enddo
      finv = 0.d0
      do i = 1,3
         do j = 1,3
            do k = 1,3
               finv(j,i) = finv(j,i)+eig(k)*f(j,k)*f(i,k)
            enddo
         enddo
      enddo
      eck = 0.d0
      do i = 1,3
         do j = 1,3
            do k = 1,3
               eck(j,i) = eck(j,i) + vec(j,k)*finv(k,i)
            enddo
         enddo
      enddo
      do i =1,3
        xn = 0.
        do k =1,3
           xn = xn + eck(k,i)*eck(k,i)
            enddo
        do k = 1,3
           eck(k,i) = eck(k,i)/sqrt(xn)
        enddo
      enddo
      x = 0.d0
      do i = 1,natom
         do j = 1,3
            do k = 1,3
               x(j,i) = x(j,i) + geom(k,i)*eck(k,j)
            enddo
         enddo
      enddo
      do i = 1,natom
         do j = 1,3
            geom(j,i) = x(j,i)
         enddo
      enddo
      return
      end
