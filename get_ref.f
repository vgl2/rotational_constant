      subroutine get_ref(ref)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      parameter (natoms = natom)
      common/masses/um(natoms),a0(3,natoms),utot
      dimension rot(3,3),brot(3),tmp(3),ref(3,natoms),com(3)
C	  Takes a reference geometry to separate the coupling in 
C	  rotations and vibrations. 
C	  Inputs:
C	  ref = reference geometry in units of bohr
C     Outputs:
C	  a0 = geometry rotated into the eckart frame. 
      com = 0.
      do i = 1,natom
         do j = 1,3
            com(j) = com(j) + um(i)*ref(j,i)
         enddo
      enddo
      do j = 1,3
         com(j) = com(j)/utot
         do i = 1,natom
            ref(j,i) = ref(j,i)-com(j)
         enddo
      enddo
      com = 0.
      do i = 1,natom
         do j = 1,3
            com(j) = com(j) + um(i)*ref(j,i)
         enddo
      enddo
      rot = 0.
      do i = 1,natom
         rot(1,1) = rot(1,1) + um(i)*(ref(2,i)**2+ref(3,i)**2)
         rot(2,2) = rot(2,2) + um(i)*(ref(1,i)**2+ref(3,i)**2)
         rot(3,3) = rot(3,3) + um(i)*(ref(2,i)**2+ref(1,i)**2)
         rot(1,2) = rot(1,2) - um(i)*(ref(2,i)*ref(1,i))
         rot(2,1) = rot(2,1) - um(i)*(ref(2,i)*ref(1,i))
         rot(1,3) = rot(1,3) - um(i)*(ref(3,i)*ref(1,i))
         rot(3,1) = rot(3,1) - um(i)*(ref(3,i)*ref(1,i))
         rot(2,3) = rot(2,3) - um(i)*(ref(2,i)*ref(3,i))
         rot(3,2) = rot(3,2) - um(i)*(ref(2,i)*ref(3,i))
      enddo
      call house(rot,3,3,brot,tmp)
      a0 = 0.
      do i = 1,natom
         do j =1,3
            do k = 1,3
               a0(j,i) = a0(j,i) + rot(k,j)*ref(k,i)
            enddo
         enddo
      enddo
      return
