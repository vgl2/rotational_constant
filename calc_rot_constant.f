        program calc_rotational_constant
        implicit real*8(a-h,o-z)
        parameter (nwavetot=200)
        parameter (nmax=70000)
        parameter (ndim=54)
        parameter (natoms=18)
        parameter (nblock=17)
        common/masses/um(natoms),a0(3,natoms),utot
        dimension x_eq(3,natoms),n(nwavetot),n0(nwavetot),time(nwavetot)
     1  ,psips(ndim,nmax),weight(nmax),coord(3,natoms,nmax),
     1  tensor(3,3),eig(3),scratch(3,3),avg_rot(3,3),
     1  tensor_diag(3,3),tensor_diag_inv(3,3),dev_rot(3,3),
     1  transform(3,3),tensor_inv(3,3),avg_tensor(3,3),
     1  tot_weight(nwavetot),rot_constant(3,3,nwavetot),
     1  rot_block(3,3,nblock,10),tot_rot(3,3,10),avg_rot_constant(3,3)

        open(unit=7,file='../combine_cage.dat',status='old',
     1  form='unformatted')
        open(unit=8,file='cage_eq.dat',status='old')
C		Calculates the rotational constant of all of the walkers buy rotating 
C		the molecule into an eckart frame and calculating the eigenvalues
C		of the moment of inertia tensor. 
C		Inputs:
C		coordinates of the walkers and their descendant weights
C		Reference geometry of the molecule.
C		Outputs:
C		Rotational Constants of the molecule in units of cm^-1. 
        do i = 1,natoms
            read(8,*) (x_eq(j,i),j=1,3)
            do j = 1,3
                x_eq(j,i) = x_eq(j,i)/0.52917721067d0
            enddo
        enddo
c       perform eckart embedding by choosing reference geometry
        utot = 0.d0
        do i = 1,natoms
            if (mod(i,3).eq.1) then
                um(i) = 15.99491461957d0
            else
                um(i) = 1.00782503223d0
            endif
            um(i) = um(i)*1822.88852962d0
            utot = utot + um(i)
        enddo
        call get_ref(x_eq)
        read(7) nwave
        do k = 1,nwave 
            tot_weight(k) = 0.d0
            read(7) n(k),n0(k),time(k)
            do i = 1,n(k)
                read(7) (psips(j,i),j=1,ndim),weight(i)
                tot_weight(k) = tot_weight(k) + weight(i)
            enddo
            do i = 1,n(k)
                ip = 0
                do j = 1,natoms
                    do l =1,3
                        ip = ip + 1
                        coord(l,j,i) = psips(ip,i)
                    enddo
                enddo
            enddo
c           perform eckart rotation
            do i = 1,3
                do j = 1,3
                    avg_tensor(i,j) = 0.d0
                enddo
            enddo
            do i = 1,n(k)
                call get_eck(coord(:,:,i))
c               calculate moment of inertia tensor
                do j = 1,3
                    do l = 1,3
                        tensor(l,j) = 0.d0
                    enddo
                enddo
                do j = 1,natoms
                    tensor(1,1)=tensor(1,1)+um(j)*(coord(2,j,i)**2+
     1              coord(3,j,i)**2)
                    tensor(2,2)=tensor(2,2)+um(j)*(coord(1,j,i)**2+
     1              coord(3,j,i)**2)
                    tensor(3,3)=tensor(3,3)+um(j)*(coord(1,j,i)**2+
     1              coord(2,j,i)**2)
                    tensor(1,2) = tensor(1,2) - um(j)*(coord(2,j,i)*
     1              coord(1,j,i))
                    tensor(2,1) = tensor(2,1) - um(j)*(coord(2,j,i)*
     1              coord(1,j,i))
                    tensor(1,3) = tensor(1,3) - um(j)*(coord(3,j,i)*
     1              coord(1,j,i))
                    tensor(3,1) = tensor(3,1) - um(j)*(coord(3,j,i)*
     1              coord(1,j,i))
                    tensor(2,3) = tensor(2,3) - um(j)*(coord(2,j,i)*
     1              coord(3,j,i))
                    tensor(3,2) = tensor(3,2) - um(j)*(coord(2,j,i)*
     1              coord(3,j,i))
                enddo
c               find the eigenvalues of the moment of inertia tensor
                transform = tensor
                call house(transform,3,3,eig,scratch)
                do j = 1,3
                    do l = 1,3
                        tensor_diag(j,l) = 0.d0
                        tensor_diag_inv(j,l) = 0.d0
                    enddo
                enddo
                do j = 1,3
                    do l =1,3
                        if (j.eq.l) then
                            tensor_diag(j,l) = eig(j)
                            tensor_diag_inv(j,l) = 1/eig(j)
                        endif
                    enddo
                enddo
                do j = 1,3
                    do l = 1,3
                        tensor_inv(j,l) = 0.d0
                    enddo
                enddo
                do j = 1,3
                    do l = 1,3
                        do m = 1,3
                            do no = 1,3
                                tensor_inv(m,no) = tensor_inv(m,no)+
     1                          transform(m,j)*tensor_diag_inv(j,l)*
     1                          transform(no,l)
                            enddo
                        enddo
                    enddo
                enddo
                do j = 1,3
                    do l = 1,3
                        avg_tensor(j,l) = avg_tensor(j,l)+weight(i)*
     1                  tensor_inv(j,l)
                    enddo
                enddo
            enddo
            do j = 1,3
                do l = 1,3
                    avg_tensor(j,l) =  avg_tensor(j,l)/tot_weight(k)
                enddo
            enddo
            do j = 1,3
                do l = 1,3
                    rot_constant(j,l,k) = avg_tensor(j,l)/2.d0
                enddo
            enddo
            do j = 1,3
                do l = 1,3
                    avg_rot_constant(j,l) = avg_rot_constant(j,l) + 
     1              rot_constant(j,l,k)
                enddo
            enddo
        enddo
c       divide into chunks of 17
        ip = 0
        is = 0
        do k = 1,nwave
            ip = ip + 1
            if (ip.eq.1) then
                is = 1
                do j = 1,3
                    do l = 1,3
                        rot_block(l,j,ip,is) = rot_constant(l,j,k)
                    enddo
                enddo
            else if (ip.lt.18) then
                do j = 1,3
                    do l = 1,3
                        rot_block(l,j,ip,is) = rot_constant(l,j,k)
                    enddo
                enddo
            else if (ip.eq.18) then
                ip = 1
                is = is + 1
                do j = 1,3
                    do l = 1,3
                        rot_block(l,j,ip,is) = rot_constant(l,j,k)
                    enddo
                enddo
            endif          
        enddo  
        do j = 1,3
            do l = 1,3
                do i = 1,9
                    tot_rot(l,j,i) = 0.d0
                enddo
                avg_rot(l,j) = 0.d0
                dev_rot(l,j) = 0.d0
            enddo
        enddo
        do j = 1,3
            do l = 1,3
                do i = 1,9
                    do k = 1,nblock
                        tot_rot(l,j,i) = tot_rot(l,j,i)+
     1                  (rot_block(l,j,k,i)/dfloat(nblock))
                    enddo
                enddo
            enddo
        enddo
        do j = 1,3
            do l = 1,3
                do i = 1,9
                    avg_rot(l,j) = avg_rot(l,j) + tot_rot(l,j,i)
                enddo
            enddo
        enddo
        do j = 1,3
            do l = 1,3
                avg_rot(l,j) = avg_rot(l,j)/9.d0
            enddo
        enddo
        do j = 1,3
            do l = 1,3
                dev_rot(l,j) = 0.d0
            enddo
        enddo
        do j = 1,3
            do l =1,3
                do i = 1,9
                    dev_rot(l,j) = dev_rot(l,j) + ((tot_rot(l,j,i)-
     1              avg_rot(l,j))**2)
                enddo
            enddo
        enddo
        do j = 1,3
            do l = 1,3
                dev_rot(l,j) = sqrt(dev_rot(l,j)/9.d0)
            enddo
        enddo
        print *, 'average'
        do j = 1,3
            print *, (avg_rot(j,l)*219474.6,l=1,3)
        enddo
        print *, 'deviations'
        do j = 1,3
            print *, (dev_rot(j,l)*219474.6,l=1,3)
        enddo
        end program
