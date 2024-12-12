
program SIMPLE
    implicit none
    integer:: N,i,j,itt,datatime
    real(8):: dx,dy,Nu,Re,alphaV,alphaP,err,maxerr,&
            & uE,uW,vN,vSs,aE,aW,aN,aS,Aeee,aee,Annn,ann,x,y,ap,u_lid
    real(8),allocatable,dimension(:,:)::u,v,us,vs,de,dn,p,pp,b,unew,vnew,pnew

    N = 101
    dx = 1.d0/(N-1.d0)
    dy = dx
    Re = 200 ; Nu = 1.d0/Re
    alphaP = 0.8d0 ; alphaV = 1.d0
    err = 10 ; maxerr = 1e-7 ; datatime = 1000
    u_lid = 1.d0

    allocate(u(N,N+1),v(N+1,N),p(N+1,N+1),us(N,N+1),vs(N+1,N),de(N,N+1),dn(N+1,N),&
            &pp(N+1,N+1),b(N+1,N+1),pnew(N+1,N+1),unew(N,N+1),vnew(N+1,N)  )

    u=0; u(:,N+1)=2.d0*u_lid ; us=0 ; de=0; v=0; vs=0; dn=0; p=1.d0; pp=0; unew=0; vnew=0; pnew=0;

    itt = 0
    do while(err>maxerr)
        ! computation of Ustar i starts from 1 to N and j starts from 1 to N+1
        do j=2,N
            do i=2,N-1
                uE  = 0.5d0*(u(i,j)+u(i+1,j))
                uW  = 0.5d0*(u(i,j)+u(i-1,j))
                vN  = 0.5d0*(v(i,j)+v(i,j+1))
                vSs = 0.5d0*(v(i,j)+v(i,j-1))

                aE  =  0.5d0*uE*dy + Nu*dy/dx
                aW  = -0.5d0*uW*dy + Nu*dy/dx
                aN  =  0.5d0*vN*dx + Nu*dx/dy
                aS  = -0.5d0*vSs*dx+ Nu*dx/dy

                aee     = aE+aW+aN+aS
                Aeee    = dy
                de(i,j) = Aeee/aee

                us(i,j) = (u(i+1,j)*aE + u(i-1,j)*aW + u(i,j+1)*aN + u(i,j-1)*aS)/aee + de(i,j)*(p(i,j)-p(i+1,j))
            end do
        end do
        !implementation of B-C fro Us
        us(:,N+1) = u_lid - us(:,N)
        us(:,1)   = -us(:,2)
        us(1,:)   = 0
        us(N,:)   = 0

        us = alphaV*us + (1.d0-alphaV)*u

        ! computation of vstar i starts from 1 to N+1 and j starts from 1 to N
        do i=2,N
            do j=2,N-1
                uE  = 0.5d0*(u(i,j)+u(i+1,j))
                uW  = 0.5d0*(u(i,j)+u(i-1,j))
                vN  = 0.5d0*(v(i,j)+v(i,j+1))
                Vss = 0.5d0*(v(i,j)+v(i,j-1))

                aE  =  0.5d0*uE*dy + Nu*dy/dx
                aW  = -0.5d0*uW*dy + Nu*dy/dx
                aN  =  0.5d0*vN*dx + Nu*dx/dy
                aS  = -0.5d0*vSs*dx+ Nu*dx/dy

                aee     = aE+aW+aN+aS
                Aeee    = dy
                dn(i,j) = Aeee/aee

                vs(i,j) = (v(i+1,j)*aE + v(i-1,j)*aW + v(i,j+1)*aN + v(i,j-1)*aS)/aee + dn(i,j)*(p(i,j)-p(i,j+1))
            end do
        end do
        !implementation of B-C for Vs
        vs(:,N)   = 0
        vs(:,1)   = 0
        vs(1,:)   = -vs(2,:)
        vs(N+1,:) = -vs(N,:)

        vs = alphaV*vs + (1.d0-alphaV)*v

        ! calculation of intermediate pressure
        pp = 0.0
        do i=2,N
            do j=2,N
                aE = de(i+1,j)*dy
                aW = de(i-1,j)*dy
                aN = dn(i,j+1)*dx
                aS = dn(i,j-1)*dx
                b(i,j) = ((us(i,j)+us(i+1,j)*0.5d0)-(us(i,j)+us(i-1,j)*0.5d0))*dy+&
                        &((vs(i,j)+vs(i,j+1)*0.5d0)-(vs(i,j)+vs(i,j-1)*0.5d0))*dx
                ap = aE+aW+aN+aS
                pp(i,j) = (aE*pp(i+1,j)+aW*pp(i-1,j)+aN*pp(i,j+1)+aS*pp(i,j-1))/ap
            end do
        end do

        ! pressure correction
        do i=2,N
            do j=2,N
                pnew(i,j) = p(i,j) + alphaP*pp(i,j)
            end do
        end do
        pnew(1,:)   = pnew(2,:)
        pnew(N+1,:) = pnew(N,:)
        pnew(:,N+1) = pnew(:,N)
        pnew(:,1)   = pnew(:,2)

        ! velocity correction
        do j=2,N
            do i=2,N-1
                unew(i,j) = us(i,j) + de(i,j)*(pp(i,j)-pp(i+1,j))
            end do
        end do
        unew(:,N+1) = u_lid - unew(:,N)
        unew(:,1)   = -unew(:,2)
        unew(1,:)   = 0
        unew(N,:)   = 0

        do i=2,N
            do j=2,N-1
                vnew(i,j) = vs(i,j) + dn(i,j)*(pp(i,j)-pp(i,j+1))
            end do
        end do
        vnew(:,N)   = 0
        vnew(:,1)   = 0
        vnew(1,:)   = -vnew(2,:)
        vnew(N+1,:) = -vnew(N,:)

        !error calculation continuty error
        err = 0
        do i=2,N
            do j=2,N
                err = err + abs(b(i,j))
            end do
        end do

        u = unew
        v = vnew
        p = pnew

        itt = itt + 1

        if(mod(itt,datatime)==0) then
            print*, 'time step = ' , itt
            print*, 'error' , err , 'time step : ',itt
        end if

    end do

    do i=1,N
        do j=1,N
            x = (i-1)*dx ; y = (j-1)*dy
            print*,x,y,0.5d0*(u(i+1,j)+u(i,j)),0.5d0*(v(i,j)+v(i,j+1)),0.25d0*&
            &(p(i,j)+p(i,j+1)+p(i+1,j)+p(i+1,j+1))
        end do
    end do

    ! deallocation of the arrays

end program

