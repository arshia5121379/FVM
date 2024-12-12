
module Variables
    implicit none
    ! definition of headers
    integer, parameter :: c_b = 0
    integer, parameter :: b_n = 1
    integer, parameter :: b_s = 2
    integer, parameter :: b_w = 4
    integer, parameter :: b_e = 8
    integer, parameter :: b_nw = b_n + b_w
    integer, parameter :: b_sw = b_s + b_w
    integer, parameter :: b_ne = b_n + b_e
    integer, parameter :: b_se = b_s + b_e

    integer, parameter :: c_f = 16
    integer, parameter :: c_x = c_f - 1
    integer, parameter :: c_a = c_f + b_n + b_s + b_e + b_w

    integer, parameter :: c_o = 256
    integer, parameter :: c_w = 512
    integer, parameter :: c_s = 1024
    integer, parameter :: c_n = 2048
    integer, parameter :: c_e = 4096

    integer, parameter :: c_wo = c_w + c_o
    integer, parameter :: c_ns = c_n + c_s
    integer, parameter :: c_sw = c_s + c_w
    integer, parameter :: c_nw = c_n + c_w
    integer, parameter :: c_no = c_n + c_o
    integer, parameter :: c_so = c_s + c_o

    integer, parameter :: c_swo = c_s + c_w + c_o
    integer, parameter :: c_nsw = c_n + c_s + c_w
    integer, parameter :: c_nwo = c_n + c_w + c_o
    integer, parameter :: c_nso = c_n + c_s + c_o

    integer, parameter :: c_nswo = c_n + c_s + c_w + c_o

    ! for flag initialization
    integer :: ihi, ilo, low, up
    real :: mx, my, rad1, x, y

    ! for computing f and g
    real :: du2dx, duvdx, duvdy, dv2dy, laplu, laplv

    ! for Poisson Equation:
    real :: add, beta_2, beta_mod, p0, rdx2, rdy2
    integer :: eps_e, eps_n, eps_s, eps_w, iter

    ! for main
    real :: delt, delx, dely, eps, gamma, gx, gy, omega, re, res, t, t_end, ui, vi
    integer :: i, j, ibound, ifull, imax, itermax, itersor, iwrite, jmax
    integer :: outcount, printable, we, wn, ws, ww, cycles
    real :: xlength, ylength
    real, dimension(:,:), allocatable :: u, v, p, rhs, uout, vout, pout, f, g, psi, zeta
    integer, dimension(:,:), allocatable :: flag
    character(len = 30) :: problem
    character(len = 2) :: outnum
end module

program Navier_Stocks_chorin
    USE Variables
    implicit none

    call Problem_Setup
    call Allocation

    ! setting initial values
    u(0:imax+1, 0:jmax+1) = ui
    v(0:imax+1, 0:jmax+1) = vi
    p(0:imax+1, 0:jmax+1) = 1.0

    if ( problem == "backwardstep" ) then
        u(0:imax+1, 0:jmax/2) = 0.0
    end if

    ! initialize flags.
    ! initialize the boundary cells.
    do i = 0, imax+1
        flag(i,0) = c_b
        flag(i, jmax+1) = c_b
    end do

    do j = 1, jmax
        flag(0, j) = c_b
        flag(imax+1, j) = c_b
    end do

    ! initialize the fluid cells.
    do i = 1, imax
        do j = 1, jmax
            flag(i, j) = c_f
        end do
    end do

    ! specialize to the problem
    if ( problem == "backwardstep" ) then
        do i = 1, jmax
            do j = 1, jmax/2
                flag(i, j) = c_b
            end do
        end do
        write ( *, '(a)' ) ' problem = ', problem
    else if ( problem == "cylinder" ) then
        mx = 20.0 / 41.0 * jmax * dely
        my = mx
        rad1 = 5.0 / 41.0 * jmax * dely

        do i = 1, imax
            do j = 1, jmax
                x = (i-0.5) * delx
                y = (j-0.5) * dely
                if ( (x-mx)*(x-mx) + (y-my)*(y-my) <= rad1 * rad1 ) then
                    flag(i, j) = c_b
                end if
            end do
        end do

    write ( *, '(a)' ) ' problem = ', problem

    else if ( problem == "cavity" ) then
        write ( *, '(a)' ) ' problem = ', problem
    end if

    ! flags for boundaries :
    ibound = 0
    do i=1 , imax
        do j=1 , jmax
            if (IAND(flag(i,j),c_f) /= c_f) then
                ibound = ibound + 1
            end if

            flag(i,j) = flag(i,j) + (IAND(flag(i-1,j),c_f) *b_w &
                                & +  IAND(flag(i+1,j),c_f) *b_e &
                                & +  IAND(flag(i,j-1),c_f) *b_s &
                                & +  IAND(flag(i,j+1),c_f) *b_n) / c_f

        end do
    end do

    ! initialize the boundary conditions.
    do j = 0, jmax+1

        ! western and eastern boundary.
        ! free slip, u = 0, d(vdn) = 0.
        if ( ww == 1 ) then
            u(0,j) = 0.0
            v(0,j) = v(1,j)

        ! no slip, u = 0, v = 0 at the boundary by averaging.
        else if ( ww == 2 ) then
            u(0,j) = 0.0
            v(0,j) = (-1.0) * v(1,j)

        ! outflow
        else if ( ww == 3 ) then
            u(0,j) = u(1,j)
            v(0,j) = v(1,j)

        ! periodic, left and right cells overlap.
        else if ( ww == 4 ) then
            u(0,j) = u(imax-1,j)
            v(0,j) = v(imax-1,j)
            u(1,j) = v(imax,j)
            v(1,j) = v(imax,j)
            p(1,j) = p(imax,j)
        end if

        ! free slip
        if ( we == 1 ) then
            u(imax,j) = 0.0
            v(imax+1,j) = v(imax,j)

        ! no slip
        else if ( we == 2 ) then
            u(imax,j) = 0.0
            v(imax+1,j) = -v(imax,j)

        ! outflow
        else if ( we == 3 ) then
            u(imax,j) = u(imax-1,j)
            v(imax+1,j) = v(imax,j)

        ! periodic
        else if ( we == 4 ) then
            u(imax,j) = u(1,j)
            v(imax+1,j) = v(2,j)
        end if

    end do
        ! northern and southern boundary
    do i = 0, imax+1

        if ( wn == 1 ) then
            v(i,jmax+1) = 0.0
            u(i,jmax+1) = u(i,jmax)

        else if ( wn == 2 ) then
            v(i,jmax) = 0.0
            u(i,jmax+1) = -u(i,jmax)

        else if ( wn == 3 ) then
            v(i,jmax) = v(i,jmax-1)
            u(i,jmax+1) = u(i,jmax)

        else if ( wn == 4 ) then
            v(i,jmax) = v(i,1)
            u(i,jmax+1) = u(i,2)
        end if

        if ( ws == 1 ) then
            v(i,0) = 0.0
            u(i,0) = u(i,1)

        else if ( ws == 2 ) then
            v(i,0) = 0.0
            u(i,0) = -u(i,1)

        else if ( ws == 3 ) then
            v(i,0) = v(i,1)
            u(i,0) = u(i,1)

        else if ( ws == 4 ) then
            v(i,0) = v(i,jmax-1)
            u(i,0) = u(i,jmax-1)
            u(i,0) = u(i,jmax)
            p(i,0) = p(i,jmax)
        end if

    end do

    ! set the boundary values at inner obstacle cells (only no-slip).
    do i = 1, imax
        do j = 1, jmax

            ! mask c_x = 000f filters the obstacle cells adjacent to fluid cells.
            if ( iand(flag(i,j), c_x) /= 0 ) then

                select case (flag(i,j))

                    case (b_n)
                        v(i,j) = 0.0
                        u(i,j) = -u(i,j+1)
                        u(i-1,j) = -u(i-1,j+1)

                    case (b_e)
                        u(i,j) = 0.0
                        v(i,j-1) = -v(i+1,j-1)

                    case (b_s)
                        v(i,j-1) = 0.0
                        u(i,j) = -u(i,j-1)
                        u(i-1,j) = -u(i-1,j-1)

                    case (b_w)
                        u(i-1,j) = 0.0
                        v(i,j-1) = -v(i-1,j-1)

                    case (b_ne)
                        v(i,j) = 0.0
                        u(i,j) = -u(i,j+1)
                        u(i-1,j) = -u(i-1,j+1)

                    case (b_se)
                        v(i,j-1) = 0.0
                        u(i,j) = -u(i,j-1)
                        u(i-1,j) = -u(i-1,j-1)

                    case (b_sw)
                        v(i,j-1) = 0.0
                        u(i-1,j) = -u(i-1,j-1)

                    case (b_nw)
                        v(i,j) = 0.0
                        u(i-1,j) = -u(i-1,j+1)
                    case default
                end select

            end if

        end do
    end do

    ! backwardstep
    ! u = 1.0 at the left boundary
    if ( problem == "backwardstep" ) then
        do j = jmax/2+1, jmax
            u(0,j) = 1.0
        end do
    ! cylinder
    ! u = 1.0 at left boundary
    else if ( problem == "cylinder" ) then
        v(0,0) = 2.0 * vi - v(1,0)
        do j = 1, jmax
            u(0,j) = ui
            v(0,j) = 2.0 * vi - v(1,j)
        end do
    ! cavity
    ! u = 1.0 at the upper boundary
    else if ( problem == "cavity" ) then
        do i = 0, imax
            u(i,jmax+1) = 2.0 - u(i,jmax)
        end do
    end if

    t = 0.0
    cycles = 0

    ! time loop.
    do while ( t < t_end )

    ! determine fluid cells
        ifull = imax * jmax - ibound

        ! compute tentative velocity field (f, g)
        ! compute flux field f
        ! only if both adjacent cells are fluid cells.
        do i = 1, imax-1
            do j = 1, jmax
                if ((( iand(flag(i,j), c_f) /= 0) .and. flag(i,j) < c_e) .and. &
                    (( iand(flag(i+1,j),c_f) /= 0) .and. flag(i+1,j) < c_e )) then

                    du2dx = ((u(i,j)+u(i+1,j))*(u(i,j)+u(i+1,j)) &
                            - gamma*abs(u(i,j)+u(i+1,j))*(u(i,j)-u(i+1,j))) &
                            / (4.0*delx)

                    duvdy = ((v(i,j)+v(i+1,j))*(u(i,j)+u(i,j+1)) &
                            + gamma*abs(v(i,j)+v(i+1,j))*(u(i,j)-u(i,j+1))) &
                            - ((v(i,j-1)+v(i+1,j-1))*(u(i,j)-u(i,j-1)) &
                            - gamma*abs(v(i,j-1)+v(i+1,j-1))*(u(i,j)-u(i,j-1))) &
                            / (4.0*dely)

                    laplu = ( u(i+1,j) - 2.0 * u(i,j) + u(i-1,j) ) / (delx * delx) &
                            + ( u(i,j+1) - 2.0 * u(i,j) + u(i,j-1) ) / (dely * dely)

                    f(i,j) = u(i,j) + delt * ( laplu / re - du2dx - duvdy + gx )
                else
                    f(i,j) = u(i,j)
                end if
            end do
        end do

        ! compute flux field g
        ! only if both adjacent cells are fluid cells
        do i = 1, imax
            do j = 1, jmax-1
                if ((( iand(flag(i,j), c_f) /= 0) .and. (flag(i,j) < c_e)) .and. &
                    (( iand(flag(i,j+1), c_f) /= 0) .and. (flag(i,j+1) < c_e))) then

                    duvdx = ((u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) &
                            + gamma*abs(u(i,j)+u(i,j+1))*(v(i,j)-v(i+1,j)) &
                            - (u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j)) &
                            - gamma*abs(u(i-1,j)+u(i-1,j+1))*(v(i,j)-v(i-1,j))) &
                            / (4.0*delx)

                    duvdy = ((v(i,j)+v(i,j+1))*(v(i,j)+v(i,j+1)) &
                            + gamma*abs(v(i,j)+v(i,j+1))*(v(i,j)-v(i,j+1)) &
                            - (v(i,j-1)+v(i,j-1))*(v(i,j)-v(i,j-1)) &
                            - gamma*abs(v(i,j-1)+v(i,j-1))*(v(i,j)-v(i,j-1))) &
                            / (4.0*dely)

                    laplv = ( v(i+1,j) - 2.0 * v(i,j) + v(i-1,j) ) / delx / delx &
                            + ( v(i,j+1) - 2.0 * v(i,j) + v(i,j-1) ) / dely / dely

                    g(i,j) = v(i,j) + delt * ( laplv / re - duvdx - duvdy + gy )
                else
                    g(i,j) = v(i,j)
                end if
            end do
        end do

        ! f and g at external boundary
        do j = 1, jmax
            f(0,j) = u(0,j)
            f(imax,j) = u(imax,j)
        end do

        do i = 1, imax
            g(i,0) = v(i,0)
            g(i,jmax) = v(i,jmax)
        end do

        !*********************************!
        ! compute right-hand side for pressure equation.
        do i = 1, imax
            do j = 1, jmax

                ! only for fluid and non-surface cells.
                if ( iand(flag(i,j), c_f) /= 0 .and. flag(i,j) < c_o ) then
                    rhs(i,j) = ( ( f(i,j) - f(i-1,j) ) / delx &
                                + ( g(i,j) - g(i,j-1) ) / dely ) / delt
                end if
            end do
        end do

        !*********************************!
        ! solve the pressure equation by successive over relaxation (SOR).
        if ( 0 < ifull ) then
            rdx2 = 1.0 / delx / delx
            rdy2 = 1.0 / dely / dely
            beta_2 = -omega / ( 2.0 * ( rdx2 + rdy2 ) )

            p0 = 0.0
            do i = 1, imax
                do j = 1, jmax
                    if ( iand( flag(i,j), c_f ) /= 0 ) then
                        p0 = p0 + p(i,j) * p(i,j)
                    end if
                end do
            end do

            p0 = sqrt( p0 / ifull )

            if ( p0 < 0.0001 ) then
                p0 = 1.0
            end if

            ! sor iteration
            do iter = 1, itermax
                    ! copy values at external boundary...
                    do i = 1, imax
                        p(i,0) = p(i,1)
                        p(i,jmax+1) = p(i,jmax)
                    end do

                    do j = 1, jmax
                        p(0,j) = p(1,j)
                        p(imax+1,j) = p(imax,j)
                    end do

                    ! and at interior boundary cells.
                    do i = 1, imax
                        do j = 1, jmax
                            if ( b_n <= flag(i,j) .and. flag(i,j) <= b_se ) then
                                if ( flag(i,j) == b_n ) then
                                    p(i,j) = p(i,j+1)
                                else if ( flag(i,j) == b_e ) then
                                    p(i,j) = p(i+1,j)
                                else if ( flag(i,j) == b_s ) then
                                    p(i,j) = p(i,j-1)
                                else if ( flag(i,j) == b_w ) then
                                    p(i,j) = p(i-1,j)
                                else if ( flag(i,j) == b_ne ) then
                                    p(i,j) = 0.5 * ( p(i,j+1) + p(i+1,j) )
                                else if ( flag(i,j) == b_se ) then
                                    p(i,j) = 0.5 * ( p(i,j-1) + p(i+1,j) )
                                else if ( flag(i,j) == b_sw ) then
                                    p(i,j) = 0.5 * ( p(i,j-1) + p(i-1,j) )
                                else if ( flag(i,j) == b_nw ) then
                                    p(i,j) = 0.5 * ( p(i,j+1) + p(i-1,j) )
                                end if
                            end if
                        end do
                    end do

                    ! relaxation method for fluid cells.
                    do i = 1, imax
                        do j = 1, jmax
                            if ( flag(i,j) < c_o ) then
                                if ( iand( flag(i,j), c_f ) /= 0 .and. flag(i,j)< c_o) then
                                    p(i,j) = ( 1.0 - omega ) * p(i,j) &
                                    - beta_2 * ( ( p(i+1,j) + p(i-1,j) ) * rdx2 &
                                    + ( p(i,j+1) + p(i,j-1) ) * rdy2 &
                                    - rhs(i,j) )
                                end if
                            end if
                        end do
                    end do
                    ! computation of the residual.
                    res = 0.0
                    do i = 1, imax
                        do j = 1, jmax

                            ! only fluid cells
                            if ( ( iand( flag(i,j), c_f ) /= 0 ) .and. flag(i,j) < c_o ) then
                                add = ( p(i+1,j) - 2.0 * p(i,j) + p(i-1,j) ) * rdx2 &
                                    + ( p(i,j+1) - 2.0 * p(i,j) + p(i,j-1) ) * rdy2 - rhs(i,j)
                                res = res + add * add
                            end if
                        end do
                    end do

                    res = sqrt( res / ifull ) / p0
                    if ( res < eps ) then
                        go to 100
                    end if
                    itersor = itersor+1
                end do

                100 write (6, '(a,g10.4,a,g10.4,a,i4,a,g12.6,a,i4,i4,i4)' ) &
                    ' t=', t+delt, ' delt=', delt, ' its=', itersor, ' res=', res

            end if
            itersor=0

            !*********************************!
            ! compute the new velocity field.
            do i = 1, imax - 1
                do j = 1, jmax
                    if ((( iand( flag(i,j), c_f ) > 0 ) .and. (flag(i,j) < c_e )) .and. &
                        (( iand( flag(i+1,j), c_f ) > 0 ) .and. (flag(i+1,j) < c_e ))) then
                        u(i,j) = f(i,j) - ( p(i+1,j) - p(i,j) ) * delt / delx
                    end if
                end do
            end do

            do i = 1, imax
                do j = 1, jmax - 1
                    if ((( iand( flag(i,j), c_f ) > 0 ) .and. (flag(i,j) < c_e )) .and. &
                        (( iand( flag(i,j+1), c_f ) > 0 ) .and. (flag(i,j+1) < c_e ))) then
                        v(i,j) = g(i,j) - ( p(i,j+1) - p(i,j) ) * delt / dely
                    end if
                end do
            end do

            ! set the boundary conditions for the next time step
            do j = 0, jmax+1

                ! western and eastern boundary.
                ! free slip, u = 0, d(vdn) = 0.
                if ( ww == 1 ) then
                    u(0,j) = 0.0
                    v(0,j) = v(1,j)

                ! no slip, u = 0, v = 0 at the boundary by averaging.
                else if ( ww == 2 ) then
                    u(0,j) = 0.0
                    v(0,j) = (-1.0) * v(1,j)

                ! outflow
                else if ( ww == 3 ) then
                    u(0,j) = u(1,j)
                    v(0,j) = v(1,j)

                ! periodic, left and right cells overlap.
                else if ( ww == 4 ) then
                    u(0,j) = u(imax-1,j)
                    v(0,j) = v(imax-1,j)
                    v(1,j) = v(imax,j)
                    p(1,j) = p(imax,j)
                end if

                ! free slip
                if ( we == 1 ) then
                    u(imax,j) = 0.0
                    v(imax+1,j) = v(imax,j)
                ! no slip
                else if ( we == 2 ) then
                    u(imax,j) = 0.0
                    v(imax+1,j) = -v(imax,j)

                ! outflow
                else if ( we == 3 ) then
                    u(imax,j) = u(imax-1,j)
                    v(imax+1,j) = v(imax,j)

                ! periodic
                else if ( we == 4 ) then
                    u(imax,j) = u(1,j)
                    v(imax+1,j) = v(2,j)
                end if

            end do

                ! northern and southern boundary
            do i = 0, imax+1

                if ( wn == 1 ) then
                    v(i,jmax) = 0.0
                    u(i,jmax+1) = u(i,jmax)
                else if ( wn == 2 ) then
                    v(i,jmax) = 0.0
                    u(i,jmax+1) = -u(i,jmax)
                else if ( wn == 3 ) then
                    v(i,jmax) = v(i,jmax-1)
                    u(i,jmax+1) = u(i,jmax)
                else if ( wn == 4 ) then
                    v(i,jmax) = v(i,1)
                    u(i,jmax+1) = u(i,2)
                end if

                if ( ws == 1 ) then
                    v(i,0) = 0.0
                    u(i,0) = u(i,1)
                else if ( ws == 2 ) then
                    v(i,0) = 0.0
                    u(i,0) = -u(i,1)
                else if ( ws == 3 ) then
                    v(i,0) = v(i,1)
                    u(i,0) = u(i,1)
                else if ( ws == 4 ) then
                    v(i,0) = v(i,jmax-1)
                    u(i,0) = u(i,jmax-1)
                    u(i,1) = u(i,jmax)
                    p(i,1) = p(i,jmax)
                end if

            end do

            ! set the boundary values at inner obstacle cells (only no-slip).
            do i = 1, imax
                do j = 1, jmax

                    ! mask c_x = 000f filters the obstacle cells adjacent to fluid cells.
                    if ( iand( flag(i,j), c_x ) /= 0 ) then
                        select case ( flag(i,j) )

                            case (b_n)
                                v(i,j) = 0.0
                                u(i,j) = -u(i,j+1)
                                u(i-1,j) = -u(i-1,j+1)

                            case (b_e)
                                u(i,j) = 0.0
                                v(i,j-1) = -v(i+1,j-1)
                                v(i,j) = -v(i+1,j)

                            case (b_s)
                                v(i,j-1) = 0.0
                                u(i,j) = -u(i,j-1)
                                u(i-1,j) = -u(i-1,j-1)

                            case (b_w)
                                u(i-1,j) = 0.0
                                v(i,j) = -v(i-1,j)
                                v(i,j-1) = -v(i-1,j-1)

                            case (b_ne)
                                v(i,j) = 0.0
                                u(i,j) = -u(i,j+1)
                                u(i-1,j) = -u(i-1,j+1)

                            case (b_se)
                                v(i,j-1) = 0.0
                                u(i,j) = -u(i,j-1)
                                u(i-1,j) = -u(i-1,j-1)

                            case (b_sw)
                                v(i,j-1) = 0.0
                                u(i-1,j) = -u(i-1,j-1)

                            case (b_nw)
                                v(i,j) = 0.0
                                u(i-1,j) = -u(i-1,j+1)

                            case default
                                ! no action for default
                        end select
                    end if
                end do
            end do

            !*********************************!
            ! setting special boundary conditions
            ! for the next time step.

            ! backwardstep
            ! u = 1.0 at the left boundary
            if ( problem == "backwardstep" ) then
                do j = jmax/2+1, jmax
                    u(0,j) = 1.0
                end do

            ! cylinder
            ! u = 1.0 at left boundary
            else if ( problem == "cylinder" ) then
                v(0,0) = 2.0 * vi - v(1,0)
                do j = 1, jmax
                    u(0,j) = ui
                    v(0,j) = 2.0 * vi - v(1,j)
                end do

            ! cavity
            ! u = 1.0 at the upper boundary
            else if ( problem == "cavity" ) then
                do i = 0, imax
                    u(i,jmax+1) = 2.0 - u(i,jmax)
                end do

            ! unrecognized problem.
            end if

            !*********************************!
            write(*,*) t
            write(*,*) delt
            write(*,*) cycles
            !*********************************!

            ! advance the time
            t = t + delt
            cycles = cycles + 1

            if (mod(cycles,outcount)==0) then
                !call digit_to_ch (cycles, outnum)
                !outnum=char(cycles)
                printable = cycles/outcount
                write (outnum,'(i2)') printable
                write (*,*) printable
                ! read(outnum,*) cycles

                ! computation of the vorticity (zeta) at the upper
                ! right corner of cell (i,j) (only if the corner is surrounded by fluid cells)
                do i = 1, imax-1
                    do j = 1, jmax-1
                        if (((iand(flag(i,j), c_f) /= 0) &
                            .and. flag(i,j) < c_e ) .and. &
                            ((iand(flag(i+1,j), c_f) /= 0) &
                            .and. flag(i+1,j) < c_e ) .and. &
                            ((iand(flag(i,j+1), c_f) /= 0) &
                            .and. flag(i,j+1) < c_e ) .and. &
                            ((iand(flag(i+1,j+1), c_f) /= 0) &
                            .and. flag(i+1,j+1) < c_e )) then
                            zeta(i,j) = ( u(i,j+1) - u(i,j) ) / dely &
                                       - ( v(i+1,j) - v(i,j) ) / delx
                        else
                            zeta(i,j) = 0.0
                        end if
                    end do
                end do
                zeta(0,0:jmax) = 0.0
                zeta(0:imax,0) = 0.0
                zeta(imax,0:jmax) = 0.0
                zeta(0:imax,jmax) = 0.0

                !*********************************!
                ! computation of the stream function at the upper
                ! right corner of cell (i,j), but only if both lower cells are fluid cells.
                do i = 0, imax
                    psi(i,0) = 0.0
                    do j = 1, jmax
                        if (((iand(flag(i,j), c_f) /= 0) &
                            .and. (flag(i,j) < c_e)) .or. &
                            ((iand(flag(i+1,j), c_f) /= 0) &
                            .and. (flag(i+1,j) < c_e))) then
                            psi(i,j) = psi(i,j-1) + u(i,j) * dely
                        else
                            psi(i,j) = psi(i,j-1)
                        end if
                    end do
                end do
                ! computing output velocity and pressure
                ! data and writing them out in separate
                ! files for post-processing

                do j = 1, jmax
                    do i = 1, imax
                        if ( iand( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
                            uout(i,j) = ( u(i,j) + u(i-1,j) ) / 2.0
                        else
                            uout(i,j) = 0.0
                        end if
                    end do
                end do
                uout(0,0:jmax) = 0.0
                uout(0:imax,0) = 0.0

                do j = 1, jmax
                    do i = 1, imax
                        if ( iand( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
                            vout(i,j) = ( v(i,j) + v(i,j-1) ) / 2.0
                        else
                            vout(i,j) = 0.0
                        end if
                    end do
                end do
                vout(0,0:jmax) = 0.0
                vout(0:imax,0) = 0.0

                do j = 1, jmax
                    do i = 1, imax
                        if ( iand( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
                            pout(i,j) = p(i,j)
                        else
                            pout(i,j) = 0.0
                        end if
                    end do
                end do
                pout(0,0:jmax) = 0.0
                pout(0:imax,0) = 0.0

                !*********************************!
                ! output u-velocity for visualization.
                open (unit=1, file='u_data_'//trim(outnum)//'.txt')
                do j = 0, jmax
                    do i = 0, imax
                        write (1, '(i5,1x,i5,1x,i5,1x,f12.6)') i, j, uout(i,j)
                    end do
                end do
                close (1)

                !*********************************!
                ! output v-velocity for visualization.
                open (2, file='v_data_'//trim(outnum)//'.txt')
                do j = 0, jmax
                    do i = 0, imax
                        write (2, '(i5,1x,i5,1x,i5,1x,f12.6)') i, j, vout(i,j)
                    end do
                end do
                close (2)

                !*********************************!
                ! output pressure field for visualization.
                open (3, file='pressure_'//trim(outnum)//'.txt')
                do j = 0, jmax
                    do i = 0, imax
                        write (3, '(i5,1x,i5,1x,i5,1x,f12.6)') i, j, pout(i,j)
                    end do
                end do
                close (3)

                !*********************************!
                ! output the stream function for visualization.
                open (4, file='psi_data_'//trim(outnum)//'.txt')
                do j = 0, jmax
                    do i = 0, imax
                        write (4, '(i5,1x,i5,1x,i5,1x,f12.6)') i, j, psi(i,j)
                    end do
                end do
                close (4)
                !*********************************!
                write(*,*) 'problem: ', problem
                write(*,*) 'results have been printed out'
                write(*,*) 'for time ', t, ' (s)'
                write(*,*) ' '

            endif

        end do
        !*********************************!
        !*********************************!
        !*********************************!
        call Deallocation

        stop

end program

subroutine Allocation
    USE Variables
    implicit none
    ! allocate arrays.
    allocate ( f(0:imax+1, 0:jmax+1) )
    allocate ( flag(0:imax+1, 0:jmax+1) )
    allocate ( g(0:imax+1, 0:jmax+1) )
    allocate ( p(0:imax+1, 0:jmax+1) )
    allocate ( psi(0:imax, 0:jmax) )
    allocate ( rhs(0:imax+1, 0:jmax+1) )
    allocate ( u(0:imax+1, 0:jmax+1) )
    allocate ( v(0:imax+1, 0:jmax+1) )
    allocate ( uout(0:imax, 0:jmax) )
    allocate ( vout(0:imax, 0:jmax) )
    allocate ( pout(0:imax, 0:jmax) )
    allocate ( zeta(0:imax, 0:jmax) )

end subroutine Allocation

subroutine Deallocation
    USE Variables
    implicit none
    ! free memory.
    deallocate (f)
    deallocate (flag)
    deallocate (g)
    deallocate (p)
    deallocate (psi)
    deallocate (rhs)
    deallocate (u)
    deallocate (v)
    deallocate (zeta)
    deallocate (uout)
    deallocate (vout)
    deallocate (pout)

end subroutine Deallocation


subroutine Problem_Setup
    USE Variables
    implicit none

    problem = 'cavity'
    xlength = 1.0
    ylength = 1.0
    imax    = 50
    jmax    = 50

    t_end = 5.0
    delt  = 0.005

    outcount = 200
    itermax  = 150

    eps = 0.001
    omega = 1.7
    gamma = 0.9

    re = 100
    gx = 0.0
    gy = 0.0
    ui = 0.0
    vi = 0.0

    ww = 2
    we = 2
    wn = 2
    ws = 2

    itersor = 0
    ifull   = 0
    ibound  = 0

end subroutine Problem_Setup
