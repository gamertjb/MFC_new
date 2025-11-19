#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_capillary.fpp'

!> @brief This module is used to compute source terms for surface tension model
module m_surface_tension

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion

    use m_weno

    use m_muscl                !< Monotonic Upstream-centered (MUSCL)
                               !! schemes for conservation laws

    use m_helper

    use m_boundary_common

    implicit none

    private; public :: s_initialize_surface_tension_module, &
 s_compute_capillary_source_flux, &
 s_get_capillary, &
 s_finalize_surface_tension_module

    !> @name color function gradient components and magnitude
    !> @{
    type(scalar_field), allocatable, dimension(:) :: c_divs
    type(scalar_field), allocatable, dimension(:) :: c2_divs
    !> @)
    $:GPU_DECLARE(create='[c_divs,c2_divs]')

    !> @name cell boundary reconstructed gradient components and magnitude
    !> @{
    real(wp), allocatable, dimension(:, :, :, :) :: gL_x, gR_x, gL_y, gR_y, gL_z, gR_z
    real(wp), allocatable, dimension(:, :, :, :) :: g2L_x, g2R_x, g2L_y, g2R_y, g2L_z, g2R_z
    !> @}
    $:GPU_DECLARE(create='[gL_x,gR_x,gL_y,gR_y,gL_z,gR_z,g2L_x,g2R_x,g2L_y,g2R_y,g2L_z,g2R_z]')

    type(int_bounds_info) :: is1, is2, is3, iv
    $:GPU_DECLARE(create='[is1,is2,is3,iv]')

    logical :: has_secondary_color_function

contains

    impure subroutine s_initialize_surface_tension_module

        integer :: j

        has_secondary_color_function = surface_tension .and. num_fluids > 2

        @:ALLOCATE(c_divs(1:num_dims + 1))

        do j = 1, num_dims + 1
            @:ALLOCATE(c_divs(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(c_divs(j))
        end do

        if (has_secondary_color_function) then
            @:ALLOCATE(c2_divs(1:num_dims + 1))
            do j = 1, num_dims + 1
                @:ALLOCATE(c2_divs(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(c2_divs(j))
            end do
        end if

        @:ALLOCATE(gL_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
        @:ALLOCATE(gR_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

        @:ALLOCATE(gL_y(idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
        @:ALLOCATE(gR_y(idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

        @:ALLOCATE(gL_z(idwbuff(3)%beg:idwbuff(3)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, num_dims + 1))
        @:ALLOCATE(gR_z(idwbuff(3)%beg:idwbuff(3)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, num_dims + 1))

        if (has_secondary_color_function) then
            @:ALLOCATE(g2L_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
            @:ALLOCATE(g2R_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

            @:ALLOCATE(g2L_y(idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
            @:ALLOCATE(g2R_y(idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

            @:ALLOCATE(g2L_z(idwbuff(3)%beg:idwbuff(3)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, num_dims + 1))
            @:ALLOCATE(g2R_z(idwbuff(3)%beg:idwbuff(3)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(1)%beg:idwbuff(1)%end, num_dims + 1))
        end if
    end subroutine s_initialize_surface_tension_module

    subroutine s_compute_capillary_source_flux( &
        vSrc_rsx_vf, vSrc_rsy_vf, vSrc_rsz_vf, &
        flux_src_vf, &
        id, isx, isy, isz)

        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsx_vf
        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsy_vf
        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsz_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_src_vf
        integer, intent(in) :: id
        type(int_bounds_info), intent(in) :: isx, isy, isz

        call s_apply_capillary_source_flux( &
            vSrc_rsx_vf, vSrc_rsy_vf, vSrc_rsz_vf, &
            flux_src_vf, id, isx, isy, isz, &
            gL_x, gR_x, gL_y, gR_y, gL_z, gR_z, c_divs, sigma)

        if (has_secondary_color_function) then
            call s_apply_capillary_source_flux( &
                vSrc_rsx_vf, vSrc_rsy_vf, vSrc_rsz_vf, &
                flux_src_vf, id, isx, isy, isz, &
                g2L_x, g2R_x, g2L_y, g2R_y, g2L_z, g2R_z, c2_divs, sigma_2)
        end if

    end subroutine s_compute_capillary_source_flux

    subroutine s_apply_capillary_source_flux( &
        vSrc_rsx_vf, vSrc_rsy_vf, vSrc_rsz_vf, &
        flux_src_vf, &
        id, isx, isy, isz, &
        gLocL_x, gLocR_x, gLocL_y, gLocR_y, gLocL_z, gLocR_z, c_fields, sigma_local)

        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsx_vf
        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsy_vf
        real(wp), dimension(-1:, 0:, 0:, 1:), intent(in) :: vSrc_rsz_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_src_vf
        integer, intent(in) :: id
        type(int_bounds_info), intent(in) :: isx, isy, isz
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(in) :: gLocL_x, gLocR_x, gLocL_y, gLocR_y
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(in) :: gLocL_z, gLocR_z
        type(scalar_field), dimension(:), intent(in) :: c_fields
        real(wp), intent(in) :: sigma_local

        real(wp), dimension(num_dims, num_dims) :: Omega
        real(wp) :: w1L, w1R, w2L, w2R, w3L, w3R, w1, w2, w3
        real(wp) :: normWL, normWR, normW
        real(wp) :: sigma_capillary
        integer :: j, k, l, i

        sigma_capillary = sigma_local

        if (id == 1) then
            #:call GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, sigma_capillary]')
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            w1L = gLocL_x(j, k, l, 1)
                            w2L = gLocL_x(j, k, l, 2)
                            w3L = 0._wp
                            if (p > 0) w3L = gLocL_x(j, k, l, 3)

                            w1R = gLocR_x(j + 1, k, l, 1)
                            w2R = gLocR_x(j + 1, k, l, 2)
                            w3R = 0._wp
                            if (p > 0) w3R = gLocR_x(j + 1, k, l, 3)

                            normWL = gLocL_x(j, k, l, num_dims + 1)
                            normWR = gLocR_x(j + 1, k, l, num_dims + 1)

                            w1 = (w1L + w1R)/2._wp
                            w2 = (w2L + w2R)/2._wp
                            w3 = (w3L + w3R)/2._wp
                            normW = (normWL + normWR)/2._wp

                            if (normW > capillary_cutoff) then
                                @:compute_capillary_stress_tensor()

                                do i = 1, num_dims

                                    flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                        flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(1, i)

                                    flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                     Omega(1, i)*vSrc_rsx_vf(j, k, l, i)

                                end do

                                flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                 sigma_local*c_fields(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 1)
                            end if
                        end do
                    end do
                end do
            #:endcall GPU_PARALLEL_LOOP
        end if

        if (id == 2) then
            #:call GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, sigma_capillary]')
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            w1L = gLocL_y(k, j, l, 1)
                            w2L = gLocL_y(k, j, l, 2)
                            w3L = 0._wp
                            if (p > 0) w3L = gLocL_y(k, j, l, 3)

                            w1R = gLocR_y(k + 1, j, l, 1)
                            w2R = gLocR_y(k + 1, j, l, 2)
                            w3R = 0._wp
                            if (p > 0) w3R = gLocR_y(k + 1, j, l, 3)

                            normWL = gLocL_y(k, j, l, num_dims + 1)
                            normWR = gLocR_y(k + 1, j, l, num_dims + 1)

                            w1 = (w1L + w1R)/2._wp
                            w2 = (w2L + w2R)/2._wp
                            w3 = (w3L + w3R)/2._wp
                            normW = (normWL + normWR)/2._wp

                            if (normW > capillary_cutoff) then
                                @:compute_capillary_stress_tensor()

                                do i = 1, num_dims

                                    flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                        flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(2, i)

                                    flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                     Omega(2, i)*vSrc_rsy_vf(k, j, l, i)

                                end do

                                flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                 sigma_local*c_fields(num_dims + 1)%sf(j, k, l)*vSrc_rsy_vf(k, j, l, 2)
                            end if
                        end do
                    end do
                end do
            #:endcall GPU_PARALLEL_LOOP
        end if

        if (id == 3) then
            if (p > 0) then
                #:call GPU_PARALLEL_LOOP(collapse=3, private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, sigma_capillary]')
                    do l = isz%beg, isz%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end

                                w1L = gLocL_z(l, k, j, 1)
                                w2L = gLocL_z(l, k, j, 2)
                                w3L = 0._wp
                                if (p > 0) w3L = gLocL_z(l, k, j, 3)

                                w1R = gLocR_z(l + 1, k, j, 1)
                                w2R = gLocR_z(l + 1, k, j, 2)
                                w3R = 0._wp
                                if (p > 0) w3R = gLocR_z(l + 1, k, j, 3)

                                normWL = gLocL_z(l, k, j, num_dims + 1)
                                normWR = gLocR_z(l + 1, k, j, num_dims + 1)

                                w1 = (w1L + w1R)/2._wp
                                w2 = (w2L + w2R)/2._wp
                                w3 = (w3L + w3R)/2._wp
                                normW = (normWL + normWR)/2._wp

                                if (normW > capillary_cutoff) then
                                    @:compute_capillary_stress_tensor()

                                    do i = 1, num_dims

                                        flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                            flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(3, i)

                                        flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                         Omega(3, i)*vSrc_rsz_vf(l, k, j, i)

                                    end do

                                    flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                                                                     sigma_local*c_fields(num_dims + 1)%sf(j, k, l)*vSrc_rsz_vf(l, k, j, 3)
                                end if
                            end do
                        end do
                    end do
                #:endcall GPU_PARALLEL_LOOP
            end if
        end if

    end subroutine s_apply_capillary_source_flux

    impure subroutine s_get_capillary(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type

        type(int_bounds_info) :: isx, isy, isz
        integer :: j, k, l, i

        isx%beg = -1; isy%beg = 0; isz%beg = 0

        if (m > 0) isy%beg = -1; if (p > 0) isz%beg = -1

        isx%end = m; isy%end = n; isz%end = p

        ! compute gradient components
        #:call GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs(1)%sf(j, k, l) = 1._wp/(x_cc(j + 1) - x_cc(j - 1))* &
                                                (f_color_primary(q_prim_vf(c_idx)%sf(j + 1, k, l)) - &
                                                 f_color_primary(q_prim_vf(c_idx)%sf(j - 1, k, l)))
                    end do
                end do
            end do
        #:endcall GPU_PARALLEL_LOOP

        if (has_secondary_color_function) then
            #:call GPU_PARALLEL_LOOP(collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            c2_divs(1)%sf(j, k, l) = 1._wp/(x_cc(j + 1) - x_cc(j - 1))* &
                                                     (f_color_secondary(q_prim_vf(c_idx)%sf(j + 1, k, l)) - &
                                                      f_color_secondary(q_prim_vf(c_idx)%sf(j - 1, k, l)))
                        end do
                    end do
                end do
            #:endcall GPU_PARALLEL_LOOP

            #:call GPU_PARALLEL_LOOP(collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            c2_divs(2)%sf(j, k, l) = 1._wp/(y_cc(k + 1) - y_cc(k - 1))* &
                                                     (f_color_secondary(q_prim_vf(c_idx)%sf(j, k + 1, l)) - &
                                                      f_color_secondary(q_prim_vf(c_idx)%sf(j, k - 1, l)))
                        end do
                    end do
                end do
            #:endcall GPU_PARALLEL_LOOP

            if (p > 0) then
                #:call GPU_PARALLEL_LOOP(collapse=3)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                c2_divs(3)%sf(j, k, l) = 1._wp/(z_cc(l + 1) - z_cc(l - 1))* &
                                                         (f_color_secondary(q_prim_vf(c_idx)%sf(j, k, l + 1)) - &
                                                          f_color_secondary(q_prim_vf(c_idx)%sf(j, k, l - 1)))
                            end do
                        end do
                    end do
                #:endcall GPU_PARALLEL_LOOP
            end if
        end if

        #:call GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs(2)%sf(j, k, l) = 1._wp/(y_cc(k + 1) - y_cc(k - 1))* &
                                                (f_color_primary(q_prim_vf(c_idx)%sf(j, k + 1, l)) - &
                                                 f_color_primary(q_prim_vf(c_idx)%sf(j, k - 1, l)))
                    end do
                end do
            end do
        #:endcall GPU_PARALLEL_LOOP

        if (p > 0) then
            #:call GPU_PARALLEL_LOOP(collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            c_divs(3)%sf(j, k, l) = 1._wp/(z_cc(l + 1) - z_cc(l - 1))* &
                                                    (f_color_primary(q_prim_vf(c_idx)%sf(j, k, l + 1)) - &
                                                     f_color_primary(q_prim_vf(c_idx)%sf(j, k, l - 1)))
                        end do
                    end do
                end do
            #:endcall GPU_PARALLEL_LOOP
        end if

        #:call GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs(num_dims + 1)%sf(j, k, l) = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            c_divs(num_dims + 1)%sf(j, k, l) = &
                                c_divs(num_dims + 1)%sf(j, k, l) + &
                                c_divs(i)%sf(j, k, l)**2._wp
                        end do
                        !c_divs(num_dims + 1)%sf(j, k, l) = &
                        !sqrt(c_divs(num_dims + 1)%sf(j, k, l))
                        c_divs(num_dims + 1)%sf(j, k, l) = &
                            sqrt(real(c_divs(num_dims + 1)%sf(j, k, l), kind=wp))
                    end do
                end do
            end do
        #:endcall GPU_PARALLEL_LOOP

        if (has_secondary_color_function) then
            #:call GPU_PARALLEL_LOOP(collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            c2_divs(num_dims + 1)%sf(j, k, l) = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                c2_divs(num_dims + 1)%sf(j, k, l) = &
                                    c2_divs(num_dims + 1)%sf(j, k, l) + &
                                    c2_divs(i)%sf(j, k, l)**2._wp
                            end do
                            c2_divs(num_dims + 1)%sf(j, k, l) = &
                                sqrt(real(c2_divs(num_dims + 1)%sf(j, k, l), kind=wp))
                        end do
                    end do
                end do
            #:endcall GPU_PARALLEL_LOOP
        end if

        call s_populate_capillary_buffers(c_divs, bc_type)

        if (has_secondary_color_function) then
            call s_populate_capillary_buffers(c2_divs, bc_type)
        end if

        iv%beg = 1; iv%end = num_dims + 1

        ! reconstruct gradient components at cell boundaries
        do i = 1, num_dims
            call s_reconstruct_cell_boundary_values_capillary(c_divs, gL_x, gL_y, gL_z, gR_x, gR_y, gR_z, i)
        end do

        if (has_secondary_color_function) then
            do i = 1, num_dims
                call s_reconstruct_cell_boundary_values_capillary(c2_divs, g2L_x, g2L_y, g2L_z, g2R_x, g2R_y, g2R_z, i)
            end do
        end if

    end subroutine s_get_capillary

    subroutine s_reconstruct_cell_boundary_values_capillary(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                            norm_dir)
        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(out) :: vL_x, vL_y, vL_z
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(out) :: vR_x, vR_y, vR_z
        integer, intent(in) :: norm_dir

        integer :: recon_dir !< Coordinate direction of the reconstruction

        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno', 'WENO_TYPE'),('muscl', 'MUSCL_TYPE')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction

                if (norm_dir == 1) then
                    is1 = idwbuff(1); is2 = idwbuff(2); is3 = idwbuff(3)
                    recon_dir = 1; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn

                elseif (norm_dir == 2) then
                    is1 = idwbuff(2); is2 = idwbuff(1); is3 = idwbuff(3)
                    recon_dir = 2; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn

                else
                    is1 = idwbuff(3); is2 = idwbuff(2); is3 = idwbuff(1)
                    recon_dir = 3; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn

                end if

                $:GPU_UPDATE(device='[is1,is2,is3,iv]')

                if (recon_dir == 1) then
                    #:call GPU_PARALLEL_LOOP(collapse=4)
                        do i = iv%beg, iv%end
                            do l = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    do j = is1%beg, is1%end
                                        vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                        vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                    end do
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                else if (recon_dir == 2) then
                    #:call GPU_PARALLEL_LOOP(collapse=4)
                        do i = iv%beg, iv%end
                            do l = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    do j = is1%beg, is1%end
                                        vL_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                        vR_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                    end do
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                else if (recon_dir == 3) then
                    #:call GPU_PARALLEL_LOOP(collapse=4)
                        do i = iv%beg, iv%end
                            do l = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    do j = is1%beg, is1%end
                                        vL_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                        vR_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                    end do
                                end do
                            end do
                        end do
                    #:endcall GPU_PARALLEL_LOOP
                end if
            end if
        #:endfor

    end subroutine s_reconstruct_cell_boundary_values_capillary

    impure subroutine s_finalize_surface_tension_module
        integer :: j

        do j = 1, num_dims
            @:DEALLOCATE(c_divs(j)%sf)
        end do

        @:DEALLOCATE(c_divs)

        @:DEALLOCATE(gL_x, gR_x)

        @:DEALLOCATE(gL_y, gR_y)
        @:DEALLOCATE(gL_z, gR_z)

        if (has_secondary_color_function) then
            do j = 1, num_dims + 1
                @:DEALLOCATE(c2_divs(j)%sf)
            end do

            @:DEALLOCATE(c2_divs)

            @:DEALLOCATE(g2L_x, g2R_x)
            @:DEALLOCATE(g2L_y, g2R_y)
            @:DEALLOCATE(g2L_z, g2R_z)
        end if

    end subroutine s_finalize_surface_tension_module

    pure elemental function f_color_primary(value) result(color_value)
        real(wp), intent(in) :: value
        real(wp) :: color_value

        color_value = min(1._wp, max(0._wp, value))
    end function f_color_primary

    pure elemental function f_color_secondary(value) result(color_value)
        real(wp), intent(in) :: value
        real(wp) :: color_value
        real(wp) :: clamped

        clamped = min(2._wp, max(0._wp, value))
        color_value = min(1._wp, max(0._wp, clamped - 1._wp))
    end function f_color_secondary

end module m_surface_tension
