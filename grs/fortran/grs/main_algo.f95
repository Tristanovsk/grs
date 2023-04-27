module grs
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use, intrinsic :: iso_fortran_env, only: real32

    integer, parameter :: dp = kind(0.d0)
    integer, parameter :: sp = kind(0.)
    integer, parameter :: rtype = sp


contains

    subroutine main_algo(nx, ny, naot, nband, nsza, nazi, nvza, &
            &                    aotlut, szalut, azilut, vzalut, &
            &                    rlut_f, rlut_c, Cext_f, Cext_c, &
            &                    vza, sza, azi, rtoa, mask, wl, &
            &                    pressure_corr, &
            &                    rg_ratio, F0, rot, aot_tot, aot_sca, aot550, fine_coef, rrs, &
            &                    rcorr, rcorrg, aot550_est, brdf_est)
        !f2py -c -m main_algo main_algo.f95

        !--------------------------------------------------------------------
        !   nband: number of satellite bands to be processed
        !   naot, nsza, nazi, nvza : dimensions of lut data arrays
        !   vza, sza, azi: satellite angle data
        !   rtoa: satellite L1C reflectance (corected from gas absorption)
        !   rg_ratio: spectral variation (referred to as epsilon in Harmel et al., RSE; 201?)
        !   TODO
        !--------------------------------------------------------------------

        implicit none

        integer, intent(in) :: nx, ny, nband, naot, nsza, nazi, nvza
        integer, dimension(nx, ny), intent(in) :: mask
        real(rtype), dimension(nx, ny), intent(in) :: sza, aot550, pressure_corr, fine_coef
        real(rtype), dimension(nband), intent(in) :: wl, rot, rg_ratio, F0
        real(rtype), dimension(nband, nx, ny), intent(in) :: vza, azi
        real(rtype), dimension(nband, nx, ny), intent(in) :: rtoa, aot_tot, aot_sca
        real(rtype), dimension(naot), intent(in) :: aotlut
        real(rtype), dimension(nsza), intent(in) :: szalut
        real(rtype), dimension(nazi), intent(in) :: azilut
        real(rtype), dimension(nvza), intent(in) :: vzalut

        real(rtype), dimension(naot, nband, nsza, nazi, nvza), intent(in) :: rlut_f, rlut_c
        real(rtype), dimension(nband), intent(in) :: Cext_f, Cext_c
        logical, intent(in) :: rrs

        real(rtype), dimension(nband, nx, ny), intent(out) :: rcorr, rcorrg
        real(rtype), dimension(nx, ny), intent(out) :: aot550_est, brdf_est

        !f2py intent(in) nx, ny,nband,naot,nsza,nazi,nvza,aotlut,szalut, razilut, vzalut
        !f2py intent(in) vza,sza,azi,rtoa,mask,wl,pressure_corr,rlut_f,rlut_c
        !f2py intent(in) Cext_f,Cext_c,rg_ratio,F0, aot_sca,rot,fine_coef,rrs
        !f2py intent(inout) aot_tot, aot550
        !f2py intent(out) rcorr, rcorrg, aot550_est, brdf_est
        !f2py depend(nx, ny) sza, aot550, mask, aot550_est, brdf_est, pressure_corr,fine_coef
        !f2py depend(nband) wl,rot,Cext_f,Cext_c,rg_ratio, F0
        !f2py depend(nband,nx, ny) vza, azi, aot_tot, aot_sca, rtoa, rcorr, rcorrg
        !f2py depend(naot) aotlut
        !f2py depend(nsza) szalut
        !f2py depend(nazi) azilut
        !f2py depend(nvza) vzalut
        !f2py depend(naot,nband,nsza,nazi,nvza) rlut_f,rlut_c

        integer :: iband, ix, iy, i, success
        real(rtype), dimension(nband) :: rsim, rsimf, rsimc, tud, brdf, aot_, rot_corr, ssa_tot
        real(rtype), dimension(nx, ny) :: mu0
        real(rtype), dimension(nband, nx, ny) :: muv, m
        real(rtype) :: rglint, tdiff_Ed, tdiff_Lu
        real(rtype) :: scale,pressure_corr_pix

        !----------------------------------------------
        ! set parameters for lut interpolation
        ! interpolation scheme (1: linear; 3 cubic interpolation)
        integer :: ier
        integer, parameter :: mint = 1, l_w = 16 * (15 + 100 + 100 + 100), & !(naot)
                &                    l_iw = 2 * (15 + 100 + 100 + 100)  !(naot)
        integer, dimension(mint), parameter :: intpol = 1
        integer, dimension(l_iw) :: iw
        real*8, dimension(l_w) :: w
        integer, parameter :: maot = 1
        real(rtype), dimension(maot) :: aotpt
        integer :: isza
        integer, dimension(nband) :: ivza, iazi
        !----------------------------------------------
        real(rtype), parameter :: pi = 4 * atan (1.0_rtype)
        real(rtype), parameter :: degrad = pi / 180._rtype

        real(real32) :: nan
        nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

        !----------------------------------------------
        ! initialize outputs
        rcorr = nan
        rcorrg = nan
        brdf_est = nan
        aot550_est = nan
        !----------------------------------------------

        scale = 0.95
        do ix = 1, nx
          do iy = 1,ny
            ! if input pixel is NaN pass to next pixel
            if (rtoa(3,ix, iy) /= rtoa(3,ix, iy)) cycle
            ! scale = (0.01/aot550(ix, iy))**(1d0/8)*1.02
            ! do not process masked pixels
            if (mask(ix, iy) .ne. 0) cycle
            pressure_corr_pix = pressure_corr(ix, iy)

            mu0(ix, iy) = cos(sza(ix, iy) * degrad)

            !---------------------------
            ! get SZA, VZA and AZI
            ! index to get data from LUT
            !---------------------------
            isza = iloc(szalut, nsza, sza(ix, iy))
            do iband = 1, nband
                ivza(iband) = iloc(vzalut, nvza, vza(iband, ix, iy))
                iazi(iband) = iloc(azilut, nazi, azi(iband, ix, iy))
            end do

            aotpt = aot550(ix, iy)
            ! print*,aotpt,fine_coef(ix, iy)
            ! correction for pressure level

            aot_ = aot_tot(:, ix, iy)
            i = 0
            success = 0
            do
                rot_corr = pressure_corr_pix * rot
                !aotpt(:) = aotpt * scale
                aotpt(:) = max(aotpt, 0.01)
                aotpt(:) = min(aotpt, 0.8)
                ! print*,'grs algo',aot_sca(:,ix, iy)
                ! print*,'grs algo aottot *rot',aot_sca(:,ix, iy),aot_sca(:,ix, iy)+rot_corr
                ssa_tot = (aot_sca(:, ix, iy) + rot_corr) / (aot_(:) + rot_corr)
                do iband = nband, 1, -1
                    !TODO generate lut for AOT 0.0001 (or 0), now lower limit is 0.01
                    !TODO  and for AOT > 0.8

                    muv(iband, ix, iy) = cos(vza(iband, ix, iy) * degrad)
                    call rgrd1(naot, aotlut, pressure_corr_pix * rlut_f(:, iband, isza, iazi(iband), ivza(iband)), &
                            &  maot, aotpt, rsimf(iband), intpol, w, l_w, iw, l_iw, ier)
                    call rgrd1(naot, aotlut, pressure_corr_pix * rlut_c(:, iband, isza, iazi(iband), ivza(iband)), &
                            &  maot, aotpt, rsimc(iband), intpol, w, l_w, iw, l_iw,ier)
                    !TODO understand why ndarray is inverted on the first dim aot
!                    if (iband == 1) then
!                        do i =1,naot
!                            print*,aotlut(i),pressure_corr(ix, iy) * rlut_c(i, iband, isza, iazi(iband), ivza(iband))
!                        end do
!                        do i =1,naot
!                            print*,'fine',aotlut(i),pressure_corr(ix, iy) * rlut_f(i, iband, isza, iazi(iband), ivza(iband))
!                        end do
!                    endif
                    if(ier/=0)then
                        print*, 'FATAL ERROR IN LUT INTERPOLATION IN SOLVER MODULE'
                        print*, 'ERROR = ', ier
                        print*, 'AOT', aotpt, max(aotpt * scale, 0.01)
                        stop
                    endif

                    rsim(iband) = fine_coef(ix, iy) * rsimf(iband) + (1 - fine_coef(ix, iy)) * rsimc(iband)
                    ! correction for absorbing aerosol
                    rsim(iband) = ssa_tot(iband) * rsim(iband)
                    rcorrg(iband, ix, iy) = rtoa(iband, ix, iy) - rsim(iband)

                    ! if negative values decrease aot
                    ! start after SWIR 2.2 microns (nband)

                    if(aotpt(1) .gt. 0.01 .and. rcorrg(iband, ix, iy) .lt. -0.0004 &
                            &  .and. i .le. 8 .and. iband .lt. nband) then !iband .le. nband - 2 .and.
                        ! print*,'aot adjustment',aotpt ,aotpt * scale,scale
                        aotpt(:) = max(aotpt * scale, 0.01)
                        !pressure_corr_pix=pressure_corr_pix*0.99
                        i = i + 1
                        exit
                    else
                        m(iband, ix, iy) = 1. / mu0(ix, iy) + 1. / muv(iband, ix, iy)
                        tud(iband) = exp(-(rot_corr(iband) + aot_(iband)) * m(iband, ix, iy))
                        brdf(iband) = max(rcorrg(iband, ix, iy) / tud(iband), 0.)
                    end if

                    if(iband==1)success = 1
                enddo

                if(success==1)then
                    ! rescale aot values
                    aot_ = aot_tot(:, ix, iy) * aotpt(1) / aot550(ix, iy)
                    aot550_est(ix, iy) = aotpt(1)
                    brdf_est(ix, iy) = brdf(nband - 1)
                    exit
                endif
            enddo

            do iband = 1, nband
                tdiff_Ed = exp(-(0.52 * rot_corr(iband) + 0.16 * aot_(iband)) * (1. / mu0(ix, iy)))
                tdiff_Lu = exp(-(0.52 * rot_corr(iband) + 0.16 * aot_(iband)) * (1. / muv(iband, ix, iy)))

                rglint = (tud(iband) * rg_ratio(iband) * 0.1*brdf(nband) + tud(iband) * rg_ratio(iband) &
                        & / rg_ratio(nband - 1) * 0.9*brdf(nband - 1))
                rcorr(iband, ix, iy) = rcorrg(iband, ix, iy) - rglint
                rcorr(iband, ix, iy) = rcorr(iband, ix, iy) / pi / tdiff_Lu / tdiff_Ed
                rcorrg(iband, ix, iy) = rcorrg(iband, ix, iy) / pi / tdiff_Lu / tdiff_Ed
                if (.not. rrs) then
                    rcorr(iband, ix, iy) = rcorr(iband, ix, iy) * F0(iband)
                    rcorrg(iband, ix, iy) = rcorrg(iband, ix, iy) * F0(iband)
                endif
                !rcorr(iband, ix, iy) = rtoa(iband, ix, iy)
                !rcorrg(iband, ix, iy)=rsim(iband)
            enddo
          enddo
        enddo

        return

    end subroutine main_algo

    function iloc(arr, Narr, value)
        ! function to find array index of the nearest value within the array arr
        ! arr :: array of real
        ! Narr :: dimension of array
        ! value :: real value to find within arr
        implicit none
        integer :: iloc
        integer :: Narr, i
        real, dimension(Narr) :: arr
        real :: value, r_, r__

        r_ = huge(r_)
        do i = 1, Narr
            r__ = abs(arr(i) - value)
            if (r__ .gt. r_) then
                iloc = i - 1
                exit
            else
                r_ = r__
                iloc = i
            end if
        end do

    end function iloc

end module grs