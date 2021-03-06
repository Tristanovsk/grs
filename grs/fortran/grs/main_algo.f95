subroutine main_algo(npix, nband, naot, &
        &                    vza, sza, azi, rtoa, mask, wl, &
        &                    pressure_corr, aotlut, rlut_f, rlut_c, Cext_f, Cext_c, &
        &                    rg_ratio, F0, rot, aot, aot550, fine_coef, nodata, rrs, &
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

    integer, parameter :: dp = kind(0.d0)
    integer, parameter :: sp = kind(0.)
    integer, parameter :: rtype = sp

    integer, intent(in) :: npix, nband, naot
    integer, dimension(npix), intent(in) :: mask
    real(rtype), dimension(npix), intent(in) :: sza, aot550, pressure_corr
    real(rtype), dimension(nband), intent(in) :: wl, aot, rot, rg_ratio, F0
    real(rtype), dimension(nband, npix), intent(in) :: vza, azi
    real(rtype), dimension(nband, npix), intent(in) :: rtoa
    real(rtype), dimension(naot), intent(in) :: aotlut
    real(rtype), dimension(naot, nband, npix), intent(in) :: rlut_f, rlut_c
    real(rtype), dimension(nband), intent(in) :: Cext_f, Cext_c
    real(rtype), intent(in) :: fine_coef
    real(rtype), intent(inout) :: nodata
    logical, intent(in) :: rrs

    real(rtype), dimension(nband, npix), intent(out) :: rcorr, rcorrg
    real(rtype), dimension(npix), intent(out) :: aot550_est, brdf_est

    !f2py intent(in) npix,nband,naot,vza,sza,azi,rtoa,mask,wl,pressure_corr,rlut_f,rlut_c,Cext_f,Cext_c,rg_ratio,F0,rot,fine_coef,rrs
    !f2py intent(inout) aot, aot550, nodata
    !f2py intent(out) rcorr, rcorrg, aot550_est, brdf_est
    !f2py depend(npix) sza, aot550, mask, aot550_est, brdf_est, pressure_corr
    !f2py depend(nband) wl,aot,rot,Cext_f,Cext_c,rg_ratio, F0
    !f2py depend(nband,npix) vza, azi, rtoa, rcorr, rcorrg
    !f2py depend(naot) aotlut
    !f2py depend(naot,nband,npix) rlut_f,rlut_c

    integer :: iband, ipix, i, success
    real(rtype), dimension(nband) :: rsim, rsimf, rsimc, tud, brdf, aot_, rot_corr
    real(rtype), dimension(npix) :: mu0
    real(rtype), dimension(nband, npix) :: muv, m
    real(rtype) :: rglint, tdiff_Ed, tdiff_Lu
    real(rtype) :: scale

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
    !----------------------------------------------
    real(rtype), parameter :: pi = 4 * atan (1.0_rtype)
    real(rtype), parameter :: degrad = pi / 180._rtype

    !----------------------------------------------
    ! initialize outputs
    rcorr = nodata
    rcorrg = nodata
    brdf_est = nodata
    aot550_est = nodata
    !----------------------------------------------

    scale = 0.975

    do ipix = 1, npix
        ! do not process masked pixels
        if (mask(ipix) .ne. 0) cycle

        ! Land filter
        !if (rtoa(ipix,4) < rtoa(ipix,8) .and. rtoa(ipix,8) > 0.15) cycle
        mu0(ipix) = cos(sza(ipix) * degrad)

        aotpt = aot550(ipix)

        ! correction for pressure level
        rot_corr = pressure_corr(ipix) * rot

        i = 0
        success = 0
        do
            !TODO generate lut for AOT 0.0001 (or 0), now lower limit is 0.01
            !TODO  and for AOT > 0.8
            aotpt(:) = max(aotpt * scale, 0.01)
            aotpt(:) = min(aotpt * scale, 0.8)
            do iband = nband, 1, -1
                muv(iband, ipix) = cos(vza(iband, ipix) * degrad)
                call rgrd1(naot, aotlut, pressure_corr(ipix) * rlut_f(:, iband, ipix),&
                        &  maot, aotpt, rsimf(iband), intpol, w, l_w, iw, l_iw, ier)
                call rgrd1(naot, aotlut, pressure_corr(ipix) * rlut_c(:, iband, ipix),&
                        &  maot, aotpt, rsimc(iband), intpol, w, l_w, iw, l_iw, ier)
                if(ier/=0)then
                    print*, 'FATAL ERROR IN LUT INTERPOLATION IN SOLVER MODULE'
                    print*, 'ERROR = ', ier
                    print*, 'AOT', aotpt,max(aotpt * scale, 0.01)
                    stop
                endif
                !write(*,*)iband,aotpt, rtoa(ipix,iband),rsimc(iband)

                rsim(iband) = fine_coef * rsimf(iband) + (1 - fine_coef) * rsimc(iband)
                rcorrg(iband, ipix) = rtoa(iband, ipix) - rsim(iband)

                ! if negative values decrease aot
                if(rcorrg(iband, ipix).lt.0. .and. iband .le. nband - 2 .and. i .le. 8)then
                    aotpt(:) = max(aotpt * scale, 0.01)

                    i = i + 1
                    exit
                else
                    m(iband, ipix) = 1. / mu0(ipix) + 1. / muv(iband, ipix)
                    tud(iband) = exp(-(rot_corr(iband) + aot(iband)) * m(iband, ipix))
                    brdf(iband) = max(rcorrg(iband, ipix) / tud(iband), 0.)
                end if

                if(iband==1)success = 1
            enddo

            if(success==1)then
                ! rescale aot values
                aot_ = aot * aotpt(1) / aot550(ipix)
                aot550_est(ipix) = aotpt(1)
                brdf_est(ipix) = brdf(nband - 2)
                exit
            endif
        enddo

        do iband = 1, nband
            tdiff_Ed = exp(-(0.52 * rot_corr(iband) + 0.16 * aot_(iband)) * (1. / mu0(ipix)))
            tdiff_Lu = exp(-(0.52 * rot_corr(iband) + 0.16 * aot_(iband)) * (1. / muv(iband, ipix)))

            rglint = 0.5 * (tud(iband) * rg_ratio(iband) * brdf(nband) + tud(iband) * rg_ratio(iband) &
                   & / rg_ratio(nband - 1) * brdf(nband - 1))
            rcorr(iband, ipix) = rcorrg(iband, ipix) - rglint
            rcorr(iband, ipix) = rcorr(iband, ipix) / pi / tdiff_Lu / tdiff_Ed
            rcorrg(iband, ipix) = rcorrg(iband, ipix) / pi / tdiff_Lu / tdiff_Ed
            if (.not. rrs) then
                rcorr(iband, ipix) = rcorr(iband, ipix) * F0(iband)
                rcorrg(iband, ipix) = rcorrg(iband, ipix) * F0(iband)
            endif
        enddo

    enddo

    return

end subroutine main_algo