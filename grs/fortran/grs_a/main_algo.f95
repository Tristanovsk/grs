subroutine main_algo(npix, nband, naot, &
&                    vza, sza, azi, rtoa, mask, wl, &
&                    aotlut, rlut_f, rlut_c, cextf, cextc, cextf550, cextc550, &
&                    rg_ratio, F0, rot, aot, aot550_in, nodata, &
&                    rcorr, rcorrg, aot550_est, brdf_est)
    !f2py -c -m main_algo main_algo.f95

    !--------------------------------------------------------------------
    !   nband: number of satellite bands to be processed
    !   naot, nsza, nazi, nvza : dimensions of lut data arrays
    !   vza, sza, azi: satellite angle data
    !   rtoa: satellite L1C reflectance (corected from gas absorption)
    !   rg_ratio: spectral variation (referred to as epsilon in Harmel et al., RSE; 2018)
    !   TODO
    !--------------------------------------------------------------------
    use common
    use lsq
    use solver

    implicit none

    integer, intent(in) :: npix, nband, naot
    integer,dimension(npix),intent(in) :: mask
    real(rtype), dimension(npix), intent(in) :: sza, aot550_in
    real(rtype), dimension(nband), intent(in) :: wl, rot, aot, rg_ratio, F0
    real(rtype), dimension(npix, nband), intent(in) :: vza, azi
    real(rtype), dimension(npix, nband), intent(in) :: rtoa
    real(rtype), dimension(naot), intent(in) :: aotlut
    real(rtype), dimension(naot, nband, npix), intent(in) :: rlut_f, rlut_c
    real(rtype), dimension(nband), intent(in) :: cextf, cextc
    real(rtype), intent(in) :: cextf550, cextc550!!
    real(rtype), dimension(npix, nband), intent(out) :: rcorr, rcorrg
    real(rtype), dimension(npix), intent(out) :: aot550_est, brdf_est
    real(rtype),intent(inout) :: nodata

    !f2py intent(in) npix,nband,naot,vza,sza,azi,rtoa,mask,wl,aotlu,rlut_f,rlut_c
    !f2py intent(in) cextf,cextc,cextf550,cextc550
    !f2py intent(in) aot,aot550_in,rg_ratio,F0,rot
    !f2py intent(inout) nodata
    !f2py intent(out) rcorr, rcorrg, aot550_est, brdf_est
    !f2py depend(npix) sza, mask, aot550_in, aot550_est, brdf_est
    !f2py depend(nband) wl,aot,rot,cextf,cextc,rg_ratio, F0
    !f2py depend(npix,nband) vza, azi, rtoa, rcorr, rcorrg
    !f2py depend(naot) aotlut
    !f2py depend(naot,nband,npix) rlut_f,rlut_c

    integer :: iband, ipix, i
    real(rtype), dimension(nband) :: rsimf, rsimc, tud, brdf, aot_est
    real(rtype), dimension(npix) :: mu0
    real(rtype), dimension(npix, nband) :: muv, m
    real(rtype) :: rglint, tdiff_Ed, tdiff_Lu


    !----------------------------------------------
    ! set parameters for lut interpolation of fitting
    ! interpolation scheme (1: linear; 3 cubic interpolation)
    integer :: ier
    integer, parameter :: mint = 1, l_w = 16 * (15 + 100 + 100 + 100), & !(naot)
            &                    l_iw = 2 * (15 + 100 + 100 + 100)  !(naot)
    integer, dimension(mint), parameter :: intpol = 1
    integer, dimension(l_iw) :: iw
    real(rtype), dimension(l_w) :: w
    integer, parameter :: maot = 1


    !----------------------------------------------
    ! set parameters for polynomial fitting
    real(rtype) :: R2
    real(rtype) :: weightaot(naot)
    real(rtype), dimension(nband, 0 : norder) :: af, ac
    real(rtype) :: aot550, beta, brdf_swir

    !----------------------------------------------
    real(rtype), parameter :: pi = 4 * atan (1.0_rtype)
    real(rtype), parameter :: degrad = pi / 180._rtype

    print*,'once again !'

    !--- set common variables
    lband = nband
    allocate(s_af(nband,0:norder), s_ac(nband,0:norder))
    allocate(s_cextf(nband), s_cextc(nband),s_rot(nband))
    allocate(s_m(nband),s_rg_ratio(nband),rsim(nband),rmes(nband),sigma(nband))

    s_cextf = cextf !(band)
    s_cextc = cextc !(band)
    s_cextf550 = cextf550
    s_cextc550 = cextc550
    s_rg_ratio = rg_ratio !(band)

    !----------------------------------------------
    ! initialize outputs
    rcorr = nodata
    rcorrg=nodata
    !----------------------------------------------
    weightaot = 1.


    do ipix = 1, npix
        ! do not process masked pixels
        if (mask(ipix) .ne. 0) cycle

        mu0(ipix)=cos(sza(ipix)*degrad)

        aot550=aot550_in(ipix)

        do iband = 1, nband
            muv(ipix, iband) = cos(vza(ipix, iband) * degrad)
            !----------------------------------------------
            !-- polynomial fitting onto aot grid
            !----------------------------------------------
            call fit_poly(aotlut, rlut_f(:, iband, ipix), &
                    &              weightaot, naot, norder, af(iband, 0 : norder), R2)
            if(R2 < 0.98)then
                print*, 'FATAL ERROR IN POLYNOMIAL FITTING (f(aot)) IN SOLVER MODULE'
                print*, 'fine mode R2 = ', R2, aotlut, rlut_f(:, iband, ipix),s_af(iband, 0 : norder)
            endif

            call fit_poly(aotlut, rlut_c(:, iband, ipix), &
                    &              weightaot, naot, norder, ac(iband, 0 : norder), R2)

            if(R2 < 0.98)then
                print*, 'FATAL ERROR IN POLYNOMIAL FITTING (f(aot)) IN SOLVER MODULE'
                print*, 'R2 = ', R2, aotlut, rlut_c(:, iband, ipix), ac(iband, 0 : norder)
                cycle
            endif
            rsimf = 0.
            rsimc = 0.
            do i = 0, norder
                rsimf = rsimf + af(iband, i) * aot550**i
                rsimc = rsimc + ac(iband, i) * aot550**i
            enddo
            !write(*,*)iband,rtoa(ipix,iband),rlut_f(iband,2,2,2,2)
            rsim(iband) = beta * rsimf(iband) + (1 - beta) * rsimc(iband)
            rcorr(ipix, iband) = rtoa(ipix, iband) - rsim(iband)
            m(ipix, iband) = 1. / mu0(ipix) + 1. / muv(ipix, iband)
            tud(iband) = exp(-(rot(iband) + aot(iband)) * m(ipix, iband))
            brdf(iband) = rcorr(ipix, iband) / tud(iband)
        enddo


!-----------------------------------------------
!       CALL SOLVER MODULE
!-----------------------------------------------
!--- set common variables
!        s_af = af(band, :)
!        s_ac = ac(band, :)
!        s_rot = rot(band)
!        s_m(1 : lband) = m(ipix, band)
!        rmes(1 : lband) = rtoa(ipix, band)
        s_af = af
        s_ac = ac
        s_rot = rot
        s_m = m(ipix, :)
        rmes = rtoa(ipix, :)
        sigma = 1._rtype
        !--- first guess
        brdf_swir = 1d-10 !max(1d-10, brdf(nband))

        beta = 0.5
        call aero_glint(aot550, beta, brdf_swir)

        aot550_est(ipix) = aot550
        brdf_est(ipix) = brdf_swir

        do iband = 1, nband
            aot_est(iband) = beta * cextf(iband) / s_cextf550 * aot550 + (1. - beta) * cextc(iband) / s_cextc550 * aot550
            tud(iband) = dexp(-(rot(iband) + aot_est(iband)) / m(ipix, iband))
            rglint = tud(iband) * rg_ratio(iband) * brdf_swir
            rsimf(iband) = 0.
            rsimc(iband) = 0.
            do i = 0, norder
                rsimf(iband) = rsimf(iband) + af(iband, i) * aot550**i
                rsimc(iband) = rsimc(iband) + ac(iband, i) * aot550**i
            enddo
            rsim(iband) = beta * rsimf(iband) + (1. - beta) * rsimc(iband)

            ! TODO improve calculation of the atmospheric transmittances
            tdiff_Ed = exp(-(0.52 * rot(iband) + 0.16 * aot_est(iband)) * (1. / mu0(ipix)))
            tdiff_Lu = exp(-(0.52 * rot(iband) + 0.16 * aot_est(iband)) * (1. / muv(ipix, iband)))

            rcorrg(ipix, iband) = rtoa(ipix, iband) - rsim(iband)
            rcorrg(ipix, iband) = rcorr(ipix, iband) * F0(iband) / pi / tdiff_Lu / tdiff_Ed

            rcorr(ipix, iband) = rtoa(ipix, iband) - rsim(iband) - rglint
            rcorr(ipix, iband) = rcorr(ipix, iband) * F0(iband) / pi / tdiff_Lu / tdiff_Ed

            !write(*,*)"band",iband,';  rcorr',rcorr(ipix,iband),tdiff_Lu,tdiff_Ed,aot_est(iband)
        enddo

    enddo

    deallocate(s_af,s_ac,s_cextf,s_cextc,s_rot,s_m,s_rg_ratio,rsim,rmes,sigma)

    return

end subroutine main_algo