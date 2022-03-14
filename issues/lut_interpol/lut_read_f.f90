subroutine read_lut(naot, nband, nsza, nazi, nvza, &
        &                    aotlut, szalut, azilut, vzalut, &
        &                    rlut)
    ! f2py -c -m read_lut lut_read_f.f90


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

    integer, intent(in) :: nband, naot, nsza, nazi, nvza
    real(rtype), dimension(naot), intent(in) :: aotlut
    real(rtype), dimension(nsza), intent(in) :: szalut
    real(rtype), dimension(nazi), intent(in) :: azilut
    real(rtype), dimension(nvza), intent(in) :: vzalut

    real(rtype), dimension(naot, nband, nsza, nazi, nvza), intent(in) :: rlut

    !f2py intent(in) nband,naot,nsza,nazi,nvza,aotlut,szalut, razilut, vzalut
    !f2py intent(in) rlut
    !f2py depend(naot) aotlut
    !f2py depend(nsza) szalut
    !f2py depend(nazi) azilut
    !f2py depend(nvza) vzalut
    !f2py depend(naot,nband,nsza,nazi,nvza) rlut

    integer ::  i, iband

    !----------------------------------------------
    real(rtype), parameter :: pi = 4 * atan (1.0_rtype)
    real(rtype), parameter :: degrad = pi / 180._rtype


    iband=1
    do i =1,naot
        print*,aotlut(i),rlut(i,iband,1,1,1)
    end do

    return

end subroutine read_lut