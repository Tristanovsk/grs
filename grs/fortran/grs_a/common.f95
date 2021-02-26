module common

integer, parameter:: rdp=kind(0.d0)
integer, parameter:: rsp=kind(0.)
integer, parameter:: rtype=rdp

! set polynomial order for lut fitting
integer,parameter :: norder=4

end module