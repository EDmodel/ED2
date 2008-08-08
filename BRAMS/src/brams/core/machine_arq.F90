! Module necessary to set the better parameters to IA32 or SX-6 machine

module machine_arq

  implicit none

  ! Where:
  ! machine=0 => Generical IA32 machine
  ! machine=1 => SX-6 Vetorial machine

#if defined(NEC_SX)
  integer, parameter :: machine=1
#endif

#if defined(PC_LINUX1)
  integer, parameter :: machine=0
#endif

#if defined(IBM)
  integer, parameter :: machine=0
#endif

end module machine_arq
