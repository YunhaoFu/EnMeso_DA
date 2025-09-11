module precision

  implicit none
  private
  public :: r8, r4  
  integer, parameter:: r8 = selected_real_kind(8)
  integer, parameter:: r4 = selected_real_kind(4)

end module precision