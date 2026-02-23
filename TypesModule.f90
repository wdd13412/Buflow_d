module TypesModule
  implicit none
  ! 统一定义64位整数和实数常量
  integer, parameter :: int64 = 8       ! 8字节整数:real64可精确存储 15-17 位十进制有效数字,指数范围是 10⁻³⁰⁸~10³⁰⁸.因此完全可存储4.06575815E-20
  integer, parameter :: real64 = 8      ! 8字节实数（双精度）
end module TypesModule
