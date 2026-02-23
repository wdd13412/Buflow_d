! 主程序：通过use模块获取main的显式接口
program RUN_MAIN_NODIFF
    use MAIN_MODULE_DIFF  ! 关键：使用包含main的模块，自动获取显式接口
    implicit none
    real(kind=8) :: input_param, input_paramd  ! 输入参数
    real(kind=8), allocatable :: cellPrimitivesout(:,:)  ! 输出结果
    real(kind=8), allocatable :: cellprimitivesoutd(:,:)  ! 输出结果
    INTRINSIC ALLOCATED
  	INTRINSIC SIZE

    ! 设置输入参数（示例值）
    input_param = -0.75d0  
    input_paramd = 1.0d0    
    ! 调用main子程序（此时编译器已知晓其接口，支持可分配数组参数）
    call MAIN_D(input_param, input_paramd, cellprimitivesout, cellprimitivesoutd)

    ! 输出结果信息（可选）
    if (allocated(cellprimitivesout)) then
        deallocate(cellprimitivesout)  ! 释放内存
    end if
    if (allocated(cellprimitivesoutd)) then
        deallocate(cellprimitivesoutd)  ! 释放内存
    end if

end program RUN_MAIN_NODIFF
