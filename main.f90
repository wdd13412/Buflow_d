module main_module
    use BuFlowModule  
    use meshdeformationn  
    implicit none    

contains  ! 模块内包含子程序，自动生成显式接口

    subroutine main(input_param, cellPrimitivesout)   
    	!$AD &input_param          ! 标记输入
        !$AD &cellPrimitives       ! 标记输出
        ! 输入参数：仅原data_4D137的第一个维度（标量）
		real(kind=8), intent(in) :: input_param  
		real(kind=8), allocatable :: cellPrimitives(:,:)		
		! 输出参数：CFD结果（二维可分配数组）
		real(kind=8), allocatable, intent(out) :: cellPrimitivesout(:,:)  		      
        ! 局部变量：构造1×4变形参数
        real(kind=8) :: data_4D137(1,4)  
        ! 传递变形网格点的变量
        integer :: nCells, meshInfo(4) 
        real(kind=8), allocatable :: point_update(:,:)  
        character(256) :: meshPath
        type(Meshh) :: mesh
 !       !$AD intent(out) :: cellPrimitives(nCells, 5)  ! 强制 Tapenade 识别维度（nCells 已定义）[-0.75d0, 0.485d0, -0.588d0, -0.912d0]
        ! 构造变形参数（仅第一个维度为输入，其余为0）
        data_4D137(1,1) = input_param  
!        data_4D137(1,2:4) = 0.0d0  ! 后三个维度固定为0
        data_4d137(1, 2) = 0.485d0
		data_4d137(1, 3) = -0.588d0
		data_4d137(1, 4) = -0.912d0  ! 后三个维度固定为0
        
        ! 读取网格
!		meshPath = "mesh/OFairfoilMesh"	
!		mesh = OpenFOAMMesh(meshPath, point_update)
		! 修改后的代码
!		meshInfo = unstructuredMeshInfo(mesh) 
!		nCells = meshInfo(1)  
!		allocate(cellPrimitives(nCells, 5)) 
        allocate(cellPrimitives(18513, 5)) 
        allocate(cellPrimitivesout(18513, 5))
        cellPrimitives = 0.0d0  ! 显式赋值，确保 Tapenade 识别为“活跃输出”
        ! 调用网格变形，获取变形后的网格点（point_update）
        call airfoil_deformation_HH(data_4D137, point_update)
        
        ! 调用CFD计算，传入变形网格点，输出结果到cellPrimitives
        call compute_CFD(point_update, cellPrimitivesout)  

        ! 释放point_update内存
        if (allocated(point_update)) deallocate(point_update)     
    end subroutine main

end module main_module
