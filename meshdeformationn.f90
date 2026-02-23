module meshdeformationn
!	use iso_fortran_env
    use TypesModule
	implicit none
    ! 2. 再禁用默认类型
    intrinsic :: verify, findloc, system, matmul, dot_product, move_alloc, new_line, huge, int
	external :: dpotrf, dpotri  ! 仅保留外部声明    
    ! 仅 Tapenade 处理时取消注释（原代码编译时注释掉）
    !include 'tapenade_includes.f90'  ! <<< 仅 Tapenade 用，原编译时注释
    !非tapenade时call move_alloc_char(originalLines, newLines)
contains 
	subroutine airfoil_deformation_HH(data_4D137, point_update)
        real*8, allocatable :: wing_update(:,:)
		real*8, intent(in) :: data_4D137(:,:)
		real*8, allocatable, intent(inout) :: point_update(:,:) 
		real*8 :: RAE2822_up(50,2), RAE2822_low(50,2)
		real*8 :: x1(50), y_up(50), y_low(50), loca(2), xx(2), fun1_result(50), funk_result(50)
		real*8 :: S(100,4), S_up(100,2), S_low(100,2)
		real*8 :: f_up(50), f_low(50), yy_up(50), yy_low(50)
		real*8 :: avg_y_51
		real*8 :: up_part(49,2), row_51(2), low_part(49,2)
		integer :: d, n, i, k
		real*8, allocatable :: wing(:,:), inoutput(:,:)
        allocate(inoutput(99, 2), source=0.0_real64)
        allocate(wing(99, 2), source=0.0_real64)
        allocate(wing_update(99, 2), source=0.0_real64)
        ! 1. 提取输入输出边界点（初始化inoutput）
		inoutput = CSV_read("boundary_inoutput_xy.csv")
		! 从文件读取翼型数据
		wing = CSV_read("wing.csv")
		
		! 提取上表面坐标
		RAE2822_up = wing(1:50,:)
		! 提取下表面坐标
		RAE2822_low = wing(50:99,:)
		
		! 提取坐标向量
		x1 = RAE2822_up(:,1)
		y_up = RAE2822_up(:,2)
		y_low = RAE2822_low(:,2)
		
		loca = (/1.0d0, 8.0d0/)
		xx = (/x1(1), x1(8)/)
		
		! 处理变形参数
		d = size(data_4D137, 2)
		S = 0.0d0
		
		do i = 1, d
			do k = 1, size(data_4D137, 1)
			    S(k,i) = 0.05d0 * cos((1.0d0 - data_4D137(k,i))/2.0d0 * 3.1415926535d0)
			end do
		end do
		
		S_up = S(:,1:2)
		S_low = S(:,3:4)
		n = size(S_up, 1)
		
		! 计算变形
		if (n > 0) then
			! 计算变形函数
			call fun1(x1, 50, fun1_result)  ! 第二个参数 50 是 x1 的长度
			call funk(x1, 50, 2, funk_result)  ! 第二个参数 50 是 X（即 x1）的长度
			f_up = fun1_result * S_up(1,1) + funk_result * S_up(1,2)
			f_low = fun1_result * S_low(1,1) + funk_result * S_low(1,2)
			
			yy_up = y_up + f_up
			yy_low = y_low + f_low
			
			! 特殊处理重叠点
			avg_y_51 = (yy_up(50) + yy_low(1)) / 2.0d0
			
			! 构建更新后的翼型
			up_part(:,1) = x1(1:49)
			up_part(:,2) = yy_up(1:49)
			
			row_51(1) = x1(50)
			row_51(2) = avg_y_51
			
			low_part(:,1) = RAE2822_low(2:50,1)
			low_part(:,2) = yy_low(2:50)
			
			! 合并三部分
			wing_update(1:49,:) = up_part
			wing_update(50,:) = row_51
			wing_update(51:99,:) = low_part
			! 输出到文件
			call CSV_write("wing_update.csv", wing_update, 99, 2)
			print *, "翼型坐标已更新并保存至wing_update.csv"
		else
			print *, "错误：无变形参数可用"
			stop
		end if
		! 释放临时分配的内存
		if (allocated(inoutput)) deallocate(inoutput)
		if (allocated(wing)) deallocate(wing)	
		allocate(point_update(18513,3))	
		
		PRINT*, "========================================"
        PRINT*, "输出wing_update的第一个点:"
        PRINT*, "wing_update(1,1) = ", wing_update(1,1)
        PRINT*, "wing_update(1,2) = ", wing_update(1,2)
        PRINT*, "wing_update(1,3) = ", wing_update(1,3)
        PRINT*, "输出wing_updated的第一个点:"
!        PRINT*, "wing_updated(1,1) = ", wing_updated(1,1)
!        PRINT*, "wing_updated(1,2) = ", wing_updated(1,2)
!        PRINT*, "wing_updated(1,3) = ", wing_updated(1,3)
        PRINT*, "========================================"
        
		call meshdeformation(wing_update, point_update)

        PRINT*, "========================================"
!第300个点变化前:(0.942606 0.0121379 0)        
		PRINT*, "输出point_update的第300个点:"
        PRINT*, "point_update(300,1) = ", point_update(300,1)
        PRINT*, "point_update(300,2) = ", point_update(300,2)
        PRINT*, "point_update(300,3) = ", point_update(300,3)
        PRINT*, "输出point_updated的第300个点:"
!        PRINT*, "point_updated(300,1) = ", point_updated(300,1)
!        PRINT*, "point_updated(300,2) = ", point_updated(300,2)
!        PRINT*, "point_updated(300,3) = ", point_updated(300,3)
        PRINT*, "========================================"
        
		if (allocated(wing_update)) deallocate(wing_update)
		print *, "aAd"
    contains
		! 变形函数fun1
		subroutine fun1(x, n, result)  ! 新增 n：x 的长度（显式维度）
			integer, intent(in) :: n  ! 数组 x 的维度（如 50）
			real*8, intent(in) :: x(n)  ! 用 n 声明 x 的维度
			real*8, intent(out) :: result(n)  ! 用 n 声明 result 的维度（与 x 一致）
			result = x**0.25d0 * (1.0d0 - x) * exp(-20.0d0 * x)
		end subroutine fun1

		! 变形函数funk
		subroutine funk(X, n, k, result)  ! 新增 n：X 的长度（显式维度）
			integer, intent(in) :: n  ! 数组 X 的维度（如 50）
			real*8, intent(in) :: X(n)  ! 用 n 声明 X 的维度
			integer, intent(in) :: k
			real*8, intent(out) :: result(n)  ! 用 n 声明 result 的维度
			real*8 :: xxx(2), ee
			xxx = [x1(1), x1(8)]
			ee = log10(0.5d0) / log10(xxx(k))
			result = (sin(3.1415926535d0 * X**ee))**3
		end subroutine funk
	end subroutine airfoil_deformation_HH
	subroutine meshdeformation(wing_update, point_update)
		implicit none
		character(256) :: meshPath
        real*8, allocatable, intent(in) :: wing_update(:,:)
        real*8, allocatable, intent(inout) :: point_update(:,:) 
		character(len=256) :: pointsFilePath
		real(kind=8), allocatable :: points(:,:), mat_inv(:,:)
		integer :: nPoints, nWing, nInoutput, nControl, i, j, k, pt, line, iostat, len_str, lineCount
		real(kind=8), allocatable :: wing_coords(:,:), wing_update_coords(:,:), inoutput_coords(:,:), r_row(:), phi_row(:)
		real(kind=8), allocatable :: control_points(:,:), d_wing_x(:), d_wing_y(:), d_inoutput_x(:), d_inoutput_y(:)
		real(kind=8), allocatable :: d_control_x(:), d_control_y(:), d_mesh_x(:), d_mesh_y(:), d_mesh_z(:)
		real(kind=8), allocatable :: phi(:,:), weights_x(:), weights_y(:), distances(:)
		real(kind=8) :: epsilonn, diff(2), r, x, y, z, d_nomal
		character(len=256), allocatable :: originalLines(:), newDataLines(:), newLines(:)
		integer :: nLine, dataStartLine
		real*8, allocatable :: wing(:,:), inoutput(:,:)
        print *, "aAa"
        meshPath = "mesh/OFairfoilMesh"
		pointsFilePath = trim(meshPath)//"/points"
		allocate(points(18513,3), source=0.0d0)  ! 空数组初始化
		points = readOFPointsFilef(pointsFilePath)
		nPoints = size(points, 1)
        if (.not. allocated(wing)) allocate(wing(99, 2), source=0.0_real64)
		if (.not. allocated(inoutput)) allocate(inoutput(99, 2), source=0.0_real64)
		inoutput = CSV_read("boundary_inoutput_xy.csv")
		wing = CSV_read("wing.csv")		
		nWing = size(wing_update, 1)
		nInoutput = size(inoutput, 1)

		allocate(wing_coords(nWing,2), wing_update_coords(nWing,2), inoutput_coords(nInoutput,2))
		do i = 1, nWing
		    wing_coords(i,:) = wing(i,:)
		    wing_update_coords(i,:) = wing_update(i,:)
		end do
		do i = 1, nInoutput
		    inoutput_coords(i,:) = inoutput(i,:)
		end do

		nControl = nWing + nInoutput
		allocate(control_points(nControl,2))
		control_points(1:nWing,:) = wing_coords
		control_points(nWing+1:nControl,:) = inoutput_coords

		allocate(d_wing_x(nWing), d_wing_y(nWing))
		d_wing_x = wing_update_coords(:,1) - wing_coords(:,1)
		d_wing_y = wing_update_coords(:,2) - wing_coords(:,2)

		allocate(d_inoutput_x(nInoutput), d_inoutput_y(nInoutput))
		d_inoutput_x = 0.0d0; d_inoutput_y = 0.0d0

		allocate(d_control_x(nControl), d_control_y(nControl))
		d_control_x(1:nWing) = d_wing_x
		d_control_x(nWing+1:nControl) = d_inoutput_x
		d_control_y(1:nWing) = d_wing_y
		d_control_y(nWing+1:nControl) = d_inoutput_y

		d_nomal = 1.0d0/4.5d0
		epsilonn = 0.9560721680067801d0

		allocate(d_mesh_x(nPoints), d_mesh_y(nPoints), d_mesh_z(nPoints))
		d_mesh_x = 0.0d0; d_mesh_y = 0.0d0; d_mesh_z = 0.0d0
      
		allocate(r_row(nControl), phi_row(nControl))
		allocate(phi(nControl, nControl))
		do j = 1, nControl
			do k = 1, nControl
				diff = control_points(j,:) - control_points(k,:)
				r_row(k) = sqrt(sum(diff**2))
			end do
			call gaussian_rbf_array(r_row, epsilonn, phi_row)
			phi(j,:) = phi_row
		end do
        allocate(weights_x(nControl), source=0.0d0)  ! 预初始化
		allocate(weights_y(nControl), source=0.0d0)
		allocate(mat_inv(nControl, nControl), source=0.0d0)
		! 2. 调用简化版 inverse（仅传 3 个参数，无 info_out）
		call inverse(nControl, phi, mat_inv)
		! 3. 矩阵乘法（用 mat_inv 替代原 inverse(phi)，维度明确）
		call matmul_real(mat_inv, d_control_x, weights_x)
		call matmul_real(mat_inv, d_control_y, weights_y)
		! 4. 释放逆矩阵（用完即释，避免内存占用）
		if (allocated(mat_inv)) deallocate(mat_inv)

		allocate(distances(nControl))
		do pt = 1, nPoints
			do j = 1, nControl
				diff = points(pt,1:2) - control_points(j,:)
				distances(j) = sqrt(sum(diff**2))
			end do
			call gaussian_rbf_array(distances, epsilonn, phi_row)
			!计算网格点的位移d_mesh_x,d_mesh_y
			call dot_product_reall(weights_x, phi_row, d_mesh_x(pt)) 
        	call dot_product_reall(weights_y, phi_row, d_mesh_y(pt)) 
		end do

		point_update = points
		do pt = 1, nPoints
		    point_update(pt,1) = point_update(pt,1) + d_mesh_x(pt)
		    point_update(pt,2) = point_update(pt,2) + d_mesh_y(pt)
		end do
        print *, "aAb"
!        if (allocated(point_update)) deallocate(point_update)
		if (allocated(points)) deallocate(points)		
		if (allocated(wing_coords)) deallocate(wing_coords)
		if (allocated(wing_update_coords)) deallocate(wing_update_coords)
		if (allocated(inoutput_coords)) deallocate(inoutput_coords)
		if (allocated(control_points)) deallocate(control_points)
		if (allocated(d_wing_x)) deallocate(d_wing_x)
		if (allocated(d_wing_y)) deallocate(d_wing_y)
		if (allocated(d_inoutput_x)) deallocate(d_inoutput_x)
		if (allocated(d_inoutput_y)) deallocate(d_inoutput_y)
		if (allocated(d_control_x)) deallocate(d_control_x)
		if (allocated(d_control_y)) deallocate(d_control_y)
		if (allocated(d_mesh_x)) deallocate(d_mesh_x)
		if (allocated(d_mesh_y)) deallocate(d_mesh_y)
		if (allocated(d_mesh_z)) deallocate(d_mesh_z)
		if (allocated(r_row)) deallocate(r_row)
		if (allocated(phi_row)) deallocate(phi_row)
		if (allocated(phi)) deallocate(phi)
		if (allocated(weights_x)) deallocate(weights_x)
		if (allocated(weights_y)) deallocate(weights_y)
		if (allocated(distances)) deallocate(distances)
		if (allocated(originalLines)) deallocate(originalLines)
		if (allocated(newDataLines)) deallocate(newDataLines)
		if (allocated(newLines)) deallocate(newLines)

		if (allocated(wing)) deallocate(wing)  ! 释放本地wing
		if (allocated(inoutput)) deallocate(inoutput)
		print *, "aAc"
	end subroutine meshdeformation
	! 矩阵×向量乘法子程序：C = A * B（A: m×k 二维矩阵, B: k 一维向量, C: m 一维向量）
	subroutine matmul_real(A, B, C)
		implicit none
		! 输入：A为二维矩阵（m行k列），B为一维向量（长度k）
		real*8, intent(in)  :: A(:,:), B(:)
		! 输出：C为一维向量（长度m，与A的行数一致）
		real*8, intent(out) :: C(:)
		! 局部变量：矩阵/向量维度、循环变量
		integer :: m, k, i, l
		
		! 初始化结果向量为0
		C = 0.0d0
		
		! 获取维度并检查合法性
		m = size(A, 1)  ! A的行数（结果C的长度必须等于m）
		k = size(A, 2)  ! A的列数（必须等于B的长度）
		
		! 检查维度匹配：A的列数 = B的长度
		if (k /= size(B)) then
		    write(*,*) "Error in matmul_real: A列数（", k, "）≠ B长度（", size(B), "）"
		    stop
		end if
		! 检查结果C的长度 = A的行数
		if (size(C) /= m) then
		    write(*,*) "Error in matmul_real: C长度（", size(C), "）≠ A行数（", m, "）"
		    stop
		end if
		
		! 矩阵×向量乘法核心逻辑（m×k 乘 k×1 = m×1）
		do i = 1, m        ! 遍历A的每一行（结果C的每个元素）
		    do l = 1, k    ! 遍历A的列/B的元素（累加求和）
		        C(i) = C(i) + A(i,l) * B(l)
		    end do
		end do
	end subroutine matmul_real
	! 向量点积子程序：result = A · B（A和B需等长）
	subroutine dot_product_reall(A, B, result)
		implicit none
		! 输入参数：两个待做点积的向量（real*8 类型）
		real*8, intent(in)  :: A(:), B(:)
		! 输出参数：点积结果（标量）
		real*8, intent(out) :: result
		! 局部变量：向量长度、循环变量
		integer :: len_a, len_b, i
		
		! 1. 获取向量长度并检查合法性（点积要求两向量等长）
		len_a = size(A)
		len_b = size(B)
		if (len_a /= len_b) then
		    write(*,*) "Error in dot_product_real: 向量长度不匹配！"
		    write(*,*) "A长度 = ", len_a, "，B长度 = ", len_b
		    stop  ! 终止程序，避免错误计算
		end if
		
		! 2. 循环计算点积（元素相乘累加）
		result = 0.0d0  ! 初始化结果
		do i = 1, len_a
		    result = result + A(i) * B(i)
		end do
	end subroutine dot_product_reall
    subroutine gaussian_rbf_array(r_array, epsilon, phi_array)
		implicit none
		real(kind=8), intent(in) :: r_array(:)
		real(kind=8), intent(in) :: epsilon
		real(kind=8), intent(out) :: phi_array(size(r_array))
		real(kind=8), parameter :: d_nomal = 1.0d0 / 4.5d0
		integer :: i

		do i = 1, size(r_array)
		    phi_array(i) = (1.0d0 - d_nomal * r_array(i))**2
		end do
	end subroutine gaussian_rbf_array
	subroutine inverse(n, mat, mat_inv)
	    integer, intent(in) :: n
		real(kind=8), intent(in) :: mat(n,n)
		real(kind=8), intent(inout) :: mat_inv(n,n)
		integer :: info = 0, i, j
		character(len=256) :: err_msg  ! 临时存储错误信息

		mat_inv = mat		
		! Cholesky分解
		call dpotrf('L', n, mat_inv, n, info)
		if (info /= 0) then
		    write(err_msg, '(a,i0,a)') "inverse: Cholesky decomposition failed (info=", info, ")"
		end if
		
		! 矩阵求逆
		call dpotri('L', n, mat_inv, n, info)
		if (info /= 0) then
		    write(err_msg, '(a,i0,a)') "inverse: Inversion failed (info=", info, ")"
		end if
		
		! 填充上三角部分（对称矩阵）
		do i = 1, n
		    do j = i+1, n
		        mat_inv(i, j) = mat_inv(j, i)
		    end do
		end do
	end subroutine inverse
	function readOFPointsFilef(filePath) result(points)
        character(len=*), intent(in) :: filePath
        real(kind=8), allocatable :: points(:,:)
        character(len=256), allocatable :: lines(:)
        integer :: startLine, pCount, i, nLines, iostat, file_unit
        character(len=256) :: line, bracketsRemoved
        real(kind=8) :: coords(3)
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Can't open points file: ", trim(filePath)
            allocate(points(0,3))
            return
        end if
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
        allocate(lines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItemss(lines, startLine, pCount)
        if (pCount <= 0) then
            print *, "Error: Invalid point count in ", trim(filePath)
            deallocate(lines)
            allocate(points(0,3))
            return
        end if
        
        allocate(points(pCount, 3))
        do i = 1, pCount
            line = trim(lines(startLine + i - 1))
            bracketsRemoved = line(2:len_trim(line)-1)
            read(bracketsRemoved, *) coords
            points(i,:) = coords
        end do
        deallocate(lines)
    end function readOFPointsFilef
    subroutine OFFile_FindNItemss(fileLines, startLine, itemCount)
        implicit none
        character(len=*), intent(in) :: fileLines(:)
        integer, intent(out) :: startLine, itemCount
        integer :: i, pos, iostat
        character(len=256) :: line
        
        itemCount = 0
        startLine = 0
        do i = 1, size(fileLines)
            line = trim(fileLines(i))
            pos = index(line, "nPoints") + index(line, "nFaces") + index(line, "size")
            if (pos == 0) cycle
            
            read(line(pos:), *, iostat=iostat) itemCount
            if (itemCount > 0) then
                startLine = i + 2
                return
            end if
        end do
        
        do i = 1, size(fileLines)
            if (isNumber(trim(fileLines(i)))) then
                read(fileLines(i), *, iostat=iostat) itemCount
                if (itemCount > 0) then
                    startLine = i + 2
                    return
                end if
            end if
        end do
        
        do i = 1, size(fileLines)
            if (index(fileLines(i), "(") > 0) then
                startLine = i + 1
                exit
            end if
        end do
        itemCount = 0
        do i = startLine, size(fileLines)
            if (index(fileLines(i), ")") > 0) exit
            itemCount = itemCount + 1
        end do
    end subroutine OFFile_FindNItemss

    ! 辅助函数：获取空闲单元号
    integer function get_free_unit()
        implicit none
        integer :: i, iostat
        logical :: opened
        do i = 10, 999
            inquire(unit=i, opened=opened, iostat=iostat)
            if (.not. opened .and. iostat == 0) then
                get_free_unit = i
                return
            end if
        end do
    end function get_free_unit
	function CSV_read(filename) result(data)
		character*(*), intent(in) :: filename
		real*8, allocatable :: data(:,:)
		integer :: unit, i, j, nrows, ncols
		
!		open(newunit=unit, file=filename, status='old')
        unit = get_free_unit()  ! 用自定义函数获取有效单元号
		open(unit=unit, file=filename, status='old')
		read(unit, *) nrows, ncols
		allocate(data(nrows, ncols))
		
		do i = 1, nrows
		    read(unit, *) (data(i,j), j=1,ncols)
		end do
		
		close(unit)
	end function CSV_read

	! CSV写入函数（实际使用空格分隔文件）
	subroutine CSV_write(filename, data, nrows, ncols)
		character*(*), intent(in) :: filename
		integer, intent(in) :: nrows, ncols
		real*8, intent(in) :: data(nrows,ncols)		
		integer :: unit, i, j
		
		! 打开文件
!		open(newunit=unit, file=filename, status='replace')
		unit = get_free_unit()
		open(unit=unit, file=filename, status='replace')
		! 写入行数和列数
		write(unit, '(i5, i5)') nrows, ncols
		
		! 写入数据
		do i = 1, nrows
			write(unit, '(99f10.6)') (data(i,j), j=1,ncols)
		end do
		
		close(unit)
	end subroutine CSV_write
	! 辅助函数：判断是否为数字
    logical function isNumber(str)
        character(len=*), intent(in) :: str
        real(kind=8) :: num
        integer :: iostat
        isNumber = .false.
        read(str, *, iostat=iostat) num
        if (iostat == 0) isNumber = .true.
    end function isNumber
end module meshdeformationn
