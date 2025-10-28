!algorithm tools

!-------------------------------------------------------------
! Subroutine: hungarian_max
! 功能：对给定的内积矩阵 P（n×n）求一对一匹配，
!      使得总内积最大。
! 思路：先计算矩阵中最大内积 P_max，然后构造成本矩阵
!      C = P_max - P，再调用标准的匈牙利算法求解最小成本匹配。
!-------------------------------------------------------------
subroutine hungarian_max(n, p, assignment)
   implicit none
   integer, intent(in) :: n
   real, intent(in) :: p(n, n)
   integer :: assignment(n)
   real, dimension(n, n) :: cost
   real :: p_max

   ! 计算内积矩阵中最大的元素
   p_max = maxval(p)

   ! 将内积矩阵转换为成本矩阵
   ! 注意：如果 p 中的数值都非负，则 cost 保证也是非负的
   cost = p_max - p

   ! 调用匈牙利算法，求解最小化成本问题，从而得到最大内积的指派关系
   call hungarian(n, cost, assignment)
end subroutine hungarian_max

!-------------------------------------------------------------
! Subroutine: hungarian
! 功能：对 n×n 的成本矩阵 cost 求解最优指派（最小总代价）。
! 输出 assignment 数组，其中 assignment(i) 表示第 i 行匹配到的列。
!-------------------------------------------------------------
subroutine hungarian(n, cost, assignment)
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: cost(n, n)
    integer :: assignment(n)
    
    integer :: i, j, k, c
    integer :: i_new, j_new, i2, j2, i3, j3, i4
    real :: minimum
    logical, dimension(n, n) :: star, prime
    logical, dimension(n) :: coveredRow, coveredCol
    logical :: found, done

    ! 初始化所有标记和 assignment 数组
    star = .false.
    prime = .false.
    coveredRow = .false.
    coveredCol = .false.
    assignment = 0

    !---------- Step 1: 行归约 ----------
    do i = 1, n
       minimum = minval(cost(i, :))
       do j = 1, n
          cost(i, j) = cost(i, j) - minimum
       end do
    end do

    !---------- Step 2: 列归约 ----------
    do j = 1, n
       minimum = minval(cost(:, j))
       do i = 1, n
          cost(i, j) = cost(i, j) - minimum
       end do
    end do

    !---------- Step 3: 星标零的初步匹配 ----------
    coveredRow = .false.
    coveredCol = .false.
    do i = 1, n
       do j = 1, n
          if ((abs(cost(i, j)) < 1.0e-6) .and. (.not. coveredRow(i)) .and. (.not. coveredCol(j))) then
             star(i, j) = .true.
             coveredRow(i) = .true.
             coveredCol(j) = .true.
          end if
       end do
    end do
    ! 取消覆盖标记，为下一步作准备
    coveredRow = .false.
    coveredCol = .false.

    !---------- Step 4: 用直线覆盖所有含有星标零的列 ----------
    do j = 1, n
       do i = 1, n
          if (star(i, j)) then
             coveredCol(j) = .true.
             exit
          end if
       end do
    end do

    !---------- Step 5: 主循环 - 不断改善匹配直至完美匹配 ----------
    do while (sum([(merge(1, 0, coveredCol(j)), j=1, n)]) < n)
        found = .false.
        do i = 1, n
            if (.not. coveredRow(i)) then
                do j = 1, n
                    if ((.not. coveredCol(j)) .and. (abs(cost(i, j)) < 1.0e-6)) then
                        prime(i, j) = .true.
                        done = .true.
                        ! 检查当前行中是否已有星标零
                        do k = 1, n
                            if (star(i, k)) then
                                done = .false.
                                c = k
                                exit
                            end if
                        end do
                        if (done) then
                        ! 找到未匹配的零：构造增广路径
                            call augmentPath(i, j, star, prime, coveredRow, coveredCol, n)
                            prime = .false.
                            coveredRow = .false.
                            coveredCol = .false.
                            ! 重新覆盖所有含有星标零的列
                            do j_new = 1, n
                            do i_new = 1, n
                                if (star(i_new, j_new)) then
                                    coveredCol(j_new) = .true.
                                    exit
                                end if
                            end do
                            end do
                            found = .true.
                            exit  ! 退出当前行循环
                        else
                            ! 如果该行已有星标零，则覆盖该行并取消覆盖星标所在的那一列
                            coveredRow(i) = .true.
                            coveredCol(c) = .false.
                        end if
                    end if
                end do
                if (found) exit
            end if
        end do

        if (.not. found) then
        ! 如果未找到未覆盖的零，则调整成本矩阵
            minimum = 1.0e9
            do i2 = 1, n
                if (.not. coveredRow(i2)) then
                    do j2 = 1, n
                        if (.not. coveredCol(j2)) then
                            if (cost(i2, j2) < minimum) then
                                minimum = cost(i2, j2)
                            end if
                        end if
                    end do
                end if
            end do
        ! 将最小值加到所有覆盖行上
            do i3 = 1, n
                if (coveredRow(i3)) then
                    do j3 = 1, n
                        cost(i3, j3) = cost(i3, j3) + minimum
                    end do
                end if
            end do
        ! 从所有未覆盖列上减去最小值
            do j3 = 1, n
                if (.not. coveredCol(j3)) then
                    do i4 = 1, n
                        cost(i4, j3) = cost(i4, j3) - minimum
                    end do
                end if
            end do
        end if
    end do


    !---------- Step 6: 从星标零中构造最终匹配 ----------
    do i = 1, n
       do j = 1, n
          if (star(i, j)) then
             assignment(i) = j
             exit
          end if
       end do
    end do
end subroutine hungarian

!-------------------------------------------------------------
! Subroutine: augmentPath
! 功能：构造一条增广路径，并沿路径交换星标与素标状态，
!      从而增加匹配的数量。
!-------------------------------------------------------------
subroutine augmentPath(row0, col0, star, prime, coveredRow, coveredCol, n)
    implicit none
    integer, intent(in) :: row0, col0, n
    logical, intent(inout) :: star(n, n), prime(n, n)
    logical, intent(inout) :: coveredRow(n), coveredCol(n)
    integer :: i, j, count, row, col
    integer, dimension(n*n, 2) :: path  ! 用于存储增广路径的序列
    logical :: done

    count = 1
    path(count,1) = row0
    path(count,2) = col0
    done = .false.
    do while (.not. done)
       row = 0
       do i = 1, n
          if (star(i, path(count,2))) then
             row = i
             exit
          end if
       end do
       if (row == 0) then
          done = .true.
       else
          count = count + 1
          path(count,1) = row
          path(count,2) = path(count-1,2)
          col = 0
          do j = 1, n
             if (prime(path(count,1), j)) then
                col = j
                exit
             end if
          end do
          count = count + 1
          path(count,1) = path(count-1,1)
          path(count,2) = col
       end if
    end do

    ! 沿增广路径交换星标和素标
    do i = 1, count
       if (star(path(i,1), path(i,2))) then
          star(path(i,1), path(i,2)) = .false.
       else
          star(path(i,1), path(i,2)) = .true.
       end if
    end do

    prime = .false.
end subroutine augmentPath
