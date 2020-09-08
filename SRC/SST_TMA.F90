!==================================================================================================
! 子程序包: SST_TMA (Tridiagonal Matrix Algorithm Package for Soil Salt Transport Process)
! 主要功能: 采用 Thomas 算法求解土壤盐分运移有限差分方程组
! 创建日期: 2020年7月8日
! 修改日志:
!==================================================================================================


SUBROUTINE SST_TMA_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('MAXITR')) THEN
    CALL RDSINT('MAXITR', MAXITR)
  ELSE
    MAXITR = 30
  END IF
  IF (RDINQR('CTOL')) THEN
    CALL RDSREA('CTOL', CTOL)
  ELSE
    CTOL = 0.0001
  END IF

  RETURN
END SUBROUTINE SST_TMA_DF


SUBROUTINE SST_TMA_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_TMA_AL


SUBROUTINE SST_TMA_AP
!**************************************************************************************************
! 求解土壤盐分运移有限差分方程组
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL :: BET
  REAL, DIMENSION(:), ALLOCATABLE :: D
  REAL, DIMENSION(:), ALLOCATABLE :: E
  REAL, DIMENSION(:), ALLOCATABLE :: F
  REAL, DIMENSION(:), ALLOCATABLE :: R
  REAL, DIMENSION(:), ALLOCATABLE :: G
  REAL, DIMENSION(:), ALLOCATABLE :: U

  ! 为局部变量分配内存
  IF (.NOT. ALLOCATED(D)) ALLOCATE(D(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(E)) ALLOCATE(E(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(F)) ALLOCATE(F(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(R)) ALLOCATE(R(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(G)) ALLOCATE(G(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(U)) ALLOCATE(U(NCEL), SOURCE = 0.0)

  ! 计算三对角线性方程组的系数和右端项
  D(:) = BBV(:, 1) + CCV(:, 1)
  E(:) = DDV(:, 1) - CCV(:, 1) + DDV(:, 2) - CCV(:, 2) + COF(:)
  F(:) = BBV(:, 2) + CCV(:, 2)
  R(:) = RHS(:)

  ! 求解三对角线性方程组
  BET = E(1)
  IF (ABS(BET) < TINY(0.0)) THEN
    WRITE(U_LOG, "('I = ', I12,' BETA = ', ES12.3)") 1, BET
    WRITE(U_LOG, "('错误 1 -- 土壤盐分运移方程组系数错误!')")
    IF (IRUNMOD == 1) WRITE(*, "('I = ', I12,' BETA = ', ES12.3)") I, BET
    IF (IRUNMOD == 1) WRITE(*, "('错误 1 -- 土壤盐分运移方程组系数错误!')")
    IERRCOD = 3
  END IF
  U(1) = R(1) / BET
  DO I = 2, NCEL
    G(I) = F(I - 1) / BET
    BET = E(I) - D(I) * G(I)
    IF (ABS(BET) < TINY(0.0)) THEN
      WRITE(U_LOG, "('I = ', I12,' BETA = ', ES12.3)") I, BET
      WRITE(U_LOG, "('错误 2 -- 土壤盐分运移方程组系数错误!')")
      IF (IRUNMOD == 1) WRITE(*, "('I = ', I12,' BETA = ', ES12.3)") I, BET
      IF (IRUNMOD == 1) WRITE(*, "('错误 2 -- 土壤盐分运移方程组系数错误!')")
      IERRCOD = 3
    END IF
    U(I) = (R(I) - D(I) * U(I - 1)) / BET
  END DO
  DO I = NCEL - 1, 1, -1
    U(I) = U(I) - G(I + 1) * U(I + 1)
  END DO

  ! 判断解是否收敛
  LCONV = .FALSE.
  IF (MAXVAL(ABS(CNEW(:) - U(:))) < CTOL) THEN
    LCONV = .TRUE.
  END IF

  ! 如果不收敛, 输出一个警告信息
  IF ((JITR >= MAXITR) .AND. (.NOT. LCONV)) THEN
    WRITE(U_LOG, "('警告 -- 土壤盐分运移方程组的解不收敛!')")
    IF (IRUNMOD == 1) WRITE(*, "('警告 -- 土壤盐分运移方程组的解不收敛!')")
    IERRCOD = 4
  END IF

  ! 保存方程组的解
  CNEW(:) = U(:)

  ! 释放局部变量的内存
  IF (ALLOCATED(D)) DEALLOCATE(D)
  IF (ALLOCATED(E)) DEALLOCATE(E)
  IF (ALLOCATED(F)) DEALLOCATE(F)
  IF (ALLOCATED(R)) DEALLOCATE(R)
  IF (ALLOCATED(G)) DEALLOCATE(G)
  IF (ALLOCATED(U)) DEALLOCATE(U)

  RETURN
END SUBROUTINE SST_TMA_AP


SUBROUTINE SST_TMA_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_TMA_DA

