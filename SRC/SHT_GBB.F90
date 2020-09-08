!==================================================================================================
! 子程序包: SHT_GBB (General Bottom Boundary Package for Soil Heat Transport Process)
! 主要功能: 模拟底部边界的热传导
! 创建日期: 2020年7月4日
! 修改日志:
!==================================================================================================


SUBROUTINE SHT_GBB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE
  
  RETURN
END SUBROUTINE SHT_GBB_DF


SUBROUTINE SHT_GBB_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(TBOT))   ALLOCATE(TBOT(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTHBOT)) ALLOCATE(WTHBOT(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SHT_GBB_AL


SUBROUTINE SHT_GBB_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('TBOT')) THEN
    CALL RDFREA('TBOT', TBOT, NPER, NPER)
  ELSE
    TBOT(:) = TINI(NCEL)
  END IF
  IF (RDINQR('WTHBOT')) THEN
    CALL RDFREA('WTHBOT', WTHBOT, NPER, NPER)
  ELSE
    WTHBOT(:) = 0.0
  END IF

  RETURN
END SUBROUTINE SHT_GBB_RP


SUBROUTINE SHT_GBB_FM
!**************************************************************************************************
! 计算表示底部边界热传导的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV2
  REAL    :: THCON2

  I = NCEL

  SELECT CASE (ISHTBOT)

  ! Dirichlet 条件
  CASE (1)
    DELV2 = DELV(I) * 0.5
    THCON2 = THCON(I)
    COF(I) = COF(I) - DELT * THCON2 / DELV2
    RHS(I) = RHS(I) - DELT * THCON2 / DELV2 * TBOT(JPER)
  
  ! Neumann 条件
  CASE (2)
    RHS(I) = RHS(I) - DELT * WTHBOT(JPER)
  
  ! 零梯度条件
  CASE (3)
    
  END SELECT

  RETURN
END SUBROUTINE SHT_GBB_FM


SUBROUTINE SHT_GBB_BD
!**************************************************************************************************
! 计算表示底部边界热传导的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV2
  REAL    :: THCON2

  I = NCEL

  SELECT CASE (ISHTBOT)

    ! Dirichlet 条件
    CASE (1)
      THCON2 = THCON(I)
      DELV2 = DELV(I) * 0.5
      WTHBOT(JPER) =  WTHBOT(JPER) + DELT * THCON2 * (TBOT(JPER) - TNEW(I)) / DELV2
      QTHV(I, 2) = QTHV(I, 2) + DELT * THCON2 * (TBOT(JPER) - TNEW(I)) / DELV2
    
    ! Neumann 条件
    CASE (2)
      QTHV(I, 2) = QTHV(I, 2) - DELT * WTHBOT(JPER)

    ! 零梯度条件
    CASE (3)

  END SELECT

  RETURN
END SUBROUTINE SHT_GBB_BD


SUBROUTINE SHT_GBB_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(TBOT))   DEALLOCATE(TBOT)
  IF (ALLOCATED(WTHBOT)) DEALLOCATE(WTHBOT)

  RETURN
END SUBROUTINE SHT_GBB_DA
