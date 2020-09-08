!==================================================================================================
! 子程序包: SHT_GTB (General Top Boundary Package for Soil Heat Transport Process)
! 主要功能: 模拟地表边界的热传导
! 创建日期: 2020年7月4日
! 修改日志:
!==================================================================================================


SUBROUTINE SHT_GTB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SHT_GTB_DF


SUBROUTINE SHT_GTB_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(TTOP))  ALLOCATE(TTOP(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTHTOP)) ALLOCATE(WTHTOP(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SHT_GTB_AL


SUBROUTINE SHT_GTB_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('TTOP')) THEN
    CALL RDFREA('TTOP', TTOP, NPER, NPER)
  ELSE
    TTOP(:) = TINI(1)
  END IF
  IF (RDINQR('WTHTOP')) THEN
    CALL RDFREA('WTHTOP', WTHTOP, NPER, NPER)
  ELSE
    WTHBOT(:) = 0.0
  END IF

  RETURN
END SUBROUTINE SHT_GTB_RP


subroutine SHT_GTB_FM
!**************************************************************************************************
! 计算表示地表边界热传导的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV1
  REAL    :: THCON1

  I = 1

  SELECT CASE (ISHTTOP)

   ! Dirichlet 条件
   CASE (1)
    DELV1 = DELV(I) * 0.5
    THCON1 = THCON(I)
    COF(I) = COF(I) - DELT * THCON1 / DELV1
    RHS(I) = RHS(I) - DELT * THCON1 / DELV1 * TTOP(JPER)
    
  ! Neumann 条件
  CASE (2)
    RHS(I) = RHS(I) - DELT * WTHTOP(JPER)

  END SELECT

  RETURN
END SUBROUTINE SHT_GTB_FM


SUBROUTINE SHT_GTB_BD
!**************************************************************************************************
! 计算表示地表边界热传导的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV1
  REAL    :: THCON1

  I = 1

  SELECT CASE (ISHTTOP)

  ! Dirichlet 条件
  CASE (1)
    DELV1 = DELV(I) * 0.5
    THCON1 = THCON(I)
    WTHTOP(JPER) = WTHTOP(JPER) + DELT * THCON1 * (TTOP(JPER) - TNEW(I)) / DELV1

  ! Neumann 条件
  CASE (2)

  END SELECT

  RETURN
END SUBROUTINE SHT_GTB_BD


SUBROUTINE SHT_GTB_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(TTOP))   DEALLOCATE(TTOP)
  IF (ALLOCATED(WTHTOP)) DEALLOCATE(WTHTOP)

  RETURN
END SUBROUTINE SHT_GTB_DA

