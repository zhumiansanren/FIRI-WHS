!==================================================================================================
! 子程序包: SWF_GTB (General Top Boundary Package for Soil Water Flow Process)
! 主要功能: 模拟地表边界的水分运动
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_GTB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('TLAGTOP')) THEN
    CALL RDSREA('TLAGTOP', TLAGTOP)
  ELSE
    TLAGTOP = 1.0
  END IF

  RETURN
END SUBROUTINE SWF_GTB_DF


SUBROUTINE SWF_GTB_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(HTOP)) ALLOCATE(HTOP(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OTOP)) ALLOCATE(OTOP(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTOP)) ALLOCATE(WTOP(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_GTB_AL


SUBROUTINE SWF_GTB_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('HTOP')) THEN
    CALL RDFREA('HTOP', HTOP, NPER, NPER)
  ELSE
    HTOP(:) = HINI(1)
  END IF
  IF (RDINQR('OTOP')) THEN
    CALL RDFREA('OTOP', OTOP, NPER, NPER)
  ELSE
    OTOP(:) = OINI(1)
  END IF
  IF (RDINQR('WTOP')) THEN
    CALL RDFREA('WTOP', WTOP, NPER, NPER)
  ELSE
    WTOP(:) = 0.0
  END IF

  RETURN
END SUBROUTINE SWF_GTB_RP


SUBROUTINE SWF_GTB_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J

  I = 1
  J = IMAT(I)

  ! 根据初始条件类型换算水势和含水率
  IF (ISWFTOP == 1) THEN
    IF (ISWFINI == 0) THEN
      CALL U_O_MVG(OTOP(JPER), HTOP(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    ELSE
      CALL U_H_MVG(HTOP(JPER), OTOP(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    END IF
  END IF

  RETURN
END SUBROUTINE SWF_GTB_ST


SUBROUTINE SWF_GTB_FM
!**************************************************************************************************
! 计算表示地表边界水分运动的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: DELV1
  REAL    :: HYCONTOP, HYCON1

  I = 1
  J = IMAT(I)

  SELECT CASE (ISWFTOP)

  ! Dirichlet 条件
  CASE (1)
    IF (ISWFINI /= 0) THEN
      CALL U_H_MVG(HTOP(JPER), OTOP(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    END IF
    CALL U_HYCON_MVG(HYCONTOP, HTOP(JPER), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
    HYCON1 = (HYCONTOP + HYCON(I)) * 0.5
    DELV1 = DELV(I) * 0.5
    COF(I) = COF(I) - DELT * HYCON1 / DELV1
    RHS(I) = RHS(I) - DELT * HYCON1 - DELT * HYCON1 * HTOP(JPER) / DELV1

  ! Neumann 条件
  CASE (2)
    RHS(I) = RHS(I) - WTOP(JPER) * DELT

  ! Cauchy 条件
  CASE (3)
    COF(I) = COF(I) - DELT / TLAGTOP
    RHS(I) = RHS(I) - DELT * (HTOP(JPER) - LOCV(I, 2)) / TLAGTOP

  END SELECT

  RETURN
END SUBROUTINE SWF_GTB_FM


SUBROUTINE SWF_GTB_BD
!**************************************************************************************************
! 计算表示地表边界水分运动的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: HYCONTOP, HYCON1
  REAL :: DELV1

  I = 1
  J = IMAT(I)

  SELECT CASE (ISWFTOP)

  ! Dirichlet 条件
  CASE (1)
    IF (ISWFINI /= 0) CALL U_H_MVG(HTOP(JPER), OTOP(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    CALL U_HYCON_MVG(HYCONTOP, HTOP(JPER), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
    DELV1 = DELV(I) * 0.5
    HYCON1 = (HYCONTOP + HYCON(I)) * 0.5
    QTOP = DELT * HYCON1 * ((HTOP(JPER) - HNEW(I)) / DELV1 + 1.0)
    WTOP(JPER) = WTOP(JPER) + QTOP

  ! Neumann 条件
  CASE (2)
    QTOP = WTOP(JPER) * DELT

  ! Cauchy 条件
  CASE (3)
    QTOP = DELT * (HTOP(JPER) - LOCV(I, 2) - HNEW(I)) / TLAGTOP
    WTOP(JPER) = WTOP(JPER) + QTOP

  END SELECT

  RETURN
END SUBROUTINE SWF_GTB_BD


SUBROUTINE SWF_GTB_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(HTOP)) DEALLOCATE(HTOP)
  IF (ALLOCATED(OTOP)) DEALLOCATE(OTOP)
  IF (ALLOCATED(WTOP)) DEALLOCATE(WTOP)

  RETURN
END SUBROUTINE SWF_GTB_DA

