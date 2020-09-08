!==================================================================================================
! 子程序包: SST_GTB (General Top Boundary Package for Soil Salt Transport Process)
! 主要功能: 模拟地表边界的盐分运移
! 创建日期: 2020年7月8日
! 修改日志: 
!==================================================================================================


SUBROUTINE SST_GTB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_GTB_DF


SUBROUTINE SST_GTB_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(CTOP))   ALLOCATE(CTOP(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSATOP)) ALLOCATE(WSATOP(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SST_GTB_AL


SUBROUTINE SST_GTB_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('CTOP')) THEN
    CALL RDFREA('CTOP', CTOP, NPER, NPER)
  ELSE
    CTOP(:) = 0.0
  END IF
  IF (RDINQR('WSATOP')) THEN
    CALL RDFREA('WSATOP', WSATOP, NPER, NPER)
  ELSE
    WSATOP(:) = 0.0
  END IF

  RETURN
END SUBROUTINE SST_GTB_RP


SUBROUTINE SST_GTB_FM
!**************************************************************************************************
! 计算地表边界盐分运移的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: DELV1
  REAL    :: O1, OSAT1, TOR1

  I = 1
  J = IMAT(I)

  ! Dirichlet 条件
  IF (ISSTTOP == 1) THEN
    IF (ISWFTOP == 1) THEN
      DELV1 = DELV(I) * 0.5
      O1 = (OTOP(JPER) + ONEW(I)) * 0.5
      OSAT1 = OSAT(J)

      SELECT CASE (ITOR)
      CASE DEFAULT
        TOR1 = 0.66 * (O1 / OSAT1) ** (8.0 / 3.0)
      CASE (2)
        TOR1 = O1 ** (7.0 / 3.0) / OSAT1 ** 2.0
      END SELECT

      DIS(I, 1) = DISL(J) * ABS(QTOP) / DELT / O1 + O1 * DIFW(J) * TOR1
      COF(I) = COF(I) + 0.5 * QTOP - DELT * DIS(I, 1) / DELV1
      RHS(I) = RHS(I) - 0.5 * QTOP * CTOP(JPER) - DELT * DIS(I, 1) / DELV1 * CTOP(JPER)
    END IF

    IF (ISWFTOP == 2) THEN
      IF (WTOP(JPER) > 0.0) THEN
        COF(I) = COF(I) + 0.5 * DELT * WTOP(JPER)
        RHS(I) = RHS(I) - 0.5 * DELT * WTOP(JPER) * CTOP(JPER)
      END IF
    END IF

    IF (ISWFTOP == 3) THEN
      COF(I) = COF(I) + 0.5 * DELT * WTOP(JPER)
      RHS(I) = RHS(I) - 0.5 * DELT * WTOP(JPER) * CTOP(JPER)
    END IF

  END IF

  ! Neumann 条件
  IF (ISSTTOP == 2) THEN
    RHS(I) = RHS(I) - DELT * WSATOP(JPER)
  END IF

  RETURN
END SUBROUTINE SST_GTB_FM


SUBROUTINE SST_GTB_BD
!**************************************************************************************************
! 计算地表边界盐分运移的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV1

  I = 1

  ! Dirichlet 条件
  IF (ISSTTOP == 1) THEN

    IF (ISWFTOP == 1) THEN
      DELV1 = DELV(I) * 0.5
      WSATOP(JPER) = WSATOP(JPER) + QTOP * (CTOP(JPER) + CNEW(I)) * 0.5 + &
                     DELT * DIS(I, 1) / DELV1 * (CTOP(JPER) - CNEW(I))
    END IF

    IF (ISWFTOP == 2) THEN
      IF (WTOP(JPER) > 0.0) THEN
        WSATOP(JPER) = WSATOP(JPER) + QTOP * CTOP(JPER)
      END IF
    END IF

  END IF

  ! Neumann 条件
  IF (ISSTTOP == 2) THEN
   
  END IF

  RETURN
END SUBROUTINE SST_GTB_BD


SUBROUTINE SST_GTB_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(CTOP))    DEALLOCATE(CTOP)
  IF (ALLOCATED(WSATOP))  DEALLOCATE(WSATOP)

  RETURN
END SUBROUTINE SST_GTB_DA

