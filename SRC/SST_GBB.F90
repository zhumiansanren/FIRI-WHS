!==================================================================================================
! 子程序包: SST_GBB (General Bottom Boundary Package for Soil Salt Transport Process)
! 主要功能: 模拟底部边界的盐分运移
! 创建日期: 2020年7月6日
! 修改日志: 
!==================================================================================================


SUBROUTINE SST_GBB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SST_GBB_DF


SUBROUTINE SST_GBB_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(CBOT))   ALLOCATE(CBOT(NPER), 	SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WSABOT)) ALLOCATE(WSABOT(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SST_GBB_AL


SUBROUTINE SST_GBB_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('CBOT')) THEN
    CALL RDFREA('CBOT', CBOT, NPER, NPER)
  ELSE
    CBOT(:) = 0.0
  END IF
  IF (RDINQR('WSABOT')) THEN
    CALL RDFREA('WSABOT', WSABOT, NPER, NPER)
  ELSE
    WSABOT(:) = 0.0
  END IF

  RETURN
END SUBROUTINE SST_GBB_RP



SUBROUTINE SST_GBB_FM
!**************************************************************************************************
! 计算底部边界盐分运移的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: DELV2
  REAL    :: O2, OSAT2, TOR2

  I = NCEL
  J = IMAT(I)

  ! Dirichlet 条件
  IF (ISSTBOT == 1) THEN

    IF (ISWFBOT == 1) THEN
      DELV2 = DELV(I) * 0.5
      O2 = (OBOT(JPER) + ONEW(I)) * 0.5
      OSAT2 = OSAT(J)
 
      SELECT CASE (ITOR)
      CASE DEFAULT
        TOR2 = 0.66 * (O2 / OSAT2) ** (8.0 / 3.0)
      CASE (2)
        TOR2 = O2 ** (7.0 / 3.0) / OSAT2 ** 2.0
      END SELECT
      
      DIS(I, 2) = DISL(J) * ABS(QV(I, 2) / DELT / O2) + O2 * DIFW(J) * TOR2
      COF(I) = COF(I) + QV(I, 2) * 0.5 - DELT * DIS(I, 2) / DELV2
      RHS(I) = RHS(I) - (QV(I, 2) * 0.5 + DELT * DIS(I, 2) / DELV2) * CBOT(JPER)
    END IF

    IF (ISWFBOT == 2) THEN
      RHS(I) = RHS(I) - WBOT(JPER) * DELT * CBOT(JPER)
    END IF

    IF (ISWFBOT == 3) THEN
      IF (WBOT(JPER) < 0.0) RHS(I) = RHS(I) - DELT * WBOT(JPER) * CBOT(JPER)
    END IF

    IF (ISWFBOT == 4) THEN
      IF (WBOT(JPER) < 0.0) RHS(I) = RHS(I) - DELT * WBOT(JPER) * CBOT(JPER)
    END IF

    IF (ISWFBOT == 5) THEN
      IF (HNEW(I) >= 0.0) THEN
        RHS(I) = RHS(I) - WBOT(JPER) * DELT * CBOT(JPER)
      END IF
    END IF

  END IF

  ! Neumann 条件
  IF (ISSTBOT == 2) THEN
    RHS(I) = RHS(I) - DELT * WSABOT(JPER)
  END IF

  ! 零梯度条件
  IF (ISSTBOT == 3) THEN
    IF (WBOT(JPER) < 0.0) THEN
      RHS(I) = RHS(I) - DELT * WBOT(JPER) * CNEW(I)
    END IF
  END IF

  RETURN
END SUBROUTINE SST_GBB_FM


SUBROUTINE SST_GBB_BD
!**************************************************************************************************
! 计算底部边界盐分运移的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV2

  I = NCEL

  ! Dirichlet 条件
  IF (ISSTBOT == 1) THEN

    IF (ISWFBOT == 1) THEN
      DELV2 = DELV(I) * 0.5
      WSABOT(JPER) = WSABOT(JPER) + WBOT(JPER) * DELT * (CBOT(JPER) + CNEW(I)) * 0.5 + &
                     DELT * DIS(I, 2) / DELV2 * (CBOT(JPER) - CNEW(I))
    END IF

    IF (ISWFBOT == 2) THEN
      WSABOT(JPER) = WSABOT(JPER) + DELT * WBOT(JPER) * CBOT(JPER)
    END IF

    IF (ISWFBOT == 3) THEN
      IF (WBOT(JPER) < 0.0) WSABOT(JPER) = WSABOT(JPER) + DELT * WBOT(JPER) * CBOT(JPER)
    END IF

    IF (ISWFBOT == 4) THEN
      IF (WBOT(JPER) < 0.0) WSABOT(JPER) = WSABOT(JPER) + DELT * WBOT(JPER) * CBOT(JPER)
    END IF

    IF (ISWFBOT == 5) THEN
      IF (HNEW(I) >= 0.0) WSABOT(JPER) = WSABOT(JPER) + DELT * WBOT(JPER) * CBOT(JPER)
    END IF

  END IF

  ! Neumann 边界条件
  IF (ISSTBOT == 2) THEN
    
  END IF

  ! 零梯度条件
  IF (ISSTBOT == 3) THEN
    IF (WBOT(JPER) < 0.0) WSABOT(JPER) = WSABOT(JPER) + DELT * WBOT(JPER) * CNEW(I)
  END IF

  RETURN
END SUBROUTINE SST_GBB_BD


SUBROUTINE SST_GBB_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(CBOT))   DEALLOCATE(CBOT)
  IF (ALLOCATED(WSABOT)) DEALLOCATE(WSABOT)

  RETURN
END SUBROUTINE SST_GBB_DA

