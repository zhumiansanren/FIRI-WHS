!==================================================================================================
! 子程序包: SWF_RWU (Root Water Uptake Pacakge for Soil Water Flow Process)
! 主要功能: 模拟根系吸水
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_RWU_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('IRWUDIST')) THEN
    CALL RDSINT('IRWUDIST', IRWUDIST)
  ELSE
    IRWUDIST = 1
  END IF
  IF (RDINQR('PRWUDIST')) THEN
    CALL RDSREA('PRWUDIST', PRWUDIST)
  ELSE
    PRWUDIST = 5.0
  END IF
  IF (RDINQR('RWUHCRT')) THEN
    CALL RDFREA('RWUHCRT', RWUHCRT, 5, 5)
  ELSE
    RWUHCRT(:) = (/0.0, -1.0, -300.0, -600.0, -15000.0/)
  END IF
   IF (RDINQR('RWUTCRT')) THEN
    CALL RDFREA('RWUTCRT', RWUTCRT, 2, 2)
  ELSE
    RWUTCRT(:) = (/0.5, 0.1/)
  END IF
  IF (RDINQR('RWUCOMPCRT')) THEN
    CALL RDSREA('RWUCOMPCRT', RWUCOMPCRT)
  ELSE
    RWUCOMPCRT = 1.0
  END IF

  RETURN
END SUBROUTINE SWF_RWU_DF


SUBROUTINE SWF_RWU_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(RDEPTH))  ALLOCATE(RDEPTH(NPER),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTPOT))   ALLOCATE(WTPOT(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WTACT))   ALLOCATE(WTACT(NPER),   SOURCE = 0.0)
  IF (.NOT. ALLOCATED(RWUDIST)) ALLOCATE(RWUDIST(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WRWUACT)) ALLOCATE(WRWUACT(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(QRWUACT)) ALLOCATE(QRWUACT(NCEL), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_RWU_AL


SUBROUTINE SWF_RWU_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('RWUDIST')) THEN
    CALL RDFREA('RWUDIST', RWUDIST, NCEL, NCEL)
  END IF

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WTPOT')) THEN
    CALL RDFREA('WTPOT', WTPOT, NPER, NPER)
  END IF
  IF (RDINQR('RDEPTH')) THEN
    CALL RDFREA('RDEPTH', RDEPTH, NPER, NPER)
  END IF

  ! 初始化根系吸水补偿系数
  RWUCOMP = 1.0

  RETURN
END SUBROUTINE SWF_RWU_RP


SUBROUTINE SWF_RWU_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL :: R1, R2, SUM1

  SELECT CASE (IRWUDIST)

  ! 采用指数函数计算根系吸水分布
  CASE (1)
    DO I = 1, NCEL
      RWUDIST(I) = 0.0
      IF (RDEPTH(JPER) >= LOCV(I, 2)) THEN
        R1 = EXP(- PRWUDIST * LOCV(I, 1) / RDEPTH(JPER)) / (1.0 - EXP(- PRWUDIST))
        R2 = EXP(- PRWUDIST * LOCV(I, 3) / RDEPTH(JPER)) / (1.0 - EXP(- PRWUDIST))
        RWUDIST(I) = R1 - R2
      END IF
    END DO

  ! 归一化根系吸水分布
  CASE (2)
    SUM1 = SUM(RWUDIST(:))
    DO I = 1, NCEL
      RWUDIST(I) = RWUDIST(I) / SUM1
    END DO

  END SELECT

  RETURN
END SUBROUTINE SWF_RWU_ST


SUBROUTINE SWF_RWU_FM
!**************************************************************************************************
! 计算表示根系吸水的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL :: WSR      ! 水分胁迫响应函数值

  QRWUACT(:) = 0.0
  DO I = 1, NCEL
    CALL U_RWUWSR_FED(WSR, HNEW(I), WTPOT(JPER) / PERLEN(JPER), RWUHCRT, RWUTCRT)
    QRWUACT(I) = DELT * WTPOT(JPER) * RWUDIST(I) * WSR / RWUCOMP
    RHS(I) = RHS(I) + QRWUACT(I)
  END DO

  RETURN
END SUBROUTINE SWF_RWU_FM


SUBROUTINE SWF_RWU_BD
!**************************************************************************************************
! 计算表示根系吸水的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I

  WRWUACT(:) = WRWUACT(:) + QRWUACT(:)
  WTACT(JPER) = WTACT(JPER) + SUM(QRWUACT(:))

  ! 更新根系吸水补偿系数
  RWUCOMP = 1.0
  IF (SUM(QRWUACT(:)) < DELT * WTPOT(JPER)) THEN
    IF (SUM(QRWUACT(:)) / (DELT * WTPOT(JPER)) > RWUCOMPCRT) THEN
      RWUCOMP = SUM(QRWUACT(:)) / SUM(DELT * WTPOT(JPER) * RWUDIST(:))
    END IF
  END IF

  RETURN
END SUBROUTINE SWF_RWU_BD


SUBROUTINE SWF_RWU_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(RDEPTH))  DEALLOCATE(RDEPTH)
  IF (ALLOCATED(WTPOT))   DEALLOCATE(WTPOT)
  IF (ALLOCATED(WTACT))   DEALLOCATE(WTACT)
  IF (ALLOCATED(RWUDIST)) DEALLOCATE(RWUDIST)
  IF (ALLOCATED(WRWUACT)) DEALLOCATE(WRWUACT)
  IF (ALLOCATED(QRWUACT)) DEALLOCATE(QRWUACT)

  RETURN
END SUBROUTINE SWF_RWU_DA

