!==================================================================================================
! 子程序包: SWF_GBB (General Bottom Boundary Package for Soil Water Flow Process)
! 主要功能: 模拟底部边界的水分运动
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_GBB_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取参数
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('TLAGBOT')) THEN
    CALL RDSREA('TLAGBOT', TLAGBOT)
  ELSE
    TLAGBOT = 1.0
  END IF

  RETURN
END SUBROUTINE SWF_GBB_DF


SUBROUTINE SWF_GBB_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(HBOT)) ALLOCATE(HBOT(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OBOT)) ALLOCATE(OBOT(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WBOT)) ALLOCATE(WBOT(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_GBB_AL


SUBROUTINE SWF_GBB_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('HBOT')) THEN
    CALL RDFREA('HBOT', HBOT, NPER, NPER)
  ELSE
    HBOT(:) = HINI(NCEL)
  END IF
  IF (RDINQR('OBOT')) THEN
    CALL RDFREA('OBOT', OBOT, NPER, NPER)
  ELSE
    OBOT(:) = OINI(NCEL)
  END IF
  IF (RDINQR('WBOT')) THEN
    CALL RDFREA('WBOT', WBOT, NPER, NPER)
  ELSE
    WBOT(:) = 0.0
  END IF

  RETURN
END SUBROUTINE SWF_GBB_RP


SUBROUTINE SWF_GBB_ST
!**************************************************************************************************
! 更新随应力期变化的数据
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J

  I = NCEL
  J = IMAT(I)

  ! 根据初始条件类型换算水势和含水率
  IF (ISWFBOT == 1) THEN
    IF (ISWFINI == 0) THEN
      CALL U_O_MVG(OBOT(JPER), HBOT(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    ELSE
      CALL U_H_MVG(HBOT(JPER), OBOT(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    END IF
  END IF

  RETURN
END SUBROUTINE SWF_GBB_ST


SUBROUTINE SWF_GBB_FM
!**************************************************************************************************
! 计算表示通用底部边界水分运动的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: DELV2
  REAL :: HYCONBOT, HYCON2

  I = NCEL
  J = IMAT(I)

  SELECT CASE (ISWFBOT)

  ! Dirichlet 条件
  CASE (1)
    IF (ISWFINI /= 0) THEN
      CALL U_H_MVG(HBOT(JPER), OBOT(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    END IF
    CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
    HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
    DELV2 = DELV(I) * 0.5
    COF(I) = COF(I) - DELT * HYCON2 / DELV2
    RHS(I) = RHS(I) - DELT * HYCON2 / DELV2 * HBOT(JPER) + DELT * HYCON2

  ! Neumann 条件
  CASE (2)
    RHS(I) = RHS(I) - DELT * WBOT(JPER)

  ! Cauchy 条件
  CASE (3)
    COF(I) = COF(I) - DELT / TLAGBOT
    RHS(I) = RHS(I) - DELT * (HBOT(JPER) - LOCV(I, 2)) / TLAGBOT

  ! 自由排水条件 (Free Drainage)
  CASE (4)
    RHS(I) = RHS(I) + DELT * HYCON(I)

  ! 渗出面条件 (Seepage Face)
  CASE (5)
    DELV2 = DELV(I) * 0.5
    HBOT(JPER) = HNEW(I) + DELV2
    IF (HBOT(JPER) >= 0.0) THEN
      HBOT(JPER) = 0.0
      CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
      HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
      COF(I) = COF(I) - DELT * HYCON2 / DELV2
      RHS(I) = RHS(I) - DELT * (HYCON2 / DELV2 * HBOT(JPER) - HYCON2)
    END IF

  END SELECT

  RETURN
END SUBROUTINE SWF_GBB_FM


SUBROUTINE SWF_GBB_BD
!**************************************************************************************************
! 计算表示底部边界水分运动的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: HYCONBOT, HYCON2
  REAL :: DELV2

  I = NCEL
  J = IMAT(I)

  ! Dirichlet 条件
  SELECT CASE (ISWFBOT)

  CASE (1)
    IF (ISWFINI /= 0) CALL U_H_MVG(HBOT(JPER), OBOT(JPER), ORES(J), OSAT(J), HYPARA(J), HYPARN(J))
    CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
    DELV2 = DELV(I) * 0.5
    HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
    QBOT = DELT * HYCON2 * ((HBOT(JPER) - HNEW(I)) / DELV2 - 1.0)
    WBOT(JPER) = WBOT(JPER) + QBOT

  ! Neumann 条件
  CASE (2)
    QBOT = DELT * WBOT(JPER)

  ! Cauchy 条件
  CASE (3)
    QBOT = DELT * (HBOT(JPER) - LOCV(I, 2) - HNEW(I)) / TLAGBOT
    WBOT(JPER) = WBOT(JPER) + QBOT

  ! 自由排水条件 (Free Drainage)
  CASE (4)
    QBOT = - DELT * HYCON(I)
    WBOT(JPER) = WBOT(JPER) + QBOT

  ! 渗出面条件 (Seepage Face)
  CASE (5)
    IF (HBOT(JPER) >= 0.0) THEN
      CALL U_HYCON_MVG(HYCONBOT, HBOT(JPER), HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
      DELV2 = DELV(I) * 0.5
      HYCON2 = (HYCONBOT + HYCON(I)) * 0.5
      QBOT = DELT * HYCON2 * ((HBOT(JPER) - HNEW(I)) / DELV2 - 1.0)
      WBOT(JPER) = WBOT(JPER) + QBOT
    END IF

  END SELECT

  RETURN
END SUBROUTINE SWF_GBB_BD


SUBROUTINE SWF_GBB_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(HBOT)) DEALLOCATE(HBOT)
  IF (ALLOCATED(OBOT)) DEALLOCATE(OBOT)
  IF (ALLOCATED(WBOT)) DEALLOCATE(WBOT)

  RETURN
end subroutine SWF_GBB_DA

