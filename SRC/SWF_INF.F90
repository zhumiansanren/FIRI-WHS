!==================================================================================================
! 子程序包: SWF_INF (Infiltration Pacakge for Soil Water Flow Process)
! 主要功能: 模拟地表水分入渗
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_INF_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('HPONDMAX')) THEN
    CALL RDSREA('HPONDMAX', HPONDMAX)
  ELSE
    HPONDMAX = 0.0
  END IF

  RETURN
END SUBROUTINE SWF_INF_DF


SUBROUTINE SWF_INF_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(WSURF)) ALLOCATE(WSURF(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WRNOF)) ALLOCATE(WRNOF(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WINF))  ALLOCATE(WINF(NPER),  SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_INF_AL


SUBROUTINE SWF_INF_RP
!**************************************************************************************************
! 为数组变量读取数据并执行准备工作
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WSURF')) THEN
    CALL RDFREA('WSURF', WSURF, NPER, NPER)
  END IF

  RETURN
END SUBROUTINE SWF_INF_RP


SUBROUTINE SWF_INF_AD
!**************************************************************************************************
! 保存前一时段结束时的地表积水深度
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  HPONDNEW = HPONDOLD

  RETURN
END SUBROUTINE SWF_INF_AD


SUBROUTINE SWF_INF_FM
!**************************************************************************************************
! 计算表示地表水分入渗的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: HYCONPOND, HYCON1
  REAL    :: DELV1
  REAL    :: QINFMAX                   ! 当前时段地表最大允许入渗水量 (cm)

  I = 1
  J = IMAT(I)

  ! 计算地表积水深度
  CALL U_HYCON_MVG(HYCONPOND, HPONDNEW, HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
  HYCON1 = (HYCONPOND + HYCON(I)) * 0.5
  DELV1 = DELV(I) * 0.5
  QINFMAX = DELT * HYCON1 * ((HPONDNEW - HNEW(I)) / DELV1 + 1.0)
  HPONDNEW = HPONDOLD + DELT * WSURF(JPER) - QINFMAX

  ! 计算地表径流水量和入渗水量
  QRNOF = 0.0
  IF (HPONDNEW > 0.0) THEN ! 地表产生积水
    IF (HPONDNEW > HPONDMAX) THEN ! 地表产生径流
      QRNOF = HPONDNEW - HPONDMAX
      HPONDNEW = HPONDMAX
    END IF
    CALL U_HYCON_MVG(HYCONPOND, HPONDNEW, HYCONSAT(J), HYPARA(J), HYPARN(J), HYPARL(J))
    HYCON1 = (HYCONPOND + HYCON(I)) * 0.5
    DELV1 = DELV(I) * 0.5
    COF(I) = COF(I) - DELT * HYCON1 / DELV1
    RHS(I) = RHS(I) - DELT * HYCON1 * (HPONDNEW / DELV1 + 1.0)
    QINF = QINFMAX
  ELSE ! 地表无积水
    RHS(I) = RHS(I) - DELT * WSURF(JPER)
    QINF = DELT * WSURF(JPER)
    HPONDNEW = 0.0
  END IF

  RETURN
END SUBROUTINE SWF_INF_FM


SUBROUTINE SWF_INF_BD
!**************************************************************************************************
! 计算表示地表水分入渗的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  WRNOF(JPER) = WRNOF(JPER) + QRNOF
  WINF(JPER) = WINF(JPER) + QINF

  RETURN
END SUBROUTINE SWF_INF_BD


SUBROUTINE SWF_INF_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(WSURF)) DEALLOCATE(WSURF)
  IF (ALLOCATED(WRNOF)) DEALLOCATE(WRNOF)
  IF (ALLOCATED(WINF))  DEALLOCATE(WINF)

  RETURN
END SUBROUTINE SWF_INF_DA

