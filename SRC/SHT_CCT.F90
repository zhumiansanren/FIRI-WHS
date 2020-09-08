!==================================================================================================
! 子程序包: SHT_CCT (Cell-Centered Transport Package for Soil Heat Transport Process)
! 主要功能: 模拟土壤单元体之间的热传导
! 创建日期: 2020年7月7日
! 修改日志:
!==================================================================================================


SUBROUTINE SHT_CCT_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SHT_CCT_DF


SUBROUTINE SHT_CCT_AL
!**************************************************************************************************
! 为数组变量分配内存并初始化
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(CCV))  ALLOCATE(CCV(NCEL, 2),  SOURCE = 0.0)
  IF (.NOT. ALLOCATED(QTHV)) ALLOCATE(QTHV(NCEL, 2), SOURCE = 0.0)

  RETURN
END SUBROUTINE SHT_CCT_AL


SUBROUTINE SHT_CCT_FM
!**************************************************************************************************
! 计算表示土壤单元体之间热传导的方程组系数
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: DELV1, DELV2
  REAL    :: THCON1, THCON2

  ! 第 2 至 NCEL - 1 土壤单元体
  DO I = 2, NCEL - 1
    DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
    DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
    THCON1 = (THCON(I - 1) * DELV(I - 1) + THCON(I) * DELV(I)) / (DELV(I - 1) + DELV(I))
    THCON2 = (THCON(I + 1) * DELV(I + 1) + THCON(I) * DELV(I)) / (DELV(I + 1) + DELV(I))
    CCV(I, 1) = DELT * THCON1 / DELV1
    CCV(I, 2) = DELT * THCON2 / DELV2
  END DO

  ! 第 1 土壤单元体
  I = 1
  DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
  THCON2 = (THCON(I + 1) * DELV(I + 1) + THCON(I) * DELV(I)) / (DELV(I + 1) + DELV(I))
  CCV(I, 2) = DELT * THCON2 / DELV2

  ! 第 NCEL 土壤单元体
  I = NCEL
  DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
  THCON1 = (THCON(I - 1) * DELV(I - 1) + THCON(I) * DELV(I)) / (DELV(I - 1) + DELV(I))
  CCV(I, 1) = DELT * THCON1 / DELV1

  RETURN
END SUBROUTINE SHT_CCT_FM


SUBROUTINE SHT_CCT_BD
!**************************************************************************************************
! 计算表示土壤单元体之间热传导的均衡项
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV1, DELV2
  REAL    :: THCON1, THCON2

  ! 第 2 至 NCEL - 1 单元体
  DO I = 2, NCEL - 1
    DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
    DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
    THCON1 = (THCON(I - 1) * DELV(I - 1) + THCON(I) * DELV(I)) / (DELV(I - 1) + DELV(I))
    THCON2 = (THCON(I + 1) * DELV(I + 1) + THCON(I) * DELV(I)) / (DELV(I + 1) + DELV(I))
    QTHV(I, 1) = 0.5 * DELT * THCON1 * ((TNEW(I - 1) - TNEW(I)) / DELV1) + &
                 0.5 * DELT * THCON1 * ((TOLD(I - 1) - TOLD(I)) / DELV1)
    QTHV(I, 2) = 0.5 * DELT * THCON2 * ((TNEW(I + 1) - TNEW(I)) / DELV2) + &
                 0.5 * DELT * THCON2 * ((TOLD(I + 1) - TOLD(I)) / DELV2)
  END DO

  ! 第 1 单元体
  I = 1
  DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
  THCON2 = (DELV(I + 1) * THCON(I + 1) + DELV(I) * THCON(I)) / (DELV(I + 1) + DELV(I))
  QTHV(I, 2) = 0.5 * DELT * THCON2 * ((TNEW(I + 1) - TNEW(I)) / DELV2) + &
               0.5 * DELT * THCON2 * ((TOLD(I + 1) - TOLD(I)) / DELV2)

  ! 第 NCEL 单元体
  I = NCEL
  DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
  THCON1 = (DELV(I - 1) * THCON(I - 1) + DELV(I) * THCON(I)) / (DELV(I - 1) + DELV(I))
  QTHV(I, 1) = 0.5 * DELT * THCON1 * ((TNEW(I - 1) - TNEW(I)) / DELV1) + &
               0.5 * DELT * THCON1 * ((TOLD(I - 1) - TOLD(I)) / DELV1)

  RETURN
END SUBROUTINE SHT_CCT_BD


SUBROUTINE SHT_CCT_DA
!**************************************************************************************************
! 释放数组变量内存
!**************************************************************************************************
  USE :: VAR
  IMPLICIT NONE

  IF (ALLOCATED(CCV))  DEALLOCATE(CCV)
  IF (ALLOCATED(QTHV)) DEALLOCATE(QTHV)

  RETURN
END SUBROUTINE SHT_CCT_DA

