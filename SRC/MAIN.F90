!==================================================================================================
! 农田水热盐运移数值模型 FIRI-WHS
! 机构: 中国农业科学院农田灌溉研究所
! 作者: 黄仲冬, 高阳, 张效先
! 邮箱: zdhuang@126.com
!
! 主程序 MAIN
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


PROGRAM MAIN

  USE :: VAR
  IMPLICIT NONE

  ! 局部变量
  REAL :: TSTRT, TEND

  CALL CPU_TIME(TSTRT)

  ! 初始化子程序包控制变量
  L_AMP_BAS = .FALSE.
  L_AMP_PCP = .FALSE.

  L_CRP_BAS = .FALSE.
  L_CRP_PET = .FALSE.
  L_CRP_INT = .FALSE.

  L_IRR_BAS = .FALSE.
  L_IRR_FWA = .FALSE.
  L_IRR_OWA = .FALSE.

  L_SWF_BAS = .FALSE.
  L_SWF_STO = .FALSE.
  L_SWF_CCF = .FALSE.
  L_SWF_GBB = .FALSE.
  L_SWF_GTB = .FALSE.
  L_SWF_INF = .FALSE.
  L_SWF_EVP = .FALSE.
  L_SWF_RWU = .FALSE.
  L_SWF_TMA = .FALSE.

  L_SHT_BAS = .FALSE.
  L_SHT_STO = .FALSE.
  L_SHT_CCT = .FALSE.
  L_SHT_GBB = .FALSE.
  L_SHT_GTB = .FALSE.
  L_SHT_TMA = .FALSE.

  L_SST_BAS = .FALSE.
  L_SST_STO = .FALSE.
  L_SST_CCD = .FALSE.
  L_SST_CCA = .FALSE.
  L_SST_GBB = .FALSE.
  L_SST_GTB = .FALSE.
  L_SST_INF = .FALSE.
  L_SST_TMA = .FALSE.

  !------------------------------------------------------------------------------------------------
  ! 1. Define (DF) - 定义模型, 读取全局参数
  !------------------------------------------------------------------------------------------------
  CALL BAS_DF

  IF (L_AMP_BAS) CALL AMP_BAS_DF
  IF (L_AMP_PCP) CALL AMP_PCP_DF

  IF (L_CRP_BAS) CALL CRP_BAS_DF
  IF (L_CRP_INT) CALL CRP_INT_DF
  IF (L_CRP_PET) CALL CRP_PET_DF

  IF (L_IRR_BAS) CALL IRR_BAS_DF
  IF (L_IRR_FWA) CALL IRR_FWA_DF
  IF (L_IRR_OWA) CALL IRR_OWA_DF

  IF (L_SWF_BAS) CALL SWF_BAS_DF
  IF (L_SWF_STO) CALL SWF_STO_DF
  IF (L_SWF_CCF) CALL SWF_CCF_DF
  IF (L_SWF_GBB) CALL SWF_GBB_DF
  IF (L_SWF_GTB) CALL SWF_GTB_DF
  IF (L_SWF_INF) CALL SWF_INF_DF
  IF (L_SWF_EVP) CALL SWF_EVP_DF
  IF (L_SWF_RWU) CALL SWF_RWU_DF
  IF (L_SWF_TMA) CALL SWF_TMA_DF

  IF (L_SHT_BAS) CALL SHT_BAS_DF
  IF (L_SHT_STO) CALL SHT_STO_DF
  IF (L_SHT_CCT) CALL SHT_CCT_DF
  IF (L_SHT_GBB) CALL SHT_GBB_DF
  IF (L_SHT_GTB) CALL SHT_GTB_DF
  IF (L_SHT_TMA) CALL SHT_TMA_DF

  IF (L_SST_BAS) CALL SST_BAS_DF
  IF (L_SST_STO) CALL SST_STO_DF
  IF (L_SST_CCD) CALL SST_CCD_DF
  IF (L_SST_CCA) CALL SST_CCA_DF
  IF (L_SST_GBB) CALL SST_GBB_DF
  IF (L_SST_GTB) CALL SST_GTB_DF
  IF (L_SST_INF) CALL SST_INF_DF
  IF (L_SST_TMA) CALL SST_TMA_DF

  !------------------------------------------------------------------------------------------------
  ! 2. Allocate (AL) - 为数组变量分配内存并初始化
  !------------------------------------------------------------------------------------------------
  CALL BAS_AL

  IF (L_AMP_BAS) CALL AMP_BAS_AL
  IF (L_AMP_PCP) CALL AMP_PCP_AL
  
  IF (L_CRP_BAS) CALL CRP_BAS_AL
  IF (L_CRP_PET) CALL CRP_PET_AL
  IF (L_CRP_INT) CALL CRP_INT_AL

  IF (L_IRR_BAS) CALL IRR_BAS_AL
  IF (L_IRR_FWA) CALL IRR_FWA_AL
  IF (L_IRR_OWA) CALL IRR_OWA_AL

  IF (L_SWF_BAS) CALL SWF_BAS_AL
  IF (L_SWF_STO) CALL SWF_STO_AL
  IF (L_SWF_CCF) CALL SWF_CCF_AL
  IF (L_SWF_GBB) CALL SWF_GBB_AL
  IF (L_SWF_GTB) CALL SWF_GTB_AL
  IF (L_SWF_INF) CALL SWF_INF_AL
  IF (L_SWF_EVP) CALL SWF_EVP_AL
  IF (L_SWF_RWU) CALL SWF_RWU_AL
  IF (L_SWF_TMA) CALL SWF_TMA_AL

  IF (L_SHT_BAS) CALL SHT_BAS_AL
  IF (L_SHT_STO) CALL SHT_STO_AL
  IF (L_SHT_CCT) CALL SHT_CCT_AL
  IF (L_SHT_GBB) CALL SHT_GBB_AL
  IF (L_SHT_GTB) CALL SHT_GTB_AL
  IF (L_SHT_TMA) CALL SHT_TMA_AL

  IF (L_SST_BAS) CALL SST_BAS_AL
  IF (L_SST_STO) CALL SST_STO_AL
  IF (L_SST_CCD) CALL SST_CCD_AL
  IF (L_SST_CCA) CALL SST_CCA_AL
  IF (L_SST_GBB) CALL SST_GBB_AL
  IF (L_SST_GTB) CALL SST_GTB_AL
  IF (L_SST_INF) CALL SST_INF_AL
  IF (L_SST_TMA) CALL SST_TMA_AL

  !------------------------------------------------------------------------------------------------
  ! 3. Read & Prepare (RP) - 为数组变量读取数据并执行准备工作
  !------------------------------------------------------------------------------------------------
  CALL BAS_RP

  IF (L_AMP_BAS) CALL AMP_BAS_RP
  IF (L_AMP_PCP) CALL AMP_PCP_RP

  IF (L_CRP_BAS) CALL CRP_BAS_RP
  IF (L_CRP_PET) CALL CRP_PET_RP
  IF (L_CRP_INT) CALL CRP_INT_RP

  IF (L_IRR_BAS) CALL IRR_BAS_RP
  IF (L_IRR_FWA) CALL IRR_FWA_RP
  IF (L_IRR_OWA) CALL IRR_OWA_RP
  
  IF (L_SWF_BAS) CALL SWF_BAS_RP
  IF (L_SWF_STO) CALL SWF_STO_RP
  IF (L_SWF_GBB) CALL SWF_GBB_RP
  IF (L_SWF_GTB) CALL SWF_GTB_RP
  IF (L_SWF_INF) CALL SWF_INF_RP
  IF (L_SWF_EVP) CALL SWF_EVP_RP
  IF (L_SWF_RWU) CALL SWF_RWU_RP

  IF (L_SHT_BAS) CALL SHT_BAS_RP
  IF (L_SHT_STO) CALL SHT_STO_RP
  IF (L_SHT_GBB) CALL SHT_GBB_RP
  IF (L_SHT_GTB) CALL SHT_GTB_RP

  IF (L_SST_BAS) CALL SST_BAS_RP
  IF (L_SST_STO) CALL SST_STO_RP
  IF (L_SST_GBB) CALL SST_GBB_RP
  IF (L_SST_GTB) CALL SST_GTB_RP

  WRITE(U_LOG, "()")
  WRITE(U_LOG, "('模拟开始...')")
  WRITE(U_LOG, "('--------------------------------------------------')")
  IF (IRUNMOD == 1) THEN
    WRITE(*, "('模拟开始...')")
    WRITE(*, "('--------------------------------------------------')")
  END IF

  DO JPER = 1, NPER

    !----------------------------------------------------------------------------------------------
    ! 4. Stress (ST) - 更新随应力期变化的数据
    !----------------------------------------------------------------------------------------------
    CALL BAS_ST

    IF (L_AMP_BAS) CALL AMP_BAS_ST
    IF (L_AMP_PCP) CALL AMP_PCP_ST

    IF (L_CRP_BAS) CALL CRP_BAS_ST
    IF (L_CRP_INT) CALL CRP_INT_ST
    IF (L_CRP_PET) CALL CRP_PET_ST
    
    IF (L_IRR_BAS) CALL IRR_BAS_ST
    IF (L_IRR_FWA) CALL IRR_FWA_ST
    IF (L_IRR_OWA) CALL IRR_OWA_ST

    IF (L_SWF_GBB) CALL SWF_GBB_ST
    IF (L_SWF_GTB) CALL SWF_GTB_ST
    IF (L_SWF_RWU) CALL SWF_RWU_ST

    DO JSTP = 1, NSTP(JPER) ! 时段循环开始

      !--------------------------------------------------------------------------------------------
      ! 5. Advance (AD) - 更新随时段变化的数据
      !--------------------------------------------------------------------------------------------
      CALL BAS_AD

      IF (L_SWF_STO) CALL SWF_STO_AD
      IF (L_SHT_STO) CALL SHT_STO_AD
      IF (L_SST_STO) CALL SST_STO_AD

      DO JITR = 1, MAXITR

        !------------------------------------------------------------------------------------------
        ! 6. Formulate (FM) - 计算土壤水分运动方程组系数
        !------------------------------------------------------------------------------------------
        
        IF (L_SWF_BAS) CALL SWF_BAS_FM
        IF (L_SWF_STO) CALL SWF_STO_FM
        IF (L_SWF_CCF) CALL SWF_CCF_FM
        IF (L_SWF_GBB) CALL SWF_GBB_FM
        IF (L_SWF_GTB) CALL SWF_GTB_FM
        IF (L_SWF_INF) CALL SWF_INF_FM
        IF (L_SWF_EVP) CALL SWF_EVP_FM
        IF (L_SWF_RWU) CALL SWF_RWU_FM

        IF (L_IRR_FWA) CALL IRR_FWA_FM
        IF (L_IRR_OWA) CALL IRR_OWA_FM

        !------------------------------------------------------------------------------------------
        ! 7. Approximate (AP) - 求解土壤水分运动方程组
        !------------------------------------------------------------------------------------------
        IF (L_SWF_TMA) CALL SWF_TMA_AP

        ! 如何解收敛，就退出迭代循环
        IF (LCONV) EXIT

      END DO

      !--------------------------------------------------------------------------------------------
      ! 8. Budget (BD) - 土壤水分运动均衡计算
      !--------------------------------------------------------------------------------------------
      IF (L_SWF_STO) CALL SWF_STO_BD
      IF (L_SWF_CCF) CALL SWF_CCF_BD
      IF (L_SWF_GBB) CALL SWF_GBB_BD
      IF (L_SWF_GTB) CALL SWF_GTB_BD
      IF (L_SWF_INF) CALL SWF_INF_BD
      IF (L_SWF_EVP) CALL SWF_EVP_BD
      IF (L_SWF_RWU) CALL SWF_RWU_BD

      !--------------------------------------------------------------------------------------------
      ! 6. Formulate (FM) - 计算土壤热传导方程组系数
      !--------------------------------------------------------------------------------------------
      IF (L_SHT_BAS) CALL SHT_BAS_FM
      IF (L_SHT_STO) CALL SHT_STO_FM
      IF (L_SHT_CCT) CALL SHT_CCT_FM
      IF (L_SHT_GBB) CALL SHT_GBB_FM
      IF (L_SHT_GTB) CALL SHT_GTB_FM

      !------------------------------------------------------------------------------------------
      ! 7. Approximate (AP) - 求解土壤热传导方程组
      !------------------------------------------------------------------------------------------
      IF (L_SHT_TMA) CALL SHT_TMA_AP

      !--------------------------------------------------------------------------------------------
      ! 8. Budget (BD) - 土壤热传导均衡计算
      !--------------------------------------------------------------------------------------------
      IF (L_SHT_STO) CALL SHT_STO_BD
      IF (L_SHT_CCT) CALL SHT_CCT_BD
      IF (L_SHT_GBB) CALL SHT_GBB_BD
      IF (L_SHT_GTB) CALL SHT_GTB_BD

      DO JITR = 1, MAXITR

        !------------------------------------------------------------------------------------------
        ! 6. Formulate (FM) - 计算土壤盐分运移方程组系数
        !------------------------------------------------------------------------------------------
        IF (L_SST_BAS) CALL SST_BAS_FM
        IF (L_SST_STO) CALL SST_STO_FM
        IF (L_SST_CCD) CALL SST_CCD_FM
        IF (L_SST_CCA) CALL SST_CCA_FM
        IF (L_SST_GBB) CALL SST_GBB_FM
        IF (L_SST_GTB) CALL SST_GTB_FM
        IF (L_SST_INF) CALL SST_INF_FM

        !------------------------------------------------------------------------------------------
        ! 7. Approximate (AP) - 求解土壤盐分运移方程组
        !------------------------------------------------------------------------------------------
        IF (L_SST_TMA) CALL SST_TMA_AP

        ! 如何解收敛，就退出迭代循环
        IF (LCONV) EXIT

      END DO

      !--------------------------------------------------------------------------------------------
      ! 8. Budget (BD) - 土壤盐分运移均衡计算
      !--------------------------------------------------------------------------------------------
      IF (L_SST_STO) CALL SST_STO_BD
      IF (L_SST_CCD) CALL SST_CCD_BD
      IF (L_SST_CCA) CALL SST_CCA_BD
      IF (L_SST_GBB) CALL SST_GBB_BD
      IF (L_SST_GTB) CALL SST_GTB_BD

    END DO

    !----------------------------------------------------------------------------------------------
    ! 9. Output (OT) - 输出模拟结果
    !----------------------------------------------------------------------------------------------
      CALL BAS_OT

    IF ((IERRCOD /= 0) .AND. (IRUNMOD /= 3)) THEN
      IF (IRUNMOD == 1) WRITE(*, "('程序停止运行!')")
      EXIT
    END IF

  END DO

  CALL CPU_TIME(TEND)
  WRITE(U_LOG, "('----------------------------------------')")
  WRITE(U_LOG, "('模拟结束，用时 ', ES10.3, ' 秒')") TEND - TSTRT
  WRITE(U_LOG, "()")
  IF (IRUNMOD == 1) THEN
    WRITE(*, "('----------------------------------------')")
    WRITE(*, "('模拟结束，用时 ', ES10.3, ' 秒')") TEND - TSTRT
    WRITE(*, *)
  END IF

  !------------------------------------------------------------------------------------------------
  ! 10. Deallocate (DA) 释放数组变量内存
  !------------------------------------------------------------------------------------------------
  CALL BAS_DA

  IF (L_AMP_BAS) CALL AMP_BAS_DA
  IF (L_AMP_PCP) CALL AMP_PCP_DA

  IF (L_CRP_BAS) CALL CRP_BAS_DA
  IF (L_CRP_INT) CALL CRP_INT_DA
  IF (L_CRP_PET) CALL CRP_PET_DA

  IF (L_IRR_BAS) CALL IRR_BAS_DA
  IF (L_IRR_FWA) CALL IRR_FWA_DA
  IF (L_IRR_OWA) CALL IRR_OWA_DA

  IF (L_SWF_BAS) CALL SWF_BAS_DA
  IF (L_SWF_STO) CALL SWF_STO_DA
  IF (L_SWF_CCF) CALL SWF_CCF_DA
  IF (L_SWF_GBB) CALL SWF_GBB_DA
  IF (L_SWF_GTB) CALL SWF_GTB_DA
  IF (L_SWF_INF) CALL SWF_INF_DA
  IF (L_SWF_EVP) CALL SWF_EVP_DA
  IF (L_SWF_RWU) CALL SWF_RWU_DA
  IF (L_SWF_TMA) CALL SWF_TMA_DA

  IF (L_SHT_BAS) CALL SHT_BAS_DA
  IF (L_SHT_STO) CALL SHT_STO_DA
  IF (L_SHT_CCT) CALL SHT_CCT_DA
  IF (L_SHT_GBB) CALL SHT_GBB_DA
  IF (L_SHT_GTB) CALL SHT_GTB_DA
  IF (L_SHT_TMA) CALL SHT_TMA_DA

  IF (L_SST_BAS) CALL SST_BAS_DA
  IF (L_SST_STO) CALL SST_STO_DA
  IF (L_SST_CCD) CALL SST_CCD_DA
  IF (L_SST_CCA) CALL SST_CCA_DA
  IF (L_SST_GBB) CALL SST_GBB_DA
  IF (L_SST_GTB) CALL SST_GTB_DA
  IF (L_SST_INF) CALL SST_INF_DA
  IF (L_SST_TMA) CALL SST_TMA_DA

  IF (IRUNMOD == 1) THEN
    WRITE(*, "('请按回车键退出...')")
    READ(*, *)
  END IF

  STOP
END PROGRAM MAIN

