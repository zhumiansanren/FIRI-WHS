# FIRI-WHS 农田水热盐数值模型 （目前已发布土壤水分运动和热传导子程序包）

FIRI-WHS 农田水热盐数值模型

是一套模拟农田水分运动、热传导和盐分运移的命令行程序包，采用模块化程序设计，使用者可根据需求选用不同的模块，也可修改现有模块或增加新的模块已满足不同模拟需求。

开发人员：黄仲冬、高阳

技术指导：张效先

机构：中国农业科学院农田灌溉研究所，英国洛桑研究所

联系方式: zdhuang@126.com

目前，我们完成了模型基本框架、土壤水分运动、土壤热传导子程序包的整理工作。

盐分运移、灌溉、作物、大气等子程序包已基本实现，但是目前还需要时间进一步整理，待相应工作完成后将逐个发布。

FIRI-WHS 土壤水分运动模块是围绕求解 Richard 方程构建的，在输入土壤水分运动参数、初始条件、边界条件、源汇项等数据后，可输出土壤剖面水势、含水率，地表入渗量、地表径流量、地表蒸发量、作物蒸腾量等数据。

您也可以将 FIRI-WHS 与其他优化程序（如PEST，UCODE）耦合，来求解逆问题。

由于目前该程序包还在缓慢开发、修改和完善中，一些模块没有经过完全测试，请在使用过程中注意。

由于时间和精力有限，目前无法提供程序使用说明，在 TEST 目录下提供了几个实例，供参考。双击 RUN.BAT 或者在 CMD 中运行RUN.BAT。

输入文件：BAS.IN, SOIL.IN, PERIOD.IN, METEO.IN, CROP.IN

输出文件: LIST.OUT, PERIOD.OUT, OBS.OUT, LOG.OUT

图形显示: Plot.R (需要安装 RStudio)

我们在程序代码中提供了大量的注释，标明了每个变量和函数的用途，以便各位老师和同学阅读和修改。非常欢迎提出批评指正，并进一步扩展该程序包。

在模型框架的搭建中，借鉴了 MODFLOW 早期版本的框架构建思路，在程序的编写中，采用了 TTUTIL FORTRAN 库，在此表示感谢！

本程序为自由软件，遵循 MIT 协议，您可以修改程序，但在发布时需要进行说明并表明出处。如有商业用途，请联系作者。

FIRI-WHS 文件结构

  -- BIN 目录，存放可执行文件。

  -- DOC 目录，存在模型文档。

  -- IVF 目录，存放 Visual Studio 2019 和 Intel Visual FORTRAN 2020 工程文件。

  -- REF 目录，存放参考文献。

  -- SRC 目录，存放程序源代码。

  -- TEST 目录，存放测试文件。

更新日志

2020年9月2日，发布 FIRI-WHS 基本框架和土壤水分运动子程序包。

2020年9月8日，发布土壤热传导子程序包。

2020年9月22日，增加了 SWF_SSK 子程序包，结合 UCODE 可用于反求土壤剖面的根系吸水速率。

