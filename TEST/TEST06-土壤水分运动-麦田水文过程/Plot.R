library(data.table)
library(ggplot2)

data <- fread(file = "PERIOD.OUT", header = T, skip = 9)
ggplot(data) +
  theme_bw() + 
  geom_line(aes(x = TIME, y = WETPOT/PERLEN), size = 1.0, colour = "red") +
  geom_line(aes(x = TIME, y = WEPOT/PERLEN),  size = 1.0, colour = "black", linetype = "dashed") +
  geom_line(aes(x = TIME, y = WEACT/PERLEN),  size = 1.0, colour = "black") +
  geom_line(aes(x = TIME, y = WTPOT/PERLEN),  size = 1.0, colour = "blue", linetype = "dashed") +
  geom_line(aes(x = TIME, y = WTACT/PERLEN),  size = 1.0, colour = "blue") +
  xlab("时间 (d)") + 
  ylab("蒸发与蒸腾速率 (cm d-1)")

ggplot(data) +
  theme_bw() + 
  geom_line(aes(x = TIME, y = C_WEPOT), size = 1.0, colour  = "black", linetype = "dashed") +
  geom_line(aes(x = TIME, y = C_WEACT), size = 1.0, colour  = "black") +
  geom_line(aes(x = TIME, y = C_WTPOT), size = 1.0, colour  = "blue", linetype = "dashed") +
  geom_line(aes(x = TIME, y = C_WTACT), size = 1.0, colour  = "blue") +
  xlab("时间 (min)") + 
  ylab("累计蒸发与蒸腾量 (cm)")

data <- fread(file = "LIST.OUT", header = T, skip = 10)
ggplot(data[TIME == 50], aes(x = MOIST, y = DEPTH, color = factor(TIME))) +
  theme_bw() + 
  geom_point() +
  scale_y_reverse() +
  xlab("含水率 (cm3/cm3)") + 
  ylab("深度 (cm)")

data <- fread(file = "OBS.OUT", header = T, skip = 9)
ggplot(data[DEPTH == 9], aes(x = TIME, y = HEAD, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("水势 (cm)")

ggplot(data[DEPTH == 9], aes(x = TIME, y = MOIST, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("含水率 (cm3/cm3)")

