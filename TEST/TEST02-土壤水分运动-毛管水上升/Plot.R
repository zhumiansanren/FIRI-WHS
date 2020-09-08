library(data.table)
library(ggplot2)

data <- fread(file = "PERIOD.OUT", header = T, skip = 9)
ggplot(data, aes(x = TIME, y = WBOT/PERLEN)) +
  theme_bw() + 
  geom_line(size = 1.0, colour  = "steelblue") +
  xlab("时间 (min)") + 
  ylab("毛管水上升速率 (cm/d)")

ggplot(data, aes(x = TIME, y = C_WBOT)) +
  theme_bw() + 
  geom_line(size = 1.0, colour  = "steelblue") +
  xlab("时间 (d)") + 
  ylab("毛管水累计上升量 (cm)")

data <- fread(file = "LIST.OUT", header = T, skip = 10)
ggplot(data, aes(x = MOIST, y = DEPTH, color = factor(TIME))) +
  theme_bw() + 
  geom_line() +
  scale_y_reverse() +
  xlab("含水率 (cm3/cm3)") + 
  ylab("深度 (cm)")

data <- fread(file = "OBS.OUT", header = T, skip = 9)
ggplot(data, aes(x = TIME, y = HEAD, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("水势 (cm)")

ggplot(data, aes(x = TIME, y = MOIST, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("含水率 (cm3/cm3)")

