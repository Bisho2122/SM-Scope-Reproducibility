library(tidyverse)
library(ggpubr)

# Load data ---------------------------------------------------------------
RMSE_Res = read.csv("../../../Data/Figure 1/Atlas_sections_RMSE_results.csv")

# Plot --------------------------------------------------------------------
RMSE_Res = RMSE_Res[order(RMSE_Res$Section),]

RMSE_Res$Section = gsub("Section_", "Plane.",RMSE_Res$Section)
RMSE_Res$section_num = sub(".*Plane.","",RMSE_Res$Section) %>% as.numeric()
RMSE_Res$section_num = RMSE_Res$section_num + 1
RMSE_Res$Section = paste0("Plane.", RMSE_Res$section_num)

RMSE_Res %>%
  ggplot(aes(x = Section, y = RMSE, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_pubclean() +
  rotate_x_text(45) +
  xlab("") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  ylim(c(200,1000))
