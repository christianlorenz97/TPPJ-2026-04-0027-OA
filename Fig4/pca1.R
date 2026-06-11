### Correlation ###  
# Packages
library(Matrix)
library(lme4)
library(lmerTest)
library(effects)
library(plotrix)
library(gridExtra)
library(piecewiseSEM)
library(nlme)
library(pgirmess)
library(splines)
library(car)
library(ISLR)
library(lattice)
library(latticeExtra)
library(HH)
library(multcomp)
library(openxlsx)
library(multcompView)
library(dplyr)
library(ggplot2)
library(ARTool)
library(sciplot)
library(MASS)
library(nortest)
library(performance)
library(ggpubr)

# Correlation scatterplot with corrplot
library(corrplot)

data_mv<-read.delim("mv_matrix.txt")
data_mv$Treatment<-ordered(data_mv$Treatment, levels=c("NW_MIX1","NW_MIXR","NW_NOI","WW_MIX1","WW_MIXR","WW_NOI"))

cor(data_mv[,2:15]) #default pearson 

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(data_mv[,2:15], method="pearson"), method="circle", col=col(200),  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label coslor and rotation
)

#remove manually correlated variables with R^2>0.8

# Correlation scatterplot with corrplot
data_mv<-read.delim("mv_matrix3.txt")

cor(data_mv[,2:9]) #default pearson

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(data_mv[,2:9], method="pearson"), method="circle", col=col(200),  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label coslor and rotation
)

### PCA ###
# Packages
library(factoextra)
library(vegan)
library(permute)
library(lattice)
library(dplyr)
library(ggplot2)

# Import
data_mv1 <- read.delim("mv_matrix3.txt")
data_mv1$Treatment <- ordered(
  data_mv1$Treatment,
  levels = c("NW_MIX1","NW_MIXR","NW_NOI","WW_MIX1","WW_MIXR","WW_NOI")
)

# PCA
pca <- prcomp(data_mv1[2:9], center = TRUE, scale. = TRUE)
biplot(pca)

# Customized ellipsoids
ind_df <- as.data.frame(pca$x[, 1:2])
colnames(ind_df) <- c("Dim.1", "Dim.2")
ind_df$Treatment <- data_mv1$Treatment

ind_df <- ind_df %>%
  mutate(
    EllType = case_when(
      grepl("_NOI$",  Treatment)  ~ "solid",
      grepl("_MIX1$", Treatment)  ~ "dashed",
      grepl("_MIXR$", Treatment)  ~ "dotted",
      TRUE ~ "solid"
    )
  )

# Color palette
pal_trt <- c("#CC6677","#AA4499","#882255","#44AA99","#117733","#332288")

# Biplot
PCAc <- fviz_pca_biplot(
  pca,
  label = "var",
  col.var = "black",
  habillage = data_mv1$Treatment,          # legend
  palette = pal_trt,
  addEllipses = FALSE
) +
  scale_shape_manual(values = c(0, 1, 2, 15, 16, 17)) +
  # filled ellipsoids, no legend
  stat_ellipse(
    data = ind_df,
    aes(
      x = Dim.1, y = Dim.2,
      color = Treatment,
      fill  = Treatment,
      linetype = EllType,
      group = Treatment
    ),
    level = 0.95,
    type = "t",
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pal_trt) +
  scale_color_manual(values = pal_trt) +
  scale_linetype_identity() +
  # legend “Groups” -> “Treatments”
  labs(color = "Treatments", shape = "Treatments")

print(PCAc)

### MANOVA One-Way ###
data_mv1 <- read.delim("mv_matrix3.txt")
var_num<-cbind(data_mv1$fRAE,
               data_mv1$fSH,
               data_mv1$fRSA,
               data_mv1$RAER,
               data_mv1$SER,
               data_mv1$SSAER,
               data_mv1$rFW,
               data_mv1$sFW
               )
pca.manova<- manova(var_num~data_mv1$Treatment) #one-way multiviariate analysis of variance
summary(pca.manova, tol=0) #if rank deficiency o rank x<x+1 error, use summary(pca.manova, tol=0)

# MANOVA - multivariate vs Treatment
fit_man <- manova(
  cbind(fRAE, fSH, fRSA, RAER, SER, SSAER, rFW, sFW) ~ Treatment,
  data = data_mv1
)

summary(fit_man, test = "Pillai")    # more robust test
summary.aov(fit_man)                 # ANOVA univariate for each variable



### MANOVA Two-Way ###
data_mv1 <- read.delim("mv_matrix3_2way.txt")
var_num<-cbind(data_mv1$fRAE,
               data_mv1$fSH,
               data_mv1$fRSA,
               data_mv1$RAER,
               data_mv1$SER,
               data_mv1$SSAER,
               data_mv1$rFW,
               data_mv1$sFW
)
pca.manova<- manova(var_num~data_mv1$Regime*data_mv1$Inoculant) #one-way multiviariate analysis of variance
summary(pca.manova, tol=0) #if rank deficiency o rank x<x+1 error, use summary(pca.manova, tol=0)

# MANOVA — multivariate vs Treatment
fit_man <- manova(
  cbind(fRAE, fSH, fRSA, RAER, SER, SSAER, rFW, sFW) ~ Regime * Inoculant,
  data = data_mv1
)

summary(fit_man, test = "Pillai")    # more robust test
summary.aov(fit_man)                 # ANOVA univariate for each variable
