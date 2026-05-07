### File to analyze outcome of all factor analyses

library(readxl)
transposed_matrix_studies_final <- read_excel("transposed matrix studies final.xlsx")
Transposed_matrix <- transposed_matrix_studies_final
colnames(Transposed_matrix)
table(Transposed_matrix[,7])
rownames(Transposed_matrix) <- Transposed_matrix$Study
colnames(Transposed_matrix)

### Look at descriptives
library(psych)
descriptive <- describe(Transposed_matrix)
descriptive
write.table(descriptive,"descrip_stats_text",sep = "\t")

###
### Look at densities for best fitting models
###

transposed_matrix_studies2 <- as.data.frame(Transposed_matrix[,7:92])
colnames(transposed_matrix_studies2)
colnames(transposed_matrix_studies2)[-1] <- paste0('x', 1:(ncol(transposed_matrix_studies2)-1))
transposed_matrix_studies2 <-subset(transposed_matrix_studies2,Type_of_scale=="uni")

library(ggplot2)
library(dplyr)
library(tidyr)


###################
### Figure 1 (unifactorial scales)
###################

## TLI

theme_set(theme_minimal())
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x8, 
                 Psych.polychor = transposed_matrix_studies2$x12,
                 Lavaan.ordered = transposed_matrix_studies2$x34,
                 Lavaan.MLR =transposed_matrix_studies2$x38)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0.4, 1.0)) +
  labs(title="Density plots for unifactorial scales", x="TLI") +
  geom_segment(x = .95,y=0, xend = .95, yend = 6,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## RMSEA
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x9, 
                 Psych.polychor = transposed_matrix_studies2$x13,
                 Lavaan.ordered = transposed_matrix_studies2$x35,
                 Lavaan.MLR =transposed_matrix_studies2$x39)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0, 0.4)) +
  labs(title="Density plots for unifactorial scales", x="RMSEA") +
  geom_segment(x = .06,y=0, xend = .06, yend = 12.5,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## SRMR
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x10, 
                 Psych.polychor = transposed_matrix_studies2$x14,
                 Lavaan.ordered = transposed_matrix_studies2$x36,
                 Lavaan.MLR =transposed_matrix_studies2$x40)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0, 0.3)) +
  labs(title="Density plots for unifactorial scales", x="SRMR") +
  geom_segment(x = .08,y=0, xend = .08, yend = 30,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## Percent variance accounted for
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x7, 
                 Psych.polychor = transposed_matrix_studies2$x11,
                 Lavaan.ordered = transposed_matrix_studies2$x33,
                 Lavaan.MLR =transposed_matrix_studies2$x37)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0.1, 1.0)) +
  labs(title="Density plots for unifactorial scales", x="Percent variance") +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))


###################
### Figure 2 (multifactorial scales)
###################

transposed_matrix_studies2 <- as.data.frame(Transposed_matrix[,7:92])
colnames(transposed_matrix_studies2)
colnames(transposed_matrix_studies2)[-1] <- paste0('x', 1:(ncol(transposed_matrix_studies2)-1))
transposed_matrix_studies2 <-subset(transposed_matrix_studies2,Type_of_scale=="multi")

library(ggplot2)
library(dplyr)
library(tidyr)

## TLI
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x16, 
                 Psych.polychor = transposed_matrix_studies2$x25,
                 Lavaan.ordered = transposed_matrix_studies2$x42,
                 Lavaan.MLR =transposed_matrix_studies2$x51)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0.4, 1.0)) +
  labs(title="Density plots for multifactorial scales", x="TLI") +
  geom_segment(x = .95,y=0, xend = .95, yend = 6,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## RMSEA
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x17, 
                 Psych.polychor = transposed_matrix_studies2$x26,
                 Lavaan.ordered = transposed_matrix_studies2$x43,
                 Lavaan.MLR =transposed_matrix_studies2$x52)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0, 0.4)) +
  labs(title="Density plots for multifactorial scales", x="RMSEA") +
  geom_segment(x = .06,y=0, xend = .06, yend = 25,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## SRMR
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x18, 
                 Psych.polychor = transposed_matrix_studies2$x27,
                 Lavaan.ordered = transposed_matrix_studies2$x44,
                 Lavaan.MLR =transposed_matrix_studies2$x53)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0, 0.3)) +
  labs(title="Density plots for multifactoria scales", x="SRMR") +
  geom_segment(x = .08,y=0, xend = .08, yend = 30,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## Percent variance accounted for
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x15, 
                 Psych.polychor = transposed_matrix_studies2$x24,
                 Lavaan.ordered = transposed_matrix_studies2$x41,
                 Lavaan.MLR =transposed_matrix_studies2$x50)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0.1, 1.0)) +
  labs(title="Density plots for multifactorial scales", x="Percent variance") +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))


###################
### Figure 3 (hierarchical scales)
###################

### Look at densities for best fitting models
### hierarchical scales
transposed_matrix_studies2 <- as.data.frame(Transposed_matrix[,7:92])
colnames(transposed_matrix_studies2)
colnames(transposed_matrix_studies2)[-1] <- paste0('x', 1:(ncol(transposed_matrix_studies2)-1))
transposed_matrix_studies2 <-subset(transposed_matrix_studies2,Type_of_scale=="hier")

library(ggplot2)
library(dplyr)
library(tidyr)

## TLI
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x16, 
                 Psych.polychor = transposed_matrix_studies2$x25,
                 Lavaan.ordered = transposed_matrix_studies2$x42,
                 Lavaan.MLR =transposed_matrix_studies2$x51)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0.4, 1.0)) +
  labs(title="Density plots for hierarchical scales", x="TLI") +
  geom_segment(x = .95,y=0, xend = .95, yend = 6,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## RMSEA
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x17, 
                 Psych.polychor = transposed_matrix_studies2$x26,
                 Lavaan.ordered = transposed_matrix_studies2$x43,
                 Lavaan.MLR =transposed_matrix_studies2$x52)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0, 0.4)) +
  labs(title="Density plots for hierarchical scales", x="RMSEA") +
  geom_segment(x = .06,y=0, xend = .06, yend = 25,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## SRMR
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x18, 
                 Psych.polychor = transposed_matrix_studies2$x27,
                 Lavaan.ordered = transposed_matrix_studies2$x44,
                 Lavaan.MLR =transposed_matrix_studies2$x53)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0, 0.3)) +
  labs(title="Density plots for hierarchical scales", x="SRMR") +
  geom_segment(x = .08,y=0, xend = .08, yend = 30,color="red",size=1.5) +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))

## Percent variance accounted for
xy <- data.frame(Psych.pearson = transposed_matrix_studies2$x15, 
                 Psych.polychor = transposed_matrix_studies2$x24,
                 Lavaan.ordered = transposed_matrix_studies2$x41,
                 Lavaan.MLR =transposed_matrix_studies2$x50)
xy <- gather(xy)
ggplot(xy, aes(x = value, fill = key)) +
  geom_density(alpha = 0.2) +
  expand_limits(x = c(0.1, 1.0)) +
  labs(title="Density plots for hierarchical scales", x="Percent variance") +
  theme(plot.title=element_text(size=25),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18))



######################################
### Code to make Figure 4 - 9
#### Show distributions three types of scales
######################################


transposed_matrix_studies2 <- as.data.frame(Transposed_matrix[,7:92])
Names = colnames(transposed_matrix_studies2)
colnames(transposed_matrix_studies2)[-1] <- paste0('c', 1:(ncol(transposed_matrix_studies2)-1))
transposed_matrix_studies2$c73 <- log(transposed_matrix_studies2$c68/transposed_matrix_studies2$c70)

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

Names
measurement_variable = 21
xy <- transposed_matrix_studies2[,c(1,measurement_variable)]
colnames(xy) <-c("Scale",'Predict')
namevar <- Names[measurement_variable]
ggdensity(xy, x = "Predict", fill = "Scale", legend="right") +
  labs(x=namevar)


#############################################################
#### Decision tree to see which variables are most important
#### For classification
#############################################################

transposed_matrix_studies2 <- as.data.frame(Transposed_matrix[,7:92])

# Give overview of types of scale
table(transposed_matrix_studies2[,1])

# take away eigenvalues
transposed_matrix_studies3 <-transposed_matrix_studies2[,c(1,6:86)]

colnames(transposed_matrix_studies3)
library(rpart)
rfit = rpart(Type_of_scale ~ ., 
             data = transposed_matrix_studies3, 
             method = "class",
             maxdepth=2,
             xval = 10)
print(rfit)
library(rpart.plot)
rpart.plot(rfit)

pred = predict(rfit, type="class")
table(pred,transposed_matrix_studies2$Type_of_scale)

library(tidyverse)
df <- data.frame(imp = rfit$variable.importance)
df2 <- df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
plot(rfit$variable.importance)
ggplot2::ggplot(df2) +
  geom_segment(aes(x = variable, y = 0, xend = variable, yend = imp), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = imp, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw()

plotcp(rfit)
printcp(rfit)



###############################################
# Use random forest to find best classification
###############################################

library(randomForest)
library("MASS")

transposed_matrix_studies3 <- as.data.frame(Transposed_matrix[,c(7,12:92)])
colnames(transposed_matrix_studies3)
names(transposed_matrix_studies3)[-1] <- paste0('x', 1:(ncol(transposed_matrix_studies3)-1))
transposed_matrix_studies3$Type_of_scale <- as.factor(transposed_matrix_studies3$Type_of_scale)

library(rsample)
data_split <- initial_split(transposed_matrix_studies3, prop = 0.75)
train_data <- training(data_split)
test_data <- testing(data_split)


fit_bag1 <- randomForest(Type_of_scale ~ ., data = train_data)
p1 <- predict(fit_bag1)
table(p1,train_data$Type_of_scale)
varImpPlot(fit_bag1)

fit_bag2 <- randomForest(Type_of_scale ~ ., data = test_data)
p1 <- predict(fit_bag2)
table(p1,test_data$Type_of_scale)


