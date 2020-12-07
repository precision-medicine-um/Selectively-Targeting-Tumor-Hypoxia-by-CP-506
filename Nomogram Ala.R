library(rms)
library(missForest)
library(pROC)
library(shiny)
library(precrec)
library(caret)
library(e1071)
library(mlr)
library(survival)
library(regplot)

Data <- read.csv("./Data4.csv", sep=",", dec=".")
Data <- Data[complete.cases(Data),]
ddist <- datadist(Data)
ddist$limits$hf[3] <- 30
ddist$limits$hf[5] <- 30
ddist$limits$hf[7] <- 30
options(datadist='ddist')
options(datadist=ddist)
model1 <- lm(ER2 ~ ic50anox + hf + dose, data = Data)
regplot(model1, plots = c("no plot", "no plot"), points=T, title= "ER2 model")

model2 <- lm(SGD2 ~ ic50anox + hf + dose, data = Data)
regplot(model2, plots = c("no plot", "no plot"), title= "SGD2 model")
