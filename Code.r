#####################################Task1#######################################################################

library(gRain)
library(Rgraphviz)
library(gRbase)
library(ggm)
library(gRim)
library(bnlearn)
library(igraph)
data(cad1, package = "gRbase")
cad1_mod<-cad1[,-c(1,8,9)]
cad.bn <- hc(cad1_mod)
net <- as(amat(cad.bn ), "graphNEL")
block <- c(3, 3, 4, 4, 4, 4,1, 1, 1, 3, 2) 
blM <- matrix(0, nrow = 11, ncol = 11)
rownames(blM) <- names(cad1_mod)
colnames(blM) <- names(cad1_mod)


for (b in 2:4){
	blM[block == b, block < b] <- 1
}

blackL <- data.frame(get.edgelist(as(blM, "igraph")))
names(blackL) <- c("from", "to")
cad.bn2 <- hc(cad1_mod, blacklist = blackL)
net.constr <- as(amat(cad.bn2), "graphNEL")
CPT<-extractCPT(cad1_mod, net.constr, smooth = 0)
CPT_grain<-grain(CPT)
compile_CPT_grain <- compile(CPT_grain)
propagate_CPT_grain <- propagate(compile_CPT_grain)
CPT.ev <- setFinding(propagate_CPT_grain, nodes = c("asia", "dysp"), states = c("yes", "yes"))
dSep(as(net.constr, "matrix"), "Smoker", "CAD", "Inherit")
dSep(as(net.constr, "matrix"), "AMI", "AngPec", "CAD")
dSep(as(net.constr, "matrix"), "QWavecode", "CAD", "Inherit")
dSep(as(net.constr, "matrix"), "QWave", "Heartfail", "CAD")
dSep(as(net.constr, "matrix"), "AMI", "Smoker", "CAD")
dSep(as(net.constr, "matrix"), "AMI", "Smoker", "QWavecode")
dSep(as(net.constr, "matrix"), "CAD", "Heartfail", "Smoker")
dSep(as(net.constr, "matrix"), "STchange", "Inherit", "CAD")
dSep(as(net.constr, "matrix"), "STcode", "Smoker", "CAD")
dSep(as(net.constr, "matrix"), "STcode", "Smoker", "CAD")
dSep(as(net.constr, "matrix"), "CAD", "Heartfail", "Hyperchol")


CPT<-extractCPT(cad1, net.constr, smooth = 0)
CPT_grain<-grain(CPT)
compile_CPT_grain <- compile(CPT_grain)
compile_CPT_grain <- propagate(compile_CPT_grain)
CAD_ev <- setFinding(compile_CPT_grain, nodes = c("Sex", "Hyperchol"), states = c("Female", "Yes"))
# probabilistic query, given evidence
with_ev <- querygrain(CAD_ev, nodes = c("Heartfail", "CAD"), type = "marginal")

# probabilistic query, given NO evidence
without_ev <- querygrain(compile_CPT_grain, nodes = c("Heartfail", "CAD"), type = "marginal")

simulate_data_with_evd<-simulate(CAD_ev, n=100)
block <- c(1, 3, 3, 4, 4, 4, 4, 1, 2, 1, 1, 1, 3, 2) 
blM <- matrix(0, nrow = 14, ncol = 14)
rownames(blM) <- names(simulate_data_with_evd)
colnames(blM) <- names(simulate_data_with_evd)


for (b in 2:4){
	blM[block == b, block < b] <- 1
}
blackL <- data.frame(get.edgelist(as(blM, "igraph")))
names(blackL) <- c("from", "to")
cad.bn2 <- hc(simulate_data_with_evd, blacklist = blackL)
net.constr <- as(amat(cad.bn2), "graphNEL")
CPT<-extractCPT(simulate_data_with_evd, net.constr, smooth = 0)
CPT_grain<-grain(CPT)
compile_CPT_grain <- compile(CPT_grain)
compile_CPT_grain <- propagate(compile_CPT_grain)
new_ev1 <- setFinding(compile_CPT_grain, nodes = c("Inherit"), states = c("Yes"))
with_ev <- querygrain(new_ev1, nodes = c("Smoker", "CAD"), type = "marginal")
without_ev <- querygrain(compile_CPT_grain, nodes = c("Smoker", "CAD"), type = "marginal")
new_ev2<-setFinding(compile_CPT_grain, nodes = c("Heartfail"), states = c("Yes"))
with_ev <- querygrain(new_ev2, nodes = c("Smoker", "CAD"), type = "marginal")
without_ev <- querygrain(compile_CPT_grain, nodes = c("Smoker", "CAD"), type = "marginal")
new_ev3<-setFinding(compile_CPT_grain, nodes = c("STchange","STcode"), states = c("Yes","Usable"))
with_ev <- querygrain(new_ev3, nodes = c("Smoker", "CAD"), type = "marginal")
without_ev <- querygrain(compile_CPT_grain, nodes = c("Smoker", "CAD"), type = "marginal")
new_ev4<-setFinding(compile_CPT_grain, nodes = c("QWave","QWavecode"), states = c("Yes","Usable"))
with_ev <- querygrain(new_ev4, nodes = c("Smoker", "CAD"), type = "marginal")
without_ev <- querygrain(compile_CPT_grain, nodes = c("Smoker", "CAD"), type = "marginal")
new_ev5<-setFinding(compile_CPT_grain, nodes = c("AMI"), states = c("Definite"))
with_ev <- querygrain(new_ev5, nodes = c("Smoker", "CAD"), type = "marginal")
without_ev <- querygrain(compile_CPT_grain, nodes = c("Smoker", "CAD"), type = "marginal")
new_ev6<-setFinding(compile_CPT_grain, nodes = c("AngPec"), states = c("Atypical"))
with_ev <- querygrain(new_ev6, nodes = c("Smoker", "CAD"), type = "marginal")
without_ev <- querygrain(compile_CPT_grain, nodes = c("Smoker", "CAD"), type = "marginal")

###################################Question No. 2 - ########################################################################
################Exploratory analysis of Titanic Data#########################################################################
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
###############Merging train & test sets for complete Exploratory Analysis#################################################
titanic_test$Survived <- NA
titanic <- rbind(titanic_train, titanic_test)
str(titanic)
summary(titanic)
####################Checking Missing Values########################################################################################
checkNA <- function(x){sum(is.na(x))/length(x)*100}
sapply(titanic,checkNA)
Title <-  gsub("^.*, (.*?)\\..*$", "\\1", titanic_train$Name)
titanic_train$Title <- as.factor(Title)
table(Title)
################Extracting Title from Name###########################################################################
Title <-  gsub("^.*, (.*?)\\..*$", "\\1", titanic_train$Name)
titanic_train$Title <- as.factor(Title)
table(Title)
#################Extracting Family Size from SibSp & Parch###########################################################
titanic_train$FamilyCount <-titanic_train$SibSp + titanic_train$Parch + 1
titanic_train$FamilySize[titanic_train$FamilyCount == 1] <- 'Single'
titanic_train$FamilySize[titanic_train$FamilyCount < 5 & titanic_train$FamilyCount >= 2] <- 'Small'
titanic_train$FamilySize[titanic_train$FamilyCount >= 5] <- 'Big'
titanic_train$FamilySize=as.factor(titanic_train$FamilySize)
table(titanic_train$FamilySize)
#######################Data Preprocessing##############################################################################
# 1.Changing names of few categorical variables for interpretability
titanic_train$Survived <- ifelse(titanic_train$Survived==1,"Yes","No")
titanic_train$Survived <- as.factor(titanic_train$Survived)

titanic_train$Embarked <- ifelse(titanic_train$Embarked=="S","Southampton",
                              ifelse(titanic_train$Embarked=="C","Cherbourg", "Queenstown"))
titanic_train$Embarked <- as.factor(titanic_train$Embarked)

# 2.Converting categorical variables from int to factor
# i) Pclass
titanic_train$Pclass <- as.factor(titanic_train$Pclass)

# ii) SibSp
titanic_train$SibSp <- as.factor(titanic_train$SibSp)

# iii) Parch
titanic_train$Parch <- as.factor(titanic_train$Parch)
############################################Missing Value Treatment####################################################

#1. Age: Replacing NA values in Age with mean
#titanic[is.na(titanic$Age),6] <- mean(titanic$Age)
titanic$Age[is.na(titanic$Age)] <- round(mean(titanic$Age, na.rm = TRUE))
titanic_train$Age[is.na(titanic_train$Age)] <- round(mean(titanic_train$Age, na.rm = TRUE))

#2. Embarked: Replacing Empty Embarked with most common value 'S'
titanic_train$Embarked <- replace(titanic_train$Embarked, which(titanic_train$Embarked==""), 'S')

#3. Cabin: Not replacing with anything as Cabin values are unique
########################Percentage of survival###########################################################
temp<-subset(titanic_train,titanic_train$Survived=="Yes")
(nrow(temp)/nrow(titanic_train))*100
#############Avergage Age Distribution ###############################################################
summary(titanic$Age)
x11()
d <- density(titanic$Age)
plot(d,main="Passenger Age Distribution",xlab="Age",ylab="Frequency",col="blue")

####################Proportion of survivors by gender###################################################
ggplot(titanic_train, aes(x=Sex,fill=Survived))+ geom_bar(position = "dodge")
 + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers")

##########################Proportion of survivors by Parch##############################################

ggplot(titanic_train, aes(x=Parch,fill=Survived))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers") + xlab("Number of Parents/Children")

##################Testing the characteristics/demographics are more likely in passengers that perished####################
##############Age distribution of Survivors & Non-Survivors

x11()
ggplot(titanic_train) + geom_freqpoly(mapping = aes(x = Age, color = Survived), binwidth = 2.5) +
ylab("Frequency")
############################distribution of Passenger Fare for Survivors & Non-Survivors##################################
ggplot(titanic_train) + geom_freqpoly(mapping = aes(x = Fare, color = Survived), binwidth = 10)
#################Passenger Class of most Non-Survivors#############################################################

ggplot(titanic_train, aes(x=Pclass,fill=Survived))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers") + xlab("Passenger Class")

########################proportion of survivors by place of Embarkment############################################
ggplot(titanic_train, aes(x=Embarked,fill=Survived))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers")
#################Characteristics of passenger Class 3 which mostly didn't survive#############################################################
ggplot(titanic_train, aes(x=SibSp,fill=Pclass))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers") + xlab("Number of Siblings")


ggplot(titanic_train, aes(x=SibSp,fill=Survived))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers")+xlab("Number of Siblings/Spouse")

ggplot(titanic_train, aes(x=Parch,fill=Pclass))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers") + xlab("Number of Parents/Children")

#####################relation between Family Size & Survival###############
ggplot(titanic_train, aes(x=FamilySize,fill=Survived))+ geom_bar(position = "dodge") + geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9),vjust=-0.2) +
ylab("Number of Passengers") + xlab("Family Size")

###########################Question 4###################################################################################################################
st <- data.frame(state.x77)
library(gRbase) 
library(gRim) 
library(gRain) 
library(glasso) 
library(graph) 
library(corrplot) 
library(ggplot2)

summary(st)
str(st)
sum(is.na(st))
library(ellipse) 
library(corrplot)
corMatrix <- cor(as.matrix(st))
col <- colorRampPalette(c("red", "yellow", "blue"))
corrplot.mixed(corMatrix, order = "AOE", lower = "number", lower.col = "black", 
               number.cex = .8, upper = "ellipse",  upper.col = col(10), 
               diag = "u", tl.pos = "lt", tl.col = "black")
M<-cor(st)
corrplot(M)
fit.pca <- prcomp(scale(st))
xlim_1 <- min(fit.pca$x[,1])-1
xlim_2 <- max(fit.pca$x[,1])+1
ylim_1 <- min(fit.pca$x[,2])-1
ylim_2 <- max(fit.pca$x[,2])+1
x11()
biplot(fit.pca, choices = c(1,2), scale = 0, xlim = c(xlim_1, xlim_2), ylim = c(ylim_1, ylim_2))
rownames(st)<-NULL
st_mod <- st[-2,]
S.body <- cov.wt(st_mod, method = "ML")
PC.body <- cov2pcor(S.body$cov)
diag(PC.body) <- 0
x11()
heatmap(PC.body)
library(ggpubr)
census_data_mod<-st_mod
x11()
ggplot(st_mod,aes(Population)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Population),
                            sd = sd(st_mod$Population)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Population before log transformation")

census_data_mod$Population<-log(census_data_mod$Population)
x11()
ggplot(census_data_mod,aes(Population)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Population after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Population),
                            sd = sd(census_data_mod$Population)),
                col = "red",
                size = 3)
x11()
ggplot(st_mod,aes(Income)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Income),
                            sd = sd(st_mod$Income)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Income before log transformation")
census_data_mod$Income<-log(census_data_mod$Income)
x11()
ggplot(census_data_mod,aes(Income)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Income after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Income),
                            sd = sd(census_data_mod$Income)),
                col = "red",
                size = 3)
x11()
ggplot(st_mod,aes(Illiteracy)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Illiteracy),
                            sd = sd(st_mod$Illiteracy)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Illiteracy before log transformation")
census_data_mod$Illiteracy<-log(census_data_mod$Illiteracy)
x11()
ggplot(census_data_mod,aes(Illiteracy)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Illiteracy after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Illiteracy),
                            sd = sd(census_data_mod$Illiteracy)),
                col = "red",
                size = 3)

x11()
ggplot(st_mod,aes(Life.Exp)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Life.Exp),
                            sd = sd(st_mod$Life.Exp)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Life.Exp before log transformation")
census_data_mod$Life.Exp<-log(census_data_mod$Life.Exp)
x11()
ggplot(census_data_mod,aes(Life.Exp)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Life.Exp after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Life.Exp),
                            sd = sd(census_data_mod$Life.Exp)),
                col = "red",
                size = 3)

x11()
ggplot(st_mod,aes(Murder)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Murder),
                            sd = sd(st_mod$Murder)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Murder before log transformation")
census_data_mod$Murder<-log(census_data_mod$Murder)
x11()
ggplot(census_data_mod,aes(Murder)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Murder after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Murder),
                            sd = sd(census_data_mod$Murder)),
                col = "red",
                size = 3)

x11()
ggplot(st_mod,aes(HS.Grad)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$HS.Grad),
                            sd = sd(st_mod$HS.Grad)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of HS Grad Students before log transformation")
census_data_mod$HS.Grad<-log(census_data_mod$HS.Grad)
x11()
ggplot(census_data_mod,aes(HS.Grad)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of HS Grad Students after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$HS.Grad),
                            sd = sd(census_data_mod$HS.Grad)),
                col = "red",
                size = 3)
Frost
x11()
ggplot(st_mod,aes(Frost)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Frost),
                            sd = sd(st_mod$Frost)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Frost before log transformation")
census_data_mod$Frost<-log(census_data_mod$Frost)
x11()
ggplot(census_data_mod,aes(Frost)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Frost after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Frost),
                            sd = sd(census_data_mod$Frost)),
                col = "red",
                size = 3)

x11()
ggplot(st_mod,aes(Area)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+stat_function(fun = dnorm,
                args = list(mean = mean(st_mod$Area),
                            sd = sd(st_mod$Area)),
                col = "red",
                size = 3)+ggtitle("Density Distribution of Area before log transformation")
census_data_mod$Area<-log(census_data_mod$Area)
x11()
ggplot(census_data_mod,aes(Area)) + geom_density( fill="#77bd89",color="#1f6e34",alpha=0.8)+ggtitle("Density Distribution of Area after log transformation")+stat_function(fun = dnorm,
                args = list(mean = mean(census_data_mod$Area),
                            sd = sd(census_data_mod$Area)),
                col = "red",
                size = 3)
st_mod1 <- stack(st_mod)
census_data_mod1<-stack(census_data_mod)
ggplot(census_data_mod1, aes(x=values, fill=ind)) + geom_density(alpha=0.5)
##############################Designing Graphical Lasso Model####
# Look at partial correlation
library("glasso")
S.body <- cov.wt(st_mod, method = "ML")
PC.body <- cov2pcor(S.body$cov)
diag(PC.body) <- 0
S <- S.body$cov
m0.lasso <- glasso(S, rho = 5) # fit the model
names(m0.lasso)
my.edges <- m0.lasso$wi != 0 # grab the non-zero edges
diag(my.edges) <- 0 # kill the diagonal
g.lasso <- as(my.edges, "graphNEL") # converting to a graphNEL object for plotting
nodes(g.lasso) <- names(st_mod)
x11()
plot(g.lasso,main="Graphical Lasso Model on Census Data with Lambda/Penalty=5" )

library(geneplotter)
graphics.off()
my_rhos <- c(2,5,10,15,25,50,75,100,125,150,175,200)
m0.lasso <- glassopath(S, rho = my_rhos)
for (i in 1:length(my_rhos)){
    my.edges <- m0.lasso$wi[ , , i] != 0 # grab the non-zero edges of the ith object
    diag(my.edges) <- 0 # kill the diagonal
    g.lasso <- as(my.edges, "graphNEL") # converting to a graphNEL object for plotting
    nodes(g.lasso) <- names(st_mod)

    x11()
    plot(g.lasso,main=paste0("Graphical Lasso Model on Census Data with Lambda/Penalty=",my_rhos[i]))
    savepdf(paste("myplot", i, sep = "_"))
 }
library(kohonen)
census_data_mod<-st_mod
census_data_mod$Population<-scale(census_data_mod$Population)
census_data_mod$Income<-scale(census_data_mod$Income)
census_data_mod$Illiteracy<-scale(census_data_mod$Illiteracy)
census_data_mod$Life.Exp<-scale(census_data_mod$Life.Exp)
census_data_mod$Murder<-scale(census_data_mod$Murder)
census_data_mod$HS.Grad<-scale(census_data_mod$HS.Grad)
census_data_mod$Frost<-scale(census_data_mod$Frost)
census_data_mod$Area<-scale(census_data_mod$Area)
census_data_mod<-as.matrix(census_data_mod)
som_grid <- somgrid(xdim = 5, ydim = 5, topo = "hexagonal")
census_data_mod.som <- som(census_data_mod, grid = som_grid, rlen = 10000)

codes <- census_data_mod.som$codes[[1]]

plot(census_data_mod.som,main="Self Organizing Maps on Census data")
census_data_mod.som$unit.classif
plot(census_data_mod.som,type="changes",main="Self Organizing Maps on Census data")