install.packages("Biodem")
library(Biodem)
library(reldist)
data("surnames") #An example from the package. It is a frequencies distribution table
hedkin <- hedrick(surnames) #Generates Hedrick (1971) kinship matrix. There are other methods (i.e. lasker, uri)
hedkin.dist <- as.dist(1-hedkin) #A very simple method. See also as.dist(), dist()
clust_hedkin <- hclust(hedkin.dist) #Hierarchical clustering of the distance matrix. There are various methods.
plot(clust_hedkin) #The tree

##Trying to get S of each population from "surnames" dataset

data_surnames <- as.data.frame(surnames)
desag_surnames <- data.frame(rep(data_surnames$Cognome,data_surnames$Freq),rep(data_surnames$Population,data_surnames$Freq)) #Disaggreagete frequencies
colnames(desag_surnames) <- c("Cognome","Population") # Put readable columns' names

fun_s <- function(x){length(unique(x$Cognome))} # a self-made function to compute "s"
ss <- by(desag_surnames, desag_surnames$Population, fun_s) #get "s" of each population
ns <- by(desag_surnames, desag_surnames$Population, nrow) #get "n" of each population

S_table <- cbind(sapply(ns,I),sapply(ss,I)) #This will make a table of n and s by population
colnames(S_table) <- c("n","s") # To change colnames

##get S
S_table <- as.data.frame(S_table)
S_table$S <- S_table$s / S_table$n
S_table

####### COMMONERS' DATA #######

##Surname's distance matrix


##get S from "commoners" data
fun_s <- function(x){length(unique(x$surname_father))} # a self-made function to compute "s"
ss <- by(commoners, commoners$community, fun_s) #get "s" of each population
ns <- by(commoners, commoners$community, nrow) #get "n" of each population
S_table <- cbind(sapply(ns,I),sapply(ss,I)) #This will make a table of n and s by population
colnames(S_table) <- c("n","s") # To change colnames
S_table <- as.data.frame(S_table)
S_table$S <- S_table$s / S_table$n
S_table


##get R from "commoners" data
fun_r <- function(x)sum(x$shares) # a self-made function to compute "r"
rr <- by(commoners, commoners$community, fun_r) #get "r" of each population
nr <- by(commoners, commoners$community, nrow) #get "n" of each population
R_table <- cbind(sapply(nr,I),sapply(rr,I)) #This will make a table of n and a by population
colnames(R_table) <- c("n","r") # To change colnames
R_table <- as.data.frame(R_table)
R_table$R <- R_table$r / R_table$n
R_table

##get G (gini) from "commoners" data
fun_g <- function(x)gini(x$shares) # a self-made function to compute "r"
gg <- by(commoners, commoners$community, fun_g) #get "r" of each population
ng <- by(commoners, commoners$community, nrow) #get "n" of each population
G_table <- cbind(sapply(ng,I),sapply(gg,I)) #This will make a table of n and a by population
colnames(G_table) <- c("n","G") # To change colnames
G_table <- as.data.frame(G_table)

##get A from "commoners" data
fun_a <- function(x)sum(x$sex == "M") # a self-made function to compute "a"
aa <- by(commoners, commoners$community, fun_a) #get "a" of each population
na <- by(commoners, commoners$community, nrow) #get "n" of each population
A_table <- cbind(sapply(na,I),sapply(aa,I)) #This will make a table of n and a by population
colnames(A_table) <- c("n","a") # To change colnames
A_table <- as.data.frame(A_table)
A_table$A <- A_table$a / A_table$n
A_table


##Compare surnames matrix with R matrix
##Matrix of R values
R_dist <- dist(R_table["R"])
#Get surnames matrix
surnames_freq <- table(commoners$surname_father, commoners$community)
hedkin <- hedrick(surnames_freq)
hedkin_dist <- as.dist(1-hedkin)
#Comparisson
library(ade4)
mantel.randtest(quasieuclid(hedkin_dist),R_dist, nrepet=9999)

##Compare surnames matrix with G matrix
##Matrix of G values
G_dist <- dist(G_table["G"])
#Get surnames matrix
surnames_freq <- table(commoners$surname_father, commoners$community)
hedkin <- hedrick(surnames_freq)
hedkin_dist <- as.dist(1-hedkin)
#Comparisson
library(ade4)
mantel.randtest(quasieuclid(hedkin_dist),G_dist, nrepet=9999)


##Define SEX
library(readr)
genero_nombres <- read_csv(../"genero_nombres.csv", 
                           col_names = FALSE)
View(genero_nombres)
colnames(genero_nombres) <- c("firstname1","sex_argent")
genero_nombres2 <- genero_nombres
colnames(genero_nombres2) <- c("firstname2","sex_argent2")
genero_nombres$firstname1 <- toupper(genero_nombres$firstname1)

commoners_sex <- join(commoners, genero_nombres, type="left", match="first")
commoners_sex <- join(commoners_sex, genero_nombres2, type="left", match="first")

commoners_sex$sex_argent3 <- coalesce(commoners_sex$sex_argent,commoners_sex$sex_argent2)
commoners_sex$sex_argent4 <- commoners_sex$sex
commoners_sex$sex_argent4[commoners_sex$sex_argent4=="A"] <- NA
commoners_sex$sex_argent5 <- coalesce(commoners_sex$sex_argent3,commoners_sex$sex_argent4)

#How to put a file with many populations in the format of the "surnames" data example?
library(readr)
example <- read_csv(example.csv)

