
######################################################################################
######################################################################################
######### SARS-CoV-2 Coefficient of Exponential Growth Alteration (CEGA) #############
######################################################################################
######################################################################################

###################################################################################################################
## 02.11.2021: DRAFT VERSION. This script WILL be annotated and cleaned further in the days to come ###############
###################################################################################################################

library(ape)
library(phangorn)
library(reshape2)
library(doParallel) # for parallelization
library(foreach) # for parallelization
library(ggplot2)
library(ggExtra)
library(bigmemory)
library(parallel)
library(phytools)

thread_number <- 4 # thread_number <- "auto" # set it to auto or to a number
deduplicated_dataset <- 0 # bolean, set to 1 (true) if dataset contains deduplicated sequences and if you have a table of re-duplication table
simulated_t_a_dataset <- 0 # set to 1 to run on simulated datasets, debug purposes # not compatible with deduplicated_dataset = 1
children_node_rule <- 1 # set to 1 to activate the filter that makes not study the nodes that have embedded homoplasic nodes 
min_offspring <- 20 # node annotated as carying an homoplasy and having less tips than "min_offspring" will be un-annotated
# they will not have their roho score computed and they will not count as embedded homoplasies neither.
tips_rule  <- 10

arbre_filtered_file <- list.files(path = ".", pattern = "^annotatedNewickTree_.*.tree")
arbre_filtered <- read.newick(arbre_filtered_file)

input <- bigmemory::read.big.matrix("input_matrix_numeric", has.row.names=T,header=T,sep="\t")

# read the node labels and melt them to have a usefull format
annotation_nodes_useful_format <- melt(strsplit(arbre_filtered$node.label,"-"))
annotation_nodes_useful_format[,1] <- as.numeric(as.character(annotation_nodes_useful_format[,1])) # added 12 may

# get the architecture of the tree 
# this function gives, for each node and tip of the tree, the name of its descendants
architecture_arbre <- Descendants(arbre_filtered,type="all")

# FIRST CONDITION
# we don't want to study nodes that were annotated as carying an homoplasy AND have less than X descendants TIPS
# architecture_arbre[[312]] from 1 to 312 are tips
nb_offspring <- unlist(lapply(architecture_arbre, length))
object <- table(nb_offspring)
names(nb_offspring) <- seq(from=1,to=length(nb_offspring),by=1)
min_offspring_tips_and_nodes <- (min_offspring * 2 )-2
enough_offsprings <- nb_offspring[nb_offspring>=min_offspring_tips_and_nodes]

# STILL FIRST CONDITION
# i delete nodes that don't have enough offsprings from list of nodes to study homoplasies from
annotation_nodes_useful_format <- annotation_nodes_useful_format[(annotation_nodes_useful_format[,2]+length(arbre_filtered$tip.label))%in%names(enough_offsprings),]

# store all homoplasic SNP positions in a vector
homoplasies <- as.numeric(as.character(unique(annotation_nodes_useful_format[,1])))

# here i can eliminate homoplasies that do not interest us
coordinates <- read.csv("correspondance2",header=T, sep="\t",dec=",", stringsAsFactors=F,check.names=FALSE)
if(! is.na(coordinates[coordinates[,1]==11083,3])){
homoplasies <-  homoplasies[!homoplasies==coordinates[coordinates[,1]==11083,3]]
}
if(! is.na(coordinates[coordinates[,1]==21575,3])){
homoplasies <-  homoplasies[!homoplasies==coordinates[coordinates[,1]==21575,3]]
}
if(! is.na(coordinates[coordinates[,1]==21987,3])){
  homoplasies <-  homoplasies[!homoplasies==coordinates[coordinates[,1]==21987,3]]
}

nodes_to_keep <- annotation_nodes_useful_format[annotation_nodes_useful_format[,1]%in%homoplasies,]
nb_tips <- length(arbre_filtered$tip.label)
nodes_to_keep[,2] <- nodes_to_keep[,2] + nb_tips

# parallel solution
arbre_filtered_tips <- arbre_filtered$t # do not delete that i need it for distance timespan calculation later
loop3 <- function(x,input){
  tips <- architecture_arbre[[x[2]]][architecture_arbre[[x[2]]]<=(nb_tips)]
  alleles <- as.character(unname(input[as.character(x[1]),arbre_filtered_tips[tips]]))
  names(alleles) <- arbre_filtered_tips[tips]
  NOT_homoplasy_count <- sum(alleles=="0")
  homoplasy_count <- sum(alleles=="1")
  
  two_lineages <- Children(arbre_filtered,x[2])
  tips_lineage_1 <- arbre_filtered_tips[architecture_arbre[[two_lineages[1]]][architecture_arbre[[two_lineages[1]]]<=(nb_tips)]]

  tips_lineage_1_alleles <- alleles[tips_lineage_1]
  lin_1_all_0 <- sum(tips_lineage_1_alleles==0)
  lin_1_all_1 <- sum(tips_lineage_1_alleles==1)
  
  tips_lineage_2_alleles <- alleles[!names(alleles)%in%tips_lineage_1]
  lin_2_all_0 <- sum(tips_lineage_2_alleles==0)
  lin_2_all_1 <- sum(tips_lineage_2_alleles==1)
  
  perfect <- length(unique(tips_lineage_1_alleles)) == 1 & length(unique(tips_lineage_2_alleles)) == 1
  
  NOT_homoplasy_count <- sum(alleles=="0")
  homoplasy_count <- sum(alleles=="1")
  return(c(NOT_homoplasy_count,homoplasy_count,perfect,lin_1_all_0,lin_1_all_1,lin_2_all_0,lin_2_all_1))
}

cl <- makeCluster(4)
datadesc <- describe(input)
clusterExport(cl, c("architecture_arbre","arbre_filtered_tips","nb_tips","loop3","datadesc","arbre_filtered"))
clusterEvalQ(cl, {
  library(bigmemory)
  library(phangorn)
})

parApply_result <-  t(parApply(cl, nodes_to_keep,1, function(x){
  input <- attach.big.matrix(datadesc)
  loop3(x,input) }))
stopCluster(cl)

# brouillon ON

# the goal now is to apply a 90% threshold on the pureness of each clade
df_90 <- df[((df$lin_1_all_0 > 0.9 * (df$lin_1_all_0 + df$lin_1_all_1)) | (df$lin_1_all_1 > 0.9 * (df$lin_1_all_0 + df$lin_1_all_1)) ) &  ((df$lin_2_all_0 > 0.9 * (df$lin_2_all_0 + df$lin_2_all_1)) | (df$lin_2_all_1 > 0.9 * (df$lin_2_all_0 + df$lin_2_all_1)) ),]

df_90$homoplasy <- as.factor(df_90$homoplasy)
df_90 <- df_90[df_90[,15]=="5_10",]
df_90 <- df_90[is.finite(df_90$CEGA2),]

df_tot <- raw_out_table2_complete_filtered
df_tot$homoplasy <- as.factor(df_tot$homoplasy)
df_tot <- df_tot[df_tot[,15]=="5_10",]
df_tot <- df_tot[is.finite(df_tot$CEGA2),]
# brouillon OFF

new_table <- cbind(nodes_to_keep,parApply_result)
new_table <- new_table[new_table[,4] >= tips_rule & new_table[,3] >= tips_rule,]
colnames(new_table) <- c("homoplasy","node","NOT_homoplasy_count","homoplasy_count","perfect","lin_1_all_0","lin_1_all_1","lin_2_all_0","lin_2_all_1")
# keep only 90% perfect nodes in new_table :
nrow(new_table)
new_table <- new_table[((new_table$lin_1_all_0 > 0.9 * (new_table$lin_1_all_0 + new_table$lin_1_all_1)) | (new_table$lin_1_all_1 > 0.9 * (new_table$lin_1_all_0 + new_table$lin_1_all_1)) ) &  ((new_table$lin_2_all_0 > 0.9 * (new_table$lin_2_all_0 + new_table$lin_2_all_1)) | (new_table$lin_2_all_1 > 0.9 * (new_table$lin_2_all_0 + new_table$lin_2_all_1)) ),]
nrow(new_table)

# for non-paralelized version, create object that will receive result
# i <- 0 ; raw_out_table <- matrix( nrow = 10000000, ncol = 4)
# debug_raw_out_table <- matrix( nrow = 10000000, ncol = 5)
# 
# debug_raw_out_table <- t(as.matrix(apply(new_table,1,function(line){
#    homoplasy <- line[1]
#    node <- line[2]
#    child_nodes <- architecture_arbre[[node]][architecture_arbre[[node]]>(length(arbre_filtered$tip.label))]
#    # si ya des child node presents dans new_table POUR LA MEME homoplasie on drop
#    count_children_node_with_homoplasy_with_rule <- sum(new_table[new_table[,2]%in%child_nodes,1]==homoplasy) 
# 
#    tip_rule_granted <- line[5] >= 10 & line[4] >= 10
#    
#    NOT_homoplasy_count <- line[5]
#    homoplasy_count <- line[4]
#    
#    return(c(homoplasy,homoplasy_count,NOT_homoplasy_count,node,count_children_node_with_homoplasy_with_rule,tip_rule_granted))
#    
#  })))
#  
 raw_out_table <- t(as.matrix(apply(new_table,1,function(line){
  homoplasy <- line[1]
  node <- line[2]
  child_nodes <- architecture_arbre[[node]][architecture_arbre[[node]]>(length(arbre_filtered$tip.label))]
    # si ya des child node presents dans new_table POUR LA MEME homoplasie on drop
    count_children_node_with_homoplasy <- sum(new_table[new_table[,2]%in%child_nodes,1]==homoplasy)
    
    if( count_children_node_with_homoplasy == 0 ){
      
      NOT_homoplasy_count <- line[3]
      homoplasy_count <- line[4]
      
      if(homoplasy_count >= tips_rule & NOT_homoplasy_count >= tips_rule){
 # i <- i+1 ; raw_out_table[i,] <- c(homoplasy,homoplasy_count,NOT_homoplasy_count,node)
        return(c(homoplasy,homoplasy_count,NOT_homoplasy_count,node))
      }else{return(c(6666666,0,0,0))}
    }else{return(c(0,0,0,0))} # fin de condition sur les children nodes
    
})))
raw_out_table

# it is normal that raw_out_table contains empty lines, they were created in order not to add lines to an existing matrix because this is not memory efficient in R. I will delete them now :
nrow(raw_out_table)
raw_out_table <- raw_out_table[raw_out_table[,1]!=0,]
nrow(raw_out_table)

#   raw_out_table[raw_out_table[,1]==1574,] # delme 11.04.2021
#   raw_out_table[raw_out_table[,4]==884881,] # delme 11.04.2021
# setwd("/media/damien/Elements/2021.03.30/04.homoplasyfinder") # delme 11.04.2021

#########################################################################
#########################  Raw data saving  #############################
#########################################################################
save.image(file = paste("Workspace_02_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))
# load(paste("Workspace_02_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))
colnames(raw_out_table) <- c("position","homoplasy_count","NOT_homoplasy_count","node")
class(raw_out_table) <- "numeric"

correspondance_file <- list.files(path = "../02.VCFs/", pattern = "^.*filtered.vcf$")
correspondance <- read.csv(paste("../02.VCFs/",correspondance_file,sep=""),header=F, sep="\t",dec=",", stringsAsFactors=F,check.names=FALSE)

correspondance$mut_type <- correspondance[,3]
correspondance$mut_type[correspondance$mut_type=="X"] <- "deletion"
correspondance$mut_type[correspondance$mut_type!="deletion"] <- "SNP"

coordinates$type <- correspondance[match(coordinates[,2],correspondance[,1]),5]
coordinates$type[grep("combo",coordinates$Gappy_aln)] <- "combo"
coordinates_restricted <- coordinates[!is.na(coordinates$type),]
correspondance$fake_position <- seq(from=1,to=nrow(correspondance),by=1)
 
# check if first line went okay 
raw_out_table[1,]
raw_out_table[,1] <- coordinates_restricted[match(raw_out_table[,1],coordinates_restricted[,3]),1]
class(raw_out_table) <- "numeric"
raw_out_table[1,]

#   raw_out_table[raw_out_table[,1]==28281,] # delme 11.04.2021
#   raw_out_table[raw_out_table[,4]==884881,] # delme 11.04.2021

coordinates_restricted$alt_count <- correspondance[match(as.numeric(coordinates_restricted[,2]),correspondance[,1]),4]

# all homoplasy positions studied :
raw_out_table[,1]
# all homoplasy combos studied :
raw_out_table[raw_out_table[,1]>30000,]

# colnames(debug_raw_out_table) <- c("position","line","node","homoplasy_count","NOT_homoplasy_count") # for debug
# debug_raw_out_table <- as.data.frame(debug_raw_out_table[!is.na(debug_raw_out_table[,1]),]) # for debug
raw_out_table <- as.data.frame(raw_out_table[!is.na(raw_out_table[,1]),])

# for debug
# write.table(debug_raw_out_table, file = paste("RAW_DEBUG_data_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")
# debug_raw_out_table$node[debug_raw_out_table$position==3037]%in%debug_raw_out_table$node[debug_raw_out_table$position==14408]

write.table(raw_out_table, file = paste("RAW_01_data_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

# calculate distances to root for all branches
library(castor)
all_distances <- get_all_distances_to_root(arbre_filtered, as_edge_count=FALSE)
# A numeric vector of size Ntips+Nnodes, with the i-th element being the distance (cumulative branch length) of the i-th tip or node to the root. Tips are indexed 1,..,Ntips and nodes are indexed (Ntips+1),..,(Ntips+Nnodes). 
# for each line of raw_out_table i need to calculate the mean of all the distances from node to tip
# 
# disssst <- t(apply(raw_out_table,1,function(x){
#   node_to_root <- all_distances[as.numeric(x[4])]
#   tips <- unlist(Descendants(arbre_filtered,as.numeric(x[4]),type="tips"))
#   # tips_to_root <- unlist(lapply(tips,function(x){ all_distances[x] }))
#   tips_to_root <- all_distances[tips]
#   tips_to_node <- tips_to_root - node_to_root
#   return(c(mean(tips_to_node),median(tips_to_node)))
# }))
# 
# # test parallelization du calcul de distance
# # optimization to check : https://www.r-bloggers.com/2016/01/strategies-to-speedup-r-code/
# cl <- makeCluster(8)
# clusterExport(cl, c("all_distances","arbre_filtered"))
# clusterEvalQ(cl, {
#   library(phangorn)
# })
# disssst <- t(parApply(cl,raw_out_table,1,function(x){
#   node_to_root <- all_distances[as.numeric(x[4])]
#   # tips <- unlist(Descendants(arbre_filtered,as.numeric(x[4]),type="tips"))
#   # tips_to_root <- all_distances[tips]
#   tips_to_node <- all_distances[unlist(Descendants(arbre_filtered,as.numeric(x[4]),type="tips"))] - node_to_root
#   return(c(mean(tips_to_node),median(tips_to_node)))
# }))
# stopCluster(cl)
# # end test parallelization du calcul de distance
# 
# # raw_out_table_fake <- raw_out_table[1:20,]
# # 
# # 
# # 
# # library(microbenchmark)
# # microbenchmark(result<-f1(raw_out_table_fake),
# #                result<-f3(raw_out_table_fake),
# #                times=1)
# # 
# # raw_out_table_fake
# # 
# # x <- raw_out_table[1,]

raw_out_table_with_dist <- raw_out_table
# raw_out_table_with_dist$mean <- disssst[,1]
# raw_out_table_with_dist$median <- disssst[,2]

raw_out_table_with_dist$mean <- NA
raw_out_table_with_dist$median <- NA

write.table(raw_out_table_with_dist, file = paste("RAW_02_distance_data_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

save.image(file = paste("Workspace_03_afterdistance_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))
# load(paste("Workspace_03_afterdistance_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))

#########################################################################
######################### Exploratory plots #############################
#########################################################################

out_table_long <- melt(raw_out_table[,-4], id=c("position"), measured=c("homoplasy_count","NOT_homoplasy_count"))
effectif <- nrow(raw_out_table)

# Density plot of values of homoplasy descendant counts and non-homoplasy descendant counts
jpeg(filename = paste("fig_1.non_filtered_count_histogram_min_offspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_long, mapping = aes (fill = variable, x = log(value))) + geom_density (alpha = .5) + ggtitle(paste("Homoplasy allele count vs non-homoplasy allele count\n non filtered dataset, n= ",effectif,sep=""))
dev.off()

out_table_ratios <- as.data.frame(cbind("position"=raw_out_table[,1],"homo/not_homo"=raw_out_table[,2]/raw_out_table[,3],"not_homo/homo"=raw_out_table[,3]/raw_out_table[,2]))
out_table_ratios_long <- melt(out_table_ratios, id=c("position"), measured=c("homo/not_homo","not_homo/homo"))

effectif <- sum(out_table_ratios_long[,2]=="homo/not_homo")
nb_unique_nodes <- length(unique(out_table_ratios_long[,1]))

# boxplot of all log10(RoHO) values along the genome. Unfiltered.
jpeg(filename = paste("fig_2.homo_vs_nothomo_ratio_along_genome_min_offspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo",],
       aes(x=position, y=log10(value), group=position))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(geom="point", fun.y=median)+
  geom_hline(yintercept = median(log10(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo","value"])), linetype='dashed')+
  ggtitle(paste("log of the ratio Homoplasy allele count / non-homoplasy allele count\nnb of studied homoplasy positions = ",nb_unique_nodes,", nb of studied paired count values= ",effectif,"\ndashed line = mediane; These are NON-FILTERED based on nreplicate",sep=""))
dev.off()

# Density plot of values of log10(RoHO)
jpeg(filename = paste("fig_3.non_filtered_ratio_histogram_min_offspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo",], mapping = aes (x = log(value))) + geom_density (alpha = .5) + ggtitle(paste("log of the ratio Homoplasy allele count / non-homoplasy allele count\n non filtered dataset, n= ",effectif,sep=""))
dev.off()

#####################################################################################
############ same plots with filtering of the data : nb of replicates ###############
############ eg. 5 nodes in the phylogeny detected as origin of a new homoplasy #####
#####################################################################################
nb_rep <- 5

out_table_restricted <- raw_out_table_with_dist[raw_out_table_with_dist[,1]%in%names(table(raw_out_table_with_dist[,1])[table(raw_out_table_with_dist[,1])>=nb_rep]),]

# implementation of the correction by the timespan of the tips rather than the branch lengths
metadata_file <- read.csv("../00.inputs/metadata_complete.tsv", row.names=NULL,sep="\t",dec=",", header=T, stringsAsFactors=F,check.names=FALSE)
colnames(metadata_file)

# vectorization to save time
nodes_to_evaluate <- unique(out_table_restricted[,4])
# nodes_to_evaluate2 <- nodes_to_evaluate[1:10] # for tests
subset_metadata_file <- metadata_file[metadata_file[,3]%in%arbre_filtered_tips,]

colnames(subset_metadata_file)
# subset_metadata_file <- subset_metadata_file[,-4]

library(HDInterval)
# test parallelization du calcul de distance
# optimization to check : https://www.r-bloggers.com/2016/01/strategies-to-speedup-r-code/
cl <- makeCluster(4)
clusterExport(cl, c("architecture_arbre","arbre_filtered_tips","subset_metadata_file"))
clusterEvalQ(cl, {
  # library(phangorn)
  library(HDInterval)
})
nb_generations_tmp <- t(parSapply(cl,nodes_to_evaluate,function(x){
  tips <- architecture_arbre[[x]][architecture_arbre[[x]]<=(length(arbre_filtered_tips))]
   homoplasy_dates <- as.Date(subset_metadata_file[match(arbre_filtered_tips[tips],subset_metadata_file[,3]),5],"%Y-%m-%d")
  homoplasy_dates <- homoplasy_dates[!is.na(homoplasy_dates)]
  #hdi <- hdi(as.numeric(homoplasy_dates),credMass = 0.95)
  # homoplasy_interv <- as.numeric(difftime(as.Date(as.numeric(hdi[2]), origin = "1970-01-01"), as.Date(as.numeric(hdi[1]), origin = "1970-01-01"), "%Y-%m-%d"), units = "days")
  homoplasy_interv_old <- as.numeric(difftime(max(homoplasy_dates), min(homoplasy_dates), "%Y-%m-%d"), units = "days")
  homoplasy_interv_gen <- homoplasy_interv_old/5
  if(homoplasy_interv_gen > 1 && !is.na(homoplasy_interv_gen)){
    # return(c(x,homoplasy_interv_old,homoplasy_interv_gen))
    return(c(x,homoplasy_interv_old,homoplasy_interv_gen,max(homoplasy_dates), min(homoplasy_dates)))
  }else{
    # return(c(x,0,0))
    return(c(x,0,0,max(homoplasy_dates), min(homoplasy_dates)))
  }
}))
stopCluster(cl)

# rm(subset_metadata_file)

nb_generations_tmp_backup <- nb_generations_tmp

# nb_generations_tmp comprises two versions of the same timespan parameter: uncorrected (2nd col) and corrected (3nd col), 1st col being node
nb_generations <- nb_generations_tmp[,c(1,3)] # for UNcorrected
# nb_generations <- nb_generations_tmp[,c(1,2)] # for corrected

plot(nb_generations_tmp[,3]~nb_generations_tmp[,5],ylab="generations",xlab="oldest_of_the_clade")

colnames(nb_generations_tmp) <- c("pos","old","new","max","min")

nb_generations_tmp <- as.data.frame(nb_generations_tmp)
class(nb_generations_tmp)

mean(nb_generations_tmp[,2])


p <- ggplot(nb_generations_tmp, aes(x=new, y=old)) + geom_point(colour = "black",alpha=0.2) +
  geom_smooth(method='lm', formula= y~x)
 ggMarginal(p, type="histogram")

homoplasy_dates <- as.Date(subset_metadata_file[match(arbre_filtered_tips,subset_metadata_file[,3]),4],"%Y-%m-%d")
homoplasy_dates <- homoplasy_dates[!is.na(homoplasy_dates)]
min(homoplasy_dates)

as.numeric(difftime(max(homoplasy_dates), min(homoplasy_dates), "%Y-%m-%d"), units = "days")

## TOKEEP :
# test si ma methode darchitecture arbre marche aussi bien que descendants pour recuperer la liste integrale des tips descendants d'un node
# architecture_arbre[[514492]][architecture_arbre[[514492]]<=(length(arbre_filtered$tip.label))]
# Descendants(arbre_filtered,514492,type="tips")
# c'est nickel !
## TOKEEP END

# TODO i need to remove the total old dist calculation and therefore the two associated columns

save.image(file = paste("Workspace_04_timespan_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))
# load(paste("Workspace_04_timespan_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))


out_table_restricted <- cbind(out_table_restricted,nb_generations[match(out_table_restricted[,4],nb_generations[,1]),2])
colnames(out_table_restricted)[length(colnames(out_table_restricted))] <- "nb_generations"

write.table(out_table_restricted, file = paste("RAW_03_timespan_data_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")


######################################################################################
############################# STATISTICAL ANALYSIS ###################################
######################################################################################
# Multiple test are performed :
# paired t test #
# Shapiro test to check wether the differences follow a normal distrib
# Sign test that do not need the diff to follow a norm. distrib.
# it only allows to say if the distrib is equal or not to zero.
# The sign test is indifferent to magnitude of difference. It only cares if it is above or below 0. If you wish to incorporate the magnitude, you should look at the Signed Rank test (wilcox.test).



library(BSDA)

multiple_t_test <- sapply(sort(unique(out_table_restricted[,1])), function(x) {
  # head(out_table_restricted)
  # max(out_table_restricted[,7])
  # x <- 23403
  # sum(out_table_restricted[,1]==23403)
  # out_table_restricted[out_table_restricted[,1]==23403,]
  homoplasy_count <- out_table_restricted[out_table_restricted[,1] == x,2]
  NOT_homoplasy_count <- out_table_restricted[out_table_restricted[,1] == x,3]
  
  res <- t.test(homoplasy_count,NOT_homoplasy_count, paired=TRUE)
  old_roho <- homoplasy_count/NOT_homoplasy_count
  
  time <- out_table_restricted[out_table_restricted[,1] == x,7] 
  CEGA2 <- (8 / time) * ((homoplasy_count / (homoplasy_count + NOT_homoplasy_count))-1/2)

  #06.06.21
  # CEGA2 <- (8 / homoplasy_interv_gen) * ((homoplasy_count/ (homoplasy_count + NOT_homoplasy_count))-1/2)
  # return(c(homoplasy,nb_offspring,homoplasy_count,NOT_homoplasy_count,node,count_children_node_with_homoplasy_with_rule,count_children_node_with_homoplasy_without_rule,tip_rule_granted,CEGA,CEGA2,homoplasy_interv_gen,roho,tips_rule,min_offspring))
  
  CEGA2 <- CEGA2[!time<=1] # fix 17.03.2021 forbid time <= 1
  mean_roho <- mean(old_roho)
  mean_CEGA <- mean(CEGA2)
  median_CEGA <- median(CEGA2)
  sd_CEGA <- sd(CEGA2)
  # if p-value < 0.05, we reject equality hypothesis, so two distributions are not equal.
  # res$estimate is the mean value of differences between members of the pair
  # if negative, it means that x < y # in our case if negative, it means that homo < not_homo
  shapiro <- shapiro.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3])$p.value
  # if p-value < 0.05, distribution do NOT follow a normal distribution
  sign_test <- SIGN.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], y = NULL, md = 0, alternative = "two.sided", conf.level = 0.95)$p.value
  sign_test_bonferroni <- p.adjust(sign_test, method = "bonferroni", n=length(unique(out_table_restricted[,1])))
  wilcox_test <- wilcox.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], paired = F, mu=0, alternative = "two.sided")$p.value
  nb_replicates <- length(out_table_restricted[out_table_restricted[,1] == x,2])
  # null hypothesis is "median = 0". We do NOT reject NULL hypothesis if p-value > 0.05
  # so no significantly different from zero if p-value > 0.05
  return(c(position=x,mean_CEGA=mean_CEGA,median_CEGA=median_CEGA,sd_CEGA=sd_CEGA,mean_roho=mean_roho,t_test=res$p.value,mean_of_diff=res$estimate,shapiro=shapiro,sign_test=sign_test,sign_test_bonferroni=sign_test_bonferroni,wilcox_test=wilcox_test,nb_replicates=nb_replicates))
})

results_2020 <- t(multiple_t_test)
Mutation_type <- coordinates_restricted[match(results_2020[,1],coordinates_restricted[,1]),7]
results_2020 <- cbind(results_2020,Mutation_type)

write.table(results_2020, file = paste("RAW_04_unfiltered_T_testMinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

out_table_restricted_for_plot <- out_table_restricted
out_table_restricted_for_plot <- out_table_restricted_for_plot[out_table_restricted_for_plot[,7]>=1,]
# 
# # activate this to check for combos only
# to_keep_for_combos <- c(23012,11288,23403,14408,23063,21765)
# out_table_restricted_for_plot <- out_table_restricted_for_plot[out_table_restricted_for_plot[,1]>30000 | out_table_restricted_for_plot[,1]%in%to_keep_for_combos,]
# 
# out_table_restricted_for_plot[out_table_restricted_for_plot[,1]==30001,1] <- 1000
# out_table_restricted_for_plot[out_table_restricted_for_plot[,1]==30002,1] <- 2000
# out_table_restricted_for_plot[out_table_restricted_for_plot[,1]==30003,1] <- 3000
# out_table_restricted_for_plot[out_table_restricted_for_plot[,1]==30009,1] <- 4000

# petite digression pour plotter : 
out_table_ratios_restricted_long <- as.data.frame(cbind(position=out_table_restricted_for_plot[,1],"homo/not_homo"=(log2(out_table_restricted_for_plot[,2]/out_table_restricted_for_plot[,3])/out_table_restricted_for_plot[,7] )))

colnames(out_table_ratios_restricted_long) <- c("position","CEGA")
write.table(out_table_ratios_restricted_long, file = paste("RAW_05_CEGA_Unconcatenated_MinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")


custom_x_titles <- c(0,sort(unique(out_table_ratios_restricted_long[,1])),29903)

table_effectifs <- table(out_table_ratios_restricted_long[,1])
table_effectifs <- as.data.frame(t(rbind(pos=as.numeric(names(table_effectifs)),effectif=as.numeric(unname(table_effectifs)))))

# will define limits for the plots in order to get the exact same plots for both filters
ylimplot_min  <- min(out_table_ratios_restricted_long[,2])
ylimplot_max <- max(out_table_ratios_restricted_long[,2])

# for figure :
svg(filename = paste("fig_5.RoHo_MinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,"_no_rep_nb.svg",sep=""), width = 10, height = 5)
ggplot(out_table_ratios_restricted_long, aes(x=position, y=out_table_ratios_restricted_long[,2], group=position))+
  ggtitle(paste("CEGA plot (ln(roho)/nb_gen) with Nb generation computed using timespan of the tips",sep="")) +
  coord_cartesian(xlim = c(0, max(out_table_ratios_restricted_long[,1])), ylim = c(ylimplot_min,ylimplot_max)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_hline(yintercept = median(out_table_ratios_restricted_long[,2]), linetype='dashed')+
  geom_hline(yintercept = mean(out_table_ratios_restricted_long[,2]))+
  geom_vline(xintercept=table_effectifs$pos, linetype="dotted", color = "grey")+ # only for nrep=10
  #geom_boxplot(width=200) +
   geom_violin(trim=T, fill='#A4A4A4', color="darkred",width=4,position=position_dodge(width = NULL)) +
  stat_summary(geom="point", fun.y=median)+
  scale_x_continuous(breaks=custom_x_titles)+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y=mean, geom="point", shape=1, size=3, color="blue")
dev.off()
# End of petite digression pour plotter

curated_deletions <- read.csv("/home/damien/Copy/42.poly_roho/curated_deletions.txt", header=F,row.names=1,sep=" ",dec=",", stringsAsFactors=F,check.names=FALSE)

# remove deletions that are not approved ;)
# those are the ones called deletions BUT that are not listed in row.names(curated_deletions)
results_2020_modified <- results_2020[results_2020[,13] != "deletion" | results_2020[,1] %in% row.names(curated_deletions) ,]

results_2020[,1]
results_2020_modified[,1]




# add mutated nucleotide here
# it will replace the perl version of the same thing
# for that i need to open the vcf filtered file
cega_to_output <- results_2020_modified

nucl_from <- coordinates_restricted$ref[match(as.numeric(cega_to_output[,1]),as.numeric(coordinates_restricted[,1]))]
nucl_to <- coordinates_restricted$alt[match(as.numeric(cega_to_output[,1]),as.numeric(coordinates_restricted[,1]))]
alt_count <- coordinates_restricted$alt_count[match(as.numeric(cega_to_output[,1]),as.numeric(coordinates_restricted[,1]))]

cega_to_output <- cbind(results_2020_modified[,1],nucl_from,nucl_to,alt_count,results_2020_modified[,2:ncol(results_2020_modified)])
colnames(cega_to_output)[1] <- "position"

cega_to_output[,1]


save.image(file = paste("Workspace_05_timespan_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))
# load(paste("Workspace_05_timespan_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))




nrow(cega_to_output)

write.table(cega_to_output, file = paste("CEGA_mut_type_fromR_T_testMinReplicates_2020_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

svg(filename = paste("fig_6.misplacement_effect_on_CEGA_",args[2],"pcent_misplaced_RoHo_MinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,"_no_rep_nb.svg",sep=""), width = 10, height = 5)

# svg(filename = paste("fig_6.misplacement_effect_on_CEGA_10pcent_misplaced_RoHo_MinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,"_no_rep_nb.svg",sep=""), width = 10, height = 5)

# hist(multiple_t_test[3,],breaks=20,xlab="Mean CEGA value",main=(paste("Hist of CEGA for 30k strains with ",args[2]," % misplaced. ",ncol(multiple_t_test)," Homoplasies",sep="")))

hist(as.numeric(results_2020_modified[,3]),breaks=20,xlab="Mean CEGA value",main=(paste("Hist of CEGA for 130,000 strains with ",args[2]," % misplaced. ",nrow(results_2020_modified)," Homoplasies",sep="")))

dev.off()


# to study how many offspring do have all the kept nodes, i have to represent only the homoplasies that were kept (after nbreplicate chosing) from out_table_restricted

unique(out_table_restricted[,1][!out_table_restricted[,1]%in%as.numeric(results_2020_modified[,1])])

out_table_restricted$homoplasy_count+out_table_restricted$NOT_homoplasy_count

?hist
hist(log10(out_table_restricted$homoplasy_count+out_table_restricted$NOT_homoplasy_count),breaks=50,xlab="log10(nb_total offspring of studied node)",main="")



