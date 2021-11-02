
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

new_table <- cbind(nodes_to_keep,parApply_result)
new_table <- new_table[new_table[,4] >= tips_rule & new_table[,3] >= tips_rule,]
colnames(new_table) <- c("homoplasy","node","NOT_homoplasy_count","homoplasy_count","perfect","lin_1_all_0","lin_1_all_1","lin_2_all_0","lin_2_all_1")
# keep only 90% perfect nodes in new_table :
new_table <- new_table[((new_table$lin_1_all_0 > 0.9 * (new_table$lin_1_all_0 + new_table$lin_1_all_1)) | (new_table$lin_1_all_1 > 0.9 * (new_table$lin_1_all_0 + new_table$lin_1_all_1)) ) &  ((new_table$lin_2_all_0 > 0.9 * (new_table$lin_2_all_0 + new_table$lin_2_all_1)) | (new_table$lin_2_all_1 > 0.9 * (new_table$lin_2_all_0 + new_table$lin_2_all_1)) ),]

 raw_out_table <- t(as.matrix(apply(new_table,1,function(line){
  homoplasy <- line[1]
  node <- line[2]
  child_nodes <- architecture_arbre[[node]][architecture_arbre[[node]]>(length(arbre_filtered$tip.label))]
    count_children_node_with_homoplasy <- sum(new_table[new_table[,2]%in%child_nodes,1]==homoplasy)
    
    if( count_children_node_with_homoplasy == 0 ){
      
      NOT_homoplasy_count <- line[3]
      homoplasy_count <- line[4]
      
      if(homoplasy_count >= tips_rule & NOT_homoplasy_count >= tips_rule){
        return(c(homoplasy,homoplasy_count,NOT_homoplasy_count,node))
      }else{return(c(6666666,0,0,0))}
    }else{return(c(0,0,0,0))} # fin de condition sur les children nodes
    
})))
raw_out_table

# it is normal that raw_out_table contains empty lines, they were created in order not to add lines to an existing matrix because this is not memory efficient in R. I will delete them now :
raw_out_table <- raw_out_table[!raw_out_table[,1]%in%c("0"),]

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


raw_out_table[,1] <- coordinates_restricted[match(raw_out_table[,1],coordinates_restricted[,3]),1]
class(raw_out_table) <- "numeric"

coordinates_restricted$alt_count <- correspondance[match(as.numeric(coordinates_restricted[,2]),correspondance[,1]),4]

# calculate distances to root for all branches
library(castor)
all_distances <- get_all_distances_to_root(arbre_filtered, as_edge_count=FALSE)

raw_out_table_with_dist <- raw_out_table

raw_out_table_with_dist$mean <- NA
raw_out_table_with_dist$median <- NA

nb_rep <- 5

out_table_restricted <- raw_out_table_with_dist[raw_out_table_with_dist[,1]%in%names(table(raw_out_table_with_dist[,1])[table(raw_out_table_with_dist[,1])>=nb_rep]),]

# implementation of the correction by the timespan of the tips rather than the branch lengths
metadata_file <- read.csv("../00.inputs/metadata_complete.tsv", row.names=NULL,sep="\t",dec=",", header=T, stringsAsFactors=F,check.names=FALSE)
colnames(metadata_file)

# vectorization to save time
nodes_to_evaluate <- unique(out_table_restricted[,4])
# nodes_to_evaluate2 <- nodes_to_evaluate[1:10] # for tests
subset_metadata_file <- metadata_file[metadata_file[,3]%in%arbre_filtered_tips,]

library(HDInterval)
cl <- makeCluster(4)
clusterExport(cl, c("architecture_arbre","arbre_filtered_tips","subset_metadata_file"))
clusterEvalQ(cl, {
  library(HDInterval)
})
nb_generations_tmp <- t(parSapply(cl,nodes_to_evaluate,function(x){
  tips <- architecture_arbre[[x]][architecture_arbre[[x]]<=(length(arbre_filtered_tips))]
   homoplasy_dates <- as.Date(subset_metadata_file[match(arbre_filtered_tips[tips],subset_metadata_file[,3]),5],"%Y-%m-%d")
  homoplasy_dates <- homoplasy_dates[!is.na(homoplasy_dates)]
  homoplasy_interv_old <- as.numeric(difftime(max(homoplasy_dates), min(homoplasy_dates), "%Y-%m-%d"), units = "days")
  homoplasy_interv_gen <- homoplasy_interv_old/5
  if(homoplasy_interv_gen > 1 && !is.na(homoplasy_interv_gen)){
    return(c(x,homoplasy_interv_old,homoplasy_interv_gen,max(homoplasy_dates), min(homoplasy_dates)))
  }else{
    return(c(x,0,0,max(homoplasy_dates), min(homoplasy_dates)))
  }
}))
stopCluster(cl)

nb_generations_tmp_backup <- nb_generations_tmp
nb_generations <- nb_generations_tmp[,c(1,3)] # for UNcorrected

colnames(nb_generations_tmp) <- c("pos","old","new","max","min")

nb_generations_tmp <- as.data.frame(nb_generations_tmp)
homoplasy_dates <- as.Date(subset_metadata_file[match(arbre_filtered_tips,subset_metadata_file[,3]),4],"%Y-%m-%d")
homoplasy_dates <- homoplasy_dates[!is.na(homoplasy_dates)]

out_table_restricted <- cbind(out_table_restricted,nb_generations[match(out_table_restricted[,4],nb_generations[,1]),2])
colnames(out_table_restricted)[length(colnames(out_table_restricted))] <- "nb_generations"

library(BSDA)

multiple_t_test <- sapply(sort(unique(out_table_restricted[,1])), function(x) {
  homoplasy_count <- out_table_restricted[out_table_restricted[,1] == x,2]
  NOT_homoplasy_count <- out_table_restricted[out_table_restricted[,1] == x,3]
  
  res <- t.test(homoplasy_count,NOT_homoplasy_count, paired=TRUE)
  old_roho <- homoplasy_count/NOT_homoplasy_count
  
  time <- out_table_restricted[out_table_restricted[,1] == x,7] 
  CEGA2 <- (8 / time) * ((homoplasy_count / (homoplasy_count + NOT_homoplasy_count))-1/2)
  CEGA2 <- CEGA2[!time<=1] # fix 17.03.2021 forbid time <= 1
  mean_roho <- mean(old_roho)
  mean_CEGA <- mean(CEGA2)
  median_CEGA <- median(CEGA2)
  sd_CEGA <- sd(CEGA2)
  shapiro <- shapiro.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3])$p.value
  sign_test <- SIGN.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], y = NULL, md = 0, alternative = "two.sided", conf.level = 0.95)$p.value
  sign_test_bonferroni <- p.adjust(sign_test, method = "bonferroni", n=length(unique(out_table_restricted[,1])))
  wilcox_test <- wilcox.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], paired = F, mu=0, alternative = "two.sided")$p.value
  nb_replicates <- length(out_table_restricted[out_table_restricted[,1] == x,2])
  return(c(position=x,mean_CEGA=mean_CEGA,median_CEGA=median_CEGA,sd_CEGA=sd_CEGA,mean_roho=mean_roho,t_test=res$p.value,mean_of_diff=res$estimate,shapiro=shapiro,sign_test=sign_test,sign_test_bonferroni=sign_test_bonferroni,wilcox_test=wilcox_test,nb_replicates=nb_replicates))
})

results_2020 <- t(multiple_t_test)
Mutation_type <- coordinates_restricted[match(results_2020[,1],coordinates_restricted[,1]),7]
results_2020 <- cbind(results_2020,Mutation_type)

curated_deletions <- read.csv("/home/damien/Copy/42.poly_roho/curated_deletions.txt", header=F,row.names=1,sep=" ",dec=",", stringsAsFactors=F,check.names=FALSE)

# remove deletions that are not approved ;)
# those are the ones called deletions BUT that are not listed in row.names(curated_deletions)
results_2020_modified <- results_2020[results_2020[,13] != "deletion" | results_2020[,1] %in% row.names(curated_deletions) ,]

cega_to_output <- results_2020_modified

nucl_from <- coordinates_restricted$ref[match(as.numeric(cega_to_output[,1]),as.numeric(coordinates_restricted[,1]))]
nucl_to <- coordinates_restricted$alt[match(as.numeric(cega_to_output[,1]),as.numeric(coordinates_restricted[,1]))]
alt_count <- coordinates_restricted$alt_count[match(as.numeric(cega_to_output[,1]),as.numeric(coordinates_restricted[,1]))]

cega_to_output <- cbind(results_2020_modified[,1],nucl_from,nucl_to,alt_count,results_2020_modified[,2:ncol(results_2020_modified)])
colnames(cega_to_output)[1] <- "position"

write.table(cega_to_output, file = paste("CEGA_mut_type_fromR_T_testMinReplicates_2020_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")
