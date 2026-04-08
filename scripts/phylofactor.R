# Load necessary libraries
library(phylofactor)
library(ape)  # For phylogenetic tree handling
library(ggtree)
library(ggplot2)

args <- commandArgs()
subcom <- args[4]

path <- "phylofactor"
out_dir <- args[5]

tree <- read.tree(args[2])  # Load your tree
all_traits <- read.csv(args[3], row.names = 1, sep = "\t")  # Load your traits matrix
traits <- all_traits[, subcom]

names(traits) <- rownames(all_traits)  # Ensure species names are set
traits <- traits[tree$tip.label]

# Perform phylogenetic factorization
phylofactor_result <- twoSampleFactor(traits,tree, method = "Fisher", nfactors = 7)

factors <- c()                                                                                           
for(i in seq_along(phylofactor_result$groups)){
    bin <- phylofactor_result$groups[[i]][[1]]
    if (length(bin)>4){ #filter out irrelevant factors                                
        node <- getMRCA(tree, bin)
        factors <- c(factors, node)
        }             
}

#write output
dir.create(paste0(out_dir, "/clades"))
for(i in seq_along(factors)){
  node <- factors[i]
  members <- extract.clade(tree, node)$tip.label
    writeLines(members, paste0(out_dir, "/clades/", node, ".txt"))
}


#make visualisation
pdf(file=paste0(out_dir,"/tree.pdf"))

n <- length(tree$tip.label)                                                                                           
lwdEdge <- rep(1.5, dim(tree$edge)[1])
layout(mat=matrix(1:3,ncol = 1),heights = rep(c(3,0.1, 1),6))
par(mar = c(0.1,4,0.1,0.1))

plot(tree, type = "phylogram", show.tip.label = F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 
for(i in seq_along(factors)){
  nodelabels("", factors[i], pch = 18, col = "red", frame = "none",cex = 2)
}
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[traits+1], pch = "|", cex = 2,xlab = '',ylab = '') 

rates <- rep(0, length(tree$tip.label))
assigned <- c()
for(i in seq_along(phylofactor_result$groups)){
    members <- phylofactor_result$groups[[i]][[1]]
    if (length(members)>4){ #filter out irrelevant factors                                
        for(j in seq_along(phylofactor_result$groups)){
            if(j!=i){
                submembers <- phylofactor_result$groups[[j]][[1]]
                if (length(submembers)>4){ #filter out irrelevant factors                               
                    if (all(submembers %in% members)){
                        members <- members[! members %in% submembers]
                        }
                }  
            }
        }  
    
        rate <- sum(traits[members])/length(members)
        for(k in seq_along(members)){
            rates[members[k]] <- rate
        }
        assigned <- c(assigned,members)
    }
}
not_assigned <- 1:n
not_assigned <- not_assigned[! not_assigned %in% assigned]
rate <- sum(traits[not_assigned])/length(not_assigned)
for(k in seq_along(not_assigned)){
        rates[not_assigned[k]] <- rate
    }
print(rates)

par(mar = c(2,4,0.1,0.1))
plot(1:n, rates,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')

axis(2, at = c(0,0.5,1),labels = c('0','0.5','1') )
mtext(side=2,'Rate',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)

dev.off()

#output rates
output <- data.frame(Genome_ID=tree$tip.label,rate=rates)
write.csv(x = output, file = paste0(out_dir, "/rates.csv"))