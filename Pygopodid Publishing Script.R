library(geomorph)
library(Morpho)
library(phytools)
library(devtools)
library(ggplot2)
library(RColorBrewer)
library(ggConvexHull)


mydata = readland.nts("PygoFinal.dta") #This is the nts file that contains the landmarks
Sliders <- rbind(define.sliders(25:34), define.sliders(35:44)) #This specifies the landmarks that are the semi-landmark curves
gecko.gpa <- gpagen(mydata[ , ,-13], curves = Sliders) #GPA-alignment (Procrustes superimposition)
info <- read.csv("InfoFull.csv", header = T)#This is for the excell file that contains all the specimen info
info <- info[-13, ] #Removing 13 because bavayia robusta is not included in overall PCA
gecko_df1 = geomorph.data.frame(gecko.gpa, Genus = info$Genus, Habitat = info$Habitat, 
                                  Species = info$Species, ID = info$ID, Geo = info$Geo, Diet = info$Diet)

Pygo.PCA <- gm.prcomp(gecko.gpa$coords)
plot(Pygo.PCA)
summary(Pygo.PCA)

Pygo.Phylomorpho <- gm.prcomp(phylo.gpa$coords, phy = edit.complete) #To create the phylomorphospace
plot(Pygo.Phylomorpho, phylo = TRUE)

Pygo.Anova <- procD.lm(gecko.gpa$coords ~ Habitat*Diet, data = gecko_df1, iter = 9999, seed = "random")
Pygo.Anova$aov.table #This is the results table that you see with the p-values

#This is how I conducted post-hoc pairwise analysis on the data
summary(pairwise(fit = Pygo.Anova, covariate = NULL, groups = gecko_df1$Habitat), test.type = "dist")

##I did not do this for the phylogenetically corrected anova because nothing was significant with that

##### Phylo Anova ###
#########################################################################################################
tree.complete <- read.tree("PygopodidTree.tre")
edit.complete <- multi2di(tree.complete) #You need to use this function to use the phylomorphospace function because it formats the tree into the type that is used by the function. This function isnt really good for polytomies because it changes the relationships within the polytomy every time the tree is created
edit.complete$tip.label
#edit.complete$edge.length[21] = 1 #The phylogeny is missing the edge length connecting the clade containing lialis to the other groups so it is assigning an edge length to that index
plot(edit.complete)
info <- read.csv("InfoFull.csv", header = T)#This is for the excell file that contains all the specimen info
Phylo.info <- info[ -c(3,13,22), ]

phylo.gpa <- gpagen(mydata[ , ,-c(3,22)], curves = Sliders)
gecko_df2 = geomorph.data.frame(phylo.gpa, Genus = Phylo.info$Genus, Habitat = Phylo.info$Habitat, 
                                Species = Phylo.info$Species, ID = Phylo.info$ID, 
                                Geo = Phylo.info$Geo, Diet = Phylo.info$Diet)
dimnames(phylo.gpa$coords) <- list( NULL, NULL, Phylo.info$Tree_ID)
phylo.anova <- procD.pgls(phylo.gpa$coords ~ Habitat*Diet, edit.complete, iter = 9999, 
                          seed = "random",data =  gecko_df2)
phylo.anova$aov.table

physignal(phy = edit.complete, phylo.gpa$coords, iter = 9999) 


#### Centroid Size Regression ####
summary(procD.lm(gecko.gpa$coords ~ log(gecko.csize))) #Regression of shape data against centroid size without phylogenetic correction

summary(procD.pgls(phylo.gpa$coords ~ log(gecko.csize[-c(3,22)]), phy = edit.complete))
#Regression of shape data aginst centroid size with phyloegentic correction

##### FULL PCA ####
brewer_test <- brewer.pal(n = 12, "Paired")
brewer_blues <- brewer.pal(n =9, "Blues")
colorfriend2 <-c(rep("black", 12), brewer_blues[c(1,7,1,7,1,7,7,1,7,7,7)], "black",
                 brewer_blues[c(7,1,7,7,7)])
Pygo_gg_shape <- c(rep(21,12), rep(23,9), 22,22,24,23,23,24,24,24)


Final.Score <- cbind(Pygo.PCA$x[ , 1]*-1, Pygo.PCA$x[ , 2]) #Multiplied the first axis by -1 so the orientation matches the phylomorphospace
#In the paper. Note, phylomorphospace is slightly different in this script because it used the depricated function PlotGMPhylomorphospace but orientation and distances are the same
dimnames(Final.Scores) <- list(info$ID, NULL) #How to carry over the names of the taxa

Pygo_gg_data <- data.frame(Final.Score)
ggplot(data = Pygo_gg_data, 
       aes(x= Pygo_gg_data[ ,1], y = Pygo_gg_data[ ,2], colour="black", fill=colorfriend2)) + 
  theme_classic() +
  #You change the shape of the data points using "shape =" in geompoint
  geom_point(shape = Pygo_gg_shape, size = 3, colour = "black", fill = colorfriend2) + 
  xlab("PC 1(53.80%) - Skull Roof Width, Parasphenoid & Occipital Length") +
  ylab("PC 2(11.93%) - Snout Elongation & Orientation" ) 
#geom_convexhull(alpha = 0.2, aes(fill = colorfriend, color = "black"))
#if you want to add taxon labels to each point hjust and vjust modify the position of the label
#geom_text(aes(label = rownames(Pygo_gg_data), hjust = 1, vjust= -1)) 



#### Aprasia and Ophidio PCA ####
fossorial.data <- mydata[ , , c(1:12,25)]
fossorial.info <- info[ c(1:12,25), ]
fossorial.gpa <- gpagen(fossorial.data, curves = Sliders)
fossorial.df1 = geomorph.data.frame(fossorial.gpa, Genus = fossorial.info$Genus, Habitat = fossorial.info$Habitat, 
                                  Species = fossorial.info$Species, ID = fossorial.info$ID, Geo = 
                                    fossorial.info$Geo, Diet = fossorial.info$Diet) 
Fossorial.PCA <- gm.prcomp(fossorial.gpa$coords)
summary(Fossorial.PCA)
Fossorial.PCA$sdev
Fossorial.ANOVA <- procD.lm(Fossorial.gpa$coords ~ Geo, data = Fossorial.df1, RRPP = TRUE, iter = 9999, seed = "random")

cbind(Fossorial.PCA$x[ ,1], Fossorial.PCA$x[ ,2])

Fossorial_gg_data <- data.frame(cbind(Fossorial.PCA$x[ ,1], Fossorial.PCA$x[ ,2]))
Fossorial_gg_data

Fossorial_gg_shape <- c(1,2,2,3,2,1,2,3,2,2,2,3,3)
ggplot(data = Fossorial_gg_data, 
       aes(x= Fossorial_gg_data[ ,1], y = Fossorial_gg_data[ ,2], colour="black")) + 
  theme_classic() +
  #You change the shape of the data points using "shape =" in geompoint
  geom_point(shape = Fossorial_gg_shape-1, size = 3, colour = "black") + 
  xlab("PC 1(27.65%)") +
  ylab("PC 2(18.21%)" ) +
  #geom_convexhull(alpha = 0.2, aes(fill = colorfriend, color = "black"))
  #if you want to add taxon labels to each point hjust and vjust modify the position of the label
  geom_text(aes(label = rownames(Fossorial_gg_data), hjust = 0, vjust= 0))

#### How to create wireframes ####
##I did not know how to write a loop at the time so I just ran this chunk over all the iterations that cover specimen number
i <- 1 #index for specimen number
spec.count <- dim(mydata) #This is how you specify the amount of iterations for the loop
mir.array <- array(data = NA, dim = c(spec.count[1], spec.count[2], spec.count[3])) #This is the array that will store the mirrored data
Full.array <- array(data = NA, dim = c(spec.count[1]*2, spec.count[2], spec.count[3]))#This is the array that will contain everything

mir.array[ , ,i] <- mirror2plane(x = mydata[ , ,i], v1 = mydata[1 , ,i], v2 = mydata[7, ,i],
                                 v3 = mydata[16, ,i]) 
Full.array[ , ,i] <- rbind(mydata[ , ,i], mir.array[ , ,i])

i = i + 1

#How to view wireframes
#Code not mine: taken from
# http://www.randigriffin.com/2017/11/10/plotting-shape-changes-geomorph.html
#With this you need to give values for points.col, points.cex,lines.col,lines.wd,bg.col,main, and legend arguments
plot.coords <- function( A, W, points.col = NULL, points.cex = NULL, lines.col = NULL, lines.wd = NULL, bg.col = NULL,
                         main = NULL, main.line = 2, main.cex = 2, legend = NULL, legend.pos = "topright"
                         , legend.title = NULL, legend.col = NULL, legend.cex = 1.2, legend.lwd = 2, legend.bty = "n",
                         params = NULL, add = FALSE) {
  if (!is.null(params)) {par3d(params)}
  points.col <- rep(points.col, length.out = nrow(A))
  points.cex <- rep(points.cex, length.out = nrow(A))
  lines.col <- rep(lines.col, length.out = nrow(W))
  lines.lwd <- rep(lines.wd, length.out = nrow(W))
  
  if (!is.null(bg.col))  rgl.bg(sphere = TRUE, color = bg.col, lit = FALSE, back = "fill")
  plot3d(A, type = "s", col = points.col, xlab = "", ylab = "", zlab = "", size = points.cex,
         aspect = FALSE, box = FALSE, axes = FALSE, add=add)
  
  if (!is.null(main) | !is.null(legend)) {
    if (!is.null(legend) & is.null(legend.col)) stop("must supply legend colors")
    bgplot3d({plot.new()
      if (!is.null(main)) title(main = maain, line = main.line, cex.main = main.cex)
      if (!is.null(legend)) legend(legend.pos, title = legend.title, legend = legend, col = legend.col,
                                   lwd = legend.lwd, cex = legend.cex, bty = legend.bty)})}
  for (i in 1:nrow(W)) {
    segments3d(rbind(A[W[i, 1], ] , A[W[i, 2], ]), lwd = 3, col= lines.col[i])
  }
}

#### Vector Plots ####
average <- mshape(gecko.gpa$coords)
average.new <-cbind(average[ , 1], average[ ,2], average[ , 3])
Average.Matrix <- as.matrix(average.new)
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,12], method = "vector", label = F)#aprasia striolata
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,6], method = "vector", label = F)#aprasia parapulchella
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,13], method = "vector", label = F)#Delma Borea
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,22], method = "vector", label = F)#Lialis burtonis
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,24], method = "vector", label = F)#ophidio
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,25], method = "vector", label = F)#Paradelma
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,26], method = "vector", label = F)#Pletholax
plotRefToTarget(M1 = Average.Matrix.m, M2 = gecko.gpa.mirrored$coords[ , ,27], method = "vector", label = F)#pygopus lepido

#######       Supplemtary Plot including bavaya robusta    ########
##Note: Plots look slightly different due to the updated geomorph functions

gecko.gpa.sup <- gpagen(mydata, curves = Sliders) #GPA-alignment (Procrustes superimposition)
info.sup <- read.csv("InfoFull.csv", header = T)
gecko_dfsup = geomorph.data.frame(gecko.gpa.sup,
                                  Species = info.sup$Species, Genus = info.sup$Genus)
PCA.sup <- gm.prcomp(gecko.gpa.sup$coords)
plot(PCA.sup)
Final.Scores.sup <- cbind(PCA.sup$x[ , 1], PCA.sup$x[ , 2])
xrange.sup <- range(PCA.sup$x[ , 1])
yrange.sup <- range(PCA.sup$x[ , 2])
plot(Final.Scores.sup[1:12 , ], cex = 1, 
     xlab = "PC 1(49.90%) - Skull Roof shape, Parasphenoid & Snout Length", ylab = "PC 2(11.03%) - Snout Elongation & Braincase Width",
     data = Final.Scores.sup, type = "p", xlim = xrange.sup, ylim = yrange.sup ,pch = 8, col = "black") #Only aprasia is plotted here
points(Final.Scores.sup[23:24, ], pch = 24, cex = 1, bg = "black") #Adding lialis ontop of the aprasia plot
points(Final.Scores.sup[14:22, ], pch = 22, cex = 1, bg = "black") #Adding Delma ontop of aprasia plot
points(Final.Scores.sup[28:30, ], pch = 21, type = "p", cex = 1, bg = "black") #Pygopus
points(Final.Scores.sup[c(27,27), ], pch = 23, cex = 1, bg = "black")#Pletholax
points(Final.Scores.sup[c(26,26), ], pch = 4, cex = 1, col = "black", lwd = 2)#Paradelma
points(Final.Scores.sup[c(25,25), ], pch = 25, cex = 1, bg = "black")#Ophidiocephalus
points(Final.Scores.sup[c(13,13), ], pch = 3, cex = 1, bg = "black")#bAVAYA ROBUSTA

