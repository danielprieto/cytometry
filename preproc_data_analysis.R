#library(GEOmap)
library(flowCore)
library(flowDensity)
library(flowViz)
library(flowStats)
library(ggcyto)
#Abro archivos
files.test1 = list.files("~/Documentos/citom/LPL-6-6-17/data/", all.files = T, full.names = TRUE)#Lee todos los archivos del directorio
fs <- read.flowSet(files.test1[3:6], transformation = F, alter.names =T)#[3:25])#Elimino las dos primeras filas que son /. y /..
#fs <- fs[, c(1,2,3,4,5,6)]
##
#Transformo logicle
bx <- logicleTransform()
bxlist <- transformList(c("FS.Lin", "SS.Lin", "FL.1.Log", "FL.2.Log", "FL.8.Log"), bx)
datostrans <- transform(fs, bxlist)
##
#Limpiamos el dataset de debris (fcs/ssc)
clean.dt <- datostrans
for (i in 1:4) { # Loop over the length of the flowSet
  f <- datostrans[[i]]
  # First restrict the FSC - A values :
  fsc.indices <- intersect (which (exprs (f)[, "FS.Lin"]<3.6) , which (exprs (f) [, "FS.Lin"] > 3.2))
  # Then restrict SSC - A values and intersect with FSC - A restriction above :
  ssc.indices <- intersect (which ( exprs ( f)[, "SS.Lin"] > 2.2) ,
                            which ( exprs (f)[, "SS.Lin"] < 3.8) )
  non.debris.indices <- intersect ( fsc.indices , ssc.indices )
  # Now only select the non debris cells and place the cleaned up flowFrame into the flowSet :
  f.clean <- f[non.debris.indices]
  clean.dt [[i]] <- f.clean
}
##
#Limpio por canal

##
#Muestro los plots
#par (mfrow = c (2, 2) , mar = c (3, 3, 3, 1))# , mgp =c (2 , 1, 0) )
#plotDens (clean.dt[[1]], c("FS.Lin" ,"SS.Lin"))
#plotDens (clean.dt[[1]], c (4, 6))
#abline (v=3, lwd=1, col="blue")
#abline (h=2.5, lwd=1, col="blue")
#plotDens (clean.dt[[1]], c (5, 6))
#plotDens (datostrans[[1]] , c (6 , 9) )
##
#Matriz de compensacion

comp.files1 = list.files("~/Documentos/citom/LPL-6-6-17/compdata/", all.files = T, full.names = TRUE)#Lee todos los archivos del directorio
comp.fs <- read.flowSet(comp.files1[3:6], transformation = F, alter.names = T)#Elimino las dos primeras filas que son /. y /..
##
cx <- logicleTransform()
cxlist <- transformList(c("FS.Lin", "SS.Lin", "FL.1.Log", "FL.2.Log", "FL.8.Log"), cx)
comp.fs <- transform(comp.fs, cxlist)
##
namepatt <- "FL.1.Log|FL.2.Log|FL.8.Log"
comp.mat <- spillover(comp.fs, unstained = sampleNames(comp.fs[4]), patt = namepatt, method = "mean", stain_match = "ordered", 
                      useNormFilt = T, pregate = F, plot = F, fsc = "FS.Lin", ssc = "SS.Lin")#Genero la matriz de compensacion
#write.table(comp.mat, file = "~/Documentos/citom/22-03-16/comp_mat_median.csv")#Mando la matriz a un archivo
#
#write.table(comp.mat, file = "~/Documentos/citom/s100a9_paper/CD3-CD14-CD19-S100A9/compensation_matrix.csv")#Mando la matriz a un archivo
##
#Compensacion
#comp.mat <- read.csv(file = "~/Documentos/citom/22-03-16/comp_mat_median.csv", sep = " ")#uso la matriz guardada para no calcular cada vez
compensados <- compensate(clean.dt, comp.mat)
#compensados <- fsApply(fs,function(frame){ 
  #extract compensation matrix from keywords 
  #comp <- keyword(frame)$`SPILL` 
 # new_frame <- compensate(frame,comp.mat) 
  #new_frame })

##
#Muestro los plots
#par (mfrow = c (2, 2) , mar = c (3, 3, 3, 1))# , mgp =c (2 , 1, 0) )
#plotDens (compensados[[1]], c("FSC.A" ,"SSC.A"))
#plotDens (compensados[[1]], c (4, 6))
#abline (v=3, lwd=1, col="blue")
#abline (h=2.5, lwd=1, col="blue")
#plotDens (clean.dt[[1]], c (5, 6))
#plotDens (datostrans[[1]] , c (6 , 9) )
##

# Apply debris gate to whole flow set :
# First , start with a copy of the original flowSet 'fs ' (normally you will work directly with 'fs ' itself )
#Hago el dot-plot 
linfos <- lymphGate(compensados, channels=c("FS.Lin", "SS.Lin"), preselection=NULL,scale=2, plot = F)
singfilt1 <- lymphGate(linfos$x, channels = c("FS.Lin", "Pulse.Width"), preselection = NULL, scale=2, bwFac = 1.5, 
                       filterId = "singGate", evaluate = T, plot =F)
linfossing <- Subset(linfos$x, singfilt1$n2gateResults)
#igmfilt <- rectangleGate("SS.Lin" = c(3.1, 3.9), "FL.2.Log" = c(2.2, 4.5))
#igmpos <- Subset(linfossing, igmfilt)
igmcd5filt <- rectangleGate("FL.1.Log" = c(0, 3.9), "FL.2.Log" = c(2.2, 4.5))
#igmcd5filt <- rectangleGate("FL.1.Log" = c(3.4, 3.6), "FL.2.Log" = c(3.5,4))
igmcd5pos <- Subset(linfossing, igmcd5filt)
LPLfilt <- rectangleGate("FL.2.Log" = c(0.75, 1.1), "FL.8.Log" = c(3.4, 4))
#igmnegs100filt <- rectangleGate("FL4.A" = c(3, 3.99), "FL2.A" = c(0.5, 2.3))
qgate1 <- quadrantGate(linfossing, stains = c("FL.2.Log","FL.8.Log"), plot = F, sd = 0.8, alpha = 3)
#qgate2 <- quadrantGate(igmpos, stains = c("FL.1.Log","FL.2.Log"), plot = F, filterId = "igmS100", sd=c(-0.2, 2))
filtrado <- filter(linfossing, qgate1)
#positivasmfi <-Subset(linfossing, igmfilt)
#xyplot(`FL.2.Log` ~ `FL.1.Log`, data = linfossing,   ylab="CD19-PE", xlab="CD5-FITC", smooth=F, xlim=c(0,4), ylim=c(1,4), stat = F)#scatter pĺot
#xyplot(`PE.A` ~ `APC.A`, data = linfossing,  filter=qgate1, ylab="IgM", xlab="LPL", smooth=F,  xlim=c(0,3), ylim=c(0,4))#pData(igmpos)$name <- c("1", "2", "3"),, stat = F, layout = c(2,3))#scatter pĺot5
xyplot(`FL.1.Log` ~ `FL.2.Log`, data = linfossing, smooth=F, xlim=c(0,1), ylim=c(3,4), stats = T)#, layout = c(3,6))#scatter pĺot
xyplot(`FL.2.Log` ~ `FL.8.Log`, data = linfossing, filter=LPLfilt, smooth=F, xlim=c(0,4), ylim=c(0.4,1.1), stats = T)#scatter pĺot
##
#
#library(flowStats)
#c2f <- curv2Filter("FL2.A", "FL4.A", bwFac = 2.1)
#c2f.results <- filter(compensados, c2f)
#xyplot(`FL2.A` ~ `FL4.A`, data = compensados, filter = c2f.results, ylab="IgM-PE", xlab="S100-APC", smooth=F, xlim=c(0,4), ylim=c(0,4), names=F)
##
#
#parallel(~., clean.dt, filter = c2f.results, alpha = 0.01)
##
#Filtro linfos (en flowSet)
#lymph <- lymphGate(compensados, channels = c("SSC-A", "FSC-A"), preselection = NULL, scale=1.5, bwFac = 2, 
#                  filterId="defaultLymphGate", evaluate = TRUE, plot = T)
##
#Cargo los FMO
#fmo.files = list.files("~/Documentos/citom/emmprin14-6-16/fmo/", all.files = T, full.names = TRUE)#Lee todos los archivos del directorio
#fmo.fs <- read.flowSet(fmo.files[3:5], transformation = F, alter.names = T)#[3:25])#Elimino las dos primeras filas que son /. y /..
#fmo.fs <- fmo.fs[, 7:12]
#Los transformo logicle
#cx <- logicleTransform()
#cxlist <- transformList(c("FL.1.Log","FL.2.Log","FL.8.Log"), cx)
#fmo.trans <- transform(fmo.fs, cxlist)
#fmo.comp <- compensate(fmo.trans, comp.mat)3.437
##
#Filtro linfos en un flowFrame del set
#lymph <- flowDensity(obj = compensados[[2]], channels = c("FSC.A", "SSC.A"), position = c(F, F), debris.gate = c(F,F))#Filtro linfos F/S
#sglt <- flowDensity(obj = lymph, singlet.gate = T)#Filtro singuletes
#IgMpos <- flowDensity(obj = sglt, channels = c(3, 2), position = c(T, NA), upper = c(F, NA))#, use.percentile = T, percentile = c(0.01, NA), use.control = c(T, F), control = c(fmo.comp[[3]], NA))#Filtro las cd19+ dentro de los singuletes
#cd3pos <- flowDensity(obj = sglt, channels = c(5, 2), position = c(T, NA), use.percentile=c(T, F), percentile=c(0.35, NA))#, use.control = c(T, F), control = c(clean.dt[[2]], NA))
#S100pos <- flowDensity(obj = sglt, channels = c(6, 2), position = c(T, NA),use.percentile=c(T, F), percentile=c(0.35, NA))#, use.control = c(T, F), control = c(fmo.comp[[2]], NA))
#IgMemmprinpos <- flowDensity(obj = sglt, channels = c(13, 15), position = c(T, T), upper = c(F, F), use.percentile = c(T, F), percentile = c(0.01, NA), use.control = c(T, T), control = c(fmo.comp[[3]], fmo.comp[[2]]))
#cd3emmprinpos <- flowDensity(obj = sglt, channels = c(11, 15), position = c(T, T), upper = c(T, F), use.percentile = c(F, T), percentile = c(NA, 0.01), use.control = c(T, T), control = c(fmo.comp[[1]], fmo.comp[[2]]))
#IgMnegcd3pos <- flowDensity(obj = sglt, channels = c(11, 13), position = c(T, T), upper = c(T, F),  use.control = c(T, T), control = c(fmo.comp[[1]], fmo.comp[[3]]))
#Ki67pos <- flowDensity(obj = sglt, channels = c(18, 9), position = c(T, NA), upper = c(T, NA), use.control = c(T, F), control = c(fmo.comp[[2]], NA))

#s100pos <- flowDensity(obj = sglt, channels = c(12, 9), position = c(T, NA), upper = c(T, NA), use.percentile =T, percentile = c(0.98, NA),  use.control = c(T, T), control = c(fmo.comp[[3]], NA))

#IgmS100pos <- flowDensity(obj = sglt, channels = c(6, 4), position = c(T, T), upper=c(T,T))#, use.control = c(T,T), control = c(comp.fs[[2]], comp.fs[[3]]))#, upper = c(F, T), use.percentile = c(F,T), percentile = c(NA, 0.5), use.control = c(T, T), control = c(fmo.comp[[1]], fmo.comp[[3]]))#Filtro las cd19+ dentro de los singuletes

#s100Ki67pos <- flowDensity(obj = sglt, channels = c(18, 12), position = c(T, T), upper = c(F, T), use.percentile = c(F,T), percentile = c(NA, 0.5), use.control = c(T, T), control = c(fmo.comp[[2]], fmo.comp[[3]]))#Filtro las cd19+ dentro de los singuletes

#IgMKi67pos <- flowDensity(obj = sglt, channels = c(15, 18), position = c(T, T), upper = c(F, T), use.control = c(T, T), control = c(fmo.comp[[1]], fmo.comp[[2]]))#Filtro las cd19+ dentro de los singuletes

#plot(getflowFrame(sglt), IgmS100pos)

##
#Histogramas comparados
#comparativo <- densityplot (~., igmpos, channels = c("FL.8.Log"), filter = igmposs100filt, layout = c(1,1), xlim = c(3,4))#Compara los histogramas de datos transformados en FL1, 3 y4
#comparativo
#Hago un marker con rangegate
#library(flowStats)
#s100 <- rangeGate(igmpos, "FL.8.Log", plot = T, alpha = "min", filterId = "igmposs100filt") #hago un marker
#conMarker <- densityplot (~.,  igmpos, channels= c("FL.8.Log"), filter = igmposs100filt, layout = c(1,1), refline = s100@min, xlim = c(2,4))
#conMarker
pData(igmpos)$name <- c("236-t1", "248-t0", "296-t1", "312-t0", "351-t0", "351-US", "352-t2", "364-t2", "370-t1", "386-t1", "387-t2", "437-t0", "446-t0", "HD")
ggcyto(igmpos, aes(x = `APC.A`, fill = name)) + geom_density(alpha = 0.2)
ggplot(igmpos, aes(x = `APC.A`, fill = name )) + geom_density(alpha = 0.2) + xlim (0,3)
