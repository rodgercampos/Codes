# Empirical identification of intra-urban subcentralities: a new methodological approach with an application for a developing country 
# Rodger B. A. Campos & Carlos R. Azzoni                  


# Memory clearing
rm(list = ls())

# Setting directory (parent folder)
setwd ("d:/Users/rantunes/Documents/Rodger/BASE DE DADOS")

# Installing packages
pkgs <- c('rgdal', 'maptools', 'raster', 'rgeos', 'dismo', 'readstata13', 
          'pastecs', 'foreing', 'spdep', 'locfit')

for (i in 1:length(pkgs)) {
  if (!library(pkgs[i], character.only =  TRUE, logical.return =  TRUE)) {
    install.packages(eval(pkgs[i]))
  }
  library(pkgs[i], character.only =  TRUE)
}

rm(pkgs, i)

# Loading files
masp_dsn <- "GIS/RMSP ShapeFiles"
masp_shp <- "RMSP_dissolve"
masp.grid_shp <- "masp.grid"

# Reading shapefiles 
masp.o <- readOGR(dsn = masp_dsn, layer = masp_shp)       # RMSP
masp <- masp.o


# Gridding shapefile 
plot(masp) #plotting the map
grid <- raster(extent(masp))
res(grid)<- 1/111.32
proj4string(grid)<-proj4string(masp) 
gridpolygon <- rasterToPolygons(grid)
masp.grid <- intersect(masp, gridpolygon)
plot(masp.grid)
masp.grid$id.masp<-1:length(masp.grid)
head(masp.grid)
summary(masp.grid$id.masp)


# Loading RAIS files 
for (year in c(2002, 2008,2014)){
  rais <- "d:/Users/rantunes/Documents/Rodger/BASE DE DADOS" 
  firmas.o<- read.dta13(paste0(rais,"/RAIS_Geo_Final/rmsp_georais",year,"_firmas_final.dta",sep=""))
  firmas <- firmas.o

# From df to shp
  head(firmas)
  dim(firmas)
  firmas$x1<-as.numeric(firmas$x)
  firmas$y1<-as.numeric(firmas$y)
  coordinates(firmas)<- ~x1+y1
  proj4string(firmas)<-proj4string(masp) # Make the grid have the same coordinate reference system (CRS) as the shapefile
  firmas.df=SpatialPointsDataFrame(firmas, data.frame(id=1:length(firmas))) #creating point shp
  geofirmas.df <- cbind(firmas, firmas.df)
  

  # Aggregating firms information by grid 
  firmas.spjoin <- over(geofirmas.df,masp.grid) #sph overlapping
  dim(firmas.spjoin)
  dim(masp.grid)
  a<-data.frame(firmas) #variáveis chaves
  b<-data.frame(firmas[,1:2])
  firmas.join<-cbind(firmas.spjoin,a)
  firmas.join$firma<-1
  head(firmas.join)
  firmas.joinb.o <-cbind(firmas.spjoin,b)
  firmas.joinb <-firmas.joinb.o[,3:5]
  

  # Aggregating var by cells' grid (id.masp)
  var <- aggregate(cbind(emprego, remdezr, firma,
                         exp1, exp2, exp3, exp4, exp5, exp6, exp7, exp8, 
                         tiposal1,tiposal2,tiposal3,tiposal4,tiposal5,tiposal6,tiposal7,
                         grauinstr1,grauinstr2,grauinstr3,grauinstr4,grauinstr5,grauinstr6,grauinstr7,
                         grauinstr8,grauinstr9,
                         sm1,sm2,sm3,sm4,sm5,sm6,sm7,sm8,sm9,
                         ibge1,ibge2,ibge3,ibge4,ibge5,ibge6,ibge7,ibge8,ibge9,ibge10,ibge11,ibge12,
                         ibge13,ibge14,ibge15,ibge16,ibge17,ibge18,ibge19,ibge20,ibge21,ibge22,
                         sexo1,sexo2
  ) ~ id_masp, data = firmas.join, sum)
  
  var[is.na(var)]<-0
  
# Combining database
  geofirmas.grid <- merge(masp.grid, var, by="id_masp")
  geofirmas.grid$emprego[is.na(geofirmas.grid$emprego)]<-0 # replacing na for 0
 
   writeOGR(geofirmas.grid , dsn="RAIS_Geo_Final", layer=paste("firmas",year,".grid"), 
           driver="ESRI Shapefile", overwrite_layer=T)}


# CA APPROACH MODEL                                

# input map files
year <- 2002 #change the year
geofirmas.grid.o <- readOGR(dsn="RAIS_Geo_Final", layer = paste0("firmas", year," .grid"))
geofirmas.grid <- geofirmas.grid.o
summary(geofirmas.grid)
dim(geofirmas.grid)
spplot(geofirmas.grid, "emprego")


# Bandwidth Selection 
g.aic.gauss.po <- gwr.sel(emprego ~ 1, data=geofirmas.grid, gweight=gwr.Gauss, method = 'aic')
g.aic.bisquare.po <- gwr.sel(emprego ~ 1, data=geofirmas.grid, gweight=gwr.bisquare, method = 'aic')
g.aic.tricube.po <- gwr.sel(emprego ~ 1, data=geofirmas.grid, gweight=gwr.tricube, method = 'aic')

g.aic <-2.063 #bandwidth with lower AIC. 

# Best Model (Generalized Least Squre)
model <- gwr(emprego ~ 1, data=geofirmas.grid, gweight = gwr.tricube, bandwidth=g.aic, hatmatrix = TRUE)
summary(model$SDF)
dim(model$SDF)
spplot(model$SDF, "X.Intercept.")

result <- cbind(geofirmas.grid,model$SDF)
summary(result)
result$ttest<-result$X.Intercept./result$X.Intercept._se
summary(result$ttest)
id <- as.data.frame(result[,1])
res <- as.data.frame(result[,72])
pred <- cbind(id,res)


# Bonferroni Adjustment
gwmx <- as.data.frame(model$SDF)
n <- dim(gwmx)[1];n
m <- dim(gwmx)[2];m
nv<-length((model$lm$coefficients))
ntests <- n*nv
enp <- 2*(model$results$v1) - (model$results$v2);enp

pvals <- round(2 * (1 - pt(abs(result$ttest), enp)), 3)
summary(pvals)
bey_pvals <- round(p.adjust(pvals, "BY", n = ntests))
beh_pvals <- round(p.adjust(pvals, "BH", n = ntests))
bon_pvals <- round(p.adjust(pvals, "bonferroni", n = ntests),3)
summary(bon_pvals)

# CBD RULE 
q99 <- quantile(result$X.Intercept., c(0.99))
q95 <- quantile(result$X.Intercept., c(0.95))
q90 <- quantile(result$X.Intercept., c(0.90))
q50 <- quantile(result$X.Intercept.,c(0.50))

# 99% Statistics Significance 
result$sbd99.1p<- ifelse((result$X.Intercept.>=q99 & bon_pvals<=0.01),1,0)
result$sbd95.1p<- ifelse((result$X.Intercept.>=q95 & bon_pvals<=0.01),1,0)
result$sbd90.1p<- ifelse((result$X.Intercept.>=q90 & bon_pvals<=0.01),1,0)
result$sbd50.1p<- ifelse((result$X.Intercept.>=q50 & bon_pvals<=0.01),1,0)

cat("Number of grids identified as part of SBD:  ",sum(result$sbd99.1p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd95.1p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd90.1p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd50.1p),"\n")

# 95% Statistic Significance 
result$sbd99.5p<- ifelse((result$X.Intercept.>=q99 & bon_pvals<=0.05),1,0)
result$sbd95.5p<- ifelse((result$X.Intercept.>=q95 & bon_pvals<=0.05),1,0)
result$sbd90.5p<- ifelse((result$X.Intercept.>=q90 & bon_pvals<=0.05),1,0)
result$sbd50.5p<- ifelse((result$X.Intercept.>=q50 & bon_pvals<=0.05),1,0)

cat("Number of grids identified as part of SBD:  ",sum(result$sbd99.5p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd95.5p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd90.5p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd50.5p),"\n")

# 90% Statistic Significance
result$sbd99.10p<- ifelse((result$X.Intercept.>=q99 & bon_pvals<=0.1),1,0)
result$sbd95.10p<- ifelse((result$X.Intercept.>=q95 & bon_pvals<=0.1),1,0)
result$sbd90.10p<- ifelse((result$X.Intercept.>=q90 & bon_pvals<=0.1),1,0)
result$sbd50.10p<- ifelse((result$X.Intercept.>=q50 & bon_pvals<=0.1),1,0)

cat("Number of grids identified as part of SBD:  ",sum(result$sbd99.10p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd95.10p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd90.10p),"\n")
cat("Number of grids identified as part of SBD:  ",sum(result$sbd50.10p),"\n")

# Save 
write.csv(pred, paste0("Rais_Geo_Final/predito_C&A",year,".csv", sep=""))
writeOGR(result, dsn="RAIS_Geo_Final", layer="SBD2014.final_bonferroni", 
         driver="ESRI Shapefile", overwrite_layer=T)
writeOGR(masp.grid, dsn= masp_dsn, layer= "masp.grid", driver="ESRI Shapefile", overwrite_layer=T)


# MS APPROACH 
data <- data.frame(geofirmas.grid$emprego,geofirmas.grid$long,geofirmas.grid$lat)
names(data) <- c("emprego","long","lat")

subnp <- function(emprego,long,lat,window=.5,pval=.01) {
  
  fit <- locfit(emprego~lp(long,lat,nn=.5,deg=1),kern="tcub",ev=dat(cv=FALSE),data=data)
  mat <- predict(fit,se.fit=TRUE,band="pred")
  yhat <- mat$fit
  sehat <- mat$se.fit
  upper <- yhat - qnorm(pval/2)*sehat
  subobs <- ifelse(emprego>upper,1,0)
  
  cat("Number of tracts identified as part of subcenters:  ",sum(subobs),"\n")
  out <- list(subobs)
  names(out) <- c("subobs")
  return(out)
}

sbd.o <- as.data.frame(lapply(data, subnp))
geofirmas.grid$sbdMS_1pp <- (sbd.o[,1])
names(geofirmas.grid) 


# Campos and Azzoni Adaptation for MS Approach (2003)
subca <- function(emprego,long,lat,window=2.095, pval=.01) {
  
  fit <- locfit(emprego~lp(long,lat,h=2,063,deg=1),kern="tcub",ev=dat(cv=FALSE),data=data)
  mat <- predict(fit,se.fit=TRUE,band="pred")
  yhat <- mat$fit
  sehat <- mat$se.fit
  upper <- yhat - qnorm(pval/2)*sehat
  subobs <- ifelse(emprego>upper,1,0)
  
  cat("Number of tracts identified as part of subcenters:  ",sum(subobs),"\n")
  out <- list(subobs)
  names(out) <- c("subobs")
  return(out)
}

yhat<-as.data.frame(yhat)
write.csv2(yhat, "Rais_Geo_Final/yhat_CA2014_predict.xls")


sbd.o <- as.data.frame(lapply(data, subca))
geofirmas.grid$sbdCA_10pp <- (sbd.o[,1])
names(geofirmas.grid) 
summary(geofirmas.grid)

writeOGR(geofirmas.grid, dsn="RAIS_Geo_Final", layer="SBD2014.final_McMillen & Smith", 
         driver="ESRI Shapefile", overwrite_layer=T)


########################################################################################
#                               MS APPROACH - 2 Stage                                  #
#######################################################################################

# input files 
geofirmas.grid.o <- readOGR(dsn="RAIS_Geo_Final", layer = "SBD2002.final_McMillen & Smith")
geofirmas.grid <- geofirmas.grid.o

shpfileMS_10p <- geofirmas.grid.o[ which(geofirmas.grid.o$sbdMS_10pp==1),]
shpfileMS_10p <- shpfileMS_10p[c("id_masp", "emprego")]
shpfileMS_5p <- geofirmas.grid.o[ which(geofirmas.grid.o$sbdMS_5pp==1),]
shpfileMS_5p <- shpfileMS_5p[c("id_masp", "emprego")]
shpfileMS_1p <- geofirmas.grid.o[ which(geofirmas.grid.o$sbdMS_1pp==1),]
shpfileMS_1p <- shpfileMS_1p[c("id_masp", "emprego")]
shpfileCA_10p <- geofirmas.grid.o[ which(geofirmas.grid.o$sbdCA_10pp==1),]
shpfileCA_10p <- shpfileCA_10p[c("id_masp", "emprego")]
shpfileCA_5p <- geofirmas.grid.o[ which(geofirmas.grid.o$sbdCA_5pp==1),]
shpfileCA_5p <- shpfileCA_5p[c("id_masp", "emprego")]
shpfileCA_1p <- geofirmas.grid.o[ which(geofirmas.grid.o$sbdCA_1pp==1),]
shpfileCA_1p <- shpfileCA_1p[c("id_masp", "emprego")]

dim(shpfileCA_1p)

neighborsMS_10p <- poly2nb(shpfileMS_10p,queen=TRUE)
wmatMS_10p <- nb2mat(neighborsMS_10p,style="B",zero.policy=TRUE)
wmatMS_10p[row(wmatMS_10p)==col(wmatMS_10p)] <- 1

neighborsMS_5p <- poly2nb(shpfileMS_5p,queen=TRUE)
wmatMS_5p <- nb2mat(neighborsMS_5p,style="B",zero.policy=TRUE)
wmatMS_5p[row(wmatMS_5p)==col(wmatMS_5p)] <- 1

neighborsMS_1p <- poly2nb(shpfileMS_1p,queen=TRUE)
wmatMS_1p <- nb2mat(neighborsMS_1p,style="B",zero.policy=TRUE)
wmatMS_1p[row(wmatMS_1p)==col(wmatMS_1p)] <- 1

neighborsCA_10p <- poly2nb(shpfileCA_10p,queen=TRUE)
wmatCA_10p <- nb2mat(neighborsCA_10p,style="B",zero.policy=TRUE)
wmatCA_10p[row(wmatCA_10p)==col(wmatCA_10p)] <- 1

neighborsCA_5p <- poly2nb(shpfileCA_5p,queen=TRUE)
wmatCA_5p <- nb2mat(neighborsCA_5p,style="B",zero.policy=TRUE)
wmatCA_5p[row(wmatCA_5p)==col(wmatCA_5p)] <- 1

neighborsCA_1p <- poly2nb(shpfileCA_1p,queen=TRUE)
wmatCA_1p <- nb2mat(neighborsCA_1p,style="B",zero.policy=TRUE)
wmatCA_1p[row(wmatCA_1p)==col(wmatCA_1p)] <- 1


shpfileMS_10p$wemp_ms_10p <- wmatMS_10p%*%shpfileMS_10p$emprego
shpfileMS_5p$wemp_ms_5p <- wmatMS_5p%*%shpfileMS_5p$emprego
shpfileMS_1p$wemp_ms_1p <- wmatMS_1p%*%shpfileMS_1p$emprego

shpfileCA_10p$wemp_ca_10p <- wmatCA_10p%*%shpfileCA_10p$emprego
shpfileCA_5p$wemp_ca_5p <- wmatCA_5p%*%shpfileCA_5p$emprego
shpfileCA_1p$wemp_ca_1p <- wmatCA_1p%*%shpfileCA_1p$emprego


shpfileMS_10p$sbd_ms_10p_e2 <- ifelse(shpfileMS_10p$wemp_ms_10p>10000,1,0)
shpfileMS_5p$sbd_ms_5p_e2 <- ifelse(shpfileMS_5p$wemp_ms_5p>10000,1,0)
shpfileMS_1p$sbd_ms_1p_e2 <- ifelse(shpfileMS_1p$wemp_ms_1p>10000,1,0)

shpfileCA_10p$sbd_ca_10p_e2 <- ifelse(shpfileCA_10p$wemp_ca_10p>10000,1,0)
shpfileCA_5p$sbd_ca_5p_e2 <- ifelse(shpfileCA_5p$wemp_ca_5p>10000,1,0)
shpfileCA_1p$sbd_ca_1p_e2 <- ifelse(shpfileCA_1p$wemp_ca_1p>10000,1,0)

geofirmas.grid <- geofirmas.grid.o
t1 <- as.data.frame(shpfileMS_10p)
t2 <- as.data.frame(shpfileMS_5p)
t3 <- as.data.frame(shpfileMS_1p)
t4 <- as.data.frame(shpfileCA_1p)
t5 <- as.data.frame(shpfileCA_5p)
t6 <- as.data.frame(shpfileCA_10p)

geofirmas.grid<-merge(geofirmas.grid, t1, by="id_masp")
geofirmas.grid<-merge(geofirmas.grid, t2, by="id_masp")
geofirmas.grid<-merge(geofirmas.grid, t3, by="id_masp")
geofirmas.grid<-merge(geofirmas.grid, t4, by="id_masp")
geofirmas.grid<-merge(geofirmas.grid, t5, by="id_masp")
geofirmas.grid<-merge(geofirmas.grid, t6, by="id_masp")

names(geofirmas.grid)
summary(geofirmas.grid$sbd_ms_1p_e2)
names(shpfileMS_1p)

result<- geofirmas.grid[c("id_masp","wemp_ca_1p","sbd_ca_1p_e2","wemp_ca_5p","sbd_ca_5p_e2","wemp_ca_10p",
                          "sbd_ca_10p_e2","wemp_ms_1p","sbd_ms_1p_e2","wemp_ms_5p","sbd_ms_5p_e2",
                          "wemp_ms_10p","sbd_ms_10p_e2")]
result<-as.data.frame(result)

export(result, "Rais_Geo_Final/SBD2002.final_McMillen & Smith.csv")
