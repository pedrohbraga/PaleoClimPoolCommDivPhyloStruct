install_github('JohnsonHsieh/iNEXT')

remotes::install_github('YanHanChen/iNEXTPD2', force = T)

remotes::install_github('KaiHsiangHu/iNEXT.3D', force = T)



## import packages
library(iNextPD)
library(ggplot2)
library(ade4)

?iNextPD

data(bird)
bird.abu <- bird$abun
bird.lab <- rownames(bird$abun)
bird.phy <- ade4::newick2phylog(bird$tre)

# set a series of sample sizes (m) for R/E computation
m <- c(1, 5, 20, 50, 100, 200, 400)
iNextPD(x=bird$abun, labels=bird.lab, phy=bird.phy, 
        q=0, datatype="abundance", size=m)


ade4::newick2phylog(write.tree(Chiroptera.FaurSven.tree.ultra))


adephylo 

#

library(iNEXTPD2)

# Datatype: incidence_raw data
data(data.inc)
data <- data.inc$data
tree <- data.inc$tree
nT <- data.inc$nT
out <- iNEXTPD(data = t(Chiroptera.NW.Comm), 
               nT = setNames(c(7708, 8189), c("Neotropical", "Nearctic")), 
               datatype = "incidence_raw", 
               tree = makeNodeLabel(Chiroptera.NW.Tree), 
               type = "meanPD",
               q = c(0, 1, 2),
               nboot = 10)

MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Hemispheric sampling") %>%
  filter(ID_Realm %in% c("Neotropical", "Nearctic")) %>%
  group_by(ID_Realm) %>% 
  summarise(n = n())


t(Chiroptera.NW.Comm)[1:6, 1:6]


write.tree(Chiroptera.FaurSven.tree)

makeNodeLabel(Chiroptera.FaurSven.tree)

