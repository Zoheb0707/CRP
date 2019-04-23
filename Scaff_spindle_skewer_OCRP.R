Main = function(Initpop, a = 0.5) {
  s = ScaffSpindGen(Initpop)
  print('hello world')
  ScaffSpindDraw(s, s[[length(s)-2 + 1]], s[[length(s)-2 + 2]])
}

Spindle <- function(InitPop, a=.5) {
     # a is alpha
  x <- InitPop
  k <- 1
  while (x[k] > 0) {
    k <- k + 1
    x[k] <- 2*rbinom(1,x[k-1],.50) + sample(c(0,-1), 1,prob=c(1-a,a))
  }
  x[k] <- max(x[k],0)
  x
}


DrawSpindle <- function(spindle,Base) {
  SpindleScale <- 0.05
  #myCol = col2rgb(sample(colors(), 1))
  polygon(Base[1] + SpindleScale * c(spindle,-spindle[length(spindle):1]), y = Base[2] + c(1:length(spindle),length(spindle):1) - 1, col=rgb(runif(1),runif(1),runif(1), 0.5), border=NULL)
}

OCRP = function(minSize, massBD) {
  totPop = 0 # actual population
  ePop = 0 # effective population, only includes large enough tables (>minSize)
  tables = 0 # total number of tables
  CM = 0 # list of tables
  while(ePop < massBD) {
    totPop = totPop + 1
    J = sample(1:totPop, 1)
    if(J == totPop) { # add a new table
      tables = tables + 1
      CM[tables] = 1
    } else { # add customer to the correct table
      peopleCount = 0
      index = 1
      while(peopleCount < J) {
        peopleCount = peopleCount + CM[index]
        index = index + 1
      }
      CM[index - 1] = CM[index-1] + 1
    }
    ePop = 0
    for(tablesize in CM) { #recalculate ePop
      if(tablesize > minSize) {
        ePop = ePop + tablesize
      }
    }
  }
  CM = CM[CM>minSize] # delete small tables
  CM
}


ScaffSpindGen <- function(InitPop, Scale = 10, a=.5) {
  slope <- 2*InitPop
    # actually, this should depend on alpha
  
    ### Setting up initial config of Chinese restaurant ###
  MassBD = floor(Scale*InitPop*rexp(1))
  #CM = NULL
  #x = floor(runif(1)^(-2)*InitPop)
  #Tot = x
  #while(Tot < MassBD) {
  #  CM = c(CM,x)
  #  x = floor(runif(1)^(-2)*InitPop)
  #  Tot = Tot + x
  #}
  #CM = c(MassBD+x-Tot,CM) #add the init leftmost block
  CM = OCRP(InitPop, MassBD)
  print(MassBD)
  print(sum(CM))
  print(CM)
    ### End of initial config ###
  
  spindles <- list(NULL)
  Posns <- NULL
	# a 3xn array; 1st row is x coord
	# 2nd row is y-coord at bottom of jump
	# 3rd row is y-coord at top of jump
  curPos = c(0,0)
	# x and y coord at top of most current jump
  DX <- 0
  i <- 0
  
  for(j in 1:length(CM)) {
    i = i+1
    spindles[[i]] <- Spindle(CM[j])
    curPos[2] <- length(spindles[[i]])
    Posns <- c(Posns,curPos[1],0,curPos[2])
    
    DX <- rexp(1)
    
    while(curPos[2] > slope*DX) {
      curPos <- curPos + DX*c(1,-slope)
      
      i = i+1
      spindles[[i]] <- Spindle(InitPop)
      
      Posns <- c(Posns,curPos,curPos[2] + length(spindles[[i]])-1)
      curPos[2] <- curPos[2] + length(spindles[[i]])-1
      
      DX <- rexp(1)
    }
    curPos <- c(curPos[1]+(curPos[2]/slope) , 0)
  }
  
  Posns <- c(Posns,curPos[1],0,0)
  
  dim(Posns) = c(3,length(Posns)/3)
  print(i)
  spindles[[i+1]] <- Posns
  spindles[[i+2]] <- SpindleClrAssign(Posns)
  spindles
}


SpindleClrAssign <- function(Posns,ShadeSD = .25) {
  gen <- NULL
  Parent <- NULL
  k<-0
  
  for (i in 1:length(Posns[1,]-1)) {
    if (Posns[2,i] == 0) {
      Parent[i] <- 0
      gen[i] <- 0
    }
    else {
      k <- i-1
      while (Posns[2,k] > Posns[2,i]) k <- k-1
      Parent[i] <- k
      gen[i] <- gen[k]+1
    }
  }
  
  CClrs <- col2rgb(sample(rainbow(sum(Posns[2,]==0)+20),sum(Posns[2,]==0)-1))/255
  Clrs <- NULL
  clrDelta <- max(ShadeSD*sqrt(3/(sum(gen)/sum(gen != 0))) , 1)
    # StDev (Unif[-a,a]) = a/sqrt(3)
    # Set ClrDelta so that a spindle of average generational distance from
    # root of clade has stDev = ShadeSD, in each of R,G,B coords.
  j<-0
  for(i in 1:length(gen)) {
    if (gen[i]==0) {
      j <- j+1
      for (k in (-2):0)
        Clrs[3*i+k] <- CClrs[3*j+k]
    }
    else {
      for (k in (-2):0)
        Clrs[3*i+k] <- Clrs[3*Parent[i]+k]+runif(1,max(-clrDelta,-Clrs[3*Parent[i]+k]), min(clrDelta,1-Clrs[3*Parent[i]+k]))
    }
  }
  
  dim(Clrs) <- c(3,length(Clrs)/3)
  Clrs
}


ScaffSpindDraw <- function(spindles,Posns,Clrs,SpindleScale=.05) {
  par(mar=c(.1,.1,.1,.1))
  plot(c(0,1),c(0,0),type="l",col="white",xlim = c(-max(spindles[[1]])*SpindleScale,max(Posns[1,])),ylim=c(0,max(Posns[3,])),axes=FALSE)
  lines(c(0,0,Posns[1,length(Posns[1,])]),c(max(Posns[3,]),0,0),lwd=2)
  for (j in 1:length(spindles)) {
    polygon(Posns[1,j] + SpindleScale * c(spindles[[j]],-spindles[[j]][length(spindles[[j]]):1]), y = Posns[2,j] + c(1:length(spindles[[j]]),length(spindles[[j]]):1) - 1, col=rgb(Clrs[[1,j]],Clrs[[2,j]],Clrs[[3,j]], 0.5), border=NA)
    lines(Posns[1,j:(j+1)],c(Posns[3,j],Posns[2,j+1]))
  }
}


ScaffSpindDrawSlant <- function(spindles,Posns,Clrs,SpindleScale=.05) {
  par(mar=c(.1,.1,.1,.1))
  plot(c(0,1),c(0,0),type="l",col="white",xlim = c(-max(spindles[[1]]*(0:(length(spindles[[1]])-1)))*SpindleScale/(length(spindles[[1]])-1),max(Posns[1,])),ylim=c(0,max(Posns[3,])),axes=FALSE)
  lines(c(0,0,Posns[1,length(Posns[1,])]),c(max(Posns[3,]),0,0),lwd=2)
  for (j in 1:length(spindles)) {
    polygon(Posns[1,j] + SpindleScale * c(spindles[[j]]*((length(spindles[[j]])-1):0),-spindles[[j]][length(spindles[[j]]):1]*((length(spindles[[j]])-1):0))/(length(spindles[[j]])-1), y = Posns[2,j] + c(1:length(spindles[[j]]),length(spindles[[j]]):1) - 1, col=rgb(Clrs[[1,j]],Clrs[[2,j]],Clrs[[3,j]], 0.65), border=NA)
    lines(Posns[1,j:(j+1)],c(Posns[3,j],Posns[2,j+1]))
  }
}


SkewerDraw <- function(spindles,Posns,Clrs,skwrLvl,MaxLen = 0) {
  skwrLvl <- round(skwrLvl)
  IP = 0
  SL = 0
  i=1
  Spindex <- NULL
  for(j in 1:length(spindles)) {
    if ((Posns[2,j] <= skwrLvl) & (Posns[3,j] > skwrLvl)) {
      Spindex[i] <- j
      SL <- SL + spindles[[j]][skwrLvl-Posns[2,j]+1]
      i = i+1
      IP[i] <- SL
    }
  }
  par(mar=c(.1,.1,.1,.1))
  plot(c(0,IP[length(IP)]),c(0,0),type="l",col="white",xlim = c(0,max(IP[length(IP)],MaxLen)),ylim=c(-1.8,1.8),axes=FALSE)
  #rect(-max(IP[length(IP)],MaxLen),-1.3,max(IP[length(IP)],MaxLen)*1.1,1.3,col="white",border=TRUE)
  
  for (j in 1:length(Spindex)) {
    rect(IP[j],-1,IP[j+1],1,border=NA,col=rgb(Clrs[1,Spindex[j]],Clrs[2,Spindex[j]],Clrs[3,Spindex[j]], 0.5))
    lines(c(IP[j],IP[j]),c(-1.8,1.8),lwd=2)
  }
  lines(c(IP[length(IP)],IP[length(IP)]),c(-1.8,1.8),lwd=2)
  
  IP
}


SkewerLen <- function(spindles,Posns,skwrLvl) {
  skwrLvl <- round(skwrLvl)
  SL = 0
  for(j in 1:length(spindles)) {
    if ((Posns[2,j] <= skwrLvl) & (Posns[3,j] > skwrLvl)) SL = SL + spindles[[j]][skwrLvl-Posns[2,j]+1]
  }
  SL
}


MaxSkewerLen <- function(spindles,Posns) {
  MaxLvl <- max(Posns[3,])
  SL <- NULL
  for (i in 0:300) {
    SL[i+1] <- SkewerLen(spindles,Posns,round(i*MaxLvl/300))
  }
  max(SL)
}
