#' run simulated annealing
#'
#' @param initialSolution a starting solution
#' @param evalFunction function to evaluate each iteration of the solution
#' @param neighborFunction the function to create a new solution from the current one
#' @param initialTemperature the starting temperature (default is 300)
#' @param deltaTemperature how much to decrease the temperature at each step (default is 1)
#' @param temperatureStep how many steps to go from \code{initialTemperature} to 0 (default is NA), overrides \code{deltaTemperature}
#' @param alpha how much the current temperature affects the neighbor solution
#' @param sameTolerance how similar two iterations have to be to stop
#' @param nTry how many overall iterations to do
#' @param tryGood how many good solutions to generate before reducing temperature (default is 10)
#' @param tryBad how many bad solutions to make before reducing temperature (default is 100)
#' 
#' @export
#' @return list of results
sa <- function(initialSolution, evalFunction, neighborFunction, initialTemperature=300, deltaTemperature=1, temperatureStep=NA, alpha=1, sameTolerance=0.00001, nTry=10000, tryGood=10, tryBad=100){
  
  if (!is.na(temperatureStep)){
    deltaTemperature <- (initialTemperature - 0) / temperatureStep
  }
  
  currSolution <- initialSolution
  currEnergy <- evalFunction(initialSolution)
  
  bestSolution <- currSolution
  bestEnergy <- currEnergy
  
  currTemp <- initialTemperature
  currTol <- 1
  iTry <- 1
  
  nIter <- (initialTemperature / deltaTemperature) + 2
  allSolutions <- vector("list", nIter)
  allEnergy <- vector("double", nIter)
  
  while ((currTemp > 0) && (currTol > sameTolerance) && (iTry < nTry)){
    temperatureTry <- 1
    iGood <- 0
    iBad <- 0
    
    # keep "good" new solutions
    tmpNewSolution <- list(10)
    tmpNewEnergy <- double(10)
    while ((iGood < tryGood) && (iBad < tryBad)){
      newSolution <- neighborFunction(currentPopulation=currSolution, currentTemperature=currTemp, alpha=alpha)
      newEnergy <- evalFunction(newSolution)
      
      deltaEnergy <- newEnergy - currEnergy
      
      if (deltaEnergy < 0){
        iGood <- iGood + 1
        
        tmpNewSolution[[iGood]] <- newSolution
        tmpNewEnergy[[iGood]] <- newEnergy
        
      } else {
        energyP <- exp((-1 * deltaEnergy) / currTemp) # calculate energy difference
        #browser(expr=TRUE)
        if (runif(1) < energyP){
          iGood <- iGood + 1
          tmpNewSolution[[iGood]] <- newSolution
          tmpNewEnergy[[iGood]] <- newEnergy
        } else {
          iBad <- iBad + 1
        }
      }
    }
    
    hasEntry <- !(sapply(tmpNewSolution, is.null))
    
    
    if (sum(hasEntry) != 0){
      tmpNewSolution <- tmpNewSolution[hasEntry]
      tmpNewEnergy <- tmpNewEnergy[hasEntry]
      bestNew <- which.min(tmpNewEnergy)
      
      currTol <- sum((currSolution - tmpNewSolution[[bestNew]])^2)
      currSolution <- tmpNewSolution[[bestNew]]
      currEnergy <- tmpNewEnergy[bestNew]
    }
    
    
    if (currEnergy < bestEnergy){
      bestSolution <- currSolution
      bestEnergy <- currEnergy
    }
    
    allSolutions[[iTry]] <- currSolution
    allEnergy[iTry] <- currEnergy
    
    currTemp <- currTemp - deltaTemperature
    
    iTry <- iTry + 1
    
    print(c(iTry, currTol, currTemp, currEnergy))
    #print(currTol)
    #print(currTemp)
    #print(iTry)
    
  }
  return(list(lastSolution=currSolution,
              lastEnergy=currEnergy,
              allSolutions=allSolutions,
              allEnergy=allEnergy,
              nTry=iTry,
              bestSolution=bestSolution,
              bestEnergy=bestEnergy))
}