# ��������� ���������� �������� �� ����� �������� ���������
# �������� ������ ���������� ������� - ���� ���������
# ��������� ����������� �� ����������� �������� v1
# �������� ����� ������ ��������� ���������� �������� � ���������� ����������� ��������
# �������� ����� ��������� ����� ������� �� ��������� ��������� ��������
  # !������ ������� ����������� ������ ����� ���� ���������, � ����� ��������� ��� � histNewData ����
  # !���� ��������� ���� �� ������� � ������������� ����������� ���������
  # !�������� �� ������� ������� ����������� � ��� �� ����������� ���������
# !��������� �����������, ������ ����� ���������� ���� ������ ����������� �� ������ ������
# !��������� ������������ ���������� ����������� ���������������
# !���������, ������ ���������� ������ ��������� - ���������� ��������� �����


setwd("D:/User Files/Directory/SolarData")

data_zero <- read.csv("data_zero.csv")
data_zero$ds <- as.POSIXct(data_zero$ds)

cyclePeriod_0 <- 48
timeZero_0 <- as.POSIXct("2018-11-23 00:00")

# ������� ����������� ������� timeZero
getIndexTimeZero <- function(x, timeZero){
  
  index <- 1
  Flag <- TRUE
  
  # ����� ������� ����� ����� �������
  for (i in 1:nrow(x)){
    if (x[i,1] == timeZero){
      index = i
      Flag = FALSE
    }
  }
  
  if (Flag){
    cat("������� ������� timeZero �� �������")
    cat('\n')
  }
  
  returnValue(index)
  
}

# ��������� ����� �� ����������� ��������
getCycleSMA <- function(x, timeZero, cyclePeriod){
  
  # C�������� ������ ������ cyclePeriod � ����������
  # � �������� �������� ���������� ������ 1, 2, 3, 4, ..., cyclePeriod ����� �������
  # �� ��� ����������� ������� ���������� �������
  
  index_0 <- getIndexTimeZero(x, timeZero)
  
  cycleSMA <- c(0)
  
  for(i in 1:cyclePeriod){
    cycleSMA[i] <- 0
  }
  
  M_index_res <- cyclePeriod*3
  
  for (i in 1:cyclePeriod){
    
    res_res <- c(0)
    
    k <- 1
    
    index <- index_0 - i + 1
    
    while (index > M_index_res){
      
      res_res[k] <- data_zero$y[index]
      
      index <- index - cyclePeriod
      
      k <- k + 1
      
    }
    
    cycleSMA[i] <- mean(res_res)
    
  }
  
  output <- c(0)
  j <- cyclePeriod
  for (i in 1:cyclePeriod){
    output[i] <- cycleSMA[cyclePeriod]
    cyclePeriod <- cyclePeriod - 1
  }
  
  returnValue(output)
  
}

# ��������� ��������� ������
getCycleSMADATA <- function(x, timeZero, cyclePeriod){
  
  index_0 <- getIndexTimeZero(x, timeZero)
  data_zero_cycle <- x
  data_zero_cycle$ds <- as.POSIXct(data_zero_cycle$ds)
  cycleSMA <- getCycleSMA(x, timeZero, cyclePeriod)
  index_down <- index_0
  index_up <- index_0
  
  while (index_down > 0){
    
    for (i in cyclePeriod_0:1){
      
      if ((index_down - cyclePeriod_0 + i) == 0){
        break
      }
      
      data_zero_cycle$y[index_down - cyclePeriod_0 + i] <- cycleSMA[i]
      
    }
    
    index_down <- index_down - cyclePeriod_0
  }
  
  while (index_up < nrow(x)){
    
    for (i in 1:cyclePeriod_0){
      
      if ((index_up + i) == nrow(x)){
        break
      }
      
      data_zero_cycle$y[index_up + i] <- cycleSMA[i]
      
    }
    
    index_up <- index_up + cyclePeriod_0
  }
  
  returnValue(data_zero_cycle)
  
}

# ��������� ������ �� ������� ��������� ������
getFlawDATA <- function(x, timeZero, cyclePeriod){
  
  index_0 <- getIndexTimeZero(x, timeZero)
  
  data_zero_flaw <- x
  data_zero_flaw$ds <- as.POSIXct(data_zero_flaw$ds)
  buffer <- getCycleSMADATA(x, timeZero, cyclePeriod)
  data_zero_flaw$y <- x$y - buffer$y
  
  returnValue(data_zero_flaw)
  
}

# ������� ��� ����������� ���������� ���������� ��� ���������
numbMaxOrMin <- function(x, number, typeFunc){
  
  buffer <- x
  output <- data.frame()
  if (typeFunc == "max"){
    buffer <- buffer[order(-buffer[,2]),]
    for (i in 1:number){
      output[i,1] <- buffer[i,1]
      output[i,2] <- buffer[i,2]
    }
  }
  if (typeFunc == "min"){
    buffer <- buffer[order(buffer[,2]),]
    for (i in 1:number){
      output[i,1] <- buffer[i,1]
      output[i,2] <- buffer[i,2]
    }
  }
  returnValue(output)
}


# ������� ����������� ������������� �������
getLikenessCoefs <- function(x, index, step, M, histNewData, compareType){
  
  likeness <- data.frame()
  k <- 1
  lastNew <- histNewData[nrow(histNewData),2]
  index <- index - step
  
  while (index > M){
    
    histOldData <- x[c((index - M):(index - 1)),]
    lastOld <- histOldData[nrow(histOldData),2]
    likeness[k,1] = index
    
    for (i in 1:nrow(histOldData)){
      if (!is.na(histOldData[i,2])){
        chekOld = i
      }
    }
    
    for (i in 1:nrow(histNewData)){
      if (!is.na(histNewData[i,2])){
        chekNew = i
      }
    }
    
    if (is.na(chekOld) || is.na(chekNew)){
      likeness[k,2] = 0
    }
    else{
      
      xy <- histOldData[,2]
      yy <- histNewData[,2]
      
      # ����������� ���������� �������
      if (compareType == "pearson"){
        likeness[k,2] = abs(cor(xy, yy, method = "pearson"))
      }
      
      # ������������ �������
      if (compareType == "quadr"){
        likeness[k,2] = sqrt((sum((xy - yy)^2))/length(xy))
      }
      
      # ������������ ������� �� ���������
      if (compareType == "quadr_diff"){
        likeness[k,2] = sqrt((sum(((xy-mean(xy)) - (yy-mean(yy)))^2))/length(xy))
      }
      
      # ������������ ������� �������������
      if (compareType == "quadr_norm"){
        likeness[k,2] = sqrt((sum(((xy/mean(xy)) - (yy/mean(yy)))^2))/length(xy))
      }
      
      # ������������ ������� � ���������� ���������� � �������� ��������
      if (compareType == "quadr_last"){
        likeness[k,2] = sqrt((sum(((xy-lastOld) - (yy-lastNew))^2))/length(xy))
      }
    }
    
    k <- k + 1
    index = index - step
    
  }
  
  returnValue(likeness)
  
}

# ������� �������� ������������ ��������
getAB <- function(x, histNewData, MSP, M){
  
  # ������� ������������� �������
  MSPData <- x[c((MSP - M):(MSP - 1)),]
  
  # ����� ������������� �������� ����������
  X <- matrix(nrow = nrow(MSPData), ncol = 2)
  for (i in 1:nrow(MSPData)){
    X[i,1] <- MSPData[i,2]
    X[i,2] <- 1
  }
  
  Y <- matrix(nrow = nrow(MSPData), ncol = 1)
  for (i in 1:nrow(histNewData)){
    Y[i,1] <- histNewData[i,2]
  }
  
  E <- matrix(nrow = nrow(MSPData), ncol = 1)
  for (i in 1:nrow(MSPData)){
    E[i] <- 1
  }
  
  tX <- t(X)
  Xn <- tX%*%X
  Yn <- tX%*%Y
  invX <- solve(Xn)
  A <- invX%*%Yn
  
  returnValue(A)
  
}

# ����������� ���������
getForecast <- function(x, histNewData, M, MSPAll){
  
  forecastAll <- data.frame()
  
  for (j in 1:length(MSPAll)){
    
    MSP <- MSPAll[j]
    
    A <- getAB(x = x, histNewData = histNewData, MSP = MSP, M = M)
    
    # ������� �������
    histBaseData <- x[c((MSP):(MSP+1)),]
    
    # ��������������
    XP <- matrix(nrow = 1, ncol = 2)
    
    for (l in 1:nrow(histBaseData)){
      
      XP[1,1] <- histBaseData[l,2]
      XP[1,2] <- 1
      
      forecastAll[l,j] <- XP%*%A
      
    }
    
  }
  
  returnValue(forecastAll)
  
}

# ���������� ���������
getMeanForecast <- function(forecastAll, compareNumber){
  
  forecastX <- forecastAll[,1]
  
  for (l in 1:nrow(forecastAll)){
    
    sumbuffer <- 0
    
    for (p in 1:compareNumber){
      
      sumbuffer <- sumbuffer + forecastAll[l,p]
      
    }
    
    forecastX[l] <- sumbuffer/compareNumber
    
  }
  
  returnValue(forecastX)
  
}

# ������
EMMSPPoints <- function(x, timeZero, P, nSegmentsCor, findNSegment, step, compareType, compareNumber, cycles, cyclePeriod){

  index <- getIndexTimeZero(x, timeZero)
  forecastSum <- data.frame()

  if (cycles){
    
    buffer <- getFlawDATA(x, timeZero, cyclePeriod)
    forecastSum <- buffer[1:(index + P),]
    forecastSum$ds <- as.POSIXct(forecastSum$ds)
    
  }
  else{
    
    forecastSum <- x[1:(index + P),]
    forecastSum$ds <- as.POSIXct(forecastSum$ds)
    
  }
  
  nSegmentsCorFalse <- nSegmentsCor
  if (!findNSegment){
    nSegmentsCor = 1
  }
  
  for(indexes in (index + 1):(index + P)){
    
    cat("������� �����:", indexes - index)
    
    # ����� ������ ����� ���������� �� ���������

    if (findNSegment){
      dispersions <- c(1:nSegmentsCor)
      for (di in  1:length(dispersions)){
        dispersions[di] <- 9999999
      }
    }
    else{
      dispersions <- c(1:nSegmentsCor)
    }
    
    forecastN <- data.frame()
    
    for (nn in 1:nSegmentsCor){
      
      if (findNSegment){
        if(nn == 1){
          dispersions[nn] <- 9999999
          forecastN[1,nn] <- 0
          next
        }
      }

      if (!findNSegment){
        M <- step*nSegmentsCorFalse
      }
      else{
        M <- step*nn
      }
      
      # ������� ����� �������
      histNewData <- forecastSum[c((indexes - M):(indexes - 1)),]
      
      # ����������� �������� �������
      likeness <- getLikenessCoefs(x = forecastSum, index = indexes, step = step, M = M, 
                                   histNewData = histNewData, compareType = compareType)
      
      # ����������� ���������� �������
      if (compareNumber > nrow(likeness)){
        compareNumber = nrow(likeness)
      }
      
      likenessValues <- data.frame()
      
      if (compareType == "pearson"){
        likenessValues <- numbMaxOrMin(x = likeness, number = compareNumber, typeFunc = "max")
      }
      else{
        likenessValues <- numbMaxOrMin(x = likeness, number = compareNumber, typeFunc = "min")
      }
      
      # ����������� ���������
      MSPAll <- likenessValues[c(1:nrow(likenessValues)),1]
      
      forecastAll <- getForecast(x = forecastSum, histNewData = histNewData, M = M, MSPAll = MSPAll)
      
      # ���������� ���������
      forecastX <- getMeanForecast(forecastAll = forecastAll, compareNumber = compareNumber)
      
      # ������ ���������
      bufferDispersoin <- c(0)
      for (d in 1:length(MSPAll)){
        bufferDispersoin[d] <- forecastAll[1,d]
      }
      
      dispersions[nn] <- var(bufferDispersoin)
      
      for(jj in 1:length(forecastX)){
        forecastN[jj,nn] <- forecastX[jj]
      }
      
    }
    
    indicatorDisp <- 1
    for (kk in 1:length(dispersions)){
      if (dispersions[kk] == min(dispersions)){
        indicatorDisp = kk
      }
    }
    
    if (!findNSegment){
      cat(", findNSegment FALSE")
      cat(", ����� �������:", M)
      cat(", � ����������:", min(dispersions))
    }
    else{
      cat(", findNSegment TRUE")
      cat(", ������ ����� �������:", indicatorDisp*step)
      cat(", c ������ ����������:", min(dispersions))
    }
    
   
    forecastSum[(indexes),2] <- forecastN[1,indicatorDisp]
    
    cat('\n')
    
  }
  
  if (cycles){
    
    dataFact <- x[c((index + 1):(index +P )),]
    cycle_buffer <- getCycleSMADATA(x, timeZero, cyclePeriod)
    forecastSum[,2] <- forecastSum[,2] + cycle_buffer[c(1:(index + P)),2]
    dataForecast <- forecastSum[c((index + 1):(index + P)),]
    
  }
  else{
    
    dataFact <- x[c((index + 1):(index + P)),]
    dataForecast <- forecastSum[c((index + 1):(index + P)),]
    
  }
  
  plot(ts(data.frame(dataFact[,2], dataForecast[,2])), plot.type="single", col = 1:2)
  
  returnValue(dataForecast)
  
}

# x - �������� ������
# timeZero - ����� ����� ����� ����� �������
# P - ����� ��������
# nSegmentsCor - ������������ ����� ��������� ����������
# findNSegment - ���� FALSE, �� ����� ������ ��������� �� ������������ � ������� nSegmentsCor*step
# step - ��� ������� (����������, ����� ����� ��������)
# compareType = c("pearson", "quadr", "quadr_diff", "quadr_norm", "quadr_last")
# compareNumber - ����� ������ ����������
# cycle - ���� TRUE, �� �������� ��������� � �������
# cyclePeriod - ������ �����������

retran <- EMMSPPoints(data_zero, timeZero_0, 48, 4, FALSE, 48, "quadr_last", 100, TRUE, cyclePeriod_0)

# ���� ������:
# �������� timeZero_0 -> index_0
# ������������ ������� �� ���� �������� data_zero �� index_0 � �� ���������� - preForecast
# ������������ ������� ����� ������� �� index_0 - �� ����� preForecast
# ������������ �������� ������� likeness
# ������������ ��������� �������� ������� likenessValuess
# �������� ������� ����� index_0 + 1 ��� ������� ������� �� likenessValuess
# �������� ����������� � ������������ � forecast
# index_0 ���������� �� 1 ������
# ������������ �������� ����������������� ����� � preForecast
# ���������� ����������� �� index_0 + P
# �������� �������
# ������� ������
