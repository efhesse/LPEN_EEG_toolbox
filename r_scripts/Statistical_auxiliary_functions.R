library(ggplot2)

#determines outliers are above x sd from the mean
outlier_determination <- function(x) {
  x0 = abs(x-mean(x,na.rm=TRUE)) > 3*sd(x,na.rm=TRUE)
  
  x[which(x0==TRUE)]=NA
  x
}

plot_logistic_regression <- function(x0,y0,title,xlabel,ylabel,glm.out,file_name_to_save,save){
  
  dat <- data.frame(matrix(ncol = 2, nrow = length(x0)))
  dat[,1]<-x0
  dat[,2]<-y0
  names(dat) <-c('x0','y0')  
  
  print(ggplot(dat, aes(x=x0, y=y0)) + geom_point() +  labs(x=xlabel,y=ylabel) + ggtitle(title) +
          stat_smooth(method="glm", se=FALSE))
  
  if(save == 1)
  {
    dev.copy(png,filename=file_name_to_save);
    dev.off()
  }
}