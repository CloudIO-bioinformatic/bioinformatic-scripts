## Load deSolve package
library(deSolve)

## Create an SIR function
sir <- function(time, state, parameters) {

  with(as.list(c(state, parameters)), {

    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I

    return(list(c(dS, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.4247, gamma = 0.14286)





## Time frame
times      <- seq(0, 189, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
out$times <- seq(0, 189, by = 1)
#contestar a
totalpoblacion <- 17000000
semana21 <- 21*7
semana23 <- 23*7
mayora21 <- out[which(out$time>=semana21),]
menora23 <- mayora21[which(mayora21$time<=semana23),]
enfermosdiariospromedio <- sum(menora23$I)/length(menora23) * totalpoblacion


#contestar b
max_i <- max(out$I)
max_i_5 <- max_i * 0.05 # son los que van quedando infectados
momento <- max(out[which(out$I>=max_i_5),]$time)

#contestar d

muertes <- 0.03 * max_i

## Delete time variable
#out$time <- NULL
## Show data
##head(out, 10)
#library(ggplot2)
#png("sir_simple.png",width=1024,height=720,units="px",bg="white")
## Plot
##matplot(x = times, y = out, type = "l",
##        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
##        lwd = 1, lty = 1, bty = "l", col = 2:4)

## Add legend
##legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")
#ggplot(out)+geom_line(size=1,aes(x=times,y=S,color="Susceptibles"))+geom_line(size=1,aes(x=times,y=I,color="Infectados"))+geom_line(size=1,aes(x=times,y=R,color="Recuperados"))+labs(title='Gráfico modelo matemático SIR, beta = 1.4247',colour='',x='Tiempo (días)', y='% de población')
#dev.off()
