## Load deSolve package
### Se instala con install.packages("deSolve"), luego seleccionan un mirror y listo.
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
###Sin cuarentena
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.4247, gamma = 0.14286)
## Time frame
times      <- seq(0, 70, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
##head(out, 10)


###Con cuarentena
## beta: infection parameter; gamma: recovery parameter
parameters_with_quarantine <- c(beta = 0.8, gamma = 0.14286)
## Solve using ode (General Solver for Ordinary Differential Equations)
out_with_quarantine <- ode(y = init, times = times, func = sir, parms = parameters_with_quarantine)
## change to data frame
out_with_quarantine <- as.data.frame(out_with_quarantine)



library(ggplot2)
png("sir_cuarentena.png",width=1024,height=720,units="px",bg="white")
## ggplot
ggplot(out)+geom_line(size=1,aes(x=times,y=S,color="Susceptibles"))+geom_line(size=1,aes(x=times,y=I,color="Infectados"))+geom_line(size=1,aes(x=times,y=R,color="Recuperados"))+geom_line(data=out_with_quarantine,size=1,linetype="dotted",aes(x=times,y=S,color="Susceptibles en cuarentena"))+geom_line(data=out_with_quarantine,size=1,linetype="dotted",aes(x=times,y=I,color="Infectados en cuarentena"))+geom_line(data=out_with_quarantine,size=1,linetype="dotted",aes(x=times,y=R,color="Recuperados en cuarentena"))+labs(title='Gráfico modelo matemático SIR con cuarentena, beta = 0.8',colour='',x='Tiempo (días)', y='% de población')
dev.off()
