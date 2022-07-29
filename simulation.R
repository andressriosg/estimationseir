# simulated data 

library(deSolve)
solution_1 <- ode(y= c(4, 0, 0.1, 0), 0:101, seir, c(4, 0.2, 0.1, 0.3, 0.2))

solution_2 <- matrix(0, nrow = 102, ncol = 5)
solution_2[1, ] = solution_1[1, ]
solution_2[, 1] = solution_1[, 1]
for (j in 4:5) {
  for (i in 1:50) {
    solution_2[2*i, 2] = solution_1[2*i, 2] + 1
    solution_2[2*i, 3] = solution_1[2*i, 3] + 0.5
    solution_2[2*i, j] = solution_1[2*i, j] + 0.2
    solution_2[(2*i + 1), 2] = solution_1[(2*i + 1), 2] - 1
    solution_2[(2*i + 1), 3] = solution_1[(2*i + 1), 3] - 0.5
    solution_2[(2*i + 1), j] = solution_1[(2*i + 1), j] - 0.2
  } 
}
plot(solution_1[, 2])
lines(solution_2[, 2])

# El modelo determinista conociendo los parametros 

se <- NULL 
ee <- NULL 
ie <- NULL 
re <- NULL 

se[1] = solution_2[1, 2]
ee[1] = solution_2[1, 3]
ie[1] = solution_2[1, 4]
re[1] = solution_2[1, 5]

parametros_e <- c(4, 0.2, 0.1, 0.3, 0.2)

for(j in 1:100){
    se[j+1] = se[j] + (parametros_e[1] - parametros_e[2]*se[j]*ie[j] - parametros_e[5]*se[j])
    ee[j+1] = ee[j] + (parametros_e[2]*se[j]*ie[j] - parametros_e[3]*ee[j] - parametros_e[5]*ee[j])
    ie[j+1] = ie[j] + (parametros_e[3]*ee[j] - parametros_e[4]*ie[j] - parametros_e[5]*ie[j])
    re[j+1] = re[j] + (parametros_e[4]*ie[j] - parametros_e[5]*re[j])
  }

plot(se, pch = 16)
lines(solution_1[, 2], col = "blue", lwd = 2)

plot(ie)
lines(solution_1[, 4])

plot(ee)
lines(solution_1[, 3])

plot(re)
lines(solution_1[, 5])

# Suavizamos los datos 

library(npregfast)
time_1 <- 1:101
suavizado_suscep <- frfast(solution_2[, 2][1:101] ~ time_1, model = "np", smooth = "kernel", kbin = 101, 
                     p =3)
suavizado_suscep <- data.frame(suavizado_suscep$p)$X1
suavizado_expos <- frfast(solution_2[, 3][1:101] ~ time_1, model = "np", smooth = "kernel", kbin = 101, 
                           p =3)
suavizado_expos <- data.frame(suavizado_expos$p)$X1
suavizado_infec <- frfast(solution_2[, 4][1:101] ~ time_1, model = "np", smooth = "kernel", kbin = 101, 
                          p =3)
suavizado_infec <- data.frame(suavizado_infec$p)$X1
suavizado_recup <- frfast(solution_2[, 5][1:101] ~ time_1, model = "np", smooth = "kernel", kbin = 101, 
                          p =3)
suavizado_recup <- data.frame(suavizado_recup$p)$X1

sigma_mv = sqrt((1/(2*100))*(sum((((suavizado_suscep[2:101] - solution_1[, 2][2:101])^2))/((solution_1[, 2][1:100]*
                                                                                            solution_1[, 4][1:100])^2)) + 
           sum((((suavizado_expos[2:101] - solution_1[, 3][2:101])^2))/((solution_1[, 2][1:100]*
                                                                                 solution_1[, 4][1:100])^2))))

SN = list() 
EN = list()
IN = list()
RN = list()
z = list()

for (k in 1:5) {
  SN[[k]] = rep(0, 100)
  EN[[k]] = rep(0, 100)
  IN[[k]] = rep(0, 100)
  RN[[k]] = rep(0, 100)
  z[[k]] = as.vector(rnorm(1000, 0, 1))
}

for (k in 1:5) {
  SN[[k]][1] = solution_fin[1, 2]
  EN[[k]][1] = solution_fin[1, 3]
  IN[[k]][1] = solution_fin[1, 4]
  RN[[k]][1] = solution_fin[1, 5]
}

for (k in 1:5) {
  for(j in 1:99){
    SN[[k]][j+1] = SN[[k]][j] + (parametros_e[1] - parametros_e[2]*SN[[k]][j]*IN[[k]][j] - parametros_e[5]*SN[[k]][j]) - 
                   sigma_mv*SN[[k]][j]*IN[[k]][j]*z[[k]][j]
    EN[[k]][j+1] = EN[[k]][j] + (parametros_e[2]*SN[[k]][j]*IN[[k]][j] - parametros_e[3]*EN[[k]][j] - 
                                parametros_e[5]*EN[[k]][j]) + sigma_mv*SN[[k]]*IN[[k]]*z[[k]][j]
    IN[[k]][j+1] = IN[[k]][j] + (parametros_e[3]*EN[[k]][j] - parametros_e[4]*IN[[k]][j] - parametros_e[5]*IN[[k]][j])
    RN[[k]][j+1] = RN[[k]][j] + (parametros_e[4]*IN[[k]][j] - parametros_e[5]*RN[[k]][j])
  }}

graf_suscp <- data.frame(solution_2[2:101, 1], solution_2[1:100, 2], solution_1[1:100, 2], suavizado_suscep[1:100], 
                         SN[[1]], SN[[2]], SN[[3]], SN[[4]], SN[[5]]) 
graf_expue <- data.frame(solution_2[2:101, 1], solution_2[1:100, 3], solution_1[1:100, 3], suavizado_expos[1:100],
                         EN[[1]], EN[[2]], EN[[3]], EN[[4]], EN[[5]])
graf_infec <- data.frame(solution_2[2:101, 1], solution_2[1:100, 4], solution_1[1:100, 4], suavizado_infec[1:100],
                         IN[[1]], IN[[2]], IN[[3]], IN[[4]], IN[[5]])
graf_recup <- data.frame(solution_2[2:101, 1], solution_2[1:100, 5], solution_1[1:100, 5], suavizado_recup[1:100],
                         RN[[1]], RN[[2]], RN[[3]], RN[[4]], RN[[5]])

library(ggplot2)

colors7 <- c("Susceptible\n data" = 'blue2', "Smoothed\n susceptible\n data" = "darkblue", 
             "Deterministic\n susceptible\n solution" = "dodgerblue", "Stochastic\n paths" = 'cadetblue2')
colors8 <- c("Exposed\n data" = 'chocolate1', "Smoothed\n exposed\n data" = 'orange1',
             "Deterministic\n susceptible\n solution" = 'darkorange4', "Stochastic\n paths" = "sandybrown")
colors9 <- c("Infected\n data" = 'red3', "Smoothed\n infected\n solution" = "red4", 
             "Deterministic\n infected\n solution" = 'firebrick1', "Stochastic\n paths" = 'palevioletred1')
colors10 <- c("Recovered\n data" = 'springgreen4', "Smoothed\n recovered\n solution" = 'darkgreen',
              "Deterministic\n recovered\n solution" = 'green3', "Stochastic\n paths" = 'palegreen3')
library(ggplot2)
p111 <- ggplot(data = graf_suscp)  + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 9], color = "Stochastic\n paths"), size = 0.8) +  
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 4], color = "Smoothed\n susceptible\n data"), size = 1.5) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 3], color = "Deterministic\n susceptible\n solution"), size = 1.5, 
            linetype = "dashed") + 
  geom_point(aes(x = graf_suscp[, 1], y = graf_suscp[, 2], color = "Susceptible\n data"), size = 1.7) + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + 
  scale_color_manual(values = colors7, name = "", guide = guide_legend(override.aes = list(linetype = 
                                                  c("blank", "solid", "dashed", "solid"), shape = c(16, NA, NA, NA)))) + 
  theme(legend.position="bottom")

p222 <- ggplot(data = graf_expue) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 9], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 4], color = "Smoothed\n exposed\n data"), size = 1.5) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 3], color = "Deterministic\n susceptible\n solution"), size = 1.5, 
            linetype = "dashed") +
  geom_point(aes(x = graf_expue[, 1], y = graf_expue[, 2], color = "Exposed\n data"), size = 1.7) + 
  labs(x = 'Time (days)', y = 'No. of exposed population') + 
  scale_color_manual(values = colors8, name = "", guide = guide_legend(override.aes = list(linetype = 
                                               c("blank", "solid", "dashed", "solid"), shape = c(16, NA, NA, NA)))) + 
  theme(legend.position="bottom")

p333 <- ggplot(data = graf_infec)  + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 9], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 4], color = "Smoothed\n infected\n solution"), size = 1.5) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 3], color = "Deterministic\n infected\n solution"), size = 1.5, 
                linetype = "dashed") +
  geom_point(aes(x = graf_infec[, 1], y = graf_infec[, 2], color = "Infected\n data"), size = 1.7)+ 
  labs(x = 'Time (days)', y = 'No. of infected population') + 
  scale_color_manual(values = colors9, name = "", guide = guide_legend(override.aes = list(linetype = 
                                              c("blank", "solid", "dashed", "solid"), shape = c(16, NA, NA, NA)))) +  
  theme(legend.position="bottom")

p444 <- ggplot(data = graf_recup) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 9], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 4], color = "Smoothed\n recovered\n solution"), size = 1.5) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 3], color = "Deterministic\n recovered\n solution"), size = 1.5, 
            linetype = "dashed") + 
  geom_point(aes(x = graf_recup[, 1], y = graf_recup[, 2], color = "Recovered\n data"), size = 1.7)+
  labs(x = 'Time (days)', y = 'No. of recovered population') + 
  scale_color_manual(values = colors10, name = "", guide = guide_legend(override.aes = list(linetype = 
                                               c("blank", "solid", "dashed", "solid"), shape = c(16, NA, NA, NA)))) +  
  theme(legend.position="bottom")

# install.packages("gridExtra")
library(gridExtra)
grid.arrange(p111, p222, p333, p444, ncol = 2)
