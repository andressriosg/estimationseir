# Aplicacion para datos reales

interalo_anuj <- list(1:152, 152:205, 205:265, 265:306, 306:363, 363:385)

cuadrados_ordinarios <- list(NULL, NULL, NULL, NULL, NULL, NULL)
x0_1 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
x0_1[[1]] <- c(proyeccion_bogota_final[1], 0, suavizado_infectados$X1[1], suavizado_recuperados$X1[1])

cuadrados_ordinarios[[1]] <- function(par, data) {
  sum((ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2}

x <- list(NULL, NULL, NULL, NULL, NULL, NULL)
y <- list(NULL, NULL, NULL, NULL, NULL, NULL)
w <- list(NULL, NULL, NULL, NULL, NULL, NULL)

x[[1]] <- seq(5e-06 , 7e-06 , length = 50)
y[[1]] <- seq(5e-05 , 7e-05, length = 50)
w[[1]] <- seq(0.5, 0.55, length = 50)
for (j in 2:length(x)) {
  x[[j]] <- seq(min(beta_no[interalo_anuj[[j]]]), max(beta_no[interalo_anuj[[j]]]), length = 50)
  y[[j]] <- seq(min(upsilon_no[interalo_anuj[[j]]]), max(upsilon_no[interalo_anuj[[j]]]), length = 50)
  w[[j]] <- seq(min(upsilon_no[interalo_anuj[[j]]]), max(upsilon_no[interalo_anuj[[j]]]), length = 50)
}

ord_1 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
ord_2 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
ord_3 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
for (j in 1:length(x)) {
  ord_1[[j]] <- data.frame(rep(x[[j]], rep(length(y[[j]]), length(x[[j]]))), rep(y[[j]], length(x[[j]])))
  ord_2[[j]] <- data.frame(rep(y[[j]], rep(length(w[[j]]), length(y[[j]]))), rep(w[[j]], length(y[[j]])))
  ord_3[[j]] <- data.frame(rep(x[[j]], rep(length(w[[j]]), length(x[[j]]))), rep(w[[j]], length(x[[j]])))
}

dim(ord_1[[1]])
dim(ord_2[[1]])
dim(ord_3[[1]])

z1 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
z2 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
z3 <- list(NULL, NULL, NULL, NULL, NULL, NULL)

library(deSolve)
for (k in 1:dim(ord_1[[1]])[1]) {
  z1[[1]][k] <- cuadrados_ordinarios[[1]](c(ord_1[[1]][k, 1], ord_1[[1]][k, 2], mean(gama[interalo_anuj[[1]]])), 
                                         data = data.frame(suavizado_infectados$X1[interalo_anuj[[1]]])) 
} 

library("plot3D")
par(mfrow=c(1,1))
?scatter3D
scatter3D(ord_1[[1]][, 1], ord_1[[1]][, 2], z1[[1]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "     β", ylab = "υ", zlab = "  SS", nticks = c(2, 2, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z1[[1]] == min(z1[[1]]))
ord_1[[1]][1494, ]

plot(ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, 6.183673e-06, 6.755102e-05, 
                                               mean(gama[interalo_anuj[[1]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[1]]])

parametros_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)
parametros_anuj[[1]] <- c(ord_1[[1]][which(z1[[1]] == min(z1[[1]])), ], mean(gama[interalo_anuj[[1]]]))

optimizacion[[1]] <- optim(par = parametros_anuj[[1]], fn = cuadrados_ordinarios[[1]], data = 
                             data.frame(suavizado_infectados$X1[interalo_anuj[[1]]]), 
                             lower = c(6.18e-06, 6.7e-05, 0.02310177), 
                             upper = c(6.183673e-06, 6.755102e-05, 0.25), method = "L-BFGS-B")

plot(ode(y= x0_1[[1]], intervalo[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[1]]])

r0_1 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
r0_1[[1]] <- (73660.35*optimizacion[[1]]$par[1]*optimizacion[[1]]$par[2])/(0.004*(0.004+optimizacion[[1]]$par[2])*(0.004+optimizacion[[1]]$par[3]))
(0.004*(0.004+optimizacion[[1]]$par[2])*(0.004+optimizacion[[1]]$par[3]))/(73660.35*optimizacion[[1]]$par[2])

## Segundo intervalo 

x0_1[[2]] <- c(ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 2][length(interalo_anuj[[1]])], 
               ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 3][length(interalo_anuj[[1]])], 
               suavizado_infectados$X1[length(interalo_anuj[[1]])], suavizado_recuperados$X1[length(interalo_anuj[[1]])])

cuadrados_ordinarios[[2]] <- function(par, data) {
  sum((ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2}

x[[2]] <- seq(0.0001, 0.2, length = 50)
y[[2]] <- seq(0.0000002, 0.0000004, length = 50)

ord_1[[2]] <- data.frame(rep(x[[2]], rep(length(y[[2]]), length(x[[2]]))), rep(y[[2]], length(x[[2]])))

library(deSolve)
for (k in 1:dim(ord_1[[1]])[1]) {
  z[[2]][k] <- cuadrados_ordinarios[[2]](c(ord_1[[2]][k, 1], ord_1[[2]][k, 2], 0.012), 
                                         data = data.frame(suavizado_infectados$X1[interalo_anuj[[2]]])) 
} 

par(mfrow=c(1,1))
?scatter3D
scatter3D(ord_1[[2]][, 1], ord_1[[2]][, 2], z[[2]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "     β", ylab = "υ", zlab = "          SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[2]] == min(z[[2]]))
ord_1[[2]][1, ]

plot(ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, 0.0001, 0.0000002, 
                                               0.012, 0.004))[, 4]) # 0.012
lines(suavizado_infectados$X1[interalo_anuj[[2]]])

parametros_anuj[[2]] <- c(0.0001, 0.0000002, 0.012)

optimizacion[[2]] <- optim(par = parametros_anuj[[2]], fn = cuadrados_ordinarios[[2]], data = 
                             data.frame(suavizado_infectados$X1[interalo_anuj[[2]]]), 
                           lower = c(1e-8, 1e-7, 0), 
                           upper = c(0.0002, 0.01, 0.04), method = "L-BFGS-B")

plot(ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[interalo_anuj[[2]]])

r0_1[[2]] <- (73660.35*optimizacion[[2]]$par[1]*optimizacion[[2]]$par[2])/(0.004*(0.004+optimizacion[[2]]$par[2])*(0.004+optimizacion[[2]]$par[3]))

## Tercer intervalo 

x0_1[[3]] <- c(ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 2][length(interalo_anuj[[2]])], 
               ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 3][length(interalo_anuj[[2]])], 
               suavizado_infectados$X1[205], suavizado_recuperados$X1[205])

cuadrados_ordinarios[[3]] <- function(par, data) {
  sum((ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2}

x[[3]] <- seq(0.0001, 0.0009, length = 50)
y[[3]] <- seq(1.2e-06, 4e-06, length = 50)

ord_1[[3]] <- data.frame(rep(x[[3]], rep(length(y[[3]]), length(x[[3]]))), rep(y[[3]], length(x[[3]])))

library(deSolve)
for (k in 1:dim(ord_1[[1]])[1]) {
  z[[3]][k] <- cuadrados_ordinarios[[3]](c(ord_1[[3]][k, 1], ord_1[[3]][k, 2], 2.600562e-03), 
                                         data = data.frame(suavizado_infectados$X1[interalo_anuj[[3]]])) 
} 

par(mfrow=c(1,1))
?scatter3D
scatter3D(ord_1[[3]][, 1], ord_1[[3]][, 2], z[[3]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "  SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[3]] == min(z[[3]]))
ord_1[[3]][1, ]

plot(suavizado_infectados$X1[interalo_anuj[[3]]])
lines(ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, 0.0001, 1.5e-6, 2.600562e-03, 0.004))[, 4])

parametros_anuj[[3]] <- c(0.0001, 1.5e-6, 2.600562e-03)
optimizacion[[3]] <- optim(par = parametros_anuj[[3]], fn = cuadrados_ordinarios[[3]], data = 
                             data.frame(suavizado_infectados$X1[interalo_anuj[[3]]]), 
                           lower = c(1e-18,  0.0000015, 0), 
                           upper = c(0.9, 0.01, 1), method = "L-BFGS-B")

plot(suavizado_infectados$X1[interalo_anuj[[3]]])
lines(ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 4])

r0_1[[3]] <- (73660.35*optimizacion[[3]]$par[1]*optimizacion[[3]]$par[2])/(0.004*(0.004+optimizacion[[3]]$par[2])*(0.004+optimizacion[[3]]$par[3]))

## Cuarto intervalo 

x0_1[[4]] <- c(ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 2][length(interalo_anuj[[3]])], 
               ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 3][length(interalo_anuj[[3]])], 
               suavizado_infectados$X1[265], 
               suavizado_recuperados$X1[265])

cuadrados_ordinarios[[4]] <- function(par, data) {
  sum((ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2}

x[[4]] <- seq(1e-7, 9e-7, length = 50)
y[[4]] <- seq(0e-05, 3e-05, length = 50)

ord_1[[4]] <- data.frame(rep(x[[4]], rep(length(y[[4]]), length(x[[4]]))), rep(y[[4]], length(x[[4]])))

library(deSolve)
for (k in 1:dim(ord_1[[1]])[1]) {
  z[[4]][k] <- cuadrados_ordinarios[[4]](c(ord_1[[4]][k, 1], ord_1[[4]][k, 2], mean(gama[interalo_anuj[[4]]])), 
                                         data = data.frame(suavizado_infectados$X1[interalo_anuj[[4]]])) 
} 

scatter3D(ord_1[[4]][, 1], ord_1[[4]][, 2], z[[4]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "  SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[4]] == min(z[[4]]))
ord_1[[4]][2468, ]

plot(suavizado_infectados$X1[interalo_anuj[[4]]])
lines(ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35,  9e-07, 1.040816e-05, mean(gama[interalo_anuj[[4]]]), 0.004))[, 4])

parametros_anuj[[4]] <- c(1e-05, 1.22449e-05, mean(gama[interalo_anuj[[4]]]))
optimizacion[[4]] <- optim(par = parametros_anuj[[4]], fn = cuadrados_ordinarios[[4]], data = 
                             data.frame(suavizado_infectados$X1[interalo_anuj[[4]]]), 
                           lower = c(0,  0.000002, 0), 
                           upper = c(0.9, 0.01, 1), method = "L-BFGS-B")

plot(suavizado_infectados$X1[interalo_anuj[[4]]])
lines(ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 4])

r0_1[[4]] <- (73660.35*optimizacion[[4]]$par[1]*optimizacion[[4]]$par[2])/(0.004*(0.004+optimizacion[[4]]$par[2])*(0.004+optimizacion[[4]]$par[3]))

## Quinto intervalo 

x0_1[[5]] <- c(ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 2][length(interalo_anuj[[4]])], 
               ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 3][length(interalo_anuj[[4]])], 
               suavizado_infectados$X1[306], 
               suavizado_recuperados$X1[306])

cuadrados_ordinarios[[5]] <- function(par, data) {
  sum((ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2}

x[[5]] <- seq(0.00001, 9e-7, length = 50)
y[[5]] <- seq(0e-05, 1.3e-05, length = 50)

ord_1[[5]] <- data.frame(rep(x[[5]], rep(length(y[[5]]), length(x[[5]]))), rep(y[[5]], length(x[[5]])))

library(deSolve)
for (k in 1:dim(ord_1[[1]])[1]) {
  z[[5]][k] <- cuadrados_ordinarios[[5]](c(ord_1[[5]][k, 1], ord_1[[5]][k, 2], 0.029), 
                                         data = data.frame(suavizado_infectados$X1[interalo_anuj[[5]]])) 
} 

scatter3D(ord_1[[5]][, 1], ord_1[[5]][, 2], z[[5]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "     β", ylab = "υ", zlab = "  SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[5]] == min(z[[5]]))
ord_1[[5]][3, ]
1e-05 == 0.00001

plot(suavizado_infectados$X1[interalo_anuj[[5]]])
lines(ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35,  0.00001, 5.306122e-07, 0.029, 0.004))[, 4])

parametros_anuj[[5]] <- c(0.00001, 5.306122e-07, 0.029)
optimizacion[[5]] <- optim(par = parametros_anuj[[5]], fn = cuadrados_ordinarios[[5]], data = 
                             data.frame(suavizado_infectados$X1[interalo_anuj[[5]]]), 
                           lower = c(0,  0.0000002, 0), 
                           upper = c(0.9, 0.01, 1), method = "L-BFGS-B")

plot(suavizado_infectados$X1[interalo_anuj[[5]]])
lines(ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 4])

r0_1[[5]] <- (73660.35*optimizacion[[5]]$par[1]*optimizacion[[5]]$par[2])/(0.004*(0.004+optimizacion[[5]]$par[2])*(0.004+optimizacion[[5]]$par[3]))

## Sexto intervalo 

x0_1[[6]] <- c(ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 2][length(interalo_anuj[[5]])], 
               ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 3][length(interalo_anuj[[5]])], 
               suavizado_infectados$X1[363], 
               suavizado_recuperados$X1[363])

cuadrados_ordinarios[[6]] <- function(par, data) {
  sum((ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2}

x[[6]] <- seq(0.0001, 0.0002, length = 50)
y[[6]] <- seq(1.8e-06, 3e-06, length = 50)

ord_1[[6]] <- data.frame(rep(x[[6]], rep(length(y[[6]]), length(x[[6]]))), rep(y[[6]], length(x[[6]])))

library(deSolve)
for (k in 1:dim(ord_1[[1]])[1]) {
  z[[6]][k] <- cuadrados_ordinarios[[6]](c(ord_1[[6]][k, 1], ord_1[[6]][k, 2], mean(gama[interalo_anuj[[6]]])), 
                                         data = data.frame(suavizado_infectados$X1[interalo_anuj[[6]]])) 
} 

scatter3D(ord_1[[6]][, 1], ord_1[[6]][, 2], z[[6]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "     β", ylab = "υ", zlab = "          SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[6]] == min(z[[6]]))
ord_1[[6]][2072, ]

plot(suavizado_infectados$X1[interalo_anuj[[6]]])
lines(ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, 0.000183673, 2.314286e-06, mean(gama[interalo_anuj[[6]]]), 0.004))[, 4])

parametros_anuj[[6]] <- c(0.000183673, 2.314286e-06, mean(gama[interalo_anuj[[6]]]))
optimizacion[[6]] <- optim(par = parametros_anuj[[6]], fn = cuadrados_ordinarios[[6]], data = 
                             data.frame(suavizado_infectados$X1[interalo_anuj[[6]]]), 
                           lower = c(0,  0.000002, 0), 
                           upper = c(0.9, 0.01, 1), method = "L-BFGS-B")

plot(suavizado_infectados$X1[interalo_anuj[[6]]])
lines(ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, optimizacion[[6]]$par, 0.004))[, 4])

r0_1[[6]] <- (73660.35*optimizacion[[6]]$par[1]*optimizacion[[6]]$par[2])/(0.004*(0.004+optimizacion[[6]]$par[2])*(0.004+optimizacion[[6]]$par[3]))

(0.004*(0.004+optimizacion[[6]]$par[2])*(0.004+optimizacion[[6]]$par[3]))/(73660.35*optimizacion[[6]]$par[2])
r0_1

plot(interalo_anuj[[1]], ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 2], xlim = c(0, 385), 
     ylim = c(0, max(ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, optimizacion[[6]]$par, 0.004))[, 3])))
lines(interalo_anuj[[2]], ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 2], type = "p")
lines(interalo_anuj[[3]], ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 2], type = "p")
lines(interalo_anuj[[4]], ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 2], type = "p")
lines(interalo_anuj[[5]], ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 2], type = "p")
lines(interalo_anuj[[6]], ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, optimizacion[[6]]$par, 0.004))[, 2], type = "p")
lines(1:153, euler_anuj_suscep[[1]], col = "blue")
lines(152:206, euler_anuj_suscep[[2]], col = "blue")
lines(205:266, euler_anuj_suscep[[3]], col = "blue")
lines(265:307, euler_anuj_suscep[[4]], col = "blue")
lines(306:364, euler_anuj_suscep[[5]], col = "blue")
lines(363:386, euler_anuj_suscep[[6]], col = "blue")

plot(interalo_anuj[[1]], ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 3], xlim = c(0, 385), 
     ylim = c(0, max(euler_anuj_expue[[6]])))
lines(interalo_anuj[[2]], ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 3], type = "p")
lines(interalo_anuj[[3]], ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 3], type = "p")
lines(interalo_anuj[[4]], ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 3], type = "p")
lines(interalo_anuj[[5]], ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 3], type = "p")
lines(interalo_anuj[[6]], ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, optimizacion[[6]]$par, 0.004))[, 3], type = "p")
lines(1:153, euler_anuj_expue[[1]], col = "blue")
lines(152:206, euler_anuj_expue[[2]], col = "blue")
lines(205:266, euler_anuj_expue[[3]], col = "blue")
lines(265:307, euler_anuj_expue[[4]], col = "blue")
lines(306:364, euler_anuj_expue[[5]], col = "blue")
lines(363:386, euler_anuj_expue[[6]], col = "blue")

euler_anuj_suscep <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(euler_anuj_suscep)) { # Los objetos euler_anuj_ son las soluciones tomando I(t) como los datos 
  euler_anuj_suscep[[k]][1] <- ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 2]
  euler_anuj_suscep[[k]][2:(length(interalo_anuj[[k]])+1)] <- 
    ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1:length(interalo_anuj[[k]]), 2] - 
    (73660.35 - optimizacion[[k]]$par[1]*ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par,0.004))[1:length(interalo_anuj[[k]]), 2]*
    suavizado_infectados$X1[interalo_anuj[[k]]] - 0.004*ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 
                                                                                                      0.004))[1:length(interalo_anuj[[k]]), 2]) 
}

euler_anuj_expue <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(euler_anuj_expue)) {
  euler_anuj_expue[[k]][1] <- ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 3]
  euler_anuj_expue[[k]][2:(length(interalo_anuj[[k]])+1)] <- 
    ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1:length(interalo_anuj[[k]]), 3] + 
    (optimizacion[[k]]$par[1]*ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par,0.004))[1:length(interalo_anuj[[k]]), 2]*
       suavizado_infectados$X1[interalo_anuj[[k]]] - (optimizacion[[k]]$par[2] + 0.004)*ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, 
                                                       optimizacion[[k]]$par, 0.004))[1:length(interalo_anuj[[k]]), 3]) 
}

euler_anuj_infec <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(euler_anuj_expue)) {
  euler_anuj_infec[[k]][1] <- suavizado_infectados$X1[interalo_anuj[[k]]][1]
  euler_anuj_infec[[k]][2:(length(interalo_anuj[[k]])+1)] <- 
    suavizado_infectados$X1[interalo_anuj[[k]]] + (optimizacion[[k]]$par[2]*ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, 
    optimizacion[[k]]$par, 0.004))[1:length(interalo_anuj[[k]]), 3] - (0.004 + optimizacion[[k]]$par[2])*suavizado_infectados$X1[interalo_anuj[[k]]])
}

euler_euler_suscep <- list(NULL, NULL, NULL, NULL, NULL, NULL) 
euler_euler_expue <- list(NULL, NULL, NULL, NULL, NULL, NULL) 
euler_euler_infec <- list(NULL, NULL, NULL, NULL, NULL, NULL)
euler_euler_recup <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:6) {
  euler_euler_suscep[[k]][1] <- ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 2]
  euler_euler_expue[[k]][1] <- ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 3]
  euler_euler_infec[[k]][1] <- ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 4]
  euler_euler_recup[[k]][1] <- ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 5]
}

for (k in 1:6) { # mu_s, mu_e Soluciones generadas del modelo y NO USANDO DATOS 
  for (j in 1:length(interalo_anuj[[k]])) {
    euler_euler_suscep[[k]][j+1] <- euler_euler_suscep[[k]][j] + 73660.35 - optimizacion[[k]]$par[1]*euler_euler_suscep[[k]][j]*euler_euler_infec[[k]][j] - 
      0.004*euler_euler_suscep[[k]][j]
    euler_euler_expue[[k]][j+1] <- euler_euler_expue[[k]][j] + optimizacion[[k]]$par[1]*euler_euler_suscep[[k]][j]*euler_euler_infec[[k]][j] - (
      optimizacion[[k]]$par[2] + 0.004)*euler_euler_expue[[k]][j]
    euler_euler_infec[[k]][j+1] <- euler_euler_infec[[k]][j] + optimizacion[[k]]$par[2]*euler_euler_expue[[k]][j] - (optimizacion[[k]]$par[3] + 
                                                                                                                       0.004)*euler_euler_infec[[k]][j]
    euler_euler_recup[[k]][j+1] <- euler_euler_recup[[k]][j] + optimizacion[[k]]$par[3]*euler_euler_infec[[k]][j] - 0.004*euler_euler_recup[[k]][j] 
  }
}

plot(suavizado_infectados$X1[1:385])
lines(1:153, euler_anuj_infec[[1]], col = "blue")
lines(152:206, euler_anuj_infec[[2]], col = "blue")
lines(205:266, euler_anuj_infec[[3]], col = "blue")
lines(265:307, euler_anuj_infec[[4]], col = "blue")
lines(306:364, euler_anuj_infec[[5]], col = "blue")
lines(363:386, euler_anuj_infec[[6]], col = "blue")

sigma_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(sigma_anuj )) {
sigma_anuj[[k]] <- (1/(2*length(euler_anuj_suscep[[k]])))*(sum(((ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, 
                     optimizacion[[k]]$par, 0.004))[2:length(interalo_anuj[[k]]), 2] - 
                     euler_anuj_suscep[[k]][2:length(interalo_anuj[[k]])])^2)/((euler_anuj_suscep[[k]][2:length(interalo_anuj[[k]])])^2*
                     (euler_anuj_infec[[k]][2:length(interalo_anuj[[k]])])^2)) + 
                    sum(((ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[2:length(interalo_anuj[[k]]), 3] - 
                    euler_anuj_expue[[k]][2:length(interalo_anuj[[k]])])^2)/((euler_anuj_suscep[[k]][2:length(interalo_anuj[[k]])])^2*
                    (euler_anuj_infec[[k]][2:length(interalo_anuj[[k]])])^2)))
}

for (k in 1:length(sigma_anuj )) {
  sigma_anuj[[k]] <- (1/(2*(length(euler_anuj_suscep[[k]])-1)))*(sum(((euler_euler_suscep[[k]][2:length(interalo_anuj[[k]])] - 
                      euler_anuj_suscep[[k]][2:length(interalo_anuj[[k]])])^2)/(((euler_euler_suscep[[k]][1:(length(interalo_anuj[[k]])-1)])^2)*
                      ((euler_euler_infec[[k]][1:(length(interalo_anuj[[k]])-1)])^2))) + 
                      sum(((euler_euler_expue[[k]][2:length(interalo_anuj[[k]])] - 
                      euler_anuj_expue[[k]][2:length(interalo_anuj[[k]])])^2)/(((euler_euler_suscep[[k]][1:(length(interalo_anuj[[k]])-1)])^2)*
                      ((euler_euler_infec[[k]][1:(length(interalo_anuj[[k]])-1)])^2))))
}

extincion <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(extincion)) {
  extincion[[k]] <- (optimizacion[[k]]$par[1])^2/(2*((sigma_anuj[[k]])^2)) - (optimizacion[[k]]$par[2] + 0.004)
}

r0_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(extincion)) {
  r0_anuj[[k]] <- r0_1[[k]] - ((73660.35*sigma_anuj[[k]]*optimizacion[[k]]$par[2])^2)/(2*((0.004^2)*(0.004 + optimizacion[[k]]$par[2])*
                                                                                              (0.004 + optimizacion[[k]]$par[3])^2))
}

r0_anuj_1 <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(extincion)) {
  r0_anuj_1[[k]] <- r0_1[[k]] - ((73660.35*sigma_anuj[[k]]*optimizacion[[k]]$par[2])^2)/(2*(0.004^2)*(0.004 + optimizacion[[k]]$par[2])*
                                                                                                   (0.004 + optimizacion[[k]]$par[3]))
}

sd_varianza <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(extincion)) {
  sd_varianza[[k]] <- sqrt(((73660.35*sigma_anuj[[k]]*optimizacion[[k]]$par[2])^2)/((2*(0.004^2)*((0.004 + optimizacion[[k]]$par[2])^3)*
                                                                                                   (0.004 + optimizacion[[k]]$par[3])^3)))
}

SN5 <- list(matrix(0, nrow = (length(interalo_anuj[[1]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[2]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[3]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[4]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[5]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[6]]) + 1), ncol = 100))
EN5 <- list(matrix(0, nrow = (length(interalo_anuj[[1]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[2]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[3]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[4]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[5]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[6]]) + 1), ncol = 100))
IN5 <- list(matrix(0, nrow = (length(interalo_anuj[[1]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[2]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[3]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[4]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[5]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[6]]) + 1), ncol = 100))
RN5 <- list(matrix(0, nrow = (length(interalo_anuj[[1]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[2]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[3]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[4]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[5]]) + 1), ncol = 100), 
            matrix(0, nrow = (length(interalo_anuj[[6]]) + 1), ncol = 100))
z5 <- list(matrix(rnorm((length(interalo_anuj[[1]]) + 1)*100, 0, 1), ncol = 100), 
           matrix(rnorm((length(interalo_anuj[[2]]) + 1)*100, 0, 1), ncol = 100),
           matrix(rnorm((length(interalo_anuj[[3]]) + 1)*100, 0, 1), ncol = 100),
           matrix(rnorm((length(interalo_anuj[[4]]) + 1)*100, 0, 1), ncol = 100),
           matrix(rnorm((length(interalo_anuj[[5]]) + 1)*100, 0, 1), ncol = 100),
           matrix(rnorm((length(interalo_anuj[[6]]) + 1)*100, 0, 1), ncol = 100))

library(deSolve)
for (k in 1:6) {
  SN5[[k]][1, 1:100] = ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 2]
  EN5[[k]][1, 1:100] = ode(y= x0_1[[k]], interalo_anuj[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[1, 3]
  IN5[[k]][1, 1:100] = suavizado_infectados$X1[interalo_anuj[[k]]][1]
  RN5[[k]][1, 1:100] = suavizado_recuperados$X1[interalo_anuj[[k]]][1]
}

for (k in 1:6) {
  for (j in 2:(length(interalo_anuj[[k]]) + 1)) {
    SN5[[k]][j, 1:100] = SN5[[k]][j-1, 1:10] + (73660.35 - optimizacion[[k]]$par[1]*SN5[[k]][j-1, 1:10]*IN5[[k]][j-1, 1:10] - 
                        0.004*SN5[[k]][j-1, 1:10]) - sigma_anuj[[k]]*z5[[k]][j-1, 1:10]*SN5[[k]][j-1, 1:10]*IN5[[k]][j-1, 1:10]
    EN5[[k]][j, 1:100] = EN5[[k]][j-1, 1:10] + (optimizacion[[k]]$par[1]*SN5[[k]][j-1, 1:10]*IN5[[k]][j-1, 1:10] - 
                  (optimizacion[[k]]$par[2] + 0.004)*EN5[[k]][j-1, 1:10]) + 
                   sigma_anuj[[k]]*z5[[k]][j-1, 1:10]*SN5[[k]][j-1, 1:10]*IN5[[k]][j-1, 1:10]
    IN5[[k]][j, 1:100] = IN5[[k]][j-1, 1:10] + (optimizacion[[k]]$par[2]*EN5[[k]][j-1, 1:10] - 
                        (optimizacion[[k]]$par[3] + 0.004)*IN5[[k]][j-1, 1:10])
    RN5[[k]][j, 1:100] = RN5[[k]][j-1, 1:10] + (optimizacion[[k]]$par[3]*IN5[[k]][j-1, 1:10] - 0.004*RN5[[k]][j-1, 1:10])
  }
}

tiempo <- seq(from = as.Date("2020-03-06"), to = as.Date("2021-03-29"), by=1)
time = c(1:152, 152:205, 205:265, 265:306, 306:363, 363:385)
susceptibles_esto_anuj <- matrix(0, ncol = 100, nrow = 6*length(time))
expuestos_esto_anuj <- matrix(0, ncol = 100, nrow = 6*length(time))
infectados_esto_anuj <- matrix(0, ncol = 100, nrow = 6*length(time))

susceptibles_esto_anuj <- rbind(SN5[[1]][1:152, ], SN5[[2]][1:length(152:205), ], SN5[[3]][1:length(205:265), ], 
                          SN5[[4]][1:length(265:306), ], SN5[[5]][1:length(306:363), ], SN5[[6]][1:length(363:385), ])
expuestos_esto_anuj <- rbind(EN5[[1]][1:152, ], EN5[[2]][1:length(152:205), ], EN5[[3]][1:length(205:265), ], 
                             EN5[[4]][1:length(265:306), ], EN5[[5]][1:length(306:363), ], EN5[[6]][1:length(363:385), ])
infectados_esto_anuj <- rbind(IN5[[1]][1:152, ], IN5[[2]][1:length(152:205), ], IN5[[3]][1:length(205:265), ], 
                              IN5[[4]][1:length(265:306), ], IN5[[5]][1:length(306:363), ], IN5[[6]][1:length(363:385), ])

a <- data.frame(c(tiempo[1:152], tiempo[152:205], tiempo[205:265], tiempo[265:306], tiempo[306:363], tiempo[363:385]), 
                c(euler_anuj_suscep[[1]][1:152], euler_anuj_suscep[[2]][1:length(152:205)], 
                  euler_anuj_suscep[[3]][1:length(205:265)], euler_anuj_suscep[[4]][1:length(265:306)], 
                  euler_anuj_suscep[[5]][1:length(306:363)], euler_anuj_suscep[[6]][1:length(363:385)]), 
                c(euler_euler_suscep[[1]][1:152], euler_euler_suscep[[2]][1:length(152:205)], 
                  euler_euler_suscep[[3]][1:length(205:265)], euler_euler_suscep[[4]][1:length(265:306)],
                  euler_euler_suscep[[5]][1:length(306:363)], euler_euler_suscep[[6]][1:length(363:385)]), 
                susceptibles_esto_anuj)
b <- data.frame(c(tiempo[1:152], tiempo[152:205], tiempo[205:265], tiempo[265:306], tiempo[306:363], tiempo[363:385]), 
                c(euler_anuj_expue[[1]][1:152], euler_anuj_expue[[2]][1:length(152:205)], 
                  euler_anuj_expue[[3]][1:length(205:265)], euler_anuj_expue[[4]][1:length(265:306)], 
                  euler_anuj_expue[[5]][1:length(306:363)], euler_anuj_expue[[6]][1:length(363:385)]), 
                c(euler_euler_expue[[1]][1:152], euler_euler_expue[[2]][1:length(152:205)], 
                  euler_euler_expue[[3]][1:length(205:265)], euler_euler_expue[[4]][1:length(265:306)],
                  euler_euler_expue[[5]][1:length(306:363)], euler_euler_expue[[6]][1:length(363:385)]), 
                  expuestos_esto_anuj)
c <- data.frame(c(tiempo[1:152], tiempo[152:205], tiempo[205:265], tiempo[265:306], tiempo[306:363], tiempo[363:385]), 
                c(casos$Infectados[1:152], casos$Infectados[152:205], 
                  casos$Infectados[205:265], casos$Infectados[265:306], 
                  casos$Infectados[306:363], casos$Infectados[363:385]), 
                c(suavizado_infectados$X1[1:152], suavizado_infectados$X1[152:205], 
                  suavizado_infectados$X1[205:265], suavizado_infectados$X1[265:306], 
                  suavizado_infectados$X1[306:363], suavizado_infectados$X1[363:385]), 
                c(ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 4][1:152], 
                  ode(y= x0_1[[2]], interalo_anuj[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 4][1:length(152:205)], 
                  ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 4][1:length(205:265)], 
                  ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 4][1:length(265:306)], 
                  ode(y= x0_1[[5]], interalo_anuj[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 4][1:length(306:363)], 
                  ode(y= x0_1[[6]], interalo_anuj[[6]], seir, c(73660.35, optimizacion[[6]]$par, 0.004))[, 4][1:length(363:385)]), 
                c(euler_anuj_infec[[1]][1:152], euler_anuj_infec[[2]][1:length(152:205)], 
                  euler_anuj_infec[[3]][1:length(205:265)], euler_anuj_infec[[4]][1:length(265:306)], 
                  euler_anuj_infec[[5]][1:length(306:363)], euler_anuj_infec[[6]][1:length(363:385)]), 
                c(euler_euler_infec[[1]][1:152], euler_euler_infec[[2]][1:length(152:205)], 
                  euler_euler_infec[[3]][1:length(205:265)], euler_euler_infec[[4]][1:length(265:306)],
                  euler_euler_infec[[5]][1:length(306:363)], euler_euler_infec[[6]][1:length(363:385)]), 
                  infectados_esto_anuj)
View(c)

colorsa <- c('Stochastic \npaths' = "gray50", 'Estimated \nsusceptible' = "royalblue1", 
             'Euler \nsusceptible \napproximation' = "midnightblue")
colorsb <- c('Stochastic \npaths' = "gray50", 'Estimated \nexposed' = "lightsalmon1", 
             'Euler \nexposed \napproximation' = "chocolate3")
colorsc <- c('Stochastic \npaths' = "gray50", 'Estimated \ninfected' = "maroon1", 
             'Euler \ninfected \napproximation' = "purple1", "Solution \nby using ode" = "mediumvioletred", 
             "Smoothed \ninfected \ndata"= "red4", "Infected \ndata" = "palevioletred3")

for(i in 5:ncol(a)) { 
p_a <- print(ggplot(data = a) + 
  geom_line(aes(x = a[, 1],  y = a[, i], color = 'Stochastic \npaths'), size = 3.3) + 
  geom_point(aes(x = a[, 1], y = a[, 3], color = 'Euler \nsusceptible \napproximation'), size = 1.2) +
  geom_point(aes(x = a[, 1],  y = a[, 2], color = 'Estimated \nsusceptible'), size = 0.7) + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colorsa, name = "", 
                                                                                      guide = guide_legend(override.aes = list(
                                                                                        linetype = c("solid", "blank", "blank"),
                                                                                        shape = c(NA, 16, 16), size = c(1,1,1)))) + 
    geom_vline(xintercept = tiempo[152]) +
    geom_vline(xintercept = tiempo[205]) +
    geom_vline(xintercept = tiempo[265]) +
    geom_vline(xintercept = tiempo[306]) +
    geom_vline(xintercept = tiempo[363]) +
  theme(legend.position="right") + scale_x_date(breaks = "3 months", date_labels = "%m/%Y"))
}

for(i in 5:ncol(b)) { 
p_b <- print(ggplot(data = b) + 
               geom_line(aes(x = b[, 1],  y = b[, i], color = 'Stochastic \npaths'), size = 3.3) + 
               geom_point(aes(x = b[, 1], y = b[, 3], color = 'Euler \nexposed \napproximation'), size = 1.2) +
               geom_point(aes(x = b[, 1],  y = b[, 2], color = 'Estimated \nexposed'), size = 0.7) + 
               labs(x = 'Time (days)', y = 'No. of exposed population') + scale_color_manual(values = colorsb, name = "", 
                                                                              guide = guide_legend(override.aes = list(
                                                                              linetype = c("solid", "blank", "blank"),
                                                                              shape = c(NA, 16, 16), size = c(1,1,1)))) + 
               geom_vline(xintercept = tiempo[152]) +
               geom_vline(xintercept = tiempo[205]) +
               geom_vline(xintercept = tiempo[265]) +
               geom_vline(xintercept = tiempo[306]) +
               geom_vline(xintercept = tiempo[363]) +
               theme(legend.position="right") + scale_x_date(breaks = "3 months", date_labels = "%m/%Y"))
}
library(gridExtra)
grid.arrange(p_a, p_b, nrow = 2)

library(ggplot2)
for(i in 7:ncol(c)) { 
  p_c <- print(ggplot(data = c) + 
                 geom_point(aes(x = c[, 1],  y = c[, 2], color = "Infected \ndata"), size = 1) +
                 geom_point(aes(x = c[, 1],  y = c[, i], color = 'Stochastic \npaths'), size = 2.5) + 
                 geom_point(aes(x = c[, 1], y = c[, 6], color = 'Euler \ninfected \napproximation'), size = 1.2) +
                 geom_line(aes(x = c[, 1],  y = c[, 4], color = "Solution \nby using ode"), size = 1.2) + 
                 geom_line(aes(x = c[, 1],  y = c[, 3], color = "Smoothed \ninfected \ndata"), size = 1.5) + 
                 geom_point(aes(x = c[, 1],  y = c[, 5], color = 'Estimated \ninfected'), size = 0.6) +
                 labs(x = 'Time (days)', y = 'No. of infected population') + scale_color_manual(values = colorsc, name = "", 
                                                                             guide = guide_legend(override.aes = list(
                                                                             linetype = c("solid", "blank", "blank", 
                                                                                          "solid", "solid", "blank"),
                                                                             shape = c(NA, 16, 16, NA, NA, 16), 
                                                                             size = c(1, 1, 1, 1, 1, 1)))) + 
                 geom_vline(xintercept = tiempo[152]) +
                 geom_vline(xintercept = tiempo[205]) +
                 geom_vline(xintercept = tiempo[265]) +
                 geom_vline(xintercept = tiempo[306]) +
                 geom_vline(xintercept = tiempo[363]) +
                 theme(legend.position="right") + scale_x_date(breaks = "3 months", date_labels = "%m/%Y"))
}
p_c
warnings()

#############################################
anuj_sus <- data.frame(tiempo[1:385], casos$Infectados[1:385], euler_anuj_suscep[1:385], SN5[[k]])

plot(c(interalo_anuj[[1]], interalo_anuj[[1]][length(interalo_anuj[[1]])] + 1) , EN5[[1]][ , 4],  
     type = "l")

plot(c(interalo_anuj[[1]], interalo_anuj[[1]][length(interalo_anuj[[1]])] + 1) , euler_anuj_suscep[[1]], xlim = c(0, 385), 
     ylim = c(min(euler_anuj_suscep[[1]]), max(euler_anuj_suscep[[1]])), type = "l")
for (k in 2:6) {
  lines(c(interalo_anuj[[k]], interalo_anuj[[k]][length(interalo_anuj[[k]])] + 1), euler_anuj_suscep[[k]], type = "l")
}
for (s in 1:6) {
  for (j in 1:100) {
    lines(c(interalo_anuj[[s]], interalo_anuj[[s]][length(interalo_anuj[[s]])] + 1), SN5[[s]][, j], col = 2:11)
  }  
}

for (i in 1:5) { 
  print(ggplot(df,aes(x,y))+geom_point()) 
}

plot(c(interalo_anuj[[1]], interalo_anuj[[1]][length(interalo_anuj[[1]])] + 1) , euler_anuj_expue[[1]], xlim = c(0, 385), 
     ylim = c(0, max(euler_anuj_expue[[6]])), type = "p", pch = 1)
for (k in 2:6) {
  lines(c(interalo_anuj[[k]], interalo_anuj[[k]][length(interalo_anuj[[k]])] + 1), euler_anuj_expue[[k]], type = "p", pch = 1)
}
for (s in 1:6) {
  for (j in 1:100) {
    lines(c(interalo_anuj[[s]], interalo_anuj[[s]][length(interalo_anuj[[s]])] + 1), EN5[[s]][ , j], col = "blue", lwd = 2)
  }  
}

## Generar los valores de mu_S a partir de I(t) generado del modelo (Euler) y S(tj) estimado de I(tj) (los datos)

plot(suavizado_infectados$X1[interalo_anuj[[4]]])
lines(ode(y= x0_1[[4]], interalo_anuj[[4]], seir, c(73660.35, ord_1[[4]][which(z[[4]] == min(z[[4]])), ], 2.600562e-03, 0.004))[, 4])

parametros_anuj[[3]] <- c(0.0001, 1.2e-6, 2.600562e-03)
optimizacion[[3]] <- optim(par = parametros_anuj[[3]], fn = cuadrados_ordinarios[[3]], data = 
                             data.frame(suavizado_infectados$X1[intervalo[[3]]]), 
                           lower = c(1e-18,  0.000002, 0), 
                           upper = c(0.9, 0.01, 1), method = "L-BFGS-B")

plot(suavizado_infectados$X1[interalo_anuj[[3]]])
lines(ode(y= x0_1[[3]], interalo_anuj[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 4])

for (j in 2:length(z)) {
  for (k in 1:dim(ord_1[[j]])[1]) {
    z[[j]][k] <- cuadrados_ordinarios[[j]](c(ord_1[[j]][k, 1], ord_1[[j]][k, 2], mean(gamma_no_para[intervalo[[j]]])), 
                                           data = data.frame(suavizado_infectados$X1[intervalo[[j]]])) 
  } 
}

for (k in 1:dim(ord_1[[j]])[1]) {
  z[[1]][k] <- cuadrados_ordinarios[[j]](c(ord_1[[j]][k, 1], ord_1[[j]][k, 2], mean(gamma_no_para[intervalo[[j]]])), 
                                         data = data.frame(suavizado_infectados$X1[intervalo[[j]]])) 
} 

library(deSolve)
for (j in 1:length(intervalo)) {
  cuadrados_ordinarios[[j]] <- function(par, data) {
    with(data, 
         sum((ode(y= x0[[j]], intervalo[[j]], seir, 
                  c(73660.35, par[1], par[2], par[3], 0.004))[, 4] - data))^2)
  }}

intervalo <- list(1:152, 152:205, 205:265, 265:306, 306:363, 363:385)
for (j in 2:length(intervalo)) {
  
  x0[[1]] <- c(susceptibles_bogota[1], expuestos_bogota[1], suavizado_infectados$X1[1], 
               suavizado_recuperados$X1[1])
  
  x0[[j]] <- c(susceptibles_bogota[intervalo[[j-1]][length(intervalo[[j-1]])]], 
               expuestos_bogota[intervalo[[j-1]][length(intervalo[[j-1]])]], 
               suavizado_infectados$X1[intervalo[[j-1]][length(intervalo[[j-1]])]], 
               suavizado_recuperados$X1[intervalo[[j-1]][length(intervalo[[j-1]])]])
}


library(deSolve)
for (j in 1:length(z)) {
  for (k in 1:dim(ord_1[[j]])[1]) {
    z[[j]][k] <- cuadrados_ordinarios[[j]](c(ord_1[[j]][k, 1], ord_1[[j]][k, 2], mean(gamma_no_para[intervalo[[j]]])), 
                                      data = data.frame(suavizado_infectados$X1[intervalo[[j]]])) 
  } 
}

library("plot3D")
par(mfrow=c(1,1))
?scatter3D
scatter3D(ord_1[[1]][, 1], ord_1[[1]][, 2], z[[1]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "  SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[1]] == min(z[[1]]))
ord_1[[1]][32, ]

plot(ode(y= x0[[1]], intervalo[[1]], seir, c(73660.35, 1.7e-06, 5.701418e-05, 
                                             mean(gama[intervalo[[1]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[1]]])

library(carData)
library(car)
scatter3D(ord_1[[2]][, 1], ord_1[[2]][, 2], z[[2]], colkey = FALSE,
          phi = 20, theta = 60, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "        SS", nticks = c(2, 2, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[2]] == min(z[[2]]))
ord_1[[2]][413, ]

plot(ode(y= x0[[2]], intervalo[[2]], seir, c(73660.35, 3.524548e-08, -0.000238623, 
                                             mean(gama[intervalo[[2]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[2]]])

scatter3D(ord_1[[3]][, 1], ord_1[[3]][, 2], z[[3]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[3]] == min(z[[3]]))
ord_1[[3]][18, ]

plot(ode(y= x0[[3]], intervalo[[3]], seir, c(73660.35, 5.806371e-08, 6.636308e-05, 
                                             mean(gama[intervalo[[3]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[3]]])

scatter3D(ord_1[[4]][, 1], ord_1[[4]][, 2], z[[4]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[4]] == min(z[[4]]))
ord_1[[4]][2223, ]

plot(ode(y= x0[[4]], intervalo[[4]], seir, c(73660.35,  5.45072e-08, 0.0006364755, 
                                             mean(gama[intervalo[[4]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[4]]])

scatter3D(ord_1[[5]][, 1], ord_1[[5]][, 2], z[[5]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "         SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[5]] == min(z[[5]]))
ord_1[[5]][13, ]

plot(ode(y= x0[[5]], intervalo[[5]], seir, c(73660.35, 2.108976e-07, -0.0007008635, 
                                             mean(gama[intervalo[[5]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[5]]])

scatter3D(ord_1[[6]][, 1], ord_1[[6]][, 2], z[[6]], colkey = FALSE,
          phi = 20, theta = 50, type = "p", ticktype = "detailed", bty = "f", pch = 20, 
          xlab = "β", ylab = "υ", zlab = "SS", nticks = c(3, 3, 2), expand = 0.9, 
          scale = 0.9, cex.axis = 0.7)
which(z[[6]] == min(z[[6]]))
ord_1[[6]][13, ]

plot(ode(y= x0[[6]], intervalo[[6]], seir, c(73660.35, 8.50562e-08, 0.0001660646, 
                                             mean(gama[intervalo[[6]]]), 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[6]]])

parametros_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)
parametros_anuj[[1]] <- c(1.7e-06, 5.701418e-05, mean(gama[intervalo[[1]]]))
parametros_anuj[[2]] <- c(3.524548e-08, -0.000238623, mean(gama[intervalo[[2]]])) 
parametros_anuj[[3]] <- c(5.806371e-08, 6.636308e-05, mean(gama[intervalo[[3]]]))
parametros_anuj[[4]] <- c(5.45072e-08, 0.0006364755, mean(gama[intervalo[[4]]]))
parametros_anuj[[5]] <- c(2.108976e-07, -0.0007008635, mean(gama[intervalo[[5]]]))
parametros_anuj[[6]] <- c(8.50562e-08, 0.0001660646, mean(gama[intervalo[[6]]]))

optimizacion <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (j in 1:length(intervalo)) {
  optimizacion[[j]] <- optim(par = parametros_anuj[[j]], fn = cuadrados_ordinarios[[j]], data = 
                               data.frame(suavizado_infectados$X1[intervalo[[j]]]), lower = 1e-16, upper = 1, 
                             method = "L-BFGS-B")
}

x0_1 <- list(NULL, NULL, NULL, NULL, NULL, NULL)
x0_1[[1]] <- c(proyeccion_bogota_final[1], 0, suavizado_infectados$X1[1], suavizado_recuperados$X1[1])

plot(ode(y= x0_1[[1]], intervalo[[1]], seir, c(73660.35, optimizacion[[1]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[1]]])

plot(ode(y= x0[[2]], intervalo[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[2]]])

plot(ode(y= x0[[3]], intervalo[[3]], seir, c(73660.35, optimizacion[[3]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[3]]])

plot(ode(y= x0[[4]], intervalo[[4]], seir, c(73660.35, optimizacion[[4]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[4]]])

plot(ode(y= x0[[5]], intervalo[[5]], seir, c(73660.35, optimizacion[[5]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[5]]])

plot(ode(y= x0[[6]], intervalo[[6]], seir, c(73660.35, optimizacion[[6]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[6]]])

plot(ode(y= x0[[2]], intervalo[[2]], seir, c(73660.35, optimizacion[[2]]$par, 0.004))[, 4])
plot(suavizado_infectados$X1)

sigma_datos <- list(NULL, NULL, NULL, NULL, NULL, NULL)

susceptibles_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)
expuestos_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)

for (k in 1:length(susceptibles_anuj)) {
  susceptibles_anuj[[k]][1] <- ode(y= x0[[k]], intervalo[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[, 2][1]
  expuestos_anuj[[k]][1] <- ode(y= x0[[k]], intervalo[[k]], seir, c(73660.35, optimizacion[[k]]$par, 0.004))[, 3][1]
}

library(deSolve)
for (i in 1:length(intervalo)) {
  sigma_datos[[i]] <- (1/2*length(intervalo[[i]]))*
    (sum((susceptibles_bogota[intervalo[[i]]] - ode(y= x0[[i]], intervalo[[i]], seir, 
    c(73660.35, optimizacion[[i]]$par, 0.004))[, 2])^2/((ode(y= x0[[i]], intervalo[[i]], seir, 
    c(73660.35, optimizacion[[i]]$par, 0.004))[, 2]*ode(y= x0[[i]], intervalo[[i]], seir, 
    c(73660.35, optimizacion[[i]]$par, 0.004))[, 4])^2)) + 
    sum((expuestos_bogota[intervalo[[i]]] - ode(y= x0[[i]], intervalo[[i]], seir, 
    c(73660.35, optimizacion[[i]]$par, 0.004))[, 3])^2/((ode(y= x0[[i]], intervalo[[i]], seir, 
    c(73660.35, optimizacion[[i]]$par, 0.004))[, 2]*ode(y= x0[[i]], intervalo[[i]], seir, 
    c(73660.35, optimizacion[[i]]$par, 0.004))[, 4])^2)))  
}

r0_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)
for (i in 1:length(intervalo)) {
  r0_anuj[[i]] <- ((optimizacion[[i]]$par[2]*optimizacion[[i]]$par[1]*73660.35)/(0.004*(0.004 + optimizacion[[i]]$par[3])*(0.004 + 
    optimizacion[[i]]$par[2])))
}

r0e_anuj <- list(NULL, NULL, NULL, NULL, NULL, NULL)
for (i in 1:length(intervalo)) {
  r0e_anuj[[i]] <- ((optimizacion[[i]]$par[2]*optimizacion[[i]]$par[1]*73660.35)/(0.004*(0.004 + optimizacion[[i]]$par[3])*(0.004 + 
  optimizacion[[i]]$par[2]))) - ((sigma_datos[[i]])^2*(73660.35)^2*(optimizacion[[i]]$par[2])^2/
  (2*((0.004)^2)*(0.004 + optimizacion[[i]]$par[3])*(0.004 + optimizacion[[i]]$par[2]))) 
}


tiempo <- seq(from = as.Date("2020-03-06"), to = as.Date("2021-03-29"), by=1)
grafica_intervalos <- data.frame(tiempo[1:385], suavizado_infectados$X1[1:385], suavizado_recuperados$X1[1:385], casos$Infectados[1:385], 
                                 casos$recuperados[1:385])

colors <- c("Smoothed infected" = "brown", 
            "Smoothed recovered" = "darkolivegreen", 
            "Infected data" = "indianred1" 
            , "Recovered data" = "olivedrab2")

library(ggplot2)
p_inter <- ggplot(data = grafica_intervalos) + 
  geom_point(aes(x = grafica_intervalos[, 1],  y = grafica_intervalos[, 5], color = 'Recovered data'), size = 1.5) +
  geom_point(aes(x = grafica_intervalos[, 1],  y = grafica_intervalos[, 4], color = 'Infected data'), size = 1.5)  + 
  geom_line(aes(x = grafica_intervalos[, 1],  y = grafica_intervalos[, 3], color = 'Smoothed recovered'), size = 1.1)  +
  geom_line(aes(x = grafica_intervalos[, 1], y = grafica_intervalos[, 2], color = 'Smoothed infected'), size = 1.1) +
  geom_vline(xintercept = grafica_intervalos[152, 1]) + 
  geom_vline(xintercept = grafica_intervalos[205, 1]) + 
  geom_vline(xintercept = grafica_intervalos[265, 1]) + 
  geom_vline(xintercept = grafica_intervalos[306, 1]) + 
  geom_vline(xintercept = grafica_intervalos[363, 1]) + 
  labs(x = 'Time (days)', y = 'Number of population', 
       subtitle = bquote("Intervals:        " ~ I[1] ~ "                         " ~ I[2] ~ "            " 
       ~ I[3] ~ "          " ~ I[4] ~ "           " ~ I[5] ~ "       " ~ I[6])) + 
  scale_color_manual(values = colors, name = "", 
       guide = guide_legend(override.aes = list(linetype = c("solid", "solid", "blank", "blank"),
       shape = c(NA, NA, 16, 16)))) + theme(legend.position="right") + 
  scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y") 

library(gridExtra)
grid.arrange(p_inter, ncol = 1) 

# Distribución de la tasa de infectividad 

histograma <- list(NULL, NULL, NULL, NULL, NULL, NULL)
normal <- list(NULL, NULL, NULL, NULL, NULL, NULL)
for (i in 1:6) {
  histograma[[i]] <- rnorm(10000, mean = optimizacion[[i]]$par[2], sd = sigma_anuj[[i]])
  normal[[i]] <- seq(optimizacion[[i]]$par[2] - 4*sigma_anuj[[i]], 
                     optimizacion[[i]]$par[2] + 4*sigma_anuj[[i]], length = 1000)
}

minimo = min(normal[[1]], normal[[2]], normal[[3]], normal[[4]], normal[[5]], normal[[6]])
maximo = max(normal[[1]], normal[[2]], normal[[3]], normal[[4]], normal[[5]], normal[[6]])

normal_final = seq(min(normal[[1]], normal[[2]], normal[[3]], normal[[4]], normal[[5]], normal[[6]]), 
                   max(normal[[1]], normal[[2]], normal[[3]], normal[[4]], normal[[5]], normal[[6]]), 
                   length = 10000)

plot(density(rnorm(10000, mean = 6.1837e-6, sd = 9.0410e-8)), type = "l", col = "cyan", main = "", 
     xlab = "Infection rate", xlim = c(minimo, maximo))
lines(density(rnorm(10000, mean = 1.0064e-4, sd = 9.3173e-09)), col = "blue")
lines(density(histograma[[3]]), col = "black")
lines(density(histograma[[4]]), col = "red")
lines(density(histograma[[5]]), col = "green")
lines(density(histograma[[6]]), col = "orange")

plot(density(histograma[[2]]))
lines(density(histograma[[3]]))
?hist

minimo_anuj = min(c(6.1837e-6 - 4*9.0410e-8, 1.0064e-4 - 4*9.3173e-09, 1.0064e-4 - 4*9.3173e-09))
maximo_anuj = max(c(6.1837e-6 - 4*9.0410e-8, 6.1837e-6 + 4*9.0410e-8, 
                    1.0064e-4 - 4*9.3173e-09, 1.0064e-4 - 4*9.3173e-09))

df <- data.frame(rnorm(10000, mean = 6.1837e-6, sd = 9.0410e-8), 
                 rnorm(10000, mean = 1.0064e-4, sd = 9.3173e-09), 
                 rnorm(10000, mean = 1.12e-4, sd = 1.9808e-11), 
                 rnorm(10000, mean = 9.9201e-6, sd = 1.1060e-09), 
                 rnorm(10000, mean = 7.6100e-5, sd = 2.7970e-09), 
                 rnorm(10000, mean = 1.8367e-4, sd = 2.6853e-09))
boxplot(df, col = 2:7, xlab = "Infection rate")
boxplot(rnorm(10000, mean = 6.1837e-6, sd = 9.0410e-8), horizontal = T, ylim = c(minimo_anuj, maximo_anuj))
boxplot(rnorm(10000, mean = 1.0064e-4, sd = 9.3173e-09), ylim = c(minimo_anuj, maximo_anuj), add = T, 
        ylim = c(minimo_anuj, maximo_anuj))
boxplot(rnorm(10000, mean = 1.0064e-4, sd = 9.3173e-09), ylim = c(minimo_anuj, maximo_anuj), add = T)


hist(rnorm(10000, mean = 6.1837e-6, sd = 9.0410e-8), col = "cyan", freq = F, breaks = 50, 
     xlab = "Infection rate", main = "")
lines(density(rnorm(10000, mean = 6.1837e-6, sd = 9.0410e-8)), col = "darkblue", lwd = 2)
hist(rnorm(10000, mean = 1.0064e-4, sd = 9.3173e-09), col = "cyan", freq = F, breaks = 50, 
     xlab = "Infection rate", main = "")
lines(density(rnorm(10000, mean = 1.0064e-4, sd = 9.3173e-09)), col = "darkblue", lwd = 2)

boxplot()

hist(histograma[[1]], col = "cyan", freq = F, breaks = 100, xlim = c(min(histograma[[1]]), 
     max(histograma[[1]])), xlab = "Transmission coefficient", main = "")
lines(density(histograma[[1]]), col = "darkblue", lwd = 2)

hist(histograma[[2]], col = "cyan", freq = F, breaks = 100, xlim = c(min(histograma[[2]]), 
     max(histograma[[2]])), xlab = "Transmission coefficient", main = "")
lines(density(histograma[[2]]), col = "darkblue", lwd = 2)

hist(histograma[[3]], col = "cyan", freq = F, breaks = 100, xlim = c(min(histograma[[3]]), 
     max(histograma[[3]])), xlab = "Transmission coefficient", main = "")
lines(density(histograma[[3]]), col = "darkblue", lwd = 2)

hist(histograma[[4]], col = "cyan", freq = F, breaks = 100, xlim = c(min(histograma[[4]]), 
     max(histograma[[4]])), xlab = "Transmission coefficient", main = "")
lines(density(histograma[[4]]), col = "darkblue", lwd = 2)

hist(histograma[[5]], col = "cyan", freq = F, breaks = 100, xlim = c(min(histograma[[5]]), 
     max(histograma[[5]])), xlab = "Transmission coefficient", main = "")
lines(density(histograma[[5]]), col = "darkblue", lwd = 2)

hist(histograma[[6]], col = "cyan", freq = F, breaks = 100, xlim = c(min(histograma[[6]]), 
     max(histograma[[6]])), xlab = "Transmission coefficient", main = "")
lines(density(histograma[[6]]), col = "darkblue", lwd = 2)
  
# Enfoque por actualización de datos 

expuestos <- (proyeccion_bogota_final[2:386] - suavizado_infectados$X1[2:386] - 
                suavizado_recuperados$X1[2:386] - 
  (proyeccion_bogota_final[1:385] - suavizado_infectados$X1[1:385] - suavizado_recuperados$X1[1:385]) + 
    73660.35 - 0.004*(proyeccion_bogota_final[1:385] - suavizado_infectados$X1[1:385] - 
                        suavizado_recuperados$X1[1:385]))
upsilon_no_param <- c(rep(optimizacion[[1]]$par[2], length(intervalo[[1]])), 
                      rep(optimizacion[[2]]$par[2], length(intervalo[[2]])), 
                      rep(optimizacion[[3]]$par[2], length(intervalo[[3]])), 
                      rep(optimizacion[[4]]$par[2], length(intervalo[[4]])), 
                      rep(optimizacion[[5]]$par[2], length(intervalo[[5]])), 
                      rep(optimizacion[[6]]$par[2], length(intervalo[[6]])))
plot(v_expuestos)

gamma_no_para <- (suavizado_recuperados$X1[2:386] - suavizado_recuperados$X1[1:385] + 
                 0.004*suavizado_recuperados$X1[1:385])/suavizado_infectados$X1[1:385]

#### Parámetros para ajustarse al número reproductivo básico 

beta_ajuste <- NULL
for (i in 1:6) {
  beta_ajuste[i] <- (0.004*(optimizacion[[i]]$par[2]+0.004)*(optimizacion[[i]]$par[3]+0.004))/
                    (73660.35*optimizacion[[i]]$par[2])
}

optimizacion123 <- list(NULL, NULL, NULL, NULL, NULL, NULL)

optimizacion123[[1]] <- optim(par = c(beta_ajuste[1], parametros_anuj[[1]][2], parametros_anuj[[1]][3]),
                              fn = cuadrados_ordinarios[[1]], 
                              data = data.frame(suavizado_infectados$X1[interalo_anuj[[i]]]), 
                              lower = c(beta_ajuste[1], optimizacion[[1]]$par[2], optimizacion[[1]]$par[3]), 
                              upper = c(1, 1, 0.25), 
                              method = "L-BFGS-B")

plot(ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion123[[1]]$par, 0.004))[, 4], 
     type = "l", col = "darkblue")
lines(suavizado_infectados$X1[intervalo[[1]]], type = "p", col = "cyan")

plot(ode(y= x0_1[[1]], interalo_anuj[[1]], seir, c(73660.35, optimizacion123[[1]]$par, 0.004))[, 4])
lines(suavizado_infectados$X1[intervalo[[1]]])
