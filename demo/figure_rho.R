data("Numden1")
zzz_merge <- rhofrommc(rbind(Numden[[1]], Numden[[2]], Numden[[3]], Numden[[4]], Numden[[5]],
                             Numden[[6]], Numden[[7]], Numden[[8]], Numden[[9]], Numden[[10]]))
# Nwna <- nrow(zzz_merge)
zzz_merge <- na.omit(zzz_merge)
# Nwona <- nrow(zzz_merge)
# Nwona / Nwna
zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho
yexp1 <- c(outer(c(1, 2, 5), 10 ^ {-3:3}, "*"))
yexp1 <- c(-yexp1, yexp1)
yexp2 <- c(0.1, 0.5, 1, 2, 10, 100)
yexp2 <- c(-yexp2, yexp2)

ggplot(data = subset(zzz_merge, .Beta == 1 & is.element(y, yexp2)),
       aes(x = .Gamma, y = rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans='log10') + 
  geom_ribbon()+facet_wrap(~factor(y))


ggplot(data = subset(zzz_merge, .Gamma == 1 & is.element(y, yexp2)),
       aes(x = .Beta, y=rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = 'log10') +
  geom_ribbon() + facet_wrap(~factor(y))


ggplot(data = subset(zzz_merge, .Gamma == 1 & .Beta == 1),
       aes(x = y, y = rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = 'log10') +
  geom_ribbon()

