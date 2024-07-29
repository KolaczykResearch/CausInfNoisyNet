# set a destination for output
outdir <- "../Data/"

set.seed(1)
norm_sd=sqrt(.5)
# School 1
O.c.school1.4 <- rnorm(153,2,norm_sd)
O.c.school1 <- rbind(O.c.school1.4*4,O.c.school1.4*3,O.c.school1.4*2.5,O.c.school1.4)

# School 2
O.c.school2.4 <- rnorm(147,2,norm_sd)
O.c.school2 <- rbind(O.c.school2.4*4,O.c.school2.4*3,O.c.school2.4*2.5,O.c.school2.4)

# School 3
O.c.school3.4 <- rnorm(159,2,norm_sd)
O.c.school3 <- rbind(O.c.school3.4*4,O.c.school3.4*3,O.c.school3.4*2.5,O.c.school3.4)

# School 4
O.c.school4.4 <- rnorm(76,2,norm_sd)
O.c.school4 <- rbind(O.c.school4.4*4,O.c.school4.4*3,O.c.school4.4*2.5,O.c.school4.4)

save(O.c.school1,O.c.school2,O.c.school3,O.c.school4, file = paste0(outdir,"outcome.RData"))
