boundcalls.data <- data.frame(
  cubes = bd_cube,
  simplices = bd_simplex,
  cross = bd_cross,
  Dimension = seq(from=5, to =100, by=5)
)

boundcalls.data %>%
gather(polytopes,OracleCalls,cubes, simplices,cross) %>%
    ggplot(aes(x=Dimension, y=OracleCalls, colour=polytopes)) +
    geom_line() + geom_point() + scale_x_continuous(breaks = seq(from=5,to=100,by=5)) +
  scale_y_continuous(breaks = c(0,5e+05,10e+05,15e+05,20e+05,22e+05), label = c(0, 5, 10, 15, 20, "x 10^5")) +
  theme(legend.position="top",text = element_text(size=20), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

TotSteps.data <- data.frame(
  cubes = sqrt(step_cube),
  simplices = sqrt(step_simplex),
  Dimension = seq(from=5, to =100, by=5)
)

TotSteps.data %>%
  gather(polytopes,sqrt_steps,cubes, simplices) %>%
  ggplot(aes(x=Dimension, y=sqrt_steps, colour=polytopes)) +
  geom_line() + geom_point() + scale_x_continuous(breaks = seq(from=5,to=100,by=5)) +
  theme(legend.position="top",text = element_text(size=20), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
