



# Figure 1 - this is demultiplexing quality as a function of coverage

reads_per_cell <- 40 * 0.01 * c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75)
percentage_correctly_assigned <- (1 / 4800) * c(1972, 3108, 3707, 4041, 4275, 4410, 4485, 4563, 4599, 4628, 4661, 4678, 4689, 4708, 4722)


demuxlet_data <- data.frame(coverage = reads_per_cell, accuracy = percentage_correctly_assigned)


p3 <- ggplot(demuxlet_data, aes(x=coverage, y=accuracy)) +
  geom_line(size=2) +
  geom_point(shape=21, fill="white", size=3, stroke=2) + 
  scale_x_continuous(breaks = round(seq(0, 30, by = 5),1)) +
  ylab(expression("Demultiplexing accuracy")) +
  xlab("Coverage (thousand reads per cell)") +
  theme(legend.position = "none",
        axis.text = element_text(size=30),
        axis.text.x = element_text(size=30),
        axis.title = element_text(size=30),
        plot.title = element_text(size=30))
p3





