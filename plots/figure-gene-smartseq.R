

# gene ENSG00000285096 - SLC14A2 
# from file final_data.down_100.ct_0.csv

# lowly expressed gene

gene1_100 <- c(0.171186235663,0.0,0.336935970084,0.0,0.500937622474,
               0.103572064785,0.120632844157,0.0587805934166,0.150226759978,0.10863779452)

gene1_50 <- c(0.12987861532,0.0,0.296654573777,0.0,0.482384655967,0.111538378984,
              0.118066066555,0.0121270411783,0.0951871695345,0.103493062837)

gene1_10 <- c(0.128948281104,0.0,0.390179296677,0.0,0.493185166975,0.0969590620861,
              0.113733689111,0.0535598335227,0.136109658866,0.094511625227)



# gene ENSG00000115263 - GCG

gene2_100 <- c(11.0385667462,0.0,11.7559233594,0.0,11.1530485348,11.7699957794,11.9357756019,
               11.6095334546,12.0229436529,11.844197698)

gene2_50 <- c(11.0213463142,0.0,11.738549851,0.0,11.1203593506,11.7423397494,11.9135362436,
              11.5816271098,11.9962678238,11.8260592202)

gene2_10 <- c(11.0336341795,0.0,11.7544524056,0.0,11.1498999466,11.7665343304,11.9322709145,
              11.6050559447,12.0194985264,11.8422469617)

df1_50 <- data.frame(gene1_100, gene1_50)
m1_50 <- lm(df1_50$gene1_100 ~ df1_50$gene1_50, df1_50)
rsq1_50 = format(summary(m1_50)$r.squared, digits = 2)
df1_50 <- data.frame(gene1_100, gene1_50)
p1 <- ggplot(df1_50, aes(x=gene1_100, y=gene1_50)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("50% downsampling") +
  ggtitle(sprintf("SLC14A2, R^2 = %s", rsq1_50)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))


df1_10 <- data.frame(gene1_100, gene1_10)
m1_10 <- lm(df1_10$gene1_100 ~ df1_10$gene1_10, df1_10)
rsq1_10 = format(summary(m1_10)$r.squared, digits = 2)
df1_10 <- data.frame(gene1_100, gene1_10)
p2 <- ggplot(df1_50, aes(x=gene1_100, y=gene1_10)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("10% downsampling") +
  ggtitle(sprintf("SLC14A2, R^2 = %s", rsq1_10)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))







df2_50 <- data.frame(gene2_100, gene2_50)
m2_50 <- lm(df2_50$gene2_100 ~ df2_50$gene2_50, df2_50)
rsq2_50 = format(summary(m2_50)$r.squared, digits = 2)
df2_50 <- data.frame(gene2_100, gene2_50)
p3 <- ggplot(df2_50, aes(x=gene2_100, y=gene2_50)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("50% downsampling") +
  ggtitle(sprintf("GCG, R^2 = %s", rsq2_50)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))


df2_10 <- data.frame(gene2_100, gene2_10)
m2_10 <- lm(df2_10$gene2_100 ~ df2_10$gene2_10, df2_10)
rsq2_10 = format(summary(m2_10)$r.squared, digits = 2)
df2_10 <- data.frame(gene2_100, gene2_10)
p4 <- ggplot(df2_50, aes(x=gene2_100, y=gene2_10)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("10% downsampling") +
  ggtitle(sprintf("GCG, R^2 = %s", rsq2_10)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))

plot_grid(p1, p2, p3, p4, scale=0.9, labels=c("A", "B", "C", "D"), label_size = 30)


