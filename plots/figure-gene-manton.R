





# gene ENSG00000116771 - AGMAT 
# from file final_data.down_100.ct_0.csv (MANTON)

# lowly expressed gene


gene1_100 <- c(0.499492341137929,0.785602232587075,0,0.306488876060667,0.329025233227684,
               0.641453029428056,0.583220767024746,0.840653295076293)

gene1_50 <- c(0.936260426421858,0,0,0.451916147436599,0.615811389462018,0.648702867682386,
              0.980766757308285,1.27339084215413)

gene1_10 <- c(0,0,0,0,1.45426993228685,1.41685772203676,0,1.39793150871455)



# gene ENSG00000082074 - FYB1

gene2_100 <- c(5.17741775699942,5.50901539876429,4.641276687476,5.34411988935911,
               4.63926429829538,5.287929786637,5.4745371595259,4.81675041129246)

gene2_50 <- c(4.93438701626118,5.51677296588069,4.58982071440857,5.30595504768129,
              4.29522110066624,5.14066152107492,5.39049310358382,4.67608046069784)

gene2_10 <- c(4.98318813886037,5.6675023945639,4.5364896924198,5.20855550049124,
              4.41924915206925,5.31845893630901,5.51465155195293,4.59003883934377)


df1_50 <- data.frame(gene1_100, gene1_50)
m1_50 <- lm(df1_50$gene1_100 ~ df1_50$gene1_50, df1_50)
rsq1_50 = format(summary(m1_50)$r.squared, digits = 2)
p1 <- ggplot(df1_50, aes(x=gene1_100, y=gene1_50)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("50% downsampling") +
  ggtitle(sprintf("AGMAT, R^2 = %s", rsq1_50)) +
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
  ggtitle(sprintf("AGMAT, R^2 = %s", rsq1_10)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))







df2_50 <- data.frame(gene2_100, gene2_50)
m2_50 <- lm(df2_50$gene2_100 ~ df2_50$gene2_50, df2_50)
rsq2_50 = format(summary(m2_50)$r.squared, digits = 2)
p3 <- ggplot(df2_50, aes(x=gene2_100, y=gene2_50)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("50% downsampling") +
  ggtitle(sprintf("FYB1, R^2 = %s", rsq2_50)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))


df2_10 <- data.frame(gene2_100, gene2_10)
m2_10 <- lm(df2_10$gene2_100 ~ df2_10$gene2_10, df2_10)
rsq2_10 = format(summary(m2_10)$r.squared, digits = 2)
p4 <- ggplot(df2_50, aes(x=gene2_100, y=gene2_10)) + geom_point(size=4) +
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Original expression") +
  ylab("10% downsampling") +
  ggtitle(sprintf("FYB1, R^2 = %s", rsq2_10)) +
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(size=20))

plot_grid(p1, p2, p3, p4, scale=0.9, labels=c("A", "B", "C", "D"), label_size = 30)




