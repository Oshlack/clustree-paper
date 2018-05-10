set.seed(1)

sim0 <- cbind(runif(1000, 0, 100), runif(1000, 0, 100))
sim1 <- cbind(rnorm(1000, mean = 50, sd = 10), rnorm(1000, mean = 50, sd = 10))
sim2 <- rbind(
    cbind(rnorm(500, mean = 25, sd = 8), rnorm(500, mean = 50, sd = 8)),
    cbind(rnorm(500, mean = 75, sd = 8), rnorm(500, mean = 50, sd = 8))
)
sim3 <- rbind(
    cbind(rnorm(333, mean = 35, sd = 8), rnorm(333, mean = 75, sd = 8)),
    cbind(rnorm(333, mean = 65, sd = 8), rnorm(333, mean = 75, sd = 8)),
    cbind(rnorm(334, mean = 50, sd = 8), rnorm(334, mean = 25, sd = 8))
)
sim4 <- rbind(
    cbind(rnorm(250, mean = 50, sd = 8), rnorm(250, mean = 75, sd = 8)),
    cbind(rnorm(250, mean = 80, sd = 8), rnorm(250, mean = 75, sd = 8)),
    cbind(rnorm(250, mean = 50, sd = 8), rnorm(250, mean = 25, sd = 8)),
    cbind(rnorm(250, mean = 20, sd = 8), rnorm(250, mean = 25, sd = 8))
)

restrict_values <- function(x, min, max) {
    x[x > max] <- max
    x[x < min] <- min

    return(x)
}

sim0 <- restrict_values(sim0, 0, 100)
sim1 <- restrict_values(sim1, 0, 100)
sim2 <- restrict_values(sim2, 0, 100)
sim3 <- restrict_values(sim3, 0, 100)
sim4 <- restrict_values(sim4, 0, 100)

add_noise <- function(sim) {
    cbind(sim,
          runif(1000, 0, 50), runif(1000, 0, 50),
          runif(1000, 0, 50), runif(1000, 0, 50),
          runif(1000, 0, 50), runif(1000, 0, 50),
          runif(1000, 0, 50), runif(1000, 0, 50))
}

sim0 <- add_noise(sim0)
sim1 <- add_noise(sim1)
sim2 <- add_noise(sim2)
sim3 <- add_noise(sim3)
sim4 <- add_noise(sim4)

cluster_sim <- function(sim, k) {
    clusterings <- sapply(1:k, function(x) {
        km <- kmeans(sim, centers = x, iter.max = 100, nstart = 10)
        km$cluster
    })
    colnames(clusterings) <- paste0("K", 1:k)

    sim_clusts <- data.frame(x = sim[, 1],
                             y = sim[, 2])

    sim_clusts <- cbind(sim_clusts, clusterings)

    return(sim_clusts)
}

clusts0 <- cluster_sim(sim0, 8)
clusts1 <- cluster_sim(sim1, 8)
clusts2 <- cluster_sim(sim2, 8)
clusts3 <- cluster_sim(sim3, 8)
clusts4 <- cluster_sim(sim4, 8)

p0 <- ggplot(clusts0, aes(x = x, y = y)) + geom_point() +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_cowplot()
p1 <- ggplot(clusts1, aes(x = x, y = y)) + geom_point() +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_cowplot()
p2 <- ggplot(clusts2, aes(x = x, y = y)) + geom_point() +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_cowplot()
p3 <- ggplot(clusts3, aes(x = x, y = y)) + geom_point() +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_cowplot()
p4 <- ggplot(clusts4, aes(x = x, y = y)) + geom_point() +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_cowplot()

t0 <- clustree(clusts0, prefix = "K")
t1 <- clustree(clusts1, prefix = "K")
t2 <- clustree(clusts2, prefix = "K")
t3 <- clustree(clusts3, prefix = "K")
t4 <- clustree(clusts4, prefix = "K")

# ss <- cbind(sim4,
#             runif(1000, 0, 50), runif(1000, 0, 50),
#             runif(1000, 0, 50), runif(1000, 0, 50),
#             runif(1000, 0, 50), runif(1000, 0, 50),
#             runif(1000, 0, 50), runif(1000, 0, 50))
# cc <- cluster_sim(ss, 8)
# ggplot(cc, aes(x = x, y = y, colour = factor(K4))) + geom_point()
# ggplot(cc, aes(x = x, y = y, colour = factor(K8))) + geom_point()
#
# clustree(cc, prefix = "K", prop_filter = 0.1)

panel <- plot_grid(p0, t0, p1, t1, p2, t2, p3, t3, p4, t4, ncol = 2,
                   rel_widths = c(1, 1.5))
save_plot("output/sim_panel.png", panel, nrow = 5, ncol = 2,
          base_width = 6, base_height = 4)


ndim <- 100
nsamples <- 1000

set.seed(1)
centre1 <- rnorm(ndim, sd = 10)
shift2  <- rnorm(ndim, sd = 2)
centre2 <- centre1 + shift2
centre3 <- (centre1 + centre2) / 2 + rnorm(ndim, sd = 5)
centre4 <- centre3 + shift2 * 0.5 + rnorm(ndim, sd = 2)

points1 <- matrix(rnorm(ndim * nsamples, mean = centre1, sd = 5),
                  nrow = nsamples, ncol = ndim, byrow = TRUE)
points2 <- matrix(rnorm(ndim * nsamples, mean = centre2, sd = 5),
                  nrow = nsamples, ncol = ndim, byrow = TRUE)
points3 <- matrix(rnorm(ndim * nsamples, mean = centre3, sd = 5),
                  nrow = nsamples, ncol = ndim, byrow = TRUE)
points4 <- matrix(rnorm(ndim * nsamples, mean = centre4, sd = 5),
                  nrow = nsamples, ncol = ndim, byrow = TRUE)

sim0 <- matrix(runif(ndim * nsamples, 0, 10),
               nrow = nsamples, ncol = ndim, byrow = TRUE)
sim1 <- points1
sim2 <- rbind(points1, points2)
sim3 <- rbind(points1, points2, points3)
sim4 <- rbind(points1, points2, points3, points4)

clusts0 <- cluster_sim(sim0, 8)
clusts1 <- cluster_sim(sim1, 8)
clusts2 <- cluster_sim(sim2, 8)
clusts3 <- cluster_sim(sim3, 8)
clusts4 <- cluster_sim(sim4, 8)

plot_sim_pca <- function(data, ngroups = 0, group_size = NA) {
    pca <- as.data.frame(prcomp(data)$x)


    if (ngroups > 0) {
        pca$Group <- factor(rep(seq_len(ngroups), each = group_size))

        gg <- ggplot(pca, aes(x = PC1, y = PC2)) +
            geom_point(aes(colour = Group), alpha = 0.5) +
            scale_colour_manual(values = c("#EC008C", "#00ADEF", "#8DC63F",
                                           "#F47920"))
    } else {
        gg <- ggplot(pca, aes(x = PC1, y = PC2)) +
            geom_point(alpha = 0.5)
    }

    gg <- gg +
        theme_cowplot() +
        theme(legend.position = "none")

    return(gg)
}

p0 <- plot_sim_pca(sim0) +
    ggtitle("Uniform noise") +
    theme(plot.title = element_text(size = 30, hjust = -0.2))
p1 <- plot_sim_pca(sim1, 1, nsamples) +
    ggtitle("Single cluster") +
    theme(plot.title = element_text(size = 30, hjust = -0.2))
p2 <- plot_sim_pca(sim2, 2, nsamples) +
    ggtitle("Two clusters") +
    theme(plot.title = element_text(size = 30, hjust = -0.2))
p3 <- plot_sim_pca(sim3, 3, nsamples) +
    ggtitle("Three clusters") +
    theme(plot.title = element_text(size = 30, hjust = -0.2))
p4 <- plot_sim_pca(sim4, 4, nsamples) +
    ggtitle("Four clusters") +
    theme(plot.title = element_text(size = 30, hjust = -0.2))

t0 <- clustree(clusts0, prefix = "K") +
    theme(legend.position = "none",
          plot.margin = unit(c(1,1,1.5,1.2),"cm"))
t1 <- clustree(clusts1, prefix = "K") +
    theme(legend.position = "none",
          plot.margin = unit(c(1,1,1.5,1.2),"cm"))
t2 <- clustree(clusts2, prefix = "K") +
    theme(legend.position = "none",
          plot.margin = unit(c(1,1,1.5,1.2),"cm"))
t3 <- clustree(clusts3, prefix = "K") +
    theme(legend.position = "none",
          plot.margin = unit(c(1,1,1.5,1.2),"cm"))
t4 <- clustree(clusts4, prefix = "K") +
    theme(legend.position = "none",
          plot.margin = unit(c(1,1,1.5,1.2),"cm"))

p4_legend <- p4 +
    guides(color = guide_legend(title.position = "top",
                                title.hjust = 0.5,
                                label.position = "top",
                                label.hjust = 0.5,
                                override.aes = list(size = 10))) +
    theme(legend.position = "bottom",
          legend.justification = "center",
          legend.title = element_text(size = 20))
t0_legend <- t0 +
    guides(size = guide_legend(title = "Cluster size",
                               title.position = "top",
                               title.hjust = 0.5,
                               label.position = "top",
                               label.hjust = 0.5,
                               order = 1),
           color = guide_legend(title.position = "top",
                                title.hjust = 0.5,
                                override.aes = list(size = 8),
                                order = 2),
           edge_colour = guide_edge_colourbar(title = "Sample count",
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              barwidth = 6,
                                              barheight = 2.2,
                                              draw.ulim = TRUE,
                                              draw.llim = TRUE,
                                              order = 3),
           edge_alpha = guide_legend(title = "In-proportion",
                                     title.position = "top",
                                     title.hjust = 0.5,
                                     label.position = "top",
                                     label.hjust = 0.5,
                                     override.aes = list(size = 10),
                                     order = 4)) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20))

l1 <- get_legend(p4_legend)
l2 <- get_legend(t0_legend)

panel <- plot_grid(p0, t0, p1, t1, p2, t2, p3, t3, p4, t4, ncol = 2,
                   rel_widths = c(1, 1.6), rel_heights = c(1, 1, 1, 1, 1))
legend <- plot_grid(l1, l2, ncol = 2, rel_widths = c(1, 1.6))
panel_legend <- plot_grid(panel, legend, nrow = 2, rel_heights = c(1, 0.05))

save_plot("output/sim_panel.png", panel_legend, nrow = 6, ncol = 2,
          base_width = 8, base_height = 6)
