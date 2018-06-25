#' Get simulations
#'
#' Get a list of demonstration simulations
#'
#' @param ndim number of dimensions to simulate in each cluster
#' @param nsamples number of samples to simulate
#'
#' @details The returned synthetic datasets are:
#'
#' 1. Random uniform noise
#' 2. A single normally distributed cluster
#' 3. Two clusters
#' 4. Three clusters
#' 5. Four clusters
#'
#' Each cluster contains `nsamples` points. Datasets with multiple clusters
#' contain the previous dataset plus a new set of points, for example dataset 3
#' contains all the points in dataset 2 plus another `nsamples` points in a new
#' cluster.
#'
#' @return list of matrices
get_sims <- function(ndim, nsamples) {
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

    sims <- list(sim0, sim1, sim2, sim3, sim4)

    return(sims)
}


#' Cluster a simulation
#'
#' Cluster a simulation using kmeans with a range of ks
#'
#' @param sim simulated data matrix to cluster
#' @param ks vector of k values to use for clustering
#'
#' @return matrix assiging samples to clusters at different k values
cluster_sim <- function(sim, ks) {
    clusterings <- sapply(ks, function(k) {
        km <- kmeans(sim, centers = k, iter.max = 100, nstart = 10)
        km$cluster
    })
    colnames(clusterings) <- paste0("K", ks)

    return(clusterings)
}


#' Plot simulation PCA
#'
#' Produce a PCA plot of a simulation
#'
#' @param sim simulated data matrix
#' @param ngroups number of clusters in the simulation
#' @param group_size number of samples in each cluster
#'
#' @return ggplot object
plot_sim_pca <- function(sim, ngroups = 0, group_size = NA) {
    pca <- as.data.frame(prcomp(sim)$x)

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
        cowplot::theme_cowplot() +
        theme(legend.position = "none",
              plot.title = element_text(size = 30, hjust = -0.2))

    return(gg)
}


#' Make simulation panel
#'
#' Produce the simulations panel shown in the paper
#'
#' @param sims list of simulated datasets
#' @param clusts list of clustering matrices
#'
#' @return ggplot panel object
make_sim_panel <- function(sims, clusts) {

    titles <- c("A - Uniform noise",
                "B - Single cluster",
                "C - Two clusters",
                "D - Three clusters",
                "E - Four clusters")

    pca_plots <- lapply(seq_along(sims), function(x) {
        plot_sim_pca(sims[[x]], ngroups = x - 1, group_size = nrow(sims[[1]])) +
            ggtitle(titles[x])
    })

    tree_plots <- lapply(clusts, function(x) {
        clustree::clustree(x, prefix = "K", node_size_range = c(10, 20),
                           node_text_size = 8) +
            expand_limits(y = 7.2) +
            theme(legend.position = "none",
                  plot.margin = unit(c(1, 1, 1.5, 1.2), "cm"))
    })

    sc3_plots <- lapply(clusts, function(x) {
        clustree::clustree(x, prefix = "K", node_size_range = c(10, 20),
                           node_text_size = 8,
                           node_colour = "sc3_stability") +
            scale_colour_viridis(option = "plasma", begin = 0.3,
                                 limits = c(0, 0.7)) +
            expand_limits(y = 7.2) +
            theme(legend.position = "none",
                  plot.margin = unit(c(1, 1, 1.5, 1.2), "cm"))
    })

    pca_legend <- cowplot::get_legend(
        pca_plots[[5]] +
        guides(color = guide_legend(title.position = "top",
                                    title.hjust = 0.5,
                                    label.position = "top",
                                    label.hjust = 0.5,
                                    override.aes = list(size = 10))) +
        theme(legend.position = "bottom",
              legend.justification = "center",
              legend.title = element_text(size = 20))
    )

    tree_legend <- cowplot::get_legend(
        tree_plots[[1]] +
            guides(size = guide_legend(title = "Cluster size",
                                       title.position = "top",
                                       title.hjust = 0.5,
                                       label.position = "top",
                                       label.hjust = 0.5,
                                       order = 1),
                   color = guide_legend(title = "k",
                                        title.position = "top",
                                        title.hjust = 0.5,
                                        override.aes = list(size = 8),
                                        order = 2),
                   edge_colour = guide_edge_colourbar(title = "Sample count",
                                                      title.position = "top",
                                                      title.hjust = 0.5,
                                                      barwidth = 8,
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
                  legend.justification = "center",
                  legend.title = element_text(size = 20))
    )

    sc3_legend <- cowplot::get_legend(
        sc3_plots[[1]] +
            guides(size = FALSE,
                   color = guide_colourbar(title = "SC3 stability",
                                           title.position = "top",
                                           title.hjust = 0.5,
                                           barwidth = 30,
                                           barheight = 2.2,
                                           draw.ulim = TRUE,
                                           draw.llim = TRUE,
                                           order = 2),
                   edge_colour = FALSE,
                   edge_alpha = FALSE) +
            theme(legend.position = "bottom",
                  legend.justification = "center",
                  legend.title = element_text(size = 20))
    )

    panel <- cowplot::plot_grid(pca_plots[[1]], tree_plots[[1]], sc3_plots[[1]],
                                pca_plots[[2]], tree_plots[[2]], sc3_plots[[2]],
                                pca_plots[[3]], tree_plots[[3]], sc3_plots[[3]],
                                pca_plots[[4]], tree_plots[[4]], sc3_plots[[4]],
                                pca_plots[[5]], tree_plots[[5]], sc3_plots[[5]],
                                ncol = 3,
                                rel_widths = c(1, 1.2, 1.2),
                                rel_heights = c(1, 1, 1, 1, 1))

    legend <- cowplot::plot_grid(pca_legend, tree_legend, sc3_legend,
                                 ncol = 3, rel_widths = c(1, 1.2, 1.2))

    panel_legend <- cowplot::plot_grid(panel, legend,
                                       nrow = 2, rel_heights = c(1, 0.05))

    return(panel_legend)
}
