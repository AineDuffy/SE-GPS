# ===================================================================
# Description: Contains plotting functions used to generate forest plots 
#              for the manuscript.
# ===================================================================


#Forest plot with OR
create_forest_plot <- function(df, yvar, groupvar, shapevar, colorvar, ci_lower, ci_upper, label_above=NULL, label_below=NULL ,y_levels,
                               dodge_width = 0.5, shape_vals = NULL, maxor = NULL, coord_limit = NULL) {
  
  
  dodge <- position_dodge2(width = dodge_width)
  
  p <- ggplot(data = df, aes(y = .data[[yvar]], x = OR, color = .data[[colorvar]])) +
    geom_point(size = 2, aes(shape = .data[[shapevar]]), stroke = 1.5, position = position_dodge2(width =0.5)) +
    scale_shape_manual(values = shape_vals) + guides(shape = "none") +
    geom_linerange(aes(xmin = .data[[ci_lower]], xmax = .data[[ci_upper]]), color = 'black', position = dodge) +
    xlab("Odds ratio") + ylab("") +
    scale_y_discrete(limits = rev(levels(y_levels))) +
    geom_vline(xintercept = 1, linetype = 'longdash', color = 'red') +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=14),strip.text = element_text(size = 13))  +
    theme_classic() +
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_line(colour = 'black', size = 0.5),
          axis.ticks = element_line(colour = "black", size = 0.5),
          axis.ticks.length = unit(0.25, "cm"))
  
  if (!is.null(coord_limit)) {
    p <- p + coord_cartesian(xlim = c(0, coord_limit))
  }
  
  if (!is.null(label_above)) {
    p <- p + geom_text(aes(x = 0.1, y = order + 0.2), label = label_above, size = 2.7, color = '#D41159')
  }
  
  if (!is.null(label_below)) {
    p <- p + geom_text(aes(x = 0.1, y = order - 0.2), label = label_below, size = 2.7, color = '#1A85FF')
  }
  

  return(p)
}

#If CI are too long add arrows at upper cut off 
add_arrows <- function(fp, maxor, yvar){
  
    fp <- fp + scale_x_continuous(limits = c(0, maxor)) +
      geom_segment(aes(x = lowerCI, xend = upperCI_cut,
                       y = .data[[yvar]], yend = .data[[yvar]]),
                   arrow = arrow(length = unit(0.25, "cm")))

}

# Include OR and CI
make_side_labels <- function(df, yvar, label, y_levels, size = 3.1, hjust_size) {
  ggplot(df) +
    geom_text(aes(x = 0, y = .data[[yvar]], label = .data[[label]]), size = size, hjust = 0) +
    scale_y_discrete(limits = rev(levels(y_levels))) +
    ggtitle('Odds ratio (95% CI)') +
    theme_void() +
    theme(plot.title = element_text(size = 11, hjust = hjust_size))
}


combine_plots <- function(fp, labels_plot, widths, legend = "bottom") {
  combined_plot <- fp + plot_spacer() + labels_plot + plot_layout(widths = widths, guides = "collect") & 
    theme(legend.position = legend, legend.box = "horizontal",legend.justification = "left")
}


# 0.3 increment plot with table underneath 
plot_binned_or_with_table <- function(binned_gps2, dataset, drugs=NULL, plots_dir, title_PLT, score_plt, shape_vals, max_uci, filename) {
  
  # Format OR and CI values
  format_values <- function(x) {
    if (abs(x) > 30 || abs(x) < 0.1) {
      return(formatC(x, format = "e", digits = 1))
    } else {
      return(round(x, 1))
    }
  }
  
  binned_gps2$upperCI2 <- sapply(binned_gps2$upperCI, format_values)
  binned_gps2$lowerCI2 <- sapply(binned_gps2$lowerCI, format_values)
  binned_gps2$CI2 <- paste0(binned_gps2$lowerCI2, " - ", binned_gps2$upperCI2)
  binned_gps2$OR2 <- sapply(binned_gps2$OR, format_values)
  
  # Plot
  plot1 <- ggplot(data = binned_gps2, aes(x = genescoresum, y = OR)) + 
    geom_point(aes(color = OR, shape = or_sig), position = position_dodge(width = 0.2), size = 2.5, stroke = 1.5) +
    xlab('SE-GPS bins') + ylab('Odds ratio (95% CI)') +
    scale_shape_manual(values = shape_vals) + guides(shape = "none") +
    ggtitle(title_PLT) +
    scale_x_discrete(limits = levels(binned_gps2$genescoresum)) +
    geom_linerange(aes(ymin = lowerCI, ymax = upperCI, color = OR), position = position_dodge(width = 0.2)) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 6), breaks = c(0, 2, 4, 6), labels = c(0, 2, 4, 6), oob = scales::squish) +
    geom_hline(yintercept = 1, color = "grey") +
    coord_cartesian(ylim = c(0, max_uci), clip = 'off') +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black", size = 13),
      axis.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      axis.line = element_line(colour = 'black', size = 0.5),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    scale_y_continuous(limits = c(-0.2, max_uci))
  
  # Optional arrows for capped CIs
  if (length(binned_gps2$upperCI_cut[!is.na(binned_gps2$upperCI_cut)]) >= 1) {
    plot1 <- plot1 +
      geom_segment(aes(x = genescoresum, xend = genescoresum, y = lowerCI, yend = upperCI_cut, color = OR),
                   arrow = arrow(length = unit(0.25, "cm")))
  }
  
  # Table for display
  table_thres <- binned_gps2 %>%
    distinct(genescoresum, Percentile, genes, parentterms, OR2, CI2, P.val) %>%
    mutate(
      Percentile = round(Percentile, 2),
      `P-value` = sprintf("%.2e", P.val)
    ) %>%
    dplyr::rename(
      `SE-GPS` = genescoresum,
      OR = OR2,
      CI = CI2,
      Genes = genes,
      Phenotypes = parentterms
    )
  
  g <- tableGrob(table_thres, rows = NULL, theme = ttheme_minimal(core = list(bg_params = list(fill = "white")),
                                                                  colhead = list(bg_params = list(fill = "white"))))
  
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))
  
  
  # Save
  png(file = filename, width = 500, height = 600)
  grid.arrange(plot1, g, nrow = 2, heights = c(2, 1))
  dev.off()
  
  return(list(plot = plot1, table = table_thres, file = filename))
}