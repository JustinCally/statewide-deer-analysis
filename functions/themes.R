delwp_cols <- c(Teal = "#00B2A9",
                Navy = "#201547",
                Environment = "#CEDC00",
                `Climate Change` = "#FDDA24",
                Energy = "#0072CE",
                Water = "#71C5E8",
                Planning = "#642667",
                Infrastructure = "#AF272F",
                FFR = "#E57200",
                Corporate = "#201547")

delwp_cols_seq <- lapply(delwp_cols, function(x) {
  tinter::tinter(x, direction = "tints", steps = 10)
})

delwp_cols_shades <- lapply(delwp_cols, function(x) {
  tinter::tinter(x, direction = "shades", steps = 10)
})

delwp_theme <- function (base_size = 11,
                         base_family = "",
                         base_line_size = base_size/22,
                         base_rect_size = base_size/22)
{
  theme_grey(base_size = base_size, base_family = base_family,
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(panel.background = element_rect(fill = "white",
                                          colour = NA),
          panel.border = element_rect(fill = NA,
                                      colour = "grey20"),
          panel.grid = element_line(colour = "grey92"),
          panel.grid.minor = element_line(size = rel(0.5)),
          strip.background = element_rect(fill = "grey85",
                                          colour = "grey20"),
          legend.key = element_rect(fill = "white",
                                    colour = NA),
          complete = TRUE)
}

shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]

  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"

  # now we just need a simple call to reposition the legend
  lemon::reposition_legend(p, 'center', panel=names)
}

delwp_palettes = function(name, n, all_palettes = list(delwp_cols), type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,
               continuous = grDevices::colorRampPalette(palette)(n),
               discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}

"#F5F8CC"
