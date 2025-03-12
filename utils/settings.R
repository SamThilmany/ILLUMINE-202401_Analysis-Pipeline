library(showtext, quietly = TRUE)

# Read the .Renviron
readRenviron("../.Renviron")

# Parallelization
detected_cores <- parallel::detectCores()
default_cores <- if (detected_cores > 1) detected_cores - 1 else 1

# Graphics
default_dpi <- 600
default_width <- 111
default_height <- 90
default_unit <- "mm"
default_font_size <- 16
default_color <- "#231F20"
default_font_family <- "Roboto"
showtext_opts(dpi = default_dpi)
showtext_auto(enable = TRUE)

# Randomization
set.seed(12345)
random_subset <- 1000

# Statistics
p_value_threshold <- 0.05
log2_fc_threshold <- log2(1.5)

# Diseases of Interest
#
#
# MeSH Tree F03.600 – mood disorders; MeSH Unique ID: D019964
#
# MeSH F03.600.150 – affective disorders, psychotic
# MeSH F03.600.150.150 – bipolar disorder
# MeSH F03.600.150.150.300 – cyclothymic disorder
# MeSH F03.600.300 – depressive disorder
# MeSH F03.600.300.350 – depression, postpartum
# MeSH F03.600.300.375 – depressive disorder, major
# MeSH F03.600.300.400 – dysthymic disorder
# MeSH F03.600.300.700 – seasonal affective disorder
#
#
# MeSH Tree F03.870 – sleep disorders; MeSH Unique ID: D012893
#
# MeSH F03.870.400 – dyssomnias
# MeSH F03.870.400.099 – sleep deprivation
# MeSH F03.870.400.200 – sleep disorders, circadian rhythm
# MeSH F03.870.400.200.500 – jet lag syndrome
# MeSH F03.870.400.800 – sleep disorders, intrinsic
# MeSH F03.870.400.800.200 – disorders of excessive somnolence
# MeSH F03.870.400.800.200.400 – hypersomnolence, idiopathic
# MeSH F03.870.400.800.200.500 – kleine-levin syndrome
# MeSH F03.870.400.800.200.750 – narcolepsy
# MeSH F03.870.400.800.200.750.500 – cataplexy
# MeSH F03.870.400.800.700 – restless legs syndrome
# MeSH F03.870.400.800.800 – sleep initiation and maintenance disorders
# MeSH F03.870.664 – parasomnias
# MeSH F03.870.664.627 – nocturnal paroxysmal dystonia
# MeSH F03.870.664.633 – rem sleep parasomnias
# MeSH F03.870.664.633.700 – rem sleep behavior disorder
# MeSH F03.870.664.633.800 – sleep paralysis
# MeSH F03.870.664.634 – restless legs syndrome
# MeSH F03.870.664.635 – sleep arousal disorders
# MeSH F03.870.664.635.600 – night terrors
# MeSH F03.870.664.635.700 – somnambulism
# MeSH F03.870.664.637 – sleep bruxism
# MeSH F03.870.664.700 – sleep-wake transition disorders
#
#
# MeSH Tree F03.080 – anxiety disorders; MeSH Unique ID: D001008
#
# MeSH F03.080.100 – agoraphobia
# MeSH F03.080.500 – neurocirculatory asthenia
# MeSH F03.080.600 – obsessive-compulsive disorder
# MeSH F03.080.700 – panic disorder
# MeSH F03.080.725 – phobic disorders
# MeSH F03.080.931 – stress disorders, traumatic
# MeSH F03.080.931.249 – combat disorders
# MeSH F03.080.931.374 – stress disorders, traumatic, acute
# MeSH F03.080.931.500 – stress disorders, post-traumatic
diseases_of_interest <- c("D019964", "D012893", "D001008")

# ggplot Theme
theme_set(
  theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      text = element_text(size = default_font_size, family = default_font_family, color = default_color),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
)
update_geom_defaults("text", list(size = 3))