Necrosis Project Figures - Histopathology
================
Derek Lamb
2025-04-14

### Load Packages

This code chunk loads the relevant packages for creating figures

Eventually, I will have my colorblind palette as a package, but for now,
here it is manually:

``` r
manual_cb <- c("#000000", "#4E3F92", "#D5005E", "#28A2D2", "#999999", "#009B77", "#CC49A7", "#00739E", "#862AAF", "#666666")
```

### Define Functions

The following code chunk defines functions for loading ulceration data
and creating histograms.

``` r
load_ulcer_data <- function(path, sheet, range){
  df <- read_xlsx(path = path, sheet = sheet, range = range) |> 
    janitor::clean_names() |> 
    select(-c(1:3)) |> 
    mutate(localization_code = str_extract(localization_code, "(?<=Irradiation Site )[[:upper:]]"),
           localization_code = ifelse(is.na(localization_code),
                                      str_extract(comments, "(?<=Back-up site )\\w{1}(?= used)"),
                                      localization_code),
           wound_id = paste0(subject_code, "_", localization_code),
           time = as.numeric(str_extract(visit_code, "(?<=Day )\\d+"))
           ) |> 
    rename(ulcer_area = x2d_ulcer_area)
  
  
  return(df)
}

make_histo_hist <- function(var, bin_width = 0.1, x_lab = NULL){
  plot <- df_histopath |> 
    ggplot(aes(x = !!sym(var))) + 
    geom_histogram(binwidth = bin_width, color = "black", alpha = 0.6) +
    theme_bw() + 
    labs(x = x_lab, y = "Frequency")
  
  ggsave(filename = paste0("images/", var, "_histogram.png"), plot = plot)
}
```

### Load data

Here we load in the data, combining the ulceration data across sheets.

``` r
# 2D Ulcer area
df_ulcer <- rbind(
  load_ulcer_data("data/ulcer_area_and_histopath.xlsx", 
                  sheet = "2D Ulcer Area 987", range = "A1:J219"),
  load_ulcer_data("data/ulcer_area_and_histopath.xlsx", 
                  sheet = "2D Ulcer Area 903", range = "A1:J234"),
  load_ulcer_data("data/ulcer_area_and_histopath.xlsx", 
                  sheet = "2 D Ulcer Area 745", range = "A1:J267"),
  load_ulcer_data("data/ulcer_area_and_histopath.xlsx", 
                  sheet = "2D Ulcer Area 384", range = "A1:J232")
) |> 
  group_by(time, wound_id) |> 
  summarise(ulcer_area = sum(ulcer_area), 
            animal_id = mean(subject_code), 
            .groups="drop")

# Histopath data
df_histopath <- read_xlsx("data/ulcer_area_and_histopath.xlsx", 
                          sheet = "Histopath meas.", range = "A1:L36") |> 
  janitor::clean_names() |> 
  drop_na(animal_number) |> 
  mutate(wound_id = paste0(animal_number, site)) 
```

### Wound area over time

I will present the ulcer area over time in two ways, first summarised
into a line graph, and additionally as a spaghetti plot.

``` r
# Summarise into df
df_summary <- df_ulcer |> 
  group_by(time) |> 
  summarise(avg_area = mean(ulcer_area), sd_area = sd(ulcer_area),
            .groups="drop")

# Line graph summary
ulcer_plot_line <- df_summary |> 
  ggplot(aes(x = time, y = avg_area)) + 
  geom_line() +
  geom_ribbon(aes(ymin = pmax(avg_area - sd_area, 0), ymax = avg_area + sd_area),
              alpha = 0.4, color = "#000000AA", linetype = "dashed")+
  theme_bw() + 
  labs(x = "Time (days)", y = parse(text = "Average~Ulcer~Area~(mm^2)")) +
  scale_x_continuous(breaks = 0:6*20)

ggsave("images/ulcer_vs_time_line_graph.png", plot = ulcer_plot_line)
```

    ## Saving 7 x 5 in image

``` r
# Overlayed lines
overlay_ulcer_plot_line <- df_ulcer |> 
  mutate(animal_id = as.factor(animal_id)) |> 
  group_by(animal_id, time) |> 
  summarise(avg_area = mean(ulcer_area), sd_area = sd(ulcer_area),
            .groups="drop") |> 
  ggplot(aes(x = time, y = avg_area, color = animal_id)) + 
  geom_line(alpha = 0.8, linetype = "dashed") +
  geom_line(data = df_summary, aes(y = avg_area, color = NULL), size = 0.8) +
  theme_bw() + 
  scale_color_manual(values = manual_cb[2:5]) +
  labs(x = "Time (days)", 
       y = parse(text = "Average~Ulcer~Area~(mm^2)"), 
       color = "Animal ID") +
  scale_x_continuous(breaks = 0:6*20)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
ggsave("images/overlayed_ulcer_vs_time_line_graph.png", plot = overlay_ulcer_plot_line)
```

    ## Saving 7 x 5 in image

``` r
# Spaghetti Plot (WIP)
df_ulcer |> 
  ggplot(aes(x = time, y = ulcer_area, group = wound_id, color = wound_id)) + 
  geom_point() + 
  geom_line() +
  scale_color_manual(values = rep(manual_cb, 4)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Time (days)", y = "Ulcer Area (units??)")
```

![](umb_cri_model_histopathology_figures_files/figure-gfm/wound%20area%20over%20times-1.png)<!-- -->

``` r
library(ggridges)

df_ulcer |> 
  ggplot(aes(y = as.factor(time), x = ulcer_area)) +
  geom_density_ridges(stat = "binline") +
  coord_flip() +
  theme_bw() +
  labs(x = "Time (days)", y = "Ulcer Area")
```

    ## `stat_binline()` using `bins = 30`. Pick better value with `binwidth`.

![](umb_cri_model_histopathology_figures_files/figure-gfm/attempt%20ridge%20plot-1.png)<!-- -->

### Histopathology histograms

In the following code chunk, I apply the previously defined
`make_histo_hist()` function to generate and save histograms of the
relevant histopathology variables.

``` r
make_histo_hist("fibrosis_granulation_tissue_area_mm2", 50,
                parse(text = "Fibrosis/Granulation~Tissue~Area~(mm^2)"))
```

    ## Saving 7 x 5 in image

``` r
make_histo_hist("deep_fibrosis_wound_depth_mm", 2,
                "Deep Fibrosis Wound Depth (mm)")
```

    ## Saving 7 x 5 in image

``` r
make_histo_hist("deep_fibrosis_wound_length_mm", 4,
                "Deep Fibrosis Wound Length (mm)")
```

    ## Saving 7 x 5 in image

``` r
make_histo_hist("epidermal_wound_length_mm", 4,
                "Epidermal Wound Length (mm)")
```

    ## Saving 7 x 5 in image

``` r
make_histo_hist("percent_re_epithelialization_percent", 5,
                "Re-epithelialization (%)")
```

    ## Saving 7 x 5 in image
