Necrosis Project Figures
================
Derek Lamb
2025-04-08

### Load Packages

This code chunk loads the relevant packages for creating figures

Eventually, I will have my colorblind palette as a package, but for now,
here it is manually:

``` r
manual_cb <- c("#000000", "#4E3F92", "#D5005E", "#28A2D2", "#999999", "#009B77", "#CC49A7",   "#00739E", "#862AAF", "#666666")
```

### Load Data

In this code chunk, I read in the data, convert it to long format, and
standardize the nomenclature.

``` r
df_injury <- read_xlsx("data/Control Group Data_Model Ver.xlsx",
                    sheet="Clinical Onset of Injury Pheno", range="A1:K36") |> 
  janitor::clean_names() |> 
  drop_na(animal_id) |> 
  rename(erythema_start = erythema_day_7, desquamation_start = desquamation, 
         desquamation_end = desquam_end, ulceration_start = ulceration, ulceration_end = ulc_end,
         necrosis_start = necrosis, necrosis_end = necro_end) |> 
  mutate(wound_id = paste0(animal_id, "_", skin_site)) |> 
  select(wound_id, animal_id, skin_site, sex, ends_with("_start"), ends_with("_end")) |> 
  pivot_longer(erythema_start:necrosis_start,
               names_to = "wound_type", names_pattern = "(.*)_start", values_to = "wound_start") |> 
  pivot_longer(erythema_end:necrosis_end,
               names_to = "wound_type2", names_pattern = "(.*)_end", values_to = "wound_end") |> 
  filter(wound_type == wound_type2) |> 
  select(-wound_type2) |> 
  mutate(
    wound_resolved = ifelse(wound_end == 120, 0, 1),
    wound_type = str_to_sentence(wound_type),
    wound_type = fct(wound_type, levels = c("Necrosis", "Ulceration", 
                                            "Desquamation", "Erythema")))
```

## Making figures

The code below creates potential figures.

### Injury Duration

#### Floating Bar Chart

First, we look at a floating bar chart where the start and end of each
bar are the mean time, with error bars denoting standard deviation
(**Note:** Previous versions of this plot showed SEM rather thand
standard deviation).

``` r
floating_bar_wound_type <- df_injury |> 
  group_by(wound_type) |> 
  summarize(avg_start = mean(wound_start, na.rm=T), 
            sem_start = sd(wound_start, na.rm=T),
            avg_end = mean(wound_end, na.rm=T),
            sem_end = sd(wound_end, na.rm=T)
            ) |> 
  ggplot() +
  geom_rect(aes(xmin = avg_start, xmax = avg_end, 
                ymin = as.numeric(wound_type) - 0.2, 
                ymax = as.numeric(wound_type) + 0.2),
                fill = "#2A45FF", 
            color = "black", alpha=0.8) +
  geom_errorbar(aes(xmin = avg_start - sem_start, xmax = avg_start, y = wound_type), width = 0.2)+
  geom_errorbar(aes(xmin = avg_end, xmax = avg_end + sem_end, y = wound_type), width = 0.2) +
  theme_bw() +
  labs(
    x = "Days after radiation",
    y = "Injury Phenotype"
  ) +
  theme(legend.position = "none")

ggsave(plot = floating_bar_wound_type, filename = "images/injury_pheno_floating_bar_start_and_end.png")
```

    ## Saving 7 x 5 in image

#### Multi boxplot

Alternatively, we can represent the wound start and end time as separate
boxplots, organized by phenotype.

``` r
boxplot_wound_type <- df_injury |> 
  group_by(wound_type) |> 
  ggplot() +
  geom_boxplot(aes(x = wound_start, y = wound_type), fill = "#2A45FF", alpha = 0.8) +
  geom_boxplot(aes(x = wound_end, y = wound_type), fill = "#FF2A45", alpha = 0.8) +
  theme_bw() +
  labs(
    x = "Days after radiation",
    y = "Injury Phenotype"
  ) +
  theme(legend.position = "none") 

ggsave(plot = boxplot_wound_type, filename = "images/injury_pheno_boxplot_start_and_end.png")
```

    ## Saving 7 x 5 in image

### Survival Analysis

Additionally, I will look at the time to wound healing. Animals were
right-censored at 120 days.

#### KM Curve for Ulceration

Below I will generate the KM curve for ulceration with a confidence
interval and without (complementary log-log).

``` r
# With CI
ulceration_km_with_ci <- df_injury |> 
  filter(wound_type == "Ulceration") |> 
  mutate(days_to_closure = wound_end - wound_start) |> 
  survfit(Surv(days_to_closure, wound_resolved) ~ 1, data = _, conf.type = "log-log") |> 
  ggsurvfit(type = "risk") +
  add_confidence_interval(type = "lines", alpha = 0.6, size = 0.4, linetype = "dashed") +
  add_confidence_interval(alpha = 0.1) + 
  labs(
    x = "Days from first ulceration",
    y = "Cumulative incidence of wound closure"
  ) +
  ylim(0, 1)

ggsave(plot = ulceration_km_with_ci, filename = "images/ulceration_kaplan_meier_with_ci.png")
```

    ## Saving 7 x 5 in image

``` r
# No CI
ulceration_km_no_ci <- df_injury |> 
  filter(wound_type == "Ulceration") |> 
  mutate(days_to_closure = wound_end - wound_start) |> 
  survfit(Surv(days_to_closure, wound_resolved) ~ 1, data = _, conf.type = "log-log") |> 
  ggsurvfit(type = "risk") +
  labs(
    x = "Days from first ulceration",
    y = "Cumulative incidence of wound closure"
  ) +
  ylim(0, 1)

ggsave(plot = ulceration_km_no_ci, filename = "images/ulceration_kaplan_meier_no_ci.png")
```

    ## Saving 7 x 5 in image

This plot shows the overlayed KM curves for all phenotypes

``` r
four_pheno_km <- df_injury |> 
  mutate(days_to_closure = wound_end - wound_start) |> 
  survfit(Surv(days_to_closure, wound_resolved) ~ wound_type, data = _, conf.type = "log-log") |> 
  ggsurvfit(type = "risk") +
  scale_color_manual(values = manual_cb, labels = c("Necrosis", "Ulceration", 
                                                    "Desquamation", "Erythema")) +
  labs(
    x = "Time from first occurence (days)",
    y = "Cumulative incidence of resolution",
    color = "Injury Phenotype"
  ) +
  ylim(0, 1) +
  theme(legend.position = "right")

ggsave(plot = four_pheno_km, filename = "images/combined_kaplan_meier_four_phenos.png")
```

    ## Saving 7 x 5 in image

Additionally, the following plot is the same, but without the erythema
phenotype.

``` r
three_pheno_km <- df_injury |> 
  filter(wound_type != "Erythema") |> 
  mutate(days_to_closure = wound_end - wound_start) |> 
  survfit(Surv(days_to_closure, wound_resolved) ~ wound_type, data = _, conf.type = "log-log") |> 
  ggsurvfit(type = "risk") +
  scale_color_manual(values = manual_cb, labels = c("Necrosis", "Ulceration", 
                                                    "Desquamation")) +
  labs(
    x = "Time from first occurence (days)",
    y = "Cumulative incidence of resolution",
    color = "Injury Phenotype"
  ) +
  ylim(0, 1) +
  theme(legend.position = "right")

ggsave(plot = three_pheno_km, filename = "images/combined_kaplan_meier_three_phenos.png")
```

    ## Saving 7 x 5 in image

Oops. Doesn’t look like you can actually consider all wound sites
independently.
