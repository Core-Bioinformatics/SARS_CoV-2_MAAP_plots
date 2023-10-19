library(tidyverse)
library(ggplot2)
library(readxl)
library(gridExtra)
library(ggpubr)
library(ggExtra)
library(ggdensity)
library(reshape2)
library(patchwork)
library(ComplexHeatmap)
library(grDevices)
library(stringr)
library(ggcorrplot)
library(Hmisc)
library(psych) 
library(GGally)
library(corrplot)
library(dplyr)
library(purrr)
library(dunn.test)


text_size_ref_1 <- theme(
  text = element_text(size=19),
  axis.title = element_text(size=19),
  axis.text = element_text(size=19),
  plot.title = element_text(size=19),
  legend.text = element_text(size=18), 
  legend.title = element_text(size=18)
)

text_size_ref_2 <- theme(
  text = element_text(size=15),
  axis.title = element_text(size=15),
  axis.text = element_text(size=15),
  plot.title = element_text(size=15),
  legend.text = element_text(size=14), 
  legend.title = element_text(size=14)
)


# Set working directory.
setwd("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Covid-paper-plots/COVID_paper_plots")

# Read in excel file.
df <- read_excel("Combined Cohort Masterfile July23.xlsx")

# Rename df$Cohort 1) Vasculitis → Rituximab 2) Patient → Convalescent.
df$Cohort <- recode(df$Cohort, "Vasculitis" = "RTX", "Patient" = "Convalescent")

#### SECTION 1: Contour plot ####
get_df_for_contour <- function(columns, cohort, timepoints) {
  print(timepoints)
  df_contour<- df %>% select(columns) %>% filter(Cohort == cohort)
  
  # Reshape data and extract necessary information
  df_KD <- df_contour  %>% 
    select(Patient_ID = `Patient ID`, contains("KD")) %>%
    pivot_longer(-Patient_ID, names_to = "Info", values_to = "KD") %>%
    mutate(ID = row_number(),
           B = str_extract(Info, timepoints),
           Variant = str_extract(Info, "WT|Omicron"))
  
  df_AB <- df_contour %>%
    select(Patient_ID = `Patient ID`, contains("[AB]")) %>%
    pivot_longer(-Patient_ID, names_to = "Info", values_to = "AB") %>%
    mutate(ID = row_number())
  
  df_joined<- inner_join(df_KD, df_AB, by = c("ID", "Patient_ID"))
  df_joined$Cohort <- cohort
  
  return(df_joined)
}


# Get contour plot data for different cohorts.
df_contour_RTX_joined <- get_df_for_contour(c("Cohort", "Patient ID", "3A [AB] WT", "3A KD WT", "3A [AB] Omicron", "3A KD Omicron"), "RTX", "2A|3A")
df_contour_HCW_joined <- get_df_for_contour(c("Cohort", "Patient ID", "3A [AB] WT", "3A KD WT", "3A [AB] Omicron", "3A KD Omicron"), "HCW", "3A")
df_contour_patients_joined <- get_df_for_contour(c("Cohort", "Patient ID", "1A [AB] WT", "1A KD WT"), "Convalescent", "1A")

# Combine the two data frames
df_combined_all_cohorts <- rbind(df_contour_RTX_joined, df_contour_HCW_joined, df_contour_patients_joined)

# Rename appropriately.
df_combined_all_cohorts$Cohort <- recode(df_combined_all_cohorts$Cohort, 
                                         "RTX" = "RTX 3A",
                                         "Convalescent" = "Conv 1A",
                                         "HCW" = "HCW 3A")

# Figure: Antibody affinity and concentration of binding sites across different patient cohorts and variants. 
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/7A.pdf", width=10, height=6)
ggplot(df_combined_all_cohorts, aes(x = KD, y = AB)) +
  geom_hdr(aes(fill=Cohort)) +
  geom_line(aes(group = interaction(Patient_ID, B)), color="black", alpha=0.3,  linetype="dashed") +
  geom_point(aes(color = Cohort, shape = B, fill = Variant), size = 3, stroke = 1) +
  scale_shape_manual(values = c(21, 21), guide = "none") +
  scale_fill_manual(values = c("WT" = "white", "Omicron" = "black", "HCW 3A" = "#FFB7C3", "Conv 1A" = "#C0FDC5", "RTX 3A" = "#FFEFBF"),
                    breaks = c("WT", "Omicron"), 
                    guide = guide_legend(override.aes = list(shape = c(NA, NA), colour = c("black", "black")))) +
  scale_color_manual(values = c("HCW 3A" = "#D81B60", "Conv 1A" = "#03A900", "RTX 3A" = "#FFC107")) + 
  theme_classic() +
  labs(x = expression(K[D]~"(M)"), y = "[Ab Binding Sites] (M)",
       color = "Cohort", fill = "Variant") +
  theme(legend.position = "right") + 
  scale_x_reverse() + text_size_ref_2
dev.off()

#### SECTION 2: group comparison plots ####
get_boxplot_ready_df <- function(columns, cohort) {
  df_boxplots <- df %>% select(columns) %>% filter(Cohort == cohort)
  
  # Reshape your data
  df_boxplots <- df_boxplots %>%
    gather(key = "measurement", value = "value", -Cohort, -Sex, -`Patient ID`) %>%
    separate(measurement, into = c("Region", "Type", "Variant"), sep = " ") %>%
    spread(key = Type, value = value) %>% rename("AB" = "[AB]") %>%
    mutate(Multiplication = AB * KD)
  
  # Split Variant into separate columns for AB and KD
  df_boxplots$Variant <- ifelse(grepl("WT", df_boxplots$Variant), "WT", "Omicron")
  
  return(df_boxplots)
}

get_custom_label_from_corr <- function(cor_test, method = "pearson") {
  if(method == "spearman") {
    symbol <- "rho"
  } else {
    symbol <- "R"
  }
  
  # Format based on magnitude and significance
  estimate_fmt <- format(cor_test$estimate, digits = 2, scientific = (abs(cor_test$estimate) < 0.01))
  p_value_fmt <- format(cor_test$p.value, digits = 2, scientific = TRUE)
  
  label <- paste0("italic(", symbol, ") == ", estimate_fmt, 
                  "*',' ~ italic(p) == ", p_value_fmt)
  return(label)
}

compute_corr <- function(dataframe, x_var, y_var, method = "pearson") {
  test_result <- cor.test(dataframe[[x_var]], dataframe[[y_var]], method = method)
  label <- get_custom_label_from_corr(test_result, method = method)
  list(test = test_result, label = label)
}

df_boxplots_conv <- get_boxplot_ready_df(c("Cohort", "Sex", "Patient ID", "1A [AB] WT", "1A KD WT"), "Convalescent")
df_boxplots_RTX <- get_boxplot_ready_df(c("Cohort", "Sex", "Patient ID", "3A [AB] WT", "3A KD WT", "3A [AB] Omicron", "3A KD Omicron"), "RTX")
df_boxplots_HCW <- get_boxplot_ready_df(c("Cohort", "Sex", "Patient ID", "3A [AB] WT", "3A KD WT", "3A [AB] Omicron", "3A KD Omicron"), "HCW")

df_boxplots_joined <- rbind(df_boxplots_conv, df_boxplots_RTX, df_boxplots_HCW)
df_boxplots_joined <- left_join(df_boxplots_joined, df[c("Patient ID", "Age", "Age Bracket (> nearest 5 years)")], by = "Patient ID") %>% rename("Age_bracket" ="Age Bracket (> nearest 5 years)")

df_boxplots_joined$Cohort <- recode(df_boxplots_joined$Cohort, 
                                    "RTX" = "RTX 3A",
                                    "Convalescent" = "Conv 1A",
                                    "HCW" = "HCW 3A")

#### Figure: Sex vs AB * KD violin plot
df_boxplots_joined_sex <- df_boxplots_joined %>% filter(!is.na(Sex))

# No correlation annotation.
df_boxplots_joined_sex$combined_group_label <-paste(df_boxplots_joined_sex$Sex, df_boxplots_joined_sex$Cohort, sep = " - ")

pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/7C.pdf", width=10, height=6)
ggplot(df_boxplots_joined_sex, aes(x = interaction(Sex, Variant, sep = " - "), y = Multiplication, fill = Variant, color = Variant)) +
  geom_violin(alpha = 0.4, position = position_dodge(width = 0.5)) +
  geom_jitter(aes(shape = Cohort), size = 2, position = position_jitter(width = 0.1)) +
  theme_classic() +
  scale_color_manual(values = c("WT" = "#8400A2", "Omicron" = "#FFC107", "Male"= "#1f77b4", "Female" = "#ff7f0e"), name="Variant") +
  scale_fill_manual(values = c("WT" = "#8400A2", "Omicron" = "#FFC107"), name = "Variant") +
  labs(x = "Sex and Variant", y = expression("[Ab]"%*% K[D]), fill = "Variant", color = "Variant") +
  scale_shape_discrete(name = "Cohort") +
  scale_x_discrete(labels = function(x) gsub(" - ", "\n", x))  + text_size_ref_1
dev.off()


# Figure: Age vs AB * KD. 
get_age_mult_plot <- function(df, x, y, x_corr = "Multiplication", y_corr = "Age_bracket", x_cord = 0.8e-14, y_cords = c(9,5), y_label = "Age") {
  results <- list()
  
  for (method in c("pearson", "spearman")) {
    age_mult_omicron_corr <- compute_corr(df_boxplots_joined[df_boxplots_joined$Variant == "Omicron",], x_corr, y_corr, method = method)
    age_mult_WT_corr <- compute_corr(df_boxplots_joined[df_boxplots_joined$Variant == "WT",], x_corr, y_corr, method = method)
    
    plot_age_mult <- ggplot(df, aes( x = x, y = y, fill = Variant, color= Variant)) +
      geom_point(aes(shape=Cohort), size=3) +
      geom_smooth(method=lm) +
      scale_color_manual(values = c("WT" = "#8400A2", "Omicron" = "#FFC107")) +
      scale_fill_manual(values = c("WT" = "#8400A2", "Omicron" = "#FFC107")) +
      theme_classic() +
      labs(y = y_label, x = expression("[Ab]"%*% K[D]), fill = "Variant") +
      annotate("text", x = x_cord, y = y_cords[1], label = age_mult_omicron_corr$label, color="#FFC107", parse = TRUE, size = 5) + 
      annotate("text", x = x_cord, y = y_cords[2], label = age_mult_WT_corr$label, color="#8400A2", parse = TRUE, size = 5) + text_size_ref_1 +
      scale_x_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
    
    # Assign the plots based on the method.
    results[[paste0("plot_age_mult_", method)]] <- plot_age_mult
  }
  
  return(results)
}


pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/7B.pdf", width=10, height=6)
df_boxplots_exact_age <- df_boxplots_joined %>% filter(!is.na(Age))
get_age_mult_plot(df_boxplots_exact_age, df_boxplots_exact_age$Multiplication, df_boxplots_exact_age$Age, y_corr ="Age", x_cord = 3.4e-15)

# All approximate
get_age_mult_plot(df_boxplots_joined, df_boxplots_joined$Multiplication, df_boxplots_joined$Age_bracket, y_label = "Age (nearest 5 years)")

# Exact for HCW and RTX, approximate for Convalescent.
df_boxplots_joined_mix_exact_approximate <- df_boxplots_joined %>%
  mutate(Age = ifelse(Cohort == "Conv 1A", Age_bracket, Age))
get_age_mult_plot(df_boxplots_joined_mix_exact_approximate, df_boxplots_joined_mix_exact_approximate$Multiplication, df_boxplots_joined_mix_exact_approximate$Age, y_corr ="Age", y_label = "Age (nearest 5 years for Conv)")
dev.off()



# Figure: disease severity vs Kd x [Ab]; Only involves the convalescent Covid cohort at 1 month post infection.
df_boxplots_Vasc <- df[,c("Cohort", "Patient ID",  "Disease Severity", "1A [AB] WT", "1A KD WT")] %>% filter(Cohort == "Convalescent") 
df_boxplots_Vasc$Multiplication <- df_boxplots_Vasc$`1A [AB] WT` * df_boxplots_Vasc$`1A KD WT`

# Summarize the number of samples for each Disease Severity where Multplication is not NA
df_counts_disease_severity <- df_boxplots_Vasc %>%
  filter(!is.na(Multiplication)) %>%
  group_by(`Disease Severity`) %>%
  summarise(n = n(), .groups = 'drop')

pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/7D.pdf", width=12, height=7)
ggplot(df_boxplots_Vasc, aes(y = `Disease Severity`, x = log10(Multiplication), fill = `Disease Severity`)) +
  geom_boxplot(alpha = 0.4, outlier.size = 0, width = 0.5, position = position_dodge(1)) +
  geom_point(size = 3) +
  geom_text(data = df_counts_disease_severity, aes(y = `Disease Severity`, x = Inf, label = paste(" n =", n)), 
            hjust = "inward", vjust = 1.5, color = "black", size = 5) +
  theme_classic() +
  labs(y = "Disease Severity", x = expression(log[10]("[Ab]" %*% K[D]))) + text_size_ref_2
dev.off()



#### STATISTICS COMPARISONS ####
# Mann Whitney U test (Wilcoxon Rank Sum test)
# Compare the distributions of Multiplication (AB * KD) for WT and Omicron variant.
# p-value = 0.4032 -> not significant difference between the distribution of multiplication in the different variants.
result_wilcox_test <- wilcox.test(Multiplication ~ Variant, data = df_boxplots_joined)
result_wilcox_test

# Is there a difference in Kd x [Ab] according to sex for each variant?
model_sex_variant_product <- lm(Multiplication ~ Sex * Variant, data = df_boxplots_joined)

# Print a summary of the model
summary(model_sex_variant_product)

#  Is there a relationship between age and Kd x [Ab]?
lm_age_product <- lm(Multiplication ~ Age, data = df_boxplots_joined)
summary(lm_age_product)

#  Is there a relationship between age bracket and Kd x [Ab]?
lm_age_bracket_product <- lm(Multiplication ~ Age_bracket, data = df_boxplots_joined)
summary(lm_age_bracket_product)

#### SECTION 3: How different serum antibody assessments relate to each other ####
create_df_variant <- function(df, cohort, timepoint, variant) {
  df %>%
    filter(Cohort == cohort) %>%
    select(
      Cohort, 
      `Patient ID`,
      ND50 = paste0(timepoint, " ", variant, " ND50"), 
      AB = paste0(timepoint, " [AB] ", variant), 
      KD = paste0(timepoint, " KD ", variant),
      LMNX = paste0(timepoint, " RBD Lmnx MFI")
    ) %>%
    mutate(Variant = variant,
           Timepoint = timepoint)
}

# Filter the datasets so that we get only data in each cohort -> perform the test on those.
seperate_df_in_cohorts <- function(df) {
  convalescent <- df %>% filter(Cohort == "Conv 1A & 1B")
  HCW <- df %>% filter(Cohort == "HCW 3A")
  RTX <- df %>% filter(Cohort == "RTX 2A & 3A")
  
  return(list(Convalescent = convalescent, HCW = HCW, RTX = RTX))
}

get_together_plot <- function(df, x, y, x_label, y_label, cords, cor_test, method = "pearson") {
  
  plot_together <- ggplot(df, aes(x= {{ x }}, y= {{ y }}, color= Cohort, shape=Cohort)) +
    geom_point() + 
    theme_classic() + 
    geom_smooth(aes(group = Variant), method="lm", color="black") + 
    labs(x = x_label, y = y_label) + 
    annotate("text", x = cords[1], y = cords[2], label = cor_test$label, parse = TRUE, size = 5)
  
  return(plot_together)
}

get_per_cohort_plot <- function(df, cor_conv, cor_HCW, cor_RTX, x, y, x_label, y_label, x_cord, y_cords, method = "pearson") {
  
  plot_per_cohort <- ggplot(df, aes(x= {{ x }}, y= {{ y }}, color=Cohort, shape=Cohort)) +
    geom_point() + theme_classic() + geom_smooth(aes(group = Cohort, fill = Cohort), method="lm") + 
    labs(x = x_label, y = y_label) + 
    annotate("text", x = x_cord, y = y_cords[1], label = cor_conv$label, color="#F8766D", parse = TRUE, size = 5) + 
    annotate("text", x = x_cord, y = y_cords[2], label = cor_HCW$label, color="#00BA38", parse = TRUE, size = 5) + 
    annotate("text", x = x_cord, y = y_cords[3], label = cor_RTX$label, color="#619CFF", parse = TRUE, size = 5)  +
    guides(color = guide_legend(override.aes = list(linetype = 0)),
           shape = guide_legend(override.aes = list(linetype = 0)))
  
  return(plot_per_cohort)
}


# Get correlation labels per Cohort.
get_corr_labels_each_cohort <- function(df_conv, df_HCW, df_RTX, x_as_string, y_as_string, method) {
  corr_conv <- compute_corr(df_conv, x_as_string, y_as_string, method = method)
  corr_HCW <- compute_corr(df_HCW, x_as_string, y_as_string, method = method)
  corr_RTX <- compute_corr(df_RTX, x_as_string, y_as_string, method = method)
  
  return(list(conv = corr_conv, HCW = corr_HCW, RTX = corr_RTX))
}

together_and_per_cohort_plots <- function(df, x, y, x_as_string, y_as_string, x_label, y_label, cords_together, df_conv, df_HCW, df_RTX, x_cord_cohort, y_cords_cohort) {
  
  results <- list() 
  
  for (method in c("pearson", "spearman")) {
    corr_together <- compute_corr(df, x_as_string, y_as_string, method = method)
    
    # Together plot
    plot_together <- get_together_plot(df, {{ x }}, {{ y }}, x_label, y_label, cords_together, corr_together)
    
    # Variation - per Cohort.
    corr_labels_cohorts <- get_corr_labels_each_cohort(df_conv, df_HCW, df_RTX, x_as_string, y_as_string, method)
    plot_per_cohort <- get_per_cohort_plot(df, corr_labels_cohorts$conv, corr_labels_cohorts$HCW, corr_labels_cohorts$RTX, {{ x }}, {{ y }}, x_label, y_label, x_cord_cohort, y_cords_cohort, method = method)
    
    # Assign the plots based on the method.
    results[[paste0("plot_together_", method)]] <- plot_together
    results[[paste0("plot_per_cohort_", method)]] <- plot_per_cohort
    
  }
  
  return(results)
}


remove_outlier_points <- function(df, n, feature, x, y, x_as_string, y_as_string, x_label, y_label, x_cord_cohort, y_cords_cohort, decreasing = TRUE, method = "pearson", x_trans = identity, y_trans = identity) {
  
  results <- list()
  
  # Remove n largest/smallest outliers
  sorted_values <- sort(df[[feature]], decreasing = TRUE)
  threshold <- sorted_values[n]
  df_no_outliers <- df[df[feature] < threshold, ]
  
  # Separate dataset without the outliers in cohorts.
  df_no_outliers_cohorts <- seperate_df_in_cohorts(df_no_outliers)
  df_no_outliers_conv <- df_no_outliers_cohorts$Convalescent
  df_no_outliers_HCW <- df_no_outliers_cohorts$HCW
  df_no_outliers_RTX <- df_no_outliers_cohorts$RTX
  
  for (method in c("pearson", "spearman")) {
    corr_labels_cohorts <- get_corr_labels_each_cohort(df_no_outliers_conv, df_no_outliers_HCW, df_no_outliers_RTX, x_as_string, y_as_string, method)
    
    # Apply transformations
    df_no_outliers$x_transformed <- x_trans(df_no_outliers[[x]])
    df_no_outliers$y_transformed <- y_trans(df_no_outliers[[y]])
    df_no_outliers <- df_no_outliers %>% filter(!is.na(Cohort))
    
    plot_per_cohort_no_outlier <- ggplot(df_no_outliers, aes(x = x_transformed, y = y_transformed, color = Cohort, shape = Cohort)) +
      geom_point() + theme_classic() + geom_smooth(aes(group = Cohort, fill=Cohort), method="lm") + 
      labs(x = x_label, y = y_label) + 
      annotate("text", x = x_cord_cohort, y = y_cords_cohort[1], label = corr_labels_cohorts$conv$label, color="#F8766D", parse = TRUE) + 
      annotate("text", x = x_cord_cohort, y = y_cords_cohort[2], label = corr_labels_cohorts$HCW$label, color="#00BA38", parse = TRUE) + 
      annotate("text", x = x_cord_cohort, y = y_cords_cohort[3], label = corr_labels_cohorts$RTX$label, color="#619CFF", parse = TRUE) +
      guides(color = guide_legend(override.aes = list(linetype = 0)),
             shape = guide_legend(override.aes = list(linetype = 0)))
    
    # Assign the plots based on the method.
    results[[paste0("plot_no_outlier_", method)]] <- plot_per_cohort_no_outlier
  }
  
  # Store df without outliers
  results[["df_no_outliers"]] <- df_no_outliers
  
  return(results)
}



# Get all points (2A and 3A timepoints for both WT and Omicron variant) for RTX cohort.
combinations_RTX  <- expand.grid(timepoint = c("2A", "3A"), variant = c("WT", "Omicron"))
serum_antibody_relation_RTX <- purrr::pmap_dfr(combinations_RTX, 
                                               ~create_df_variant(df, "RTX", ..1, ..2))

# Get all points (3A timepoint for both WT and Omicron variant) for HCW cohort.
serum_antibody_relation_HCW <- purrr::pmap_dfr(list(c("3A"), c("WT", "Omicron")), 
                                               ~create_df_variant(df, "HCW", ..1, ..2))

# Get all points (1A and 1B timepoint for WT variant) for patients cohort.
serum_antibody_relation_patient_1A <- create_df_variant(df, "Convalescent", "1A", "WT")
serum_antibody_relation_patient_1B <- create_df_variant(df, "Convalescent", "1B", "WT")
serum_antibody_relation_patient <- rbind(serum_antibody_relation_patient_1A, serum_antibody_relation_patient_1B)

serum_antibody_relation_all_combined <- rbind(serum_antibody_relation_RTX, serum_antibody_relation_HCW, serum_antibody_relation_patient)

# Keep only WT variant. UPDATE - change names to showcase timepoints.
serum_antibody_relation_all_combined_only_WT <-
  serum_antibody_relation_all_combined %>% filter(Variant == "WT") %>% mutate(Multiplication = AB * KD, min_logKD = -log10(KD), Cohort = recode(Cohort,
                                                                                                                                                "Convalescent" = "Conv 1A & 1B",
                                                                                                                                                "HCW" = "HCW 3A",
                                                                                                                                                "RTX" = "RTX 2A & 3A"))

df_serum_antibody_relation_cohorts <- seperate_df_in_cohorts(serum_antibody_relation_all_combined_only_WT)
df_serum_antibody_relation_WT_only_Convalescent <- df_serum_antibody_relation_cohorts$Convalescent %>% mutate("-logKD" = -log10(KD))
df_serum_antibody_relation_WT_only_HCW <- df_serum_antibody_relation_cohorts$HCW %>% mutate("-logKD" = -log10(KD))
df_serum_antibody_relation_WT_only_RTX <- df_serum_antibody_relation_cohorts$RTX %>% mutate("-logKD" = -log10(KD))


#### Figure: AB vs ND50
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/8B.pdf", width= 10, height=6)
AB_ND50 <-
  together_and_per_cohort_plots(
    df = serum_antibody_relation_all_combined_only_WT,
    x = AB,
    y = log10(ND50),
    x_as_string = "AB",
    y_as_string = "ND50",
    x_label = "[Ab Binding Sites] (M)",
    y_label =expression(log[10] ~ "(ND50)"),
    cords_together = c(5.0e-07, 4.5),
    df_conv = df_serum_antibody_relation_WT_only_Convalescent,
    df_HCW = df_serum_antibody_relation_WT_only_HCW,
    df_RTX = df_serum_antibody_relation_WT_only_RTX,
    x_cord_cohort = 1.5e-06,
    y_cords_cohort = c(2.6, 2.4, 2.2)
  )

# Make text larger refinements.
AB_ND50 <- lapply(AB_ND50, function(plot) {
  plot + text_size_ref_1
})

AB_ND50
dev.off()




####  Figure: KD vs ND50
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/supp_5A.pdf", width=10, height=6)
KD_ND50 <-
  together_and_per_cohort_plots(
    df = serum_antibody_relation_all_combined_only_WT,
    x = min_logKD,
    y = log10(ND50),
    x_as_string = "min_logKD",
    y_as_string = "ND50",
    x_label = expression(-log[10](K[D]~"(M)")),
    y_label = expression(log[10]~"(ND50)"),
    cords_together = c(8, 1),
    df_conv = df_serum_antibody_relation_WT_only_Convalescent,
    df_HCW = df_serum_antibody_relation_WT_only_HCW,
    df_RTX = df_serum_antibody_relation_WT_only_RTX,
    x_cord_cohort = 8,
    y_cords_cohort =  c(4.5, 4.3, 4.1)
  )


KD_ND50 <- lapply(KD_ND50, function(plot) {
  plot + text_size_ref_2
})

KD_ND50

dev.off()

#### Figure: LMNX vs ND50
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/8C.pdf", width=10, height=6)
LMNX_ND50 <-
  together_and_per_cohort_plots(
    df = serum_antibody_relation_all_combined_only_WT,
    x = LMNX,
    y = log10(ND50),
    x_as_string = "LMNX",
    y_as_string = "ND50",
    x_label = "Lmnx RBD MFI",
    y_label = expression(log[10]~"(ND50)"),
    cords_together = c(10000, 4),
    df_conv = df_serum_antibody_relation_WT_only_Convalescent,
    df_HCW = df_serum_antibody_relation_WT_only_HCW,
    df_RTX = df_serum_antibody_relation_WT_only_RTX,
    x_cord_cohort = 10000,
    y_cords_cohort = c(4.5, 4.3, 4.1)
  )


# Make text larger refinements.
LMNX_ND50 <- lapply(LMNX_ND50, function(plot) {
  plot + text_size_ref_1
})

LMNX_ND50
dev.off()

#### Figure: LMNX vs AB
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/supp_5C.pdf", width=10, height=6)
LMNX_AB <-
  together_and_per_cohort_plots(
    df = serum_antibody_relation_all_combined_only_WT,
    x = LMNX,
    y = AB,
    x_as_string = "LMNX",
    y_as_string = "AB",
    x_label = "Lmnx RBD MFI",
    y_label = "[Ab Binding Sites] (M)",
    cords_together = c(10000,  2e-06),
    df_conv = df_serum_antibody_relation_WT_only_Convalescent,
    df_HCW = df_serum_antibody_relation_WT_only_HCW,
    df_RTX = df_serum_antibody_relation_WT_only_RTX,
    x_cord_cohort = 10000,
    y_cords_cohort = c(2e-06, 1.9e-06, 1.8e-06)
  )

LMNX_AB <- lapply(LMNX_AB, function(plot) {
  plot + text_size_ref_2
})

LMNX_AB

dev.off()

#### Figure: LMNX vs KD
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Refinement_plots/supp_5B.pdf", width=10, height=6)
LMNX_KD <-
  together_and_per_cohort_plots(
    df = serum_antibody_relation_all_combined_only_WT,
    x = min_logKD,
    y = LMNX,
    x_as_string = "min_logKD",
    y_as_string = "LMNX",
    x_label = expression(-log[10](K[D](M))),
    y_label = "Lmnx RBD MFI",
    cords_together = c(7.7,  8000),
    df_conv = df_serum_antibody_relation_WT_only_Convalescent,
    df_HCW = df_serum_antibody_relation_WT_only_HCW,
    df_RTX = df_serum_antibody_relation_WT_only_RTX,
    x_cord_cohort = 8.8,
    y_cords_cohort = c(14000, 12000, 10000)
  )

LMNX_KD <- lapply(LMNX_KD, function(plot) {
  plot + text_size_ref_2
})


LMNX_KD
dev.off()

#### SECTION 4: heatmap ####
data_heatmap <- df %>% select("Cohort", "2A S Lmnx MFI", "3A S Lmnx MFI", "2A WT ND50", "3A WT ND50", "2A Delta ND50", "3A Delta ND50", "2A [AB] WT", "3A [AB] WT")

create_modular_heatmap <- function(name_heatmap, df_selected, colour_first, colour_second, new_column_names, annotation_gp = NULL, annotate_left = F, column_name = name_heatmap) {
  mat = as.matrix(df_selected)
  ha_left = NULL
  
  # Create time point annotation.
  timepoint = str_extract(colnames(mat), "\\d+$")
  ha_top = HeatmapAnnotation(timepoint = timepoint, annotation_name_side = "left", col = list(timepoint = c("1" = "#add8e6", "2" = "#53A5C0", "3" = "#00008b")), annotation_name_gp = annotation_gp)
  
  if (annotate_left) {
    # Create cohort annotation.
    annotation_df <- data.frame(Cohort = data_heatmap$Cohort)
    cohort_colors <- list(Cohort = c("RTX" = "#ED1D24", "HCW" = "#1DB52B", "Convalescent" = "#005EB8"))
    ha_left = rowAnnotation(df = annotation_df, col = cohort_colors)
  }
  
  # Create colour map for heatmap.
  col_fun <- colorRampPalette(c(colour_first, colour_second))
  
  # Change the column names
  colnames(mat) <- new_column_names
  
  return(Heatmap(mat, name = name_heatmap, column_title = column_name,
                 top_annotation = ha_top, left_annotation = ha_left, cluster_rows = F, cluster_columns = F, col = col_fun(50), na_col = "#F5F5F5"))
}

data_heatmap_lmnx <- df %>% select(N_Lmnx_MFI_1 = "1A N Lmnx MFI", S_LMNX_MFI_1 = "1A S Lmnx MFI", RBD_Lmnx_MFI_1 = "1A RBD Lmnx MFI", N_Lmnx_MFI_2 = "2A N Lmnx MFI", S_LMNX_MFI_2 = "2A S Lmnx MFI", RBD_Lmnx_MFI_2 = "2A RBD Lmnx MFI",  N_Lmnx_MFI_3 = "3A N Lmnx MFI", S_LMNX_MFI_3 = "3A S Lmnx MFI", RBD_Lmnx_MFI_3 = "3A RBD Lmnx MFI") %>%
  mutate(across(everything(), ~ log10(.))) # apply log10 transformation.
data_heatmap_KD <- df %>% select(KD_WT_1 = "1A KD WT", KD_WT_2 = "2A KD WT", KD_Omicron_2 = "2A KD Omicron", KD_WT_3 = "3A KD WT", KD_Omicron_3 = "3A KD Omicron") %>%  mutate(across(everything(), ~ log10(.))) # apply log10 transformation.
data_heatmap_AB <- df %>% select(AB_WT_1 = "1A [AB] WT", AB_WT_2 = "2A [AB] WT", AB_Omicron_2 = "2A [AB] Omicron", AB_WT_3 = "3A [AB] WT", AB_Omicron_3 = "3A [AB] Omicron") %>% mutate(across(everything(), ~ log10(.))) # apply log10 transformation.
data_heatmap_ND <- df %>% select(ND_WT_1 = "1A WT ND50", ND_WT_2 = "2A WT ND50",  ND_Omicron_2 = "2A Omicron ND50", ND_WT_3 = "3A WT ND50",  ND_Omicron_3 = "3A Omicron ND50") %>%
  mutate(across(everything(), ~ log10(.))) # apply log10 transformation.


ht_lmnx <- create_modular_heatmap("log10(Luminex MFI)", data_heatmap_lmnx, "white", "#653475", new_column_names = c("N", "S", "RBD","N", "S", "RBD", "N", "S", "RBD"), annotate_left = T, column_name = expression(log[10]~"(Luminex MFI)"))
ht_KD <- create_modular_heatmap("KD", data_heatmap_KD, "#a61900", "white", new_column_names = c("WT", "WT", "Omicron", "WT", "Omicron"), gpar(fontsize = 0), column_name = expression(log[10]~"("~K[D]~"(M)"~")"))
ht_AB <- create_modular_heatmap("[Ab]", data_heatmap_AB, "white", "#00821a", new_column_names = c("WT", "WT", "Omicron", "WT", "Omicron"), gpar(fontsize = 0), column_name = expression(log[10]~"([Ab])"))
ht_ND <- create_modular_heatmap("log10(ND50)", data_heatmap_ND, "white", "#006EB3", new_column_names = c("WT", "WT", "Omicron", "WT", "Omicron"), gpar(fontsize = 0), column_name = expression(log[10]~"(ND50)"))
ht_age <- Heatmap(df$`Age Bracket (> nearest 5 years)`, name = "Age", na_col = "#F5F5F5", col = colorRampPalette(c("#f9ffab", "red"))(50) )
ht_sex  <- Heatmap(df$Sex, name = "Sex", na_col = "#F5F5F5", col = c("Male" = "black", "Female" = "yellow"))
ht_ethnicity <- Heatmap(df$Ethnicity, name = "Ethnicity", na_col = "#F5F5F5", col = c("White British" = "#AADB1A", "White Irish" = "#d95f02", "Black African" = "#00A3E1", "Other Ethnic Group" = "#e7298a", "Other Mixed Background" = "#737373", "Other Black Background" = "#e6ab02", "Other White Background" = "brown"))

# Figure: summary heatmap
pdf("/Users/rafaelkoll/Desktop/4/Uni/Research/second_paper/Covid-paper-plots/COVID_paper_plots/plots_pdf/figure12.pdf", width=14, height=7)
ht_lmnx + ht_KD + ht_AB + ht_ND + ht_age + ht_sex + ht_ethnicity
dev.off()