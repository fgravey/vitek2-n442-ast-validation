library(dplyr)
library(ggplot2)
library(openxlsx)
library(readxl)
library(stringr)
library(tidyr)
library(readr)
library(rstatix)
library(tibble)


# Fonction pour nettoyer et ajuster les CMI hors échelle selon les limites communes
clean_and_adjust_mic_common_limits <- function(vitek_values, bmd_values) {
  vitek_values <- as.character(vitek_values)
  bmd_values <- as.character(bmd_values)
  
  get_limits <- function(values, direction) {
    # Supprimer les NA pour éviter les erreurs avec grepl
    values <- na.omit(values)
    
    if (direction == "lower") {
      if (any(grepl("<", values))) {
        return(min(as.numeric(gsub("<", "", values[grepl("<", values)]))) / 2)
      } else {
        return(NA)
      }
    } else if (direction == "upper") {
      if (any(grepl(">", values))) {
        return(max(as.numeric(gsub(">", "", values[grepl(">", values)]))) * 2)
      } else {
        return(NA)
      }
    }
  }
  
  # Déterminer les limites
  vitek_lower_limit <- get_limits(vitek_values, "lower")
  vitek_upper_limit <- get_limits(vitek_values, "upper")
  bmd_lower_limit <- get_limits(bmd_values, "lower")
  bmd_upper_limit <- get_limits(bmd_values, "upper")
  
  # Déterminer les limites communes
  # Pour la limite inférieure commune
  if (all(is.na(c(vitek_lower_limit, bmd_lower_limit)))) {
    common_lower_limit <- NA
  } else {
    common_lower_limit <- max(c(vitek_lower_limit, bmd_lower_limit), na.rm = TRUE)
  }
  
  # Pour la limite supérieure commune
  if (all(is.na(c(vitek_upper_limit, bmd_upper_limit)))) {
    common_upper_limit <- NA
  } else {
    common_upper_limit <- min(c(vitek_upper_limit, bmd_upper_limit), na.rm = TRUE)
  }
  
  
  # Fonction pour ajuster les valeurs
  adjust_value <- function(x, lower_limit, upper_limit) {
    x <- as.character(x)
    if (grepl("<", x)) {
      num <- as.numeric(gsub("<", "", x))
      if (!is.na(num) && !is.na(lower_limit)) {
        return(max(num / 2, lower_limit))
      } else {
        return(num / 2)
      }
    } else if (grepl(">", x)) {
      num <- as.numeric(gsub(">", "", x))
      if (!is.na(num) && !is.na(upper_limit)) {
        return(min(num * 2, upper_limit))
      } else {
        return(num * 2)
      }
    } else if (toupper(x) %in% c("SFP", "ZIT")) {
      return(NA)
    } else {
      num <- suppressWarnings(as.numeric(x))
      if (!is.na(num) && !is.na(lower_limit) && num < lower_limit) {
        return(lower_limit)
      } else if (!is.na(num) && !is.na(upper_limit) && num > upper_limit) {
        return(upper_limit)
      } else {
        return(num)
      }
    }
  }
  
  vitek_adjusted <- sapply(vitek_values, adjust_value, lower_limit = common_lower_limit, upper_limit = common_upper_limit)
  bmd_adjusted <- sapply(bmd_values, adjust_value, lower_limit = common_lower_limit, upper_limit = common_upper_limit)
  
  return(list(vitek_adjusted = vitek_adjusted, bmd_adjusted = bmd_adjusted))
}

# Fonction pour transformer le tableau d'origine
transform_data_for_cmi <- function(data) {
  # Liste des antibiotiques
  antibiotics <- c("PTZ", "CZD", "CRO", "FEP", "ETP", "IMP", "AKN", "LVX", "NFE", "SXT", "AMX", "AMC", "FIX", "FOX", "FOS", "MEC")
  
  # Initialiser le tableau de sortie avec les colonnes de stratification
  output_data <- data %>% select(Groupe, Especes, Phénotype)
  
  # Pour chaque antibiotique, ajouter les colonnes CMI initiales, corrigées et d'interprétation
  for (antibio in antibiotics) {
    vitek_col <- paste0(antibio, ".VITEK.res")
    bmd_col <- paste0(antibio, ".BMD.res")
    vitek_inter_col <- paste0(antibio, ".VITEK.inter")
    bmd_inter_col <- paste0(antibio, ".BMD.inter")
    
    if (all(c(vitek_col, bmd_col, vitek_inter_col, bmd_inter_col) %in% colnames(data))) {
      # Ajouter les colonnes CMI initiales
      output_data <- output_data %>%
        bind_cols(
          data %>%
            select(!!sym(vitek_col), !!sym(bmd_col)) %>%
            rename(!!paste0(antibio, "_VITEK_initial") := !!sym(vitek_col),
                   !!paste0(antibio, "_BMD_initial") := !!sym(bmd_col))
        )
      
      # Ajouter les colonnes CMI corrigées
      adjusted_values <- clean_and_adjust_mic_common_limits(data[[vitek_col]], data[[bmd_col]])
      output_data <- output_data %>%
        bind_cols(
          tibble(!!paste0(antibio, "_VITEK_corrected") := adjusted_values$vitek_adjusted,
                 !!paste0(antibio, "_BMD_corrected") := adjusted_values$bmd_adjusted)
        )
      
      # Ajouter les colonnes d'interprétation
      output_data <- output_data %>%
        bind_cols(
          data %>%
            select(!!sym(vitek_inter_col), !!sym(bmd_inter_col)) %>%
            rename(!!paste0(antibio, "_VITEK_inter") := !!sym(vitek_inter_col),
                   !!paste0(antibio, "_BMD_inter") := !!sym(bmd_inter_col))
        )
    }
  }
  
  return(output_data)
}

#####################################

# Fonction principale pour analyser les données
analyze_mic_data <- function(data, test_col, reference_col, molecule_name = "Molecule", test_name = "Test", reference_name = "Reference", output_dir = "plots", threshold = 16) {
  initial_nrow <- nrow(data)
  
  # Élimination des lignes contenant le symbole "/"
  data <- data %>%
    filter(
      !grepl("/", !!sym(test_col)) &
        !grepl("/", !!sym(reference_col))
    )
  
  slash_excluded <- initial_nrow - nrow(data)
  if (slash_excluded > 0) {
    message(paste(slash_excluded, "lignes ont été exclues car elles contenaient le symbole '/' dans les colonnes ", test_col, " ou ", reference_col, "."))
  }
  
  initial_nrow_after_slash <- nrow(data)
  
  # Nettoyage et ajustement des valeurs de la colonne test
  data <- data %>%
    mutate(
      !!paste0(test_col, "_num") := clean_and_adjust_mic(!!sym(test_col)),
      !!paste0(reference_col, "_num") := clean_and_adjust_mic(!!sym(reference_col))
    )
  
  # Filtrer les lignes avec des valeurs NA
  na_excluded <- sum(is.na(data[[paste0(test_col, "_num")]]) | is.na(data[[paste0(reference_col, "_num")]]))
  if (na_excluded > 0) {
    message(paste(na_excluded, "lignes ont été exclues car elles contenaient des valeurs non convertibles en numériques dans les colonnes ", test_col, " ou ", reference_col, "."))
  }
  
  data <- data %>%
    filter(!is.na(!!sym(paste0(test_col, "_num"))) & !is.na(!!sym(paste0(reference_col, "_num"))))
  
  if (nrow(data) == 0) {
    warning(paste("Aucune donnée valide pour", test_col, "et", reference_col))
    return(NULL)
  }
  
  # Détermination des valeurs minimale et maximale de la colonne test_num
  min_test <- min(data[[paste0(test_col, "_num")]], na.rm = TRUE)
  max_test <- max(data[[paste0(test_col, "_num")]], na.rm = TRUE)
  
  # Ajustement des valeurs de la colonne reference_num en fonction des valeurs minimale et maximale de test_num
  data <- data %>%
    mutate(
      !!paste0(reference_col, "_num_adjusted") := case_when(
        !!sym(paste0(reference_col, "_num")) < min_test ~ min_test,
        !!sym(paste0(reference_col, "_num")) > max_test ~ max_test,
        TRUE ~ !!sym(paste0(reference_col, "_num"))
      )
    )
  
  # Notification des ajustements
  adjusted_low <- sum(data[[paste0(reference_col, "_num")]] < min_test, na.rm = TRUE)
  adjusted_high <- sum(data[[paste0(reference_col, "_num")]] > max_test, na.rm = TRUE)
  
  if (adjusted_low > 0) {
    message(paste(adjusted_low, "valeurs de ", reference_col, " ont été ajustées à la valeur minimale de ", min_test, " car elles étaient inférieures à cette valeur."))
  }
  
  if (adjusted_high > 0) {
    message(paste(adjusted_high, "valeurs de ", reference_col, " ont été ajustées à la valeur maximale de ", max_test, " car elles étaient supérieures à cette valeur."))
  }
  
  # Calcul de l'EA (Essential Agreement)
  data <- data %>%
    mutate(Is_EA = abs(log2(!!sym(paste0(test_col, "_num"))) - log2(!!sym(paste0(reference_col, "_num_adjusted")))) <= 1)
  
  # Calcul de la contribution au biais
  data <- data %>%
    mutate(
      bias_contribution = case_when(
        !!sym(paste0(test_col, "_num")) > !!sym(paste0(reference_col, "_num_adjusted")) ~ "higher",
        !!sym(paste0(test_col, "_num")) < !!sym(paste0(reference_col, "_num_adjusted")) ~ "lower",
        TRUE ~ "neutral"
      ),
      bias_value = !!sym(paste0(test_col, "_num")) - !!sym(paste0(reference_col, "_num_adjusted"))
    )
  
  # ---------------------------------------------------------
  # 4. CRÉATION DU SCATTERGRAM
  # ---------------------------------------------------------
  
  # Définition des limites pour le graphique
  min_val <- min(c(data[[paste0(test_col, "_num")]], data[[paste0(reference_col, "_num_adjusted")]]), na.rm = TRUE)
  max_val <- max(c(data[[paste0(test_col, "_num")]], data[[paste0(reference_col, "_num_adjusted")]]), na.rm = TRUE)
  
  
  # Création des graduations logiques
  breaks_log <- 2^seq(log2(min_val), log2(max_val))
  
  plot_scatter <- ggplot(data, aes_string(x = paste0(reference_col, "_num"), y = paste0(test_col, "_num"))) +
    # A. ZONE D'ACCORD ESSENTIEL
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) +
    geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") +
    geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") +
    # B. LES POINTS (Bubble Plot)
    geom_count(aes(color = Is_EA), alpha = 0.7) +
    # C. ÉCHELLES
    scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
    scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
    # D. COULEURS ET ESTHÉTIQUE
    scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"),
                       labels = c("Discordance (>1 dilution)", "Essential agreement (±1 dilution)")) +
    scale_size_continuous(name = "Number of strains") +
    # E. THÈME GRAPHIQUE
    theme_bw() +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    ) +
    # F. TITRES
    labs(
      title = paste("Scattergram, comparaison of MICs for ", molecule_name, " (", "VITEK", " vs ", "BMD", ")", sep = ""),
      subtitle = "Dotted lines: Area of essential agreement (± 1 dilution)",
      x = paste("MICs with ", "BMD", " (Reference method) [mg/L]"),
      y = paste("MICs with ", "VITEK", " [mg/L]")
    )
  
  # Sauvegarder le scattergram dans un fichier image
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  plot_filename <- file.path(output_dir, paste0(gsub("\\.", "_", test_col), "_scattergram.png"))
  ggsave(plot_filename, plot = plot_scatter, width = 10, height = 10, dpi = 300)
  
  # Calcul du biais global
  bias_analysis <- data %>%
    summarise(
      higher = sum(bias_contribution == "higher", na.rm = TRUE),
      lower = sum(bias_contribution == "lower", na.rm = TRUE),
      total = n()
    ) %>%
    mutate(
      bias = ifelse(total == 0, NA, (higher - lower) / total * 100)
    )
  
  # Calcul du % d'EA global
  ea_global <- mean(data$Is_EA) * 100
  
  # Calcul de la sensibilité et spécificité
  data <- data %>%
    mutate(
      test_positive = !!sym(paste0(test_col, "_num")) > threshold,
      reference_positive = !!sym(paste0(reference_col, "_num_adjusted")) > threshold
    )
  
  sens_spec <- data %>%
    summarise(
      TP = sum(test_positive & reference_positive),
      FP = sum(test_positive & !reference_positive),
      TN = sum(!test_positive & !reference_positive),
      FN = sum(!test_positive & reference_positive)
    ) %>%
    mutate(
      Sensibility = ifelse(TP + FN == 0, NA, TP / (TP + FN)),
      Specificity = ifelse(TN + FP == 0, NA, TN / (TN + FP))
    )
  
  # Calcul des intervalles de confiance
  if (nrow(data) > 0) {
    ea_values <- replicate(1000, {
      sample_data <- data[sample(nrow(data), replace = TRUE), ]
      mean(sample_data$Is_EA) * 100
    })
    ci <- quantile(ea_values, probs = c(0.025, 0.975))
  } else {
    ci <- c(NA, NA)
  }
  
  # Calcul du test de McNemar
  if (nrow(data) > 0) {
    mcnemar_result <- tryCatch(
      mcnemar.test(table(data$test_positive, data$reference_positive)),
      error = function(e) list(p.value = NA)
    )
  } else {
    mcnemar_result <- list(p.value = NA)
  }
  
  # Sélection des colonnes spécifiques pour le tableau final
  data_clean_final <- data %>%
    select(
      !!sym(test_col),
      !!sym(reference_col),
      !!sym(paste0(test_col, "_num")),
      !!sym(paste0(reference_col, "_num")),
      !!sym(paste0(reference_col, "_num_adjusted")),
      Is_EA,
      bias_contribution,
      bias_value
    )
  
  message(paste("Analyse terminée pour", test_col, "et", reference_col, ". Nombre final de lignes analysées :", nrow(data_clean_final)))
  
  return(
    list(
      data_clean = data_clean_final,
      ea_global = ea_global,
      bias = bias_analysis$bias,
      sens_spec = sens_spec,
      ci = ci,
      mcnemar_result = mcnemar_result$p.value,
      plot = plot_scatter
    )
  )
}

analyze_multiple_antibiotics <- function(data, antibiotic_pairs, molecule_name = "Molecule", output_file = "resultats_antibiotiques.xlsx", output_dir = "plots") {
  
  # Créer le répertoire de sortie pour les images s'il n'existe pas
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Créer un fichier Excel pour les résultats
  wb <- createWorkbook()
  
  for (pair in antibiotic_pairs) {
    test_col <- pair$columns[1]
    reference_col <- pair$columns[2]
    threshold <- pair$threshold
    
    message(paste("Analyse en cours pour", test_col, "et", reference_col, "avec un seuil de", threshold))
    
    # Appeler la fonction analyze_mic_data pour chaque paire
    result <- analyze_mic_data(data, test_col, reference_col, molecule_name, test_col, reference_col, output_dir, threshold)
    
    if (is.null(result)) {
      warning(paste("Aucun résultat valide pour", test_col, "et", reference_col))
      next
    }
    
    # Ajouter un onglet pour les données détaillées de chaque antibiotique
    sheet_name <- paste0(gsub("\\.", "_", test_col), "_détails")
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = result$data_clean)
    
    # Ajouter un onglet pour les résultats globaux de chaque antibiotique
    sheet_name_global <- paste0(gsub("\\.", "_", test_col), "_résultats")
    results_df <- data.frame(
      Molecule = molecule_name,
      Antibiotique = test_col,
      Seuil = threshold,
      EA_global = result$ea_global,
      Bias = result$bias,
      Sensibility = result$sens_spec$Sensibility,
      Specificity = result$sens_spec$Specificity,
      CI_lower = result$ci[1],
      CI_upper = result$ci[2],
      McNemar_p_value = result$mcnemar_result
    )
    addWorksheet(wb, sheetName = sheet_name_global)
    writeData(wb, sheet = sheet_name_global, x = results_df)
  }
  
  # Sauvegarder le fichier Excel
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message(paste("Tous les résultats ont été sauvegardés dans", output_file))
}

########################################

# Fonction pour analyser une paire d'antibiotiques par strate
analyze_mic_data_by_strata <- function(data, test_col, reference_col, strata_columns, output_dir = "plots", threshold = 16) {
  # Nettoyage et ajustement des valeurs de la colonne test
  data <- data %>%
    mutate(
      !!paste0(test_col, "_num") := clean_and_adjust_mic(!!sym(test_col)),
      !!paste0(reference_col, "_num") := clean_and_adjust_mic(!!sym(reference_col))
    ) %>%
    filter(!is.na(!!sym(paste0(test_col, "_num"))) & !is.na(!!sym(paste0(reference_col, "_num"))))
  
  if (nrow(data) == 0) {
    warning(paste("Aucune donnée valide pour", test_col, "et", reference_col))
    return(NULL)
  }
  
  # Détermination des valeurs minimale et maximale de la colonne test_num
  min_test <- min(data[[paste0(test_col, "_num")]], na.rm = TRUE)
  max_test <- max(data[[paste0(test_col, "_num")]], na.rm = TRUE)
  
  # Ajustement des valeurs de la colonne reference_num
  data <- data %>%
    mutate(
      !!paste0(reference_col, "_num_adjusted") := case_when(
        !!sym(paste0(reference_col, "_num")) < min_test ~ min_test,
        !!sym(paste0(reference_col, "_num")) > max_test ~ max_test,
        TRUE ~ !!sym(paste0(reference_col, "_num"))
      )
    )
  
  # Calcul de l'EA (Essential Agreement)
  data <- data %>%
    mutate(Is_EA = abs(log2(!!sym(paste0(test_col, "_num"))) - log2(!!sym(paste0(reference_col, "_num_adjusted")))) <= 1)
  
  # Calcul de la contribution au biais
  data <- data %>%
    mutate(
      bias_contribution = case_when(
        !!sym(paste0(test_col, "_num")) > !!sym(paste0(reference_col, "_num_adjusted")) ~ "higher",
        !!sym(paste0(test_col, "_num")) < !!sym(paste0(reference_col, "_num_adjusted")) ~ "lower",
        TRUE ~ "neutral"
      ),
      bias_value = !!sym(paste0(test_col, "_num")) - !!sym(paste0(reference_col, "_num_adjusted")),
      test_positive = !!sym(paste0(test_col, "_num")) > threshold,
      reference_positive = !!sym(paste0(reference_col, "_num_adjusted")) > threshold
    )
  
  # Calcul des statistiques par strate
  results_by_strata <- data %>%
    group_by(across(all_of(strata_columns))) %>%
    summarise(
      count = n(),
      EA_global = mean(Is_EA) * 100,
      higher = sum(bias_contribution == "higher", na.rm = TRUE),
      lower = sum(bias_contribution == "lower", na.rm = TRUE),
      bias = (higher - lower) / count * 100,
      TP = sum(test_positive & reference_positive, na.rm = TRUE),
      FP = sum(test_positive & !reference_positive, na.rm = TRUE),
      TN = sum(!test_positive & !reference_positive, na.rm = TRUE),
      FN = sum(!test_positive & reference_positive, na.rm = TRUE),
      Sensibility = ifelse(TP + FN == 0, NA, TP / (TP + FN)),
      Specificity = ifelse(TN + FP == 0, NA, TN / (TN + FP)),
      .groups = "drop_last"
    ) %>%
    mutate(McNemar_p_value = NA_real_)
  
  # Calcul du test de McNemar par strate
  for (i in 1:nrow(results_by_strata)) {
    strata_values <- as.list(results_by_strata[i, strata_columns])
    strata_data <- data
    for (col in strata_columns) {
      strata_data <- strata_data %>%
        filter(!!sym(col) == strata_values[[col]])
    }
    if (nrow(strata_data) > 0) {
      contingency_table <- table(strata_data$test_positive, strata_data$reference_positive)
      if (all(dim(contingency_table) == 2)) {
        test_result <- tryCatch(
          mcnemar.test(contingency_table),
          error = function(e) list(p.value = NA)
        )
        results_by_strata$McNemar_p_value[i] <- test_result$p.value
      }
    }
  }
  
  # Calcul des intervalles de confiance par strate
  calculate_ci <- function(sub_data) {
    if (nrow(sub_data) > 0) {
      ea_values <- replicate(1000, {
        sample_data <- sub_data[sample(nrow(sub_data), replace = TRUE), ]
        mean(sample_data$Is_EA) * 100
      })
      ci <- quantile(ea_values, probs = c(0.025, 0.975))
    } else {
      ci <- c(NA, NA)
    }
    return(ci)
  }
  
  results_by_strata <- results_by_strata %>%
    mutate(
      CI_lower = NA_real_,
      CI_upper = NA_real_
    )
  
  for (i in 1:nrow(results_by_strata)) {
    strata_values <- as.list(results_by_strata[i, strata_columns])
    strata_data <- data
    for (col in strata_columns) {
      strata_data <- strata_data %>%
        filter(!!sym(col) == strata_values[[col]])
    }
    ci <- calculate_ci(strata_data)
    results_by_strata$CI_lower[i] <- ci[1]
    results_by_strata$CI_upper[i] <- ci[2]
  }
  
  return(list(data_clean = data, results_by_strata = results_by_strata))
}

# Fonction pour analyser plusieurs antibiotiques par strate
analyze_multiple_antibiotics_by_strata <- function(data, antibiotic_pairs, strata_columns, output_file = "resultats_antibiotiques_stratifies.xlsx", output_dir = "plots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  wb <- createWorkbook()
  
  for (pair in antibiotic_pairs) {
    test_col <- pair$columns[1]
    reference_col <- pair$columns[2]
    threshold <- pair$threshold
    
    message(paste("Analyse en cours pour", test_col, "et", reference_col, "avec un seuil de", threshold))
    
    # Appeler la fonction analyze_mic_data_by_strata pour chaque paire
    result <- analyze_mic_data_by_strata(data, test_col, reference_col, strata_columns, output_dir, threshold)
    
    if (is.null(result)) {
      warning(paste("Aucun résultat valide pour", test_col, "et", reference_col))
      next
    }
    
    # Ajouter un onglet pour les données détaillées de chaque antibiotique
    sheet_name <- paste0(gsub("\\.", "_", test_col), "_détails")
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = result$data_clean)
    
    # Ajouter un onglet pour les résultats globaux de chaque antibiotique
    sheet_name_global <- paste0(gsub("\\.", "_", test_col), "_résultats")
    addWorksheet(wb, sheetName = sheet_name_global)
    writeData(wb, sheet = sheet_name_global, x = result$results_by_strata)
  }
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message(paste("Tous les résultats ont été sauvegardés dans", output_file))
}

# Fonction pour analyser plusieurs antibiotiques par type de strate et écrire dans des fichiers texte
analyze_multiple_antibiotics_by_strata_types_to_txt <- function(data, antibiotic_pairs, strata_types, output_dir = "results") {
  # Créer le répertoire de sortie s'il n'existe pas
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  for (stratum_type in strata_types) {
    # Créer un data frame vide pour stocker tous les résultats
    all_results <- data.frame()
    
    for (pair in antibiotic_pairs) {
      test_col <- paste0(pair$columns[1], "")
      reference_col <- paste0(pair$columns[2], "")
      threshold <- pair$threshold
      
      message(paste("Analyse en cours pour", test_col, "et", reference_col, "avec un seuil de", threshold, "par", stratum_type))
      
      # Appeler la fonction pour chaque paire
      result <- analyze_mic_data_by_stratum(data, test_col, reference_col, stratum_type, threshold = threshold)
      
      if (!is.null(result)) {
        # Ajouter le nom de l'antibiotique aux résultats
        result$Antibiotique <- gsub("\\..*", "", pair$columns[1])
        
        # Ajouter les résultats au data frame global
        all_results <- bind_rows(all_results, result)
      }
    }
    
    # Écrire les résultats dans un fichier texte avec séparateur de tabulation
    output_file <- file.path(output_dir, paste0(stratum_type, "_results.txt"))
    write.table(all_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  }
}

analyze_mic_data_by_stratum <- function(data, test_col, reference_col, stratum_column, output_dir = "plots", threshold = 16) {
  # Nettoyage et ajustement des valeurs de la colonne test
  data <- data %>%
    mutate(
      !!paste0(test_col, "_num") := clean_and_adjust_mic(!!sym(test_col)),
      !!paste0(reference_col, "_num") := clean_and_adjust_mic(!!sym(reference_col))
    ) %>%
    filter(!is.na(!!sym(paste0(test_col, "_num"))) & !is.na(!!sym(paste0(reference_col, "_num"))))
  
  if (nrow(data) == 0) {
    warning(paste("Aucune donnée valide pour", test_col, "et", reference_col))
    return(NULL)
  }
  
  # Détermination des valeurs minimale et maximale de la colonne test_num
  min_test <- min(data[[paste0(test_col, "_num")]], na.rm = TRUE)
  max_test <- max(data[[paste0(test_col, "_num")]], na.rm = TRUE)
  
  # Ajustement des valeurs de la colonne reference_num
  data <- data %>%
    mutate(
      !!paste0(reference_col, "_num_adjusted") := case_when(
        !!sym(paste0(reference_col, "_num")) < min_test ~ min_test,
        !!sym(paste0(reference_col, "_num")) > max_test ~ max_test,
        TRUE ~ !!sym(paste0(reference_col, "_num"))
      )
    )
  
  # Calcul de l'EA (Essential Agreement)
  data <- data %>%
    mutate(Is_EA = abs(log2(!!sym(paste0(test_col, "_num"))) - log2(!!sym(paste0(reference_col, "_num_adjusted")))) <= 1)
  
  # Calcul de la contribution au biais
  data <- data %>%
    mutate(
      bias_contribution = case_when(
        !!sym(paste0(test_col, "_num")) > !!sym(paste0(reference_col, "_num_adjusted")) ~ "higher",
        !!sym(paste0(test_col, "_num")) < !!sym(paste0(reference_col, "_num_adjusted")) ~ "lower",
        TRUE ~ "neutral"
      ),
      bias_value = !!sym(paste0(test_col, "_num")) - !!sym(paste0(reference_col, "_num_adjusted")),
      test_positive = !!sym(paste0(test_col, "_num")) > threshold,
      reference_positive = !!sym(paste0(reference_col, "_num_adjusted")) > threshold
    )
  
  # Calcul des statistiques par strate
  results_by_stratum <- data %>%
    group_by(!!sym(stratum_column)) %>%
    summarise(
      Total_valide = n(),
      EA_global = mean(Is_EA) * 100,
      higher = sum(bias_contribution == "higher", na.rm = TRUE),
      lower = sum(bias_contribution == "lower", na.rm = TRUE),
      bias = ifelse(Total_valide == 0, NA, (higher - lower) / Total_valide * 100),
      TP = sum(test_positive & reference_positive, na.rm = TRUE),
      FP = sum(test_positive & !reference_positive, na.rm = TRUE),
      TN = sum(!test_positive & !reference_positive, na.rm = TRUE),
      FN = sum(!test_positive & reference_positive, na.rm = TRUE),
      Sensibility = ifelse(TP + FN == 0, NA, TP / (TP + FN)),
      Specificity = ifelse(TN + FP == 0, NA, TN / (TN + FP)),
      .groups = "drop"
    )
  
  # Calcul du test de McNemar par strate
  results_by_stratum$McNemar_p_value <- NA_real_
  
  for (i in 1:nrow(results_by_stratum)) {
    stratum_value <- results_by_stratum[[stratum_column]][i]
    stratum_data <- data %>%
      filter(!!sym(stratum_column) == stratum_value)
    
    if (nrow(stratum_data) > 0) {
      contingency_table <- table(stratum_data$test_positive, stratum_data$reference_positive)
      if (all(dim(contingency_table) == 2)) {
        test_result <- tryCatch(
          mcnemar.test(contingency_table, correct = FALSE),
          error = function(e) list(p.value = NA)
        )
        results_by_stratum$McNemar_p_value[i] <- ifelse(is.null(test_result$p.value), NA, test_result$p.value)
      }
    }
  }
  
  # Calcul des intervalles de confiance par strate
  calculate_ci <- function(sub_data) {
    if (nrow(sub_data) > 0) {
      ea_values <- replicate(1000, {
        sample_data <- sub_data[sample(nrow(sub_data), replace = TRUE), ]
        mean(sample_data$Is_EA) * 100
      })
      ci <- quantile(ea_values, probs = c(0.025, 0.975))
    } else {
      ci <- c(NA, NA)
    }
    return(ci)
  }
  
  results_by_stratum <- results_by_stratum %>%
    mutate(
      CI_lower = NA_real_,
      CI_upper = NA_real_
    )
  
  for (i in 1:nrow(results_by_stratum)) {
    stratum_value <- results_by_stratum[[stratum_column]][i]
    stratum_data <- data %>%
      filter(!!sym(stratum_column) == stratum_value)
    ci <- calculate_ci(stratum_data)
    results_by_stratum$CI_lower[i] <- ci[1]
    results_by_stratum$CI_upper[i] <- ci[2]
  }
  
  return(results_by_stratum)
}

########################################

# Fonction pour analyser les erreurs catégorielles
analyze_categorical_agreement <- function(data, test_col_inter, reference_col_inter) {
  # Filtrer les lignes contenant "\", "ZIT" ou des valeurs manquantes
  data_filtered <- data %>%
    filter(!str_detect(toupper(!!sym(test_col_inter)), "[\\\\]|ZIT") &
             !str_detect(toupper(!!sym(reference_col_inter)), "[\\\\]|ZIT") &
             !is.na(!!sym(test_col_inter)) &
             !is.na(!!sym(reference_col_inter)))
  
  # Compter le nombre total de données avant filtrage
  total_data <- nrow(data)
  
  # Compter le nombre de données après filtrage
  filtered_data <- nrow(data_filtered)
  
  message(paste("Total de données:", total_data, "| Données valides:", filtered_data))
  
  # Calcul des erreurs catégorielles
  data_filtered <- data_filtered %>%
    mutate(
      CA = !!sym(test_col_inter) == !!sym(reference_col_inter),
      VME = (toupper(!!sym(test_col_inter)) == "S") & (toupper(!!sym(reference_col_inter)) == "R"),
      ME = (toupper(!!sym(test_col_inter)) == "R") & (toupper(!!sym(reference_col_inter)) == "S"),
      miE = ((toupper(!!sym(test_col_inter)) == "S" & toupper(!!sym(reference_col_inter)) == "SFP") |
               (toupper(!!sym(test_col_inter)) == "SFP" & toupper(!!sym(reference_col_inter)) == "S") |
               (toupper(!!sym(test_col_inter)) == "R" & toupper(!!sym(reference_col_inter)) == "SFP") |
               (toupper(!!sym(test_col_inter)) == "SFP" & toupper(!!sym(reference_col_inter)) == "R"))
    )
  
  # Calcul des statistiques
  results <- data_filtered %>%
    summarise(
      Total_valide = n(),
      CA_count = sum(CA, na.rm = TRUE),
      CA_percentage = mean(CA, na.rm = TRUE) * 100,
      VME_count = sum(VME, na.rm = TRUE),
      VME_percentage = mean(VME, na.rm = TRUE) * 100,
      ME_count = sum(ME, na.rm = TRUE),
      ME_percentage = mean(ME, na.rm = TRUE) * 100,
      miE_count = sum(miE, na.rm = TRUE),
      miE_percentage = mean(miE, na.rm = TRUE) * 100,
      Total_initial = total_data,
      Valides = filtered_data
    )
  
  return(results)
}


# Fonction pour analyser plusieurs antibiotiques et tout mettre dans un seul onglet Excel
analyze_multiple_antibiotics_inter_single_sheet <- function(data, antibiotic_pairs, output_file = "resultats_antibiotiques_interpretation.xlsx") {
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Résultats")
  
  # Créer un data frame vide pour stocker tous les résultats
  all_results <- data.frame()
  
  for (pair in antibiotic_pairs) {
    test_col <- paste0(pair$columns[1], "_inter")
    reference_col <- paste0(pair$columns[2], "_inter")
    
    message(paste("Analyse en cours pour", test_col, "et", reference_col))
    
    # Appeler la fonction pour chaque paire
    result <- analyze_categorical_agreement(data, test_col, reference_col)
    
    # Ajouter le nom de l'antibiotique aux résultats
    result$Antibiotique <- gsub("\\..*", "", pair$columns[1])
    
    # Ajouter les résultats au data frame global
    all_results <- bind_rows(all_results, result)
  }
  
  # Écrire tous les résultats dans un seul onglet
  writeData(wb, sheet = "Résultats", x = all_results)
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message(paste("Tous les résultats ont été sauvegardés dans un seul onglet du fichier", output_file))
}

############################################

#Fonction pour analyser les erreurs catégorielles par strate

# Fonction pour lire les données
read_data <- function(file_path) {
  data <- read_delim(file_path, delim = ";", show_col_types = FALSE)
  return(data)
}

# Fonction pour catégoriser les erreurs
categorize_errors <- function(vitek, bmd) {
  if (is.na(vitek) | is.na(bmd)) {
    return(NA)
  } else if (vitek == bmd) {
    return("CA")
  } else if (vitek == "S" & bmd == "R") {
    return("VME")
  } else if (vitek == "R" & bmd == "S") {
    return("ME")
  } else {
    return("mE")
  }
}

# Fonction pour calculer la répartition des erreurs
calculate_error_distribution <- function(data, category) {
  # Liste des antibiotiques
  antibiotics <- c("PTZ", "CZD", "CRO", "FEP", "ETP", "IMP", "AKN", "LVX", "NFE", "SXT", "AMX", "AMC", "FIX", "FOX", "FOS", "MEC")
  
  # Initialiser une liste pour stocker les résultats
  results <- list()
  
  for (antibio in antibiotics) {
    vitek_inter_col <- paste0(antibio, "_VITEK_inter")
    bmd_inter_col <- paste0(antibio, "_BMD_inter")
    
    if (vitek_inter_col %in% colnames(data) & bmd_inter_col %in% colnames(data)) {
      # Créer un data frame temporaire pour l'antibiotique actuel
      temp_df <- data %>%
        select(all_of(c(category, vitek_inter_col, bmd_inter_col))) %>%
        mutate(Error_Type = mapply(categorize_errors, !!sym(vitek_inter_col), !!sym(bmd_inter_col))) %>%
        filter(!is.na(Error_Type)) %>%
        group_by(!!sym(category), Error_Type) %>%
        summarise(Count = n(), .groups = "drop") %>%
        pivot_wider(names_from = Error_Type, values_from = Count, values_fill = 0) %>%
        mutate(Antibiotique = antibio)
      
      # Stocker les résultats
      results[[length(results) + 1]] <- temp_df
    }
  }
  
  # Combiner les résultats pour tous les antibiotiques
  if (length(results) > 0) {
    combined_results <- bind_rows(results) %>%
      pivot_longer(cols = -c(all_of(category), Antibiotique), names_to = "Error_Type", values_to = "Count") %>%
      pivot_wider(id_cols = all_of(category), names_from = c(Antibiotique, Error_Type), values_from = Count, values_fill = 0)
    return(combined_results)
  } else {
    return(data.frame())
  }
}

# Fonction principale pour effectuer les analyses
main_analysis <- function(file_path) {
  # Lire les données
  data <- read_data(file_path)
  
  # Calculer la répartition des erreurs pour chaque catégorie
  categories <- c("Groupe", "Especes", "Phenotype")
  all_results <- list()
  
  for (category in categories) {
    all_results[[category]] <- calculate_error_distribution(data, category)
  }
  
  # Créer un fichier Excel avec un onglet pour chaque catégorie
  wb <- createWorkbook()
  
  for (category in names(all_results)) {
    addWorksheet(wb, category)
    writeData(wb, sheet = category, x = all_results[[category]])
  }
  
  # Sauvegarder le fichier Excel
  saveWorkbook(wb, "./tables/error_distribution_combined.xlsx", overwrite = TRUE)
}


################################################################################

# --- Fonction d'analyse pour une ligne donnée (un groupe, une espèce, etc.) ---
analyze_row_variation <- function(row_data, id_name) {
  
  # 1. Restructuration des données : Passer de "Colonnes larges" à "Format long"
  # On ne garde que les colonnes qui finissent par _CA, _ME, _mE ou _VME
  long_data <- row_data %>%
    pivot_longer(
      cols = matches("_(CA|ME|mE|VME)$"), 
      names_to = "Original_Col", 
      values_to = "Count"
    ) %>%
    # Extraction propre du nom de l'antibiotique et du type d'erreur
    # Regex : Capture tout avant le dernier "_" comme nom d'antibio
    extract(
      col = Original_Col, 
      into = c("Antibiotique", "Error_Type"), 
      regex = "^(.*)_(CA|ME|mE|VME)$"
    ) %>%
    mutate(Count = ifelse(is.na(Count), 0, Count)) # Remplacer les vides par 0
  
  # 2. Création de la table de contingence (Lignes=Antibios, Colonnes=Erreurs)
  contingency_table <- long_data %>%
    group_by(Antibiotique, Error_Type) %>%
    summarise(Count = sum(Count), .groups = 'drop') %>%
    pivot_wider(names_from = Error_Type, values_from = Count, values_fill = 0) %>%
    column_to_rownames(var = "Antibiotique")
  
  # 3. Test du Chi² (si la table contient des données)
  # On vérifie qu'il y a au moins 2 colonnes et 2 lignes et que la somme n'est pas nulle
  if (nrow(contingency_table) > 1 && ncol(contingency_table) > 1 && sum(contingency_table) > 0) {
    
    # Utilisation de tryCatch pour éviter que le script ne plante si le test échoue (ex: données insuffisantes)
    test_result <- tryCatch({
      chisq.test(as.matrix(contingency_table))
    }, error = function(e) return(NULL))
    
    if (!is.null(test_result)) {
      return(data.frame(
        Sous_Categorie = id_name,
        Chi2_Stat = round(test_result$statistic, 3),
        Df = test_result$parameter,
        P_Value = test_result$p.value
      ))
    }
  }
  
  # Si le test n'est pas possible
  return(data.frame(Sous_Categorie = id_name, Chi2_Stat = NA, Df = NA, P_Value = NA))
}

# --- Fonction Principale ---
run_global_analysis <- function(input_file, output_file) {
  
  # Charger le fichier Excel source
  wb_source <- loadWorkbook(input_file)
  sheet_names <- names(wb_source)
  
  # Créer le fichier Excel de sortie
  wb_results <- createWorkbook()
  
  print(paste("Début de l'analyse du fichier :", input_file))
  
  # Boucle sur chaque feuille (Groupe, Espèce, Phénotype...)
  for (sheet in sheet_names) {
    
    cat("Traitement de la feuille :", sheet, "\n")
    
    # Lire les données de la feuille
    data <- read.xlsx(wb_source, sheet = sheet)
    
    # Initialiser une liste pour stocker les résultats de cette feuille
    sheet_results_list <- list()
    
    # Boucle sur chaque ligne de la feuille (chaque Groupe/Espèce individuellement)
    # On suppose que la 1ère colonne contient l'identifiant (ex: "1", "E. coli", "AmpC")
    id_col_name <- colnames(data)[1]
    
    for (i in 1:nrow(data)) {
      # Récupérer l'identifiant de la ligne (ex: "Groupe 1")
      row_id <- as.character(data[i, 1])
      
      # Lancer l'analyse pour cette ligne spécifique
      row_result <- analyze_row_variation(data[i, ], row_id)
      
      # Ajouter au tableau de résultats
      sheet_results_list[[i]] <- row_result
    }
    
    # Combiner tous les résultats de la feuille
    final_sheet_results <- do.call(rbind, sheet_results_list)
    
    # Ajouter l'interprétation
    final_sheet_results <- final_sheet_results %>%
      mutate(
        Interpretation = case_when(
          P_Value < 0.001 ~ "Très significatif (p < 0.001)",
          P_Value < 0.01  ~ "Significatif (p < 0.01)",
          P_Value < 0.05  ~ "Faiblement significatif (p < 0.05)",
          P_Value >= 0.05 ~ "Non significatif",
          TRUE ~ "Données insuffisantes"
        )
      )
    
    # Écrire les résultats dans une nouvelle feuille du fichier de sortie
    addWorksheet(wb_results, sheet = substr(paste0("Res_", sheet), 1, 31)) # Max 31 chars pour Excel
    writeData(wb_results, sheet = substr(paste0("Res_", sheet), 1, 31), x = final_sheet_results)
  }
  
  # Sauvegarder
  saveWorkbook(wb_results, file = output_file, overwrite = TRUE)
  print(paste("Terminé ! Fichier sauvegardé :", output_file))
}

# --- Exécution ---
input_path <- "./tables/error_distribution_combined.xlsx"
output_path <- "./tables/chi_square_results_by_row.xlsx"

# Vérification et lancement
if (file.exists(input_path)) {
  run_global_analysis(input_path, output_path)
} else {
  message("Erreur : Le fichier d'entrée est introuvable.")
}



# --- Fonction : Analyse "Un contre le Reste" pour une ligne de données ---
analyze_one_vs_rest <- function(row_data, id_name, category_type) {
  
  # 1. Mise en format long (Même logique que précédemment)
  long_data <- row_data %>%
    pivot_longer(
      cols = matches("_(CA|ME|mE|VME)$"), 
      names_to = "Original_Col", 
      values_to = "Count"
    ) %>%
    extract(
      col = Original_Col, 
      into = c("Antibiotique", "Error_Type"), 
      regex = "^(.*)_(CA|ME|mE|VME)$"
    ) %>%
    mutate(Count = ifelse(is.na(Count), 0, Count))
  
  # Liste des antibiotiques présents dans cette ligne
  antibio_list <- unique(long_data$Antibiotique)
  
  # Stockage des résultats pour cette ligne
  results_list <- list()
  
  # 2. Boucle sur chaque antibiotique (PTZ, puis CZD, etc.)
  for (target_ab in antibio_list) {
    
    # A. Données de l'antibiotique cible
    target_data <- long_data %>%
      filter(Antibiotique == target_ab) %>%
      select(Error_Type, Count) %>%
      rename(Target = Count)
    
    # B. Données du "Reste" (Tous les autres antibiotiques combinés)
    rest_data <- long_data %>%
      filter(Antibiotique != target_ab) %>%
      group_by(Error_Type) %>%
      summarise(Rest = sum(Count), .groups = 'drop')
    
    # C. Fusionner pour créer la table de contingence (2 colonnes : Target vs Rest)
    contingency_df <- full_join(target_data, rest_data, by = "Error_Type") %>%
      replace(is.na(.), 0) %>% # Remplacer les NA par 0 si une erreur manque
      column_to_rownames(var = "Error_Type")
    
    # Transposer pour avoir : Lignes = Target/Rest, Colonnes = CA/ME/mE/VME
    contingency_matrix <- t(as.matrix(contingency_df))
    
    # D. Test statistique
    # On ne fait le test que si on a des données
    res_pval <- NA
    res_stat <- NA
    
    if (sum(contingency_matrix) > 0 && ncol(contingency_matrix) > 1) {
      test_result <- tryCatch({
        chisq.test(contingency_matrix)
      }, warning = function(w) {
        # Si warning (souvent dû à des petits effectifs), on tente un Fisher simulé si besoin, 
        # ou on garde le Chi2 avec avertissement. Ici on garde le résultat standard.
        suppressWarnings(chisq.test(contingency_matrix))
      }, error = function(e) return(NULL))
      
      if (!is.null(test_result)) {
        res_pval <- test_result$p.value
        res_stat <- test_result$statistic
      }
    }
    
    # E. Enregistrement du résultat
    results_list[[target_ab]] <- data.frame(
      Type_Categorie = category_type,
      Nom_Categorie = id_name,       # Ex: "Groupe 1" ou "E. coli"
      Antibiotique = target_ab,      # Ex: "PTZ"
      Chi2_Stat = round(res_stat, 3),
      P_Value = res_pval
    )
  }
  
  # Combiner les résultats de tous les antibios de cette ligne
  return(do.call(rbind, results_list))
}

# --- Fonction Principale ---
run_detailed_analysis <- function(input_file, output_file) {
  
  wb_source <- loadWorkbook(input_file)
  sheet_names <- names(wb_source)
  
  # On va créer un seul grand tableau de résultats ou par feuille
  wb_results <- createWorkbook()
  
  for (sheet in sheet_names) {
    message(paste("Analyse détaillée de la feuille :", sheet))
    
    data <- read.xlsx(wb_source, sheet = sheet)
    
    # Liste pour accumuler tous les résultats de cette feuille
    sheet_all_results <- list()
    
    # Boucle sur chaque ligne (Chaque Groupe/Espèce)
    for (i in 1:nrow(data)) {
      row_id <- as.character(data[i, 1]) # Nom de la ligne (ex: "Groupe 1")
      
      # Analyse par antibiotique pour cette ligne
      row_res <- analyze_one_vs_rest(data[i, ], row_id, sheet)
      
      sheet_all_results[[i]] <- row_res
    }
    
    # Fusionner tout
    final_df <- do.call(rbind, sheet_all_results)
    
    # Ajout des interprétations
    final_df <- final_df %>%
      mutate(
        Interpretation = case_when(
          P_Value < 0.001 ~ "*** (p < 0.001)",
          P_Value < 0.01  ~ "** (p < 0.01)",
          P_Value < 0.05  ~ "* (p < 0.05)",
          P_Value >= 0.05 ~ "NS",
          TRUE ~ "NA"
        )
      ) %>%
      arrange(Nom_Categorie, P_Value) # Trier par catégorie puis par p-value
    
    # Écriture dans l'Excel de sortie
    sheet_out_name <- substr(paste0("Détail_", sheet), 1, 31)
    addWorksheet(wb_results, sheet = sheet_out_name)
    writeData(wb_results, sheet = sheet_out_name, x = final_df)
    
    # Colorier les cellules significatives (Bonus visuel)
    conditionalFormatting(wb_results, sheet_out_name, cols = 6, rows = 2:(nrow(final_df)+1), 
                          rule = "<0.05", style = createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE"))
  }
  
  saveWorkbook(wb_results, file = output_file, overwrite = TRUE)
  message(paste("Terminé ! Analyse par antibiotique sauvegardée :", output_file))
}

# --- Exécution ---
input_path <- "./tables/error_distribution_combined.xlsx"
output_path <- "./tables/chi_square_detailed_by_antibio.xlsx"

if (file.exists(input_path)) {
  run_detailed_analysis(input_path, output_path)
} else {
  message("Fichier introuvable.")
}

# --- Fonction : Analyse d'un antibiotique à travers différentes catégories ---
analyze_antibio_across_groups <- function(data, sheet_name) {
  
  # 1. Identifier le nom de la colonne catégorie (ex: "Groupe", "Especes", "Phenotype")
  category_col_name <- colnames(data)[1]
  
  # 2. Transformer en format long pour manipuler les données
  long_data <- data %>%
    rename(Category = all_of(category_col_name)) %>% # On renomme temporairement pour faciliter le code
    pivot_longer(
      cols = matches("_(CA|ME|mE|VME)$"), 
      names_to = "Original_Col", 
      values_to = "Count"
    ) %>%
    extract(
      col = Original_Col, 
      into = c("Antibiotique", "Error_Type"), 
      regex = "^(.*)_(CA|ME|mE|VME)$"
    ) %>%
    mutate(Count = ifelse(is.na(Count), 0, Count)) # Remplacer NA par 0
  
  # Liste des antibiotiques uniques
  antibio_list <- unique(long_data$Antibiotique)
  
  results_list <- list()
  
  # 3. Boucle sur chaque antibiotique (ex: On analyse PTZ pour tous les groupes)
  for (ab in antibio_list) {
    
    # Créer la table de contingence pour CET antibiotique
    # Lignes = Les différents Groupes/Espèces
    # Colonnes = Les types d'erreurs (CA, ME, mE, VME)
    contingency_table <- long_data %>%
      filter(Antibiotique == ab) %>%
      select(Category, Error_Type, Count) %>%
      pivot_wider(names_from = Error_Type, values_from = Count, values_fill = 0) %>%
      column_to_rownames(var = "Category")
    
    # Nettoyage : On ne garde que les lignes (Groupes) qui ont des données pour cet antibio
    # Si un groupe n'a pas été testé pour le PTZ (somme = 0), on l'enlève du test
    contingency_table <- contingency_table[rowSums(contingency_table) > 0, ]
    
    # Test statistique
    res_pval <- NA
    res_stat <- NA
    df <- NA
    
    # Il faut au moins 2 groupes comparables et des données
    if (nrow(contingency_table) > 1 && sum(contingency_table) > 0) {
      test_result <- tryCatch({
        chisq.test(as.matrix(contingency_table))
      }, warning = function(w) {
        suppressWarnings(chisq.test(as.matrix(contingency_table)))
      }, error = function(e) return(NULL))
      
      if (!is.null(test_result)) {
        res_pval <- test_result$p.value
        res_stat <- test_result$statistic
        df <- test_result$parameter
      }
    }
    
    # Enregistrement
    results_list[[ab]] <- data.frame(
      Source_Fichier = sheet_name,
      Antibiotique = ab,
      Chi2_Stat = round(res_stat, 3),
      Degres_Liberte = df,
      P_Value = res_pval
    )
  }
  
  return(do.call(rbind, results_list))
}

# --- Fonction Principale ---
run_horizontal_analysis <- function(input_file, output_file) {
  
  wb_source <- loadWorkbook(input_file)
  sheet_names <- names(wb_source)
  wb_results <- createWorkbook()
  
  for (sheet in sheet_names) {
    message("Analyse 'Horizontale' (Comparaison entre catégories) pour : ", sheet)
    
    data <- read.xlsx(wb_source, sheet = sheet)
    
    # Analyse
    sheet_res <- analyze_antibio_across_groups(data, sheet)
    
    if (!is.null(sheet_res)) {
      # Ajout Interprétation
      sheet_res <- sheet_res %>%
        mutate(
          Interpretation = case_when(
            P_Value < 0.001 ~ "*** (Différence Très Signif.)",
            P_Value < 0.01  ~ "** (Différence Signif.)",
            P_Value < 0.05  ~ "* (Différence Faible)",
            P_Value >= 0.05 ~ "NS (Pas de différence)",
            TRUE ~ "Données insuffisantes"
          )
        ) %>%
        arrange(P_Value) # Les antibios les plus impactés en premier
      
      # Sauvegarde onglet
      sheet_out <- substr(paste0("Comp_", sheet), 1, 31)
      addWorksheet(wb_results, sheet = sheet_out)
      writeData(wb_results, sheet = sheet_out, x = sheet_res)
      
      # Couleur si p < 0.05
      conditionalFormatting(wb_results, sheet_out, cols = 6, rows = 2:(nrow(sheet_res)+1),
                            rule = "<0.05", style = createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE"))
    }
  }
  
  saveWorkbook(wb_results, file = output_file, overwrite = TRUE)
  message("Terminé ! Analyse sauvegardée : ", output_file)
}

# --- Exécution ---
input_path <- "./tables/error_distribution_combined.xlsx"
output_path <- "./tables/chi_square_antibio_across_groups.xlsx"

if (file.exists(input_path)) {
  run_horizontal_analysis(input_path, output_path)
} else {
  message("Fichier introuvable.")
}

################################################################################
################################################################################

# --- Fonction pour calculer les Odds Ratios (Facteurs de risque) ---
calculate_odds_ratios <- function(data, sheet_name) {
  
  # 1. Préparation des données en format long
  category_col <- colnames(data)[1]
  
  long_data <- data %>%
    rename(Category_ID = all_of(category_col)) %>%
    pivot_longer(
      cols = matches("_(CA|ME|mE|VME)$"), 
      names_to = "Original_Col", 
      values_to = "Count"
    ) %>%
    extract(
      col = Original_Col, 
      into = c("Antibiotique", "Error_Type"), 
      regex = "^(.*)_(CA|ME|mE|VME)$"
    ) %>%
    mutate(Count = ifelse(is.na(Count), 0, Count))
  
  # Liste des éléments à itérer
  categories <- unique(long_data$Category_ID)
  antibiotics <- unique(long_data$Antibiotique)
  error_types <- c("VME", "ME", "mE") # On se concentre sur les erreurs
  
  results_list <- list()
  counter <- 1
  
  # 2. Boucles d'analyse
  for (ab in antibiotics) {
    for (err in error_types) {
      
      # Filtrer les données pour cet antibiotique
      ab_data <- long_data %>% filter(Antibiotique == ab)
      
      # Calculer le total des tests pour cet antibio (par catégorie)
      # (Pour savoir combien de "Non-Erreur" il y a)
      total_tests_per_cat <- ab_data %>%
        group_by(Category_ID) %>%
        summarise(Total = sum(Count), .groups = 'drop')
      
      # Combien d'erreurs spécifiques (ex: VME) par catégorie ?
      error_counts_per_cat <- ab_data %>%
        filter(Error_Type == err) %>%
        select(Category_ID, Count)
      
      # Fusionner pour avoir : Catégorie | N_Erreur | N_Total
      stats_df <- left_join(total_tests_per_cat, error_counts_per_cat, by = "Category_ID") %>%
        mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
        mutate(No_Error = Total - Count) # "Pas d'erreur" = Total - (VME ou ME ou mE)
      
      # Vérifier qu'il y a des données valides
      if(sum(stats_df$Total) == 0) next
      
      # 3. Analyse "Catégorie vs Reste"
      for (cat in categories) {
        
        # A. Définir la table 2x2 pour le Fisher Test
        #                 | Erreur (Oui) | Erreur (Non)
        # ----------------|--------------|-------------
        # Groupe Cible    |      a       |      b
        # Reste du Monde  |      c       |      d
        
        # Ligne "Cible"
        target_row <- stats_df %>% filter(Category_ID == cat)
        if(nrow(target_row) == 0) next
        
        a <- target_row$Count
        b <- target_row$No_Error
        
        # Ligne "Reste"
        rest_rows <- stats_df %>% filter(Category_ID != cat)
        c <- sum(rest_rows$Count)
        d <- sum(rest_rows$No_Error)
        
        # Si aucun test n'a été fait pour cette catégorie sur cet antibio
        if ((a + b) == 0) next 
        
        # Matrice
        matrix_fisher <- matrix(c(a, c, b, d), nrow = 2)
        
        # B. Test de Fisher
        # On utilise tryCatch car fisher.test peut échouer si données bizarres
        res_test <- tryCatch({
          fisher.test(matrix_fisher, conf.int = TRUE)
        }, error = function(e) return(NULL))
        
        if (!is.null(res_test)) {
          
          odds_ratio <- as.numeric(res_test$estimate)
          p_val <- res_test$p.value
          ci_low <- res_test$conf.int[1]
          ci_high <- res_test$conf.int[2]
          
          # Interprétation du risque
          risk_label <- "Neutre"
          if (p_val < 0.05) {
            if (odds_ratio > 1) risk_label <- "FACTEUR DE RISQUE"
            if (odds_ratio < 1) risk_label <- "Facteur Protecteur"
          }
          
          # On ne sauvegarde que si on a au moins une erreur quelque part (sinon inintéressant)
          if ((a + c) > 0) {
            results_list[[counter]] <- data.frame(
              Type_Donnee = sheet_name,
              Categorie = cat,
              Antibiotique = ab,
              Type_Erreur = err,
              Nb_Erreurs_Groupe = a,
              Total_Tests_Groupe = a+b,
              Nb_Erreurs_Reste = c,
              Total_Tests_Reste = c+d,
              Odds_Ratio = round(odds_ratio, 2),
              IC_Min = round(ci_low, 2),
              IC_Max = round(ci_high, 2),
              P_Value = p_val,
              Interpretation = risk_label
            )
            counter <- counter + 1
          }
        }
      }
    }
  }
  
  return(do.call(rbind, results_list))
}

# --- Fonction Principale ---
run_risk_factor_analysis <- function(input_file, output_file) {
  
  wb_source <- loadWorkbook(input_file)
  sheet_names <- names(wb_source)
  wb_results <- createWorkbook()
  
  for (sheet in sheet_names) {
    message("Calcul des facteurs de risque (OR) pour : ", sheet)
    
    data <- read.xlsx(wb_source, sheet = sheet)
    res_df <- calculate_odds_ratios(data, sheet)
    
    if (!is.null(res_df)) {
      # Filtrer pour mettre en avant les résultats significatifs
      # On trie par P-value pour voir les "vrais" risques en haut
      res_df <- res_df %>%
        arrange(P_Value) %>%
        mutate(
          Signif = case_when(
            P_Value < 0.001 ~ "***",
            P_Value < 0.01  ~ "**",
            P_Value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
        )
      
      # Sauvegarde
      sheet_out <- substr(paste0("Risk_", sheet), 1, 31)
      addWorksheet(wb_results, sheet = sheet_out)
      writeData(wb_results, sheet = sheet_out, x = res_df)
      
      # Coloration automatique des Facteurs de Risque (Rouge) et Protecteurs (Vert)
      # Colonne 13 = Interpretation
      conditionalFormatting(wb_results, sheet_out, cols = 13, rows = 2:(nrow(res_df)+1),
                            rule = '=="FACTEUR DE RISQUE"', 
                            style = createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE"))
      
      conditionalFormatting(wb_results, sheet_out, cols = 13, rows = 2:(nrow(res_df)+1),
                            rule = '=="Facteur Protecteur"', 
                            style = createStyle(fontColour = "#006100", bgFill = "#C6EFCE"))
    }
  }
  
  saveWorkbook(wb_results, file = output_file, overwrite = TRUE)
  message("Terminé ! Analyse des risques sauvegardée : ", output_file)
}

# --- Exécution ---
input_path <- "./tables/error_distribution_combined.xlsx"
output_path <- "./tables/risk_factors_odds_ratios.xlsx"

if (file.exists(input_path)) {
  run_risk_factor_analysis(input_path, output_path)
} else {
  message("Attention : Fichier introuvable.")
}

################################################################################
################################################################################
################################################################################


# Exemple d'utilisation
df <- read_excel("./data/Tableaux_final_analyses_pour_R.xlsx",
                 sheet = "data")


# Remplacer les virgules par des points dans toutes les colonnes de type caractère
df <- df %>%
  mutate(across(where(is.character), ~ gsub(",", ".", .)))

transformed_data <- transform_data_for_cmi(df)
write.table(transformed_data, "./data/transformed_cmi_data.csv",
            row.names = FALSE,
            sep = ";")

################################################################################
################################################################################

df <- read.csv("./data/transformed_cmi_data.csv", sep = ";")


antibiotics <- list(
  list(columns = c("PTZ_VITEK_corrected", "PTZ_BMD_corrected"), threshold = 8),
  list(columns = c("CZD_VITEK_corrected", "CZD_BMD_corrected"), threshold = 4),
  list(columns = c("CRO_VITEK_corrected", "CRO_BMD_corrected"), threshold = 2),
  list(columns = c("FEP_VITEK_corrected", "FEP_BMD_corrected"), threshold = 4),
  list(columns = c("ETP_VITEK_corrected", "ETP_BMD_corrected"), threshold = 0.5),
  list(columns = c("IMP_VITEK_corrected", "IMP_BMD_corrected"), threshold = 4),
  list(columns = c("AKN_VITEK_corrected", "AKN_BMD_corrected"), threshold = 8),
  list(columns = c("LVX_VITEK_corrected", "LVX_BMD_corrected"), threshold = 1),
  list(columns = c("NFE_VITEK_corrected", "NFE_BMD_corrected"), threshold = 64),
  list(columns = c("SXT_VITEK_corrected", "SXT_BMD_corrected"), threshold = 80),
  list(columns = c("AMX_VITEK_corrected", "AMX_BMD_corrected"), threshold = 8),
  list(columns = c("AMC_VITEK_corrected", "AMC_BMD_corrected"), threshold = 8),
  list(columns = c("FIX_VITEK_corrected", "FIX_BMD_corrected"), threshold = 1),
  list(columns = c("FOX_VITEK_corrected", "FOX_BMD_corrected"), threshold = 8),
  list(columns = c("FOS_VITEK_corrected", "FOS_BMD_corrected"), threshold = 32),
  list(columns = c("MEC_VITEK_corrected", "MEC_BMD_corrected"), threshold = 8)
)


# Appeler la fonction pour analyser plusieurs antibiotiques
analyze_multiple_antibiotics(df,
                             antibiotics,
                             molecule_name = "Antibiotics",
                             output_file = "./tables/resultats_antibiotiques.xlsx",
                             output_dir = "./figures/")




################################################################################
################################# Analyses stratifiées sur EA et biais #########
################################################################################

# Colonnes pour la stratification
strata_columns <- c("Groupe", "Especes", "Phenotype")

# Appeler la fonction pour analyser plusieurs antibiotiques par strate
analyze_multiple_antibiotics_by_strata_types_to_txt(df,
                                       antibiotics,
                                       strata_columns,
                                       output_dir = "./tables/")



################################################################################
################################# Analyses interprétations #####################
################################################################################


# Définition des antibiotiques
antibiotics_inter <- list(
  list(columns = c("PTZ_VITEK", "PTZ_BMD")),
  list(columns = c("CZD_VITEK", "CZD_BMD")),
  list(columns = c("CRO_VITEK", "CRO_BMD")),
  list(columns = c("FEP_VITEK", "FEP_BMD")),
  list(columns = c("ETP_VITEK", "ETP_BMD")),
  list(columns = c("IMP_VITEK", "IMP_BMD")),
  list(columns = c("AKN_VITEK", "AKN_BMD")),
  list(columns = c("LVX_VITEK", "LVX_BMD")),
  list(columns = c("NFE_VITEK", "NFE_BMD")),
  list(columns = c("SXT_VITEK", "SXT_BMD")),
  list(columns = c("AMX_VITEK", "AMX_BMD")),
  list(columns = c("AMC_VITEK", "AMC_BMD")),
  list(columns = c("FIX_VITEK", "FIX_BMD")),
  list(columns = c("FOX_VITEK", "FOX_BMD")),
  list(columns = c("FOS_VITEK", "FOS_BMD")),
  list(columns = c("MEC_VITEK", "MEC_BMD"))
)



# Appeler la fonction pour analyser plusieurs antibiotiques
analyze_multiple_antibiotics_inter_single_sheet(df,
                                   antibiotics_inter,
                                   output_file = "./tables/resultats_antibiotiques_interpretation.xlsx")



################################################################################
################################################################################


main_analysis("./data/transformed_cmi_data.csv")
