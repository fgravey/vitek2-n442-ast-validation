######## travail sur les données epidémiologie du CHU
library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(xlsx)
library(viridis)
library(tidyverse)
library(gtsummary)
library(writexl)
library(scales)


# -----------------------------------------------------------------------------
# DEFINITION DES FONCTIONS
# -----------------------------------------------------------------------------

# Fonction logique basée sur : S, I (incluant SFP), R
classer_erreur <- function(vitek, bmd) {
  
  if (vitek == bmd) {
    # Ex: S/S, R/R ou I/I (incluant SFP/SFP ou SFP/I)
    return("Concordance")
    
  } else if (bmd == "S" & vitek == "R") {
    # La référence dit Sensible, Vitek dit Résistant
    return("Erreur Majeure (ME)")
    
  } else if (bmd == "R" & vitek == "S") {
    # La référence dit Résistant, Vitek dit Sensible (DANGER)
    return("Erreur Très Majeure (VME)")
    
  } else {
    # Tous les autres cas impliquent forcément un "I" (ou ancien SFP) face à un S ou R
    return("Erreur Mineure (mE)")
  }
}

# Fonction pour créer le graphique
creer_graphique_publication <- function(data, variable_x, titre) {
  
  # On s'assure que l'ordre des erreurs est respecté pour l'empilement
  # On inverse l'ordre des niveaux pour que "Concordance" soit en bas ou en haut selon préférence
  # Ici : Concordance en bas, erreurs "empilées" dessus
  data$Type_Erreur <- factor(data$Type_Erreur, 
                             levels = c("Erreur Très Majeure (VME)", 
                                        "Erreur Majeure (ME)", 
                                        "Erreur Mineure (mE)", 
                                        "Concordance"))
  
  ggplot(data, aes(x = !!sym(variable_x), fill = Type_Erreur)) +
    
    # Barres empilées normalisées à 100% (position = "fill")
    geom_bar(position = "fill", width = 0.7) +
    
    # Axe Y en pourcentage
    scale_y_continuous(labels = scales::percent_format()) +
    
    # Application des couleurs définies
    scale_fill_manual(values = couleurs_erreurs) +
    
    # Thème épuré "Journal Scientifique"
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 10), # Italique pour espèces
      axis.title.x = element_blank(), # Pas besoin de titre en bas si les noms sont clairs
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom", # Légende en bas
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(), # Pas de grille verticale
      panel.grid.minor = element_blank()
    ) +
    
    # Titres
    labs(
      title = titre,
      y = "Proportion des souches",
      subtitle = "Comparaison VITEK vs BMD (Référence)"
    )
}


# Fonction pour extraire la valeur numérique des CMI
# On enlève les <, >, = pour pouvoir placer le point sur le graphique
clean_mic_values <- function(x) {
  # On remplace les virgules par des points (au cas où)
  x <- str_replace(x, ",", ".")
  # On extrait uniquement les chiffres
  nums <- as.numeric(str_extract(x, "[0-9]+(\\.[0-9]+)?"))
  return(nums)
}



# -----------------------------------------------------------------------------
# DEFINITION DES COULEURS
# -----------------------------------------------------------------------------


# 1. Définition des couleurs standards "Publication"
# Vert ou Gris pour Concordance, Jaune pour Mineure, 
# Orange pour Majeure, Rouge vif pour Très Majeure (Danger)
couleurs_erreurs <- c(
  "Concordance" = "#e0e0e0",          # Gris clair (neutre)
  "Erreur Mineure (mE)" = "#f1c40f",  # Jaune
  "Erreur Majeure (ME)" = "#e67e22",  # Orange
  "Erreur Très Majeure (VME)" = "#c0392b" # Rouge foncé
)


# ---------------------------------------------------------
# 1. IMPORTATION ET NETTOYAGE DES DONNÉES
# ---------------------------------------------------------
df <- read_excel("./data/Tableaux_final_analyses_pour_R.xlsx", sheet = "data")


################################################################################
########################## 1 - TAZO ################################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(PTZ.VITEK.inter != "ZIT" & PTZ.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(PTZ.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(PTZ.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau croisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}



# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité croissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_TZP.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------

# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/PTZ_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/PTZ_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------


data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(PTZ.VITEK.res),
    res_BMD_num = clean_mic_values(PTZ.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(1), log2(128)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Tazocilline)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/PTZ_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")

################################################################################
########################## 2 - Ceftazidime #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(CZD.VITEK.inter != "ZIT" & CZD.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(CZD.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(CZD.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau croisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}



# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité croissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_CZD.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------


# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/CZD_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/CZD_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------


data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(CZD.VITEK.res),
    res_BMD_num = clean_mic_values(CZD.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.0625), log2(128)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Ceftazidime)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/CZD_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 3 - CRO  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(CRO.VITEK.inter != "ZIT" & CRO.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(CRO.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(CRO.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau croisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}



# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité croissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_CRO.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------

# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/CRO_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/CRO_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------


data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(CRO.VITEK.res),
    res_BMD_num = clean_mic_values(CRO.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.03125), log2(128)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Ceftriaxone)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/CRO_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 4 - FEP  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(FEP.VITEK.inter != "ZIT" & FEP.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(FEP.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(FEP.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau FEPisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}



# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité FEPissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_FEP.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------

# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/FEP_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/FEP_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------


data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(FEP.VITEK.res),
    res_BMD_num = clean_mic_values(FEP.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.0625), log2(128)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Essential agreement (±1 dilution)")) +
  scale_size_continuous(name = "Number of strains") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram, comparaison of cefepime MICs",
    subtitle = "Dotted lines: Area of essential agreement (± 1 dilution)",
    x = "MICs with BMD (reference method) [mg/L]",
    y = "MICs with VITEK [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/scattergram_cefepime.png", plot = plot_scatter, width = 10, height = 10, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")

################################################################################
########################## 5 - ETP  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(ETP.VITEK.inter != "ZIT" & ETP.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(ETP.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(ETP.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau ETPisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}



# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité ETPissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_ETP.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------

# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/ETP_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/ETP_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------


data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(ETP.VITEK.res),
    res_BMD_num = clean_mic_values(ETP.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.0625), log2(8)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Imipenem)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/ETP_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 6 - IMP  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(IMP.VITEK.inter != "ZIT" & IMP.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(IMP.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(IMP.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau IMPisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité IMPissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_IMP.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/IMP_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/IMP_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(IMP.VITEK.res),
    res_BMD_num = clean_mic_values(IMP.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.125), log2(32)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Essential agreement (±1 dilution)")) +
  scale_size_continuous(name = "Number of strains") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram, comparaison of imipenem MICs",
    subtitle = "Dotted lines: Area of essential agreement (± 1 dilution)",
    x = "MICs with BMD (reference method) [mg/L]",
    y = "MICs with VITEK [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/IMP_Scattergram_CMI.png", plot = plot_scatter, width = 10, height = 10, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")

################################################################################
########################## 7 - AKN  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(AKN.VITEK.inter != "ZIT" & AKN.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(AKN.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(AKN.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau AKNisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité AKNissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_AKN.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/AKN_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/AKN_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(AKN.VITEK.res),
    res_BMD_num = clean_mic_values(AKN.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.125), log2(32)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très AKNortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Amikacin)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/AKN_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")




################################################################################
########################## 8 - LVX  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(LVX.VITEK.inter != "ZIT" & LVX.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(LVX.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(LVX.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau LVXisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité LVXissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_LVX.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/LVX_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/LVX_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(LVX.VITEK.res),
    res_BMD_num = clean_mic_values(LVX.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.001953125), log2(32)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très LVXortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Levofloxacin)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/LVX_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 9 - NFE  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(NFE.VITEK.inter != "ZIT" & NFE.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(NFE.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(NFE.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau NFEisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité NFEissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_NFE.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")



# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/NFE_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/NFE_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(NFE.VITEK.res),
    res_BMD_num = clean_mic_values(NFE.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.25), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Essential agreement (±1 dilution)")) +
  scale_size_continuous(name = "Number of strains") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram, comparaison of nitrofurantoins MICs",
    subtitle = "Dotted lines: Area of essential agreement (± 1 dilution)",
    x = "MICs with BMD (reference method) [mg/L]",
    y = "MICs with VITEK [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/NFE_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 10, height = 10, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 10 - SXT  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(SXT.VITEK.inter != "ZIT" & SXT.BMD.inter != "ZIT") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(SXT.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(SXT.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau SXTisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité SXTissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_SXT.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/SXT_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/SXT_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(SXT.VITEK.res),
    res_BMD_num = clean_mic_values(SXT.BMD.res),
    res_BMD_num = res_BMD_num*20,
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(10), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très SXTortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Cotrimoxazol)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/SXT_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 11 - AMX #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(AMX.VITEK.inter != "ZIT" & AMX.BMD.inter != "ZIT"
         & AMX.VITEK.inter != "\\" & AMX.BMD.inter != "\\") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(AMX.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(AMX.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau AMXisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité AMXissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_AMX.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/AMX_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/AMX_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(AMX.VITEK.res),
    res_BMD_num = clean_mic_values(AMX.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.5), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très AMXortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Amoxicillin)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/AMX_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")





################################################################################
########################## 12 - AMC  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(AMC.VITEK.inter != "ZIT" & AMC.BMD.inter != "ZIT"
         & AMC.VITEK.inter != "\\" & AMC.BMD.inter != "\\") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(AMC.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(AMC.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau AMCisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité AMCissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_AMC.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/AMC_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/AMC_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(AMC.VITEK.res),
    res_BMD_num = clean_mic_values(AMC.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.25), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très AMCortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Amoxicillin - clavulanic acid)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/AMC_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")




################################################################################
########################## 13 - FIX  #########################################
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(FIX.VITEK.inter != "ZIT" & FIX.BMD.inter != "ZIT"
         & FIX.VITEK.inter != "\\" & FIX.BMD.inter != "\\") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(FIX.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(FIX.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau FIXisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité FIXissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_FIX.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/FIX_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/FIX_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(FIX.VITEK.res),
    res_BMD_num = clean_mic_values(FIX.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.015625), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Essential agreement (±1 dilution)")) +
  scale_size_continuous(name = "Number of strains") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très important pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram, comparaison of cefixime MICs",
    subtitle = "Dotted lines: Area of essential agreement (± 1 dilution)",
    x = "MICs with BMD (reference method) [mg/L]",
    y = "MICs with VITEK [mg/L]"
  )
# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/FIX_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 10, height = 10, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 14 - FOX ######################################### 
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(FOX.VITEK.inter != "ZIT" & FOX.BMD.inter != "ZIT"
         & FOX.VITEK.inter != "\\" & FOX.BMD.inter != "\\") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(FOX.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(FOX.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau FOXisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité FOXissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_FOX.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/FOX_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/FOX_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(FOX.VITEK.res),
    res_BMD_num = clean_mic_values(FOX.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.125), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très FOXortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (Cefoxitin)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/FOX_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")



################################################################################
########################## 11 FOS ######################################### 
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(FOS.VITEK.inter != "ZIT" & FOS.BMD.inter != "ZIT"
         & FOS.VITEK.inter != "\\" & FOS.BMD.inter != "\\") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(FOS.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(FOS.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau FOSisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité FOSissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_FOS.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/FOS_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/FOS_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(FOS.VITEK.res),
    res_BMD_num = clean_mic_values(FOS.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.0625), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très FOSortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (CeFOSim)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/FOS_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")




################################################################################
########################## 15 - MEC #########################################   
################################################################################

data_clean <- df %>%
  # 1. Exclusion des résultats "ZIT" (ni l'un ni l'autre)
  filter(MEC.VITEK.inter != "ZIT" & MEC.BMD.inter != "ZIT"
         & MEC.VITEK.inter != "\\" & MEC.BMD.inter != "\\") %>%
  
  # 2. Normalisation et Conversion SFP -> I
  mutate(
    # On nettoie les espaces et on met en majuscules
    raw_vitek = str_to_upper(str_trim(MEC.VITEK.inter)),
    raw_bmd   = str_to_upper(str_trim(MEC.BMD.inter)),
    
    # Règle : Si "SFP", alors on remplace par "I", sinon on garde la valeur
    VITEK = ifelse(raw_vitek == "SFP", "I", raw_vitek),
    BMD   = ifelse(raw_bmd == "SFP", "I", raw_bmd)
  )

cat("Nombre de souches analysées :", nrow(data_clean), "\n")

# ---------------------------------------------------------
# 2. CLASSIFICATION DES ERREURS
# ---------------------------------------------------------
# Application de la classification
data_analyzed <- data_clean %>%
  rowwise() %>%
  mutate(Type_Erreur = classer_erreur(VITEK, BMD)) %>%
  ungroup() %>%
  # Création d'une variable binaire pour les stats (Discordance Oui/Non)
  mutate(Is_Discordant = ifelse(Type_Erreur == "Concordance", "Non", "Oui"))

# ---------------------------------------------------------
# 3. CALCULS DES TAUX
# ---------------------------------------------------------

# Tableau récapitulatif des erreurs (Nombre et %)
synthese_erreurs <- data_analyzed %>%
  count(Type_Erreur) %>%
  mutate(
    Total = sum(n),
    Pourcentage = round(n / Total * 100, 2)
  )

print("--- RÉPARTITION DES RÉSULTATS ---")
print(synthese_erreurs)

# Calcul précis du taux de concordance global
taux_global <- sum(data_analyzed$Type_Erreur == "Concordance") / nrow(data_analyzed) * 100
cat(paste0("\nTAUX DE CONCORDANCE GLOBAL : ", round(taux_global, 2), " %\n"))

# ---------------------------------------------------------
# 4. ANALYSE STATISTIQUE DES FACTEURS D'INFLUENCE
# ---------------------------------------------------------
# Test du Chi-2 ou Fisher pour voir si l'erreur dépend de la variable

variables_a_tester <- c("Especes", "Groupe", "Phénotype", "Origine")

print("\n--- INFLUENCE DES VARIABLES SUR LA DISCORDANCE ---")

res_stats <- list() # Pour stocker les résultats si besoin

for (var in variables_a_tester) {
  
  # Création du tableau MECisé : Variable vs Discordance (Oui/Non)
  tbl <- table(data_analyzed[[var]], data_analyzed$Is_Discordant)
  
  # Test statistique (Fisher est plus robuste si tu as peu d'échantillons dans une catégorie)
  test_stat <- tryCatch({
    chisq.test(tbl)
  }, warning = function(w) {
    fisher.test(tbl, simulate.p.value = TRUE)
  })
  
  p_val <- test_stat$p.value
  signif <- ifelse(p_val < 0.05, "OUI (Significatif)", "NON")
  
  cat(paste0("\nVariable : ", var))
  cat(paste0("\n  -> P-value : ", format.pval(p_val, digits = 3, eps = 0.001)))
  cat(paste0("\n  -> Influence significative ? ", signif, "\n"))
}

# ------------------------------------------------------------------
# 1. PRÉPARATION : ORDONNER LES ERREURS
# ------------------------------------------------------------------
# On définit un ordre logique pour les colonnes du tableau : 
# D'abord le succès (Concordance), puis les erreurs par gravité MECissante.
data_analyzed$Type_Erreur <- factor(data_analyzed$Type_Erreur, 
                                    levels = c("Concordance", 
                                               "Erreur Mineure (mE)", 
                                               "Erreur Majeure (ME)", 
                                               "Erreur Très Majeure (VME)"))

# Liste des variables que vous voulez analyser
# (Assurez-vous que les noms correspondent à vos colonnes nettoyées)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

# ---------------------------------------------------------
# PRÉPARATION DE L'EXPORT EXCEL
# ---------------------------------------------------------

# 1. On initialise la liste avec votre tableau de synthèse global
# (Cela créera le premier onglet du fichier Excel)
liste_resultats <- list("Synthese_Globale" = synthese_erreurs)

print("--- GÉNÉRATION DES DONNÉES ---")

# 2. On ajoute les autres onglets via la boucle (Détails par variable)
vars_to_analyze <- c("Especes", "Groupe", "Phénotype", "Origine")

for (var in vars_to_analyze) {
  
  col_name <- grep(var, names(data_analyzed), value = TRUE)[1]
  
  if (!is.na(col_name)) {
    
    # --- Création du tableau gtsummary ---
    mes_labels <- list()
    mes_labels[[col_name]] <- var
    
    table_gtsummary <- data_analyzed %>%
      select(all_of(col_name), Type_Erreur) %>% 
      tbl_summary(
        by = Type_Erreur,
        percent = "row",
        label = mes_labels,
        statistic = all_categorical() ~ "{n} ({p}%)"
      ) %>%
      add_n() %>%
      add_p() %>%
      modify_header(label = "**Variable**") 
    
    # --- Conversion en format brut pour Excel ---
    table_brute <- as_tibble(table_gtsummary)
    
    # --- Ajout à la liste des onglets ---
    # On coupe le nom à 30 caractères max (limite Excel)
    nom_onglet <- substr(col_name, 1, 30) 
    liste_resultats[[nom_onglet]] <- table_brute
    
    cat(paste0("Ajouté à l'export : ", nom_onglet, "\n"))
  }
}

# 3. EXPORT FINAL VERS EXCEL
# Cela va créer un fichier avec un onglet par variable
write_xlsx(liste_resultats, path = "./tables/Resultats_Antibio_MEC.xlsx")

cat("\nTERMINE ! Le fichier 'Resultats_Antibio_Complets.xlsx' a été créé.")


# ---------------------------------------------------------
# FIGURE : PROFIL D'ERREUR PAR GROUPE (STYLE AAC/CMI)
# ---------------------------------------------------------
# ---------------------------------------------------------
# 2. GÉNÉRATION DES FIGURES
# ---------------------------------------------------------

# Figure 1 : Par Espèce (Idéal pour voir si une espèce pose problème)
p1 <- creer_graphique_publication(data_analyzed, "Especes", "Performance par Espèce Bactérienne")
print(p1)

# Figure 2 : Par Phénotype (Pour voir si les Carbapénémases sont mal détectées)
p2 <- creer_graphique_publication(data_analyzed, "Phénotype", "Performance par Phénotype de Résistance")
print(p2)

# ---------------------------------------------------------
# 4. SAUVEGARDE EN HAUTE RÉSOLUTION (300 DPI pour publication)
# ---------------------------------------------------------
ggsave("./figures/MEC_Figure_Especes_300dpi.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("./figures/MEC_Figure_Phenotypes_300dpi.png", plot = p2, width = 10, height = 6, dpi = 300)

cat("Graphiques générés et sauvegardés en haute qualité (PNG 300 DPI).")


# ---------------------------------------------------------
# 1. NETTOYAGE DES CMI (CONVERSION EN NOMBRES)
# ---------------------------------------------------------

data_scatter <- data_analyzed %>%
  mutate(
    # Conversion en nombres pour le graphique
    res_VITEK_num = clean_mic_values(MEC.VITEK.res),
    res_BMD_num = clean_mic_values(MEC.BMD.res),
    
    # On crée une étiquette pour dire si c'est concordant ou non (pour la couleur)
    Is_EA = abs(log2(res_VITEK_num) - log2(res_BMD_num)) <= 1
  ) %>%
  # On filtre les éventuels ratés (NA)
  filter(!is.na(res_VITEK_num) & !is.na(res_BMD_num))

# ---------------------------------------------------------
# 2. CRÉATION DU SCATTERGRAM (BUBBLE PLOT)
# ---------------------------------------------------------

# Définition des limites pour le graphique (de la plus petite à la plus grande CMI trouvée)
min_val <- min(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))
max_val <- max(c(data_scatter$res_VITEK_num, data_scatter$res_BMD_num))

# Création des graduations logiques (2, 4, 8, 16, 32, 64...)
breaks_log <- 2^seq(log2(0.0625), log2(512)) # Génère 2, 4, 8, 16... jusqu'à 128

plot_scatter <- ggplot(data_scatter, aes(x = res_BMD_num, y = res_VITEK_num)) +
  
  # A. ZONE D'ACCORD ESSENTIEL (Les lignes diagonales)
  # Ligne centrale (Identité parfaite x=y)
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) + 
  # Lignes pointillées (+/- 1 dilution) = Zone d'acceptabilité technique
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "gray50") + # +1 dilution
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", color = "gray50") + # -1 dilution
  
  # B. LES POINTS (Bubble Plot)
  # geom_count compte combien de souches sont au même endroit et adapte la taille
  geom_count(aes(color = Is_EA), alpha = 0.7) +
  
  # C. ÉCHELLES (Logarithmique base 2 pour les antibios)
  scale_x_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  scale_y_log10(breaks = breaks_log, labels = breaks_log, limits = c(min_val, max_val)) +
  
  # D. COULEURS ET ESTHÉTIQUE
  scale_color_manual(values = c("FALSE" = "#c0392b", "TRUE" = "#27ae60"), 
                     labels = c("Discordance (>1 dilution)", "Accord Essentiel (±1 dilution)")) +
  scale_size_continuous(name = "Nombre de souches") + # Légende pour la taille des bulles
  
  # E. THEME GRAPHIQUE
  theme_bw() +
  theme(
    aspect.ratio = 1, # Force le graphique à être carré (très MECortant pour un scattergram)
    legend.position = "bottom",
    panel.grid.minor = element_blank(), # Enlève les petites grilles inutiles
    axis.title = element_text(face = "bold")
  ) +
  
  # F. TITRES
  labs(
    title = "Scattergram : Comparaison des CMI (CeMECim)",
    subtitle = "Lignes pointillées : Zone d'Accord Essentiel (± 1 dilution)",
    x = "CMI Référence (BMD) [mg/L]",
    y = "CMI Test (VITEK) [mg/L]"
  )

# Affichage
print(plot_scatter)

# ---------------------------------------------------------
# 3. SAUVEGARDE
# ---------------------------------------------------------
ggsave("./figures/MEC_Scattergram_CMI_300dpi.png", plot = plot_scatter, width = 8, height = 8, dpi = 300)
cat("Scattergram sauvegardé sous 'Scattergram_CMI_300dpi.png'")

















#read the two dataFrames
setwd('./data')
res = read.csv('vitekvscmi_maj_gatien.csv', sep = ';')

###################################################################################################################
###################################################################################################################
###################################################################################################################
# CRO
sort(unique(res$CMI_CRO))
sort(unique(res$VITEK_CRO))

ggplot(res, aes(x=VITEK_CRO, fill = factor(CMI_CRO, levels = c('<=0.008', '0.016', '0.03', '0.06', '0.125', '0.5', '1', '2',
                                                                   '4', '8', '16', '>16')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c("#001e00", "#003200", "#004600", "#005900", "#006d00", "#008000", "#009400", "#f0d600", "#cc7000", "#ff8c00", "#be0000", "#a50000"))+
  scale_fill_manual(values=c("#001e00", "#003200", "#004600", "#005900", "#006d00", "#008000", "#009400", "#f0d600", "#cc7000", "#ff8c00", "#be0000", "#a50000"), name = 'CRO BMD MICs')+
  scale_x_discrete(limits = c("<=0.25", "0.5", "1", "2", "4", "8", "16", "32", ">32"))+
  geom_vline(xintercept=3, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=3, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=4, label='Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "VITEK MICs", y = "Number of strains", 
       title = "Distribution of ceftriaxone MICs obtained with vitek and BMD",
       caption = expression("Global agreement between the two methods 91%")) +
  theme_minimal()


ggsave('../figures/ceftriaxone.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)



###################################################################################################################
###################################################################################################################
###################################################################################################################
# CAZ
sort(unique(res$CMI_CAZ))
sort(unique(res$VITEK_CAZ))

ggplot(res, aes(x=VITEK_CAZ, fill = factor(CMI_CAZ, levels = c('<=0.5', '1', '2', '4', '8', '16', '32', '>32')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c(  "#008000", "#009400", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"))+
  scale_fill_manual(values=c( "#008000", "#009400", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"), name = 'CAZ BMD MICs')+
  scale_x_discrete(limits = c("<=0.12", "0.25", "0.5", "1", "2", "4", "8", "16", "32", ">32"))+
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=6, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=4, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=6, label='Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distribution of ceftazidime MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 88%")) +
  theme_minimal()


ggsave('../figures/ceftazidime.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)


###################################################################################################################
###################################################################################################################
###################################################################################################################
# FEP
sort(unique(res$CMI_FEP))
sort(unique(res$VITEK_FEP))

ggplot(res, aes(x=VITEK_FEP, fill = factor(CMI_FEP, levels = c('<=0.5', '1', '2', '4', '8', '16', '32', '>32')))) + 
  geom_histogram(stat = 'count', bins = 0.1) +
  scale_color_manual(values=c( "#005900", "#008000", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"))+
  scale_fill_manual(values=c("#005900", "#008000", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"), name = 'FEP BMD MICs')+
  scale_x_discrete(limits = c("<=0.12", "0.25", "0.5", "1", "2", "4", "8", "16", ">16"))+
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=6, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=4, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=6, label='Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distribution of cefepime MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 87%")) +
  theme_minimal()


ggsave('../figures/cefepime.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)


###################################################################################################################
###################################################################################################################
###################################################################################################################
# PTZ

unique(res$VITEK_PTZ)
ptz = res[res$VITEK_PTZ != "\\",]
sort(unique(ptz$CMI_PTZ))
sort(unique(ptz$VITEK_PTZ))

ggplot(ptz, aes(x=VITEK_PTZ, fill = factor(CMI_PTZ, levels = c('<=2', '4', '8', '16', "32",'64', '>64')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#006d00", "#008000", "#009400",  "#ff8c00", "#cc7000", "#be0000", "#a50000"))+
  scale_fill_manual(values=c("#006d00", "#008000", "#009400",  "#ff8c00", "#cc7000", "#be0000", "#a50000"), name = 'PTZ BMD MICs')+
  scale_x_discrete(limits = c("<=4", "8", "16", "32", "64", ">64"))+
  geom_vline(xintercept=2, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=2, label='Sensitive/Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distribution of piperacillin - tazobactam MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 84%")) +
  theme_minimal()


ggsave('../figures/piperacillin_tazobactam.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)


###################################################################################################################
###################################################################################################################
###################################################################################################################
# ertapeneme

unique(res$CMI_ERT)
unique(res$VITEK_ERT)

ggplot(res, aes(x=VITEK_ERT, fill = factor(CMI_ERT, levels = c('<=0.06', '0.12', '0.25', '0.5', '1', '>1')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c(  "#005900", "#006d00", "#008000", "#009400", "#be0000", "#a50000"))+
  scale_fill_manual(values=c( "#005900", "#006d00", "#008000", "#009400", "#be0000", "#a50000"), name = 'Ertapenem BMD MICs')+
  scale_x_discrete(limits = c("<=0.12", "0.25", "0.5", "1", "2", "4", ">4"))+
  geom_vline(xintercept=3, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=3, label='Sensitive/Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distribution of ertapenem MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 93%")) +
  theme_minimal()

ggsave('../figures/ertapenem.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)



###################################################################################################################
###################################################################################################################
###################################################################################################################
# imipenem

unique(res$VITEK_IMP)
imi = res[res$VITEK_IMP != "\\",]
unique(imi$CMI_IMP)
unique(imi$VITEK_IMP)

ggplot(imi, aes(x=VITEK_IMP, fill = factor(CMI_IMP, levels = c('<=0.5', '1', '2', '4')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#006d00", "#008000", "#009400", "#f0d600"))+
  scale_fill_manual(values=c( "#006d00", "#008000", "#009400", "#f0d600"), name = 'Imipenem BMD MICs')+
  scale_x_discrete(limits = c("<=0.25", "0.5", "1", "2", "4", "8", ">8"))+
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=5, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=4, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=5, label='Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distrubution of imipenem MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 81%")) +
  theme_minimal()


ggsave('../figures/imipenem.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)


###################################################################################################################
###################################################################################################################
# levofloxacine

unique(res$CMI_LVX)
unique(res$VITEK_LVX)

ggplot(res, aes(x=VITEK_LVX, fill = factor(CMI_LVX, levels = c( "<=0.01", "0.02", "0.03", "0.06", "0.125", "0.25",  "0.5", "1", "2", "4", "8", ">8")))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#001e00", "#003200", "#004600", "#005900", "#006d00", "#008000", "#009400", "#f0d600", "#cc7000", "#ff8c00", "#be0000", "#a50000"))+
  scale_fill_manual(values=c( "#001e00", "#003200", "#004600", "#005900", "#006d00", "#008000", "#009400", "#f0d600", "#cc7000", "#ff8c00", "#be0000", "#a50000"),
                    name = 'Levofloxacine BMD MICs')+
  scale_x_discrete(limits = c("<=0.12", "0.25", "0.5", "1", "2", "4"))+
  geom_vline(xintercept=3, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=3, label='Sensitive breakpoint', y=27),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=4, label='Resistant breakpoint', y=27),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distrubution of levofloxacine MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 90%")) +
  theme_minimal()


ggsave('../figures/levofloxacine.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)

###################################################################################################################
###################################################################################################################
###################################################################################################################
# cotrimoxazole

unique(res$CMI_SXT)
unique(res$VITEK_SXT)

ggplot(res, aes(x=VITEK_SXT, fill = factor(CMI_SXT, levels = c('<=20','80', ">160")))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#006d00", "#f0d600",  "#a50000"))+
  scale_fill_manual(values=c("#006d00", "#f0d600",  "#a50000"), name = 'SXT BMD MICs')+
  scale_x_discrete(limits = c("<=20", "40", "80", "160", ">160"))+
  geom_vline(xintercept=2, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=3, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=2, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=3, label='Resistant breakpoint', y=20),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distrubution of cotrimoxazole MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 94%")) +
  theme_minimal()


ggsave('../figures/cotrimoxazole.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)

###################################################################################################################
###################################################################################################################
###################################################################################################################

# Furannes

unique(res$CMI_FURANES)
unique(res$VITEK_FURANES)

ggplot(res, aes(x=VITEK_FURANES, fill = factor(CMI_FURANES, levels = c('<=0.5', '1', '2', '4', '8', '16', '32', '64', '128', '256')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c("#000500","#001e00", "#003200", "#004600", "#005900", "#006d00", "#008000", "#009400", "#be0000", "#a50000"))+
  scale_fill_manual(values=c("#000500","#001e00", "#003200", "#004600", "#005900", "#006d00", "#008000", "#009400", "#be0000", "#a50000"),
                    name = 'Nitrofurantoin BMD MICs')+
  scale_x_discrete(limits = c("<=16", "32", "64", "128", "256"))+
  geom_vline(xintercept=3, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=3, label='Sensitive/Resistant breakpoint', y=29),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distrubution of nitrofurantoin MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 65%")) +
  theme_minimal()

ggsave('../figures/nitrofurantoin.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)


###################################################################################################################
###################################################################################################################
# amikacine

unique(res$CMI_AKN)
unique(res$VITEK_AKN)

ak = res[res$VITEK_AKN != "\\",]

unique(ak$CMI_AKN)
unique(ak$VITEK_AKN)

ggplot(ak, aes(x=VITEK_AKN, fill = factor(CMI_AKN, levels = c( "<=1", "2", "4", "8", "16")))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#005900", "#006d00", "#008000", "#009400", "#a50000"))+
  scale_fill_manual(values=c( "#005900", "#006d00", "#008000", "#009400", "#a50000"),
                    name = 'Amikacine BMD MICs')+
  scale_x_discrete(limits = c("<=1", "2", "4", "8", "16"))+
  geom_vline(xintercept=3, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=3, label='Sensitive/Resistant breakpoint', y=29),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  labs(x = "Vitek AST-N442", y = "Number of strains", 
       title = "Distrubution of amikacine MICs obtained with Vitek AST-N442 and BMD",
       caption = expression("Global agreement between the two methods 97%")) +
  theme_minimal()


ggsave('../figures/amikacine.jpeg', width = 40, height = 30, units = 'cm', dpi = 600)


###################################################################################################################
#UNIQUEMENT LES GROUPES 3
###################################################################################################################

#read the two dataFrames
g3 = read.csv('Desktop/vitek/data/greoupe3.csv', sep = ';')

###################################################################################################################
###################################################################################################################
###################################################################################################################
# FEP
unique(res$CMI_FEP)

ggplot(g3, aes(x=VITEK_FEP, fill = factor(CMI_FEP, levels = c('<0,5', '1', '2', '4', '8', '16', '32', '>32')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#005900", "#008000", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"))+
  scale_fill_manual(values=c("#005900", "#008000", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"), name = 'FEP Manual MICs')+
  scale_x_discrete(limits = c("<0,12", "0,25", "0,5", "1", "2", "4", "8", "16", ">16"))+
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=6, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=4, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=6, label='Resistant breakpoints', y=20),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "VITEK MICs", y = "Number of strains", 
       title = "Distrubution between Cefepime MICs obtained by vitek and BMD",
       caption = expression("Breakpoints according to the CASFM-EUCAST 2024")) +
  theme_minimal()



###################################################################################################################
#UNIQUEMENT LES GROUPES 1 et 2
###################################################################################################################

#read the two dataFrames
g12 = read.csv('Desktop/vitek/data/g1g2.csv', sep = ';')

###################################################################################################################
###################################################################################################################
###################################################################################################################
# FEP
unique(res$CMI_FEP)

ggplot(g12, aes(x=VITEK_FEP, fill = factor(CMI_FEP, levels = c('<0,5', '1', '2', '4', '8', '16', '32', '>32')))) + 
  geom_histogram(stat = 'count') +
  scale_color_manual(values=c( "#005900", "#008000", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"))+
  scale_fill_manual(values=c("#005900", "#008000", "#ffe724", "#f0d600",  "#ff8c00", "#cc7000", "#be0000", "#a50000"), name = 'FEP Manual MICs')+
  scale_x_discrete(limits = c("<0,12", "0,25", "0,5", "1", "2", "4", "8", "16", ">16"))+
  geom_vline(xintercept=4, color = "grey",linetype="dashed", size = 0.8) +
  geom_vline(xintercept=6, color = "grey",linetype="dashed", size = 0.8) +
  geom_text(aes(x=4, label='Sensitive breakpoint', y=20),
            colour="grey",angle=90, vjust=2, family = "Arial") +
  geom_text(aes(x=6, label='Resistant breakpoints', y=20),
            colour="grey",angle=90, vjust=2,family = "Arial") +
  labs(x = "VITEK MICs", y = "Number of strains", 
       title = "Distrubution between Cefepime MICs obtained by vitek and BMD",
       caption = expression("Breakpoints according to the CASFM-EUCAST 2024")) +
  theme_minimal()



= 'FEP  MICs')+
