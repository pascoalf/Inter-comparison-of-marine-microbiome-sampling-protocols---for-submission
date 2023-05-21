prok_16S_diversity_rarefied_metadata %>% 
  pivot_longer(cols = c(Species_richness, Shannon_index, Simpson_index), names_to = "Index", values_to = "Score") %>% 
  ggplot(aes(x = effected_volume, y = Score)) +
  geom_point()+
  facet_wrap(~Index, scales = "free")+
  scale_x_log10()
#
plot(prok_16S_diversity_rarefied_metadata$Species_richness,
prok_16S_diversity_rarefied_metadata$Shannon_index)
#
prok_16S_diversity_rarefied_metadata %>% 
  select("Species richness" = Species_richness, 
         "Shannon index" = Shannon_index, 
         "Simpson index" = Simpson_index) %>% 
  plot(main = "MetaB16SV4V5 (rarefied)")
###


grid.arrange(
prok_16S_diversity_rarefied_metadata %>%
  ggplot(aes(x = Species_richness, y = Shannon_index)) + 
  geom_point() + 
  geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
  labs(x = "Species richness",
       y = "Shannon index",
       title = paste("Pearson correlation", 
                     round(cor(prok_16S_diversity_rarefied_metadata$Species_richness,
                               prok_16S_diversity_rarefied_metadata$Shannon_index),
                           digits = 3)))+
  theme(panel.grid = element_blank()),
prok_16S_diversity_rarefied_metadata %>%
  ggplot(aes(x = Species_richness, y = Simpson_index)) + 
  geom_point() + 
  geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
  labs(x = "Species richness",
       y = "Simpson index",
       title = paste("Pearson correlation", 
                     round(cor(prok_16S_diversity_rarefied_metadata$Species_richness,
                               prok_16S_diversity_rarefied_metadata$Simpson_index),
                           digits = 3)))+
  theme(panel.grid = element_blank()),
prok_16S_diversity_rarefied_metadata %>%
  ggplot(aes(x = Simpson_index, y = Shannon_index)) + 
  geom_point() + 
  geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
  labs(x = "Simpson index",
       y = "Shannon index",
       title = paste("Pearson correlation", 
                     round(cor(prok_16S_diversity_rarefied_metadata$Simpson_index,
                               prok_16S_diversity_rarefied_metadata$Shannon_index),
                           digits = 3))) +
  theme(panel.grid = element_blank())
)


















#
shapes <- seq_along(unique(prok_16S_diversity_rarefied_metadata$device))
#
pairs(~ Species_richness + Shannon_index + Simpson_index, 
      data = prok_16S_diversity_rarefied_metadata,
      #labels = colnames(data),  # Variable names
      pch = 20,                 # Pch symbol
      col = qualitative_colors[prok_16S_diversity_rarefied_metadata$size_fraction],  # Background color of the symbol (pch 21 to 25)
      main = "MetaB16SV4V5",    # Title of the plot
      row1attop = TRUE,         # If FALSE, changes the direction of the diagonal
      gap = 1,                  # Distance between subplots
      cex.labels = NULL,        # Size of the diagonal text
      font.labels = 1)
##

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

##
alpha_cor_16S <- prok_16S_diversity_rarefied_metadata %>% 
  mutate(col_vector = case_when(device == "membrane" ~ qualitative_colors[3],
                                device == "sterivex" ~ qualitative_colors[5],
                                TRUE ~ "grey")) %>% 
  mutate(fill_col = case_when(
    effected_volume <= 10 ~ "cyan",
    effected_volume <= 120 ~ "blue2",
    effected_volume <= 1000 ~ "blue4")) %>% 
  mutate(shape_vector = case_when(size_fraction == ">0.22 µm" ~ 22,
                                  size_fraction == "0.22-3 µm" ~ 23,
                                  size_fraction == "3-20 µm" ~ 24,
                                  size_fraction == ">20 µm" ~ 25,
                                  TRUE ~ 26))

