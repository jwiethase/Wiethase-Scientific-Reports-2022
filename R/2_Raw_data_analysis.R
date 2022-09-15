library(tidyverse)
library(raster)
library(gridExtra)
library(grid)
library(car)

fileprefix = "BG200m"
theme_set(theme_gray(base_size = 11))
dir.create(file.path("output", fileprefix), showWarnings = FALSE)


# Load and prepare data
# --------------------------------------------------------------------------------------------------------
annRain <- read.csv("data/PrecipitationYearly_99-20.csv") %>% 
   mutate(year = as.factor(1999:2020)) %>% 
   filter(year != 1999)
load(paste0("data/", fileprefix, "_raw_spdf.RData"))
deg <- do.call(stack, lapply(list.files(path = "data", pattern = '_200m', full.names = T), raster))

deg_change <- subset(deg, 2:nlayers(deg)) - subset(deg, 1:(nlayers(deg)-1))

long_deg <- data.frame("deg" = unlist(raw_spdf@data[, 1:21]), "year" = rep(2000:2020, each = NROW(raw_spdf)), 
                       "CA" = raw_spdf@data[, "CA"],
                       "pop" = raw_spdf@data[, "pop"], "pop_q"= raw_spdf@data[, "pop_q"], "rain" = raw_spdf@data[, "rain"],
                       "rain_q"= raw_spdf@data[, "rain_q"], "tlu"= raw_spdf@data[, "tlu"], "tlu_q"= raw_spdf@data[, "tlu_q"],
                       "deg_q"= raw_spdf@data[, "deg_q"],
                       "deg_c" = c(rep(NA, NROW(raw_spdf)), unlist(raw_spdf@data[,grep("Change", names(raw_spdf))])))

# scale absolute values 0 - 1:
long_deg$deg <- long_deg$deg + (-1 * min(long_deg$deg))
long_deg$deg <- long_deg$deg / max(long_deg$deg)

# scale changes -1 to +1:
long_deg$deg_c <- long_deg$deg_c / max(long_deg$deg_c, na.rm = TRUE)

long_deg$year <- as.factor(long_deg$year)
long_deg$CA <- factor(long_deg$CA, levels = c(0, 1, 3, 2), labels = c("NONE", "CCRO", "WMA", "NP"))
long_deg$deg_q <- factor(long_deg$deg_q, levels = c(-1, 0, 1), labels = c("Low", "Medium", "High"))
long_deg$rain_q <- factor(long_deg$rain_q, levels = c(-1, 0, 1), labels = c("Low", "Medium", "High"))
long_deg$pop_q <- factor(long_deg$pop_q, levels = c(-1, 0, 1), labels = c("Low", "Medium", "High"))
long_deg$tlu_q <- factor(long_deg$tlu_q, levels = c(-1, 0, 1), labels = c("Low", "Medium", "High"))


deg_18_20 <- deg[[18:20]]
deg_mean <- calc(deg_18_20, fun = mean, na.rm = TRUE)

ff <- function(i) {
   p <- which(i >= 0)
   n <- which(i <= 0)
   # positive values
   if (length(p) > 0) {
      i[p] <- i[p] - min(i[p], na.rm=TRUE)
      i[p] <- i[p] / max(i[p])
   }
   # negative values
   if (length(n) > 0) {
      i[n] <- i[n] - max(i[n], na.rm=TRUE)
      i[n] <- i[n] / abs(min(i[n]))
   }
   i
}

deg_mean_norm <- calc(deg_mean, ff)

# Explore the data
# --------------------------------------------------------------------------------------------------------
plot_values <- function(x, varname, show_rain){
   p1 <- ggplot() +
      geom_line(data = annRain, 
                 aes(x = as.factor(year), y = mean/1500, group = 0),
                 col = "black", lwd = 4, alpha = 0.3) +
      geom_line(data = annRain, 
                aes(x = as.factor(year), y = mean/1500, group = 0),
                col = "black", lty = 2, lwd = 0.5) +
      # geom_blank(data = data.frame(year = as.factor(1999:2020), deg = c(0,1)), aes(x = year, y = deg)) +
      # geom_rect(data = annRain, aes(xmin = as.character(1999:2019), xmax = as.character(2000:2020),
      #                               ymin = 0, ymax = 1, fill = mean), alpha = 0.5) +
      labs(fill = "Rainfall") 
   if(show_rain == FALSE ){
      p1 <- p1 + scale_fill_continuous(guide = FALSE)
   }
   p1 <- p1 +
      ggnewscale::new_scale_fill() +
      # stat_boxplot(geom ='errorbar', data = x, aes(x=year, y=deg, fill = get(varname)),
      #              stat_params = list(width = 0.01)) +
      geom_boxplot(data = x, aes(x=year, y=deg, fill = get(varname)),
                   # outlier.alpha = 0.05, outlier.size = 0.4, 
                   outlier.shape = NA,
                   lwd = .3,
                   notch = TRUE) +
      labs(x = "Year", y = "Normalized bare ground score",
           fill = varname) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
            axis.title = element_text(size=10),
            axis.text = element_text(size=10),
            legend.key.size = unit(0.3, 'cm'),
            legend.title = element_text(size=10),
            legend.position = "bottom") +
      scale_y_continuous(name = "Normalized bare ground score",
                         sec.axis = sec_axis(~ . * 1500, name = "Annual rainfall [mm]", breaks = seq(250, 1250, by = 250)), 
                         limits = c(0, 0.85))
   return(p1)
}

CA_plot <- plot_values(long_deg, "CA", show_rain = FALSE) +
   scale_fill_manual(values = c('black', "#E69F00", "#56B4E9", "#009E73")) + 
   xlab("") + labs(fill = "Land use designation:")
deg_q_plot <- plot_values(long_deg, "deg_q", show_rain = TRUE) + labs(fill = "Bare ground percentile") +
   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00"))

CA_plot <- arrangeGrob(CA_plot, top = textGrob("A)", x = unit(0.01, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                                    gp=gpar(col="black", fontsize=13, fontfamily="Arial")))
deg_q_plot <- arrangeGrob(deg_q_plot, top = textGrob("B)", x = unit(0.01, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=13, fontfamily="Arial")))

#### FIGURE 4
CA_deg_q <- grid.arrange(CA_plot, deg_q_plot)


#### Changes and sudden shocks
## Find years with biggest increases in mean value:
deg_change <- subset(deg, 2:nlayers(deg)) - subset(deg, 1:(nlayers(deg) - 1))
annual_mean_change <- cellStats(deg_change, 'mean')
big_increases <- which(rank(annual_mean_change) > 16)

# These numbers are converted to years in the data frame by:
big_increase_years <- 2000 + big_increases

# Plot the years when degradation increased the most:
big_increase_data <- long_deg[long_deg$year %in% big_increase_years,]

#### FIGURE 5
(biggest_increase_years <- ggplot() +
   geom_blank(data = data.frame(year = as.factor(big_increase_years), deg_c= c(-1, 1, -1, 1)), aes(x=year, y=deg_c)) +
   geom_boxplot(data = big_increase_data, aes(x = year, y = deg_c, fill = deg_q), 
                # outlier.alpha = 0.5, outlier.size = 0.6,
                outlier.shape = NA,
                notch = TRUE) +
   labs(x = "Year", y = "Bare ground cover change", fill = "Bare ground\npercentile") +
   theme_minimal() +
   theme(axis.title=element_text(size=10),
         axis.text = element_text(size=10),
         legend.key.size = unit(0.3, 'cm'),
         legend.title = element_text(size=8))  +
      annotate("segment", x = 1, xend = 1.25, 
               y = 0.4, yend = 0.4,
               colour = "black") +
      annotate("segment", x = 1, xend = 1, 
               y = 0.4, yend = 0.38,
               colour = "black") +
      annotate("segment", x = 1.25, xend = 1.25, 
               y =  0.4, yend = 0.38,
               colour = "black") +
      annotate("text", x = 1.125,  y = 0.5, 
               label = "n.s.", size = 4.5) +
      ylim(-0.3, 0.6) +
   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00")))
   

# Statistical test:
mod1 <- aov(deg_c ~ year + deg_q + year*deg_q, data = big_increase_data)
hist(mod1$residuals)
summary(mod1)
car::Anova(mod1)

## Absolute recovery in years following big increase years:
recovery_years <- big_increase_years + 1
recovery_years[recovery_years == 2012] <- 2013  # 2012 was continued decline year

recov_data <- long_deg[long_deg$year %in% recovery_years, ]
(absolute_recovery <- ggplot()  + 
   geom_point(data = annRain[annRain$year %in% recovery_years, ], 
             aes(x = as.factor(recovery_years), y = mean/1500),
             col = "black", fill = "blue", cex = 4, alpha = 0.5, pch = 23) +
   geom_blank(data = data.frame(year = as.factor(recovery_years), deg_c= c(-1, 1, -1, 1)), aes(x=year, y=deg_c)) +
   geom_boxplot(data = recov_data, aes(x=year, y=deg_c, fill = deg_q), 
                # outlier.alpha = 0.5, outlier.size = 0.6,
                outlier.shape = NA,
                notch = TRUE) +
   labs(x = "Year", y = "Bare ground cover change", fill = "Bare ground percentile: ") +
   theme_minimal() +
   theme(axis.title=element_text(size=10),
         axis.text = element_text(size=10),
         legend.key.size = unit(0.5, 'cm'),
         legend.position = 'none') +
   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00")) +
   scale_x_discrete(breaks = as.factor(recovery_years)) +
   scale_y_continuous(name = "Bare ground cover change", breaks = seq(-1, 1, by = 0.5),
                      sec.axis = sec_axis(~ . * 1500, name = "Annual rainfall [mm]", breaks = seq(500, 1500, by = 250)),
                      limits = c(-0.7, 0.9)))

# Statistical test:
mod2 <- aov(deg_c ~ year + deg_q + year*deg_q, data = recov_data)
hist(mod2$residuals)
summary(mod2)
car::Anova(mod2)
pairs(emmeans::emmeans(mod2, ~ deg_q|year, type = "response"))

## Relative recovery:
## Proportion of the decline (positive score in big increase years) that returns in the recovery year
long_recov_years <- long_deg[long_deg$year %in% recovery_years, ]
long_recov_years$rel_recov <- long_deg$deg_c[long_deg$year %in% recovery_years] / (long_deg$deg_c[long_deg$year %in% big_increase_years] * -1)

(relative_recovery <- ggplot()  + 
   geom_blank(data = data.frame(year = as.factor(recovery_years), rel_recov= c(-3, 3, -3, 3)), aes(x=year, y=rel_recov)) +
   geom_boxplot(data = long_recov_years, aes(x=year, y=rel_recov, fill = deg_q), 
                # outlier.alpha = 0.5, outlier.size = 0.6,
                outlier.shape = NA,
                notch = TRUE) +
   labs(x = "Year", y = "Relative recovery", fill = "Bare ground percentile:") +
   theme_minimal() +
   theme(axis.title=element_text(size=10),
         axis.text = element_text(size=10),
         legend.key.size = unit(0.5, 'cm'),
         legend.position = 'bottom') +
   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00")) +
   scale_x_discrete(breaks = as.factor(recovery_years)) +
   scale_y_continuous(name = "Bare ground score change", breaks = seq(-3, 3, by = 1),
                      limits = c(-3, 3),
                      sec.axis = sec_axis(~ . * 300, name = "      ", breaks = seq(0, 0, by = 0), labels = "    "))
   )
# ggsave("output/relative_recovery.png", plot =relative_recovery, dpi = 150)
# 
# Statistical test:
mod3 <- aov(deg_c ~ year + deg_q + year*deg_q, data = long_recov_years)
hist(mod3$residuals)
summary(mod3)
car::Anova(mod3)

pairs(emmeans::emmeans(mod3, ~ deg_q|year, type = "response"))

absolute_recovery <- arrangeGrob(absolute_recovery, top = textGrob("A)", x = unit(0.01, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
relative_recovery <- arrangeGrob(relative_recovery, top = textGrob("B)", x = unit(0.01, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                                                     gp=gpar(col="black", fontsize=18, fontfamily="Arial")))

#### FIGURE 6
abs_rel_recovery <- grid.arrange(absolute_recovery, relative_recovery, nrow = 2)

