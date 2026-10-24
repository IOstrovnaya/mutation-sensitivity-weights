#________________________________________________________
#___ A Weighted Estimation Approach to Correct for 
#      Imperfect Mutation Detection Due to Intra-tumor 
#      Heterogeneity in Tumor Sequencing
#___ Audrey Mauguen, Irina Ostrovnaya
#___ October 2026
#________________________________________________________


# Function to represent the simulation results
# Example for N=400; can be adapted for N=200 or N=1000

#___ Packages
# general
library(dplyr)
# Graphs
library(ggplot2)
library(patchwork)

# theme for plot
theme_all <-  theme_bw(base_size = 15) + theme(strip.background = element_rect(fill = "white",  # Background color
                                                                               color = "gray",  # Border color
                                                                               linewidth = 0.5,            # Border size
                                                                               linetype = "solid"   # Border line type
),
strip.text = element_text(color = "black",  # Text color
                          face = "bold",       # Text face
                          size = 14            # Text size
),
legend.position = "bottom")




load("tres.Rdata")
tres <- as.data.frame(tres)

tres$th0 <- round(tres$th0, 2)
tres$th1 <- round(tres$th1, 2)

# Different scenarios to represent
tres$scenario <- paste0("pm=", tres$pm, "; pr=", tres$pr, "; th0=", tres$th0, "; th1=", tres$th1)

# Keep N= 400
tres <- tres %>% filter(N == 400)

# Formatting data for plots
tp <- rbind(cbind.data.frame(Method = "Weighted", pm = tres$pm, OR = tres$OR, th0 = tres$th0, th1 = tres$th1, power = tres$power.uncond, estOR = tres$or.uncond),
            cbind.data.frame(Method = "Unweighted", pm = tres$pm, OR = tres$OR, th0 = tres$th0, th1 = tres$th1, power = tres$power.nw, estOR = tres$or.nw))

tp <- tp %>% mutate(Positive = case_when(estOR >= 1 ~ "Estimated OR ≥ 1",
                                         estOR < 1 ~ "Estimated OR <1"),
                    Case = case_when(Method == "Weighted" & estOR >= 1 ~ "Weighted; estimated OR ≥ 1",
                                     Method == "Weighted" & estOR < 1 ~ "Weighted; estimated OR < 1",
                                     Method == "Unweighted" & estOR >= 1 ~ "Unweighted; estimated OR ≥ 1",
                                     Method == "Unweighted" & estOR < 1 ~ "Unweighted; estimated OR < 1"))
tp$th0 <- paste0("θ", "0", " = ", tp$th0)
tp$th1 <- paste0("θ", "1", " = ", tp$th1)


#__ Power plot
set.seed(161)

# 10% mutation rate
p10 <- tp %>% filter(pm == "0.1") %>%
  ggplot(aes(x = as.factor(OR), y = power, col = Method, group = Method)) + geom_jitter(size = 3, width=0.1, height=0) +
  geom_line(linewidth = 1.3, data = tp %>% filter(pm == "0.1", OR > 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(th1~th0) +
  xlab("True OR") + ylab("Type I error and Power") +
  theme_all +  
  ggtitle("Mutation prevalence = 10% (N=400)")

# 50% mutation rate
p50 <- tp %>% filter(pm == "0.5") %>%
  ggplot(aes(x = as.factor(OR), y = power, col = Method, group = Method)) + geom_jitter(size = 3, width=0.1, height=0) +
  geom_line(linewidth = 1.3, data = tp %>% filter(pm == "0.5", OR > 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(th1~th0) +
  xlab("True OR") + ylab("") +
  theme_all +  
  ggtitle("Mutation prevalence = 50%  (N=400)")

# Figure
p10 | p50


#_ OR plot
set.seed(159)

# Calculate the 95% confidence intervals for the Odds Ratios
tor <- tres %>%
  mutate(estOR.uncond = exp(b1.uncond),
         estOR.low.uncond = exp(b1.uncond - 1.96*b1.se.uncond),
         estOR.high.uncond = exp(b1.uncond + 1.96*b1.se.uncond),
         estOR.nw = exp(b1.nw),
         estOR.low.nw = exp(b1.nw - 1.96*b1.se.nw),
         estOR.high.nw = exp(b1.nw + 1.96*b1.se.nw))

# Labels for the facets
tor$th0 <- paste0("θ", "0", " = ", tor$th0)
tor$th1 <- paste0("θ", "1", " = ", tor$th1)

# 10% mutation rate
or10 <- tor %>% filter(pm == 0.1) %>%
  ggplot(aes(x = OR)) + 
  geom_point(aes(x = OR + 0.05, y = estOR.uncond, col = "Weighted"), alpha = 1, size = 2) +
  geom_errorbar(aes(x = OR + 0.05, ymin = estOR.low.uncond, ymax = estOR.high.uncond, col = "Weighted"), width = 0.01, alpha = 0.5,
                linewidth = 1.3) + 
  geom_point(aes(y = estOR.nw, col = "Unweighted"), alpha = 1, size = 2) +
  geom_errorbar(aes(ymin = estOR.low.nw, ymax = estOR.high.nw, col = "Unweighted"), width = 0.01, alpha = 0.5,
                linewidth = 1.3) + 
  # scale_color_manual(name = "Simulated p0", values = msk_palette("contrast")) +
  scale_color_manual(name = "Method",
                     breaks = c("Unweighted", "Weighted"),
                     values = c("Unweighted" = "#007CBA", "Weighted" = "#DF4602")) +
  scale_x_continuous(breaks = c(1, 1.5, 2)) +
  #  coord_cartesian(ylim = c(0, 7)) +
  facet_grid(th1~th0, scales = "free_y") +
  geom_errorbar(data = tor, aes(y = NULL, ymax = OR, ymin = OR), color = "black") +
  xlab("True OR") + ylab("Mean estimate with 95% confidence interval") +
  ggtitle("Mutation prevalence = 10%  (N=400)") +
  theme_all   

# 50% mutation rate
or50 <- tor %>% filter(pm == 0.5) %>%
  ggplot(aes(x = OR)) + 
  geom_point(aes(x = OR + 0.05, y = estOR.uncond, col = "Weighted"), alpha = 1, size = 2) +
  geom_errorbar(aes(x = OR + 0.05, ymin = estOR.low.uncond, ymax = estOR.high.uncond, col = "Weighted"), width = 0.01, alpha = 0.5,
                linewidth = 1.3) + 
  geom_point(aes(y = estOR.nw, col = "Unweighted"), alpha = 1, size = 2) +
  geom_errorbar(aes(ymin = estOR.low.nw, ymax = estOR.high.nw, col = "Unweighted"), width = 0.01, alpha = 0.5,
                linewidth = 1.3) + 
  # scale_color_manual(name = "Simulated p0", values = msk_palette("contrast")) +
  scale_color_manual(name = "Method",
                     breaks = c("Unweighted", "Weighted"),
                     values = c("Unweighted" = "#007CBA", "Weighted" = "#DF4602")) +
  scale_x_continuous(breaks = c(1, 1.5, 2)) +
  #  coord_cartesian(ylim = c(0, 9)) +
  facet_grid(th1~th0, scales = "free_y") +
  geom_errorbar(data = tor, aes(y = NULL, ymax = OR, ymin = OR), color = "black") +
  xlab("True OR") + ylab("Mean estimate with 95% confidence interval") +
  ggtitle("Mutation prevalence = 50%  (N=400)") +
  theme_all   

# Figure
or10 | or50




