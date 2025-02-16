######## Trossulus dissolution analysis ########

### Load required packages
library(tidyverse)
library(sf)
library(purrr)
library(ggplot2)
library(car)
library(rcompanion)

### Load dataset
tross <- read.csv("/Users/rachelcarlson/Documents/Berkeley/Research/Postdoc-2022-present/Trossulus/dissolutions_alltime.csv")

### Convert all varieties of unpainted shell treatment in californianus from previous dataset to "shell.californianus"
tross$treatment2 <- NA
tross$treatment2 <- ifelse(tross$treatment == "shell" | tross$treatment == "shell.accl" | tross$treatment == "shell.long", "shell.californianus", tross$treatment)
cali <- tross %>% filter(species == "mytilus_californianus")
trossulus <- tross %>% filter(species == "mytilus_trossulus") %>% filter(duration > 24) %>% filter(alk.sd.0 <= 10) %>% filter(alk.sd.1 <= 10)

### Filter for californianus unpainted treatments that came from a consistent sump and were comparable length of incubation (> 24 hours)
cali_unpaint <- cali %>% filter(duration > 24 & sump == 34)

### Filter out mussels with periostracum < 50% in Alisha's dataset (peri was intentionally shaved)
alisha <- read.csv("/Users/rachelcarlson/Documents/Berkeley/Research/Postdoc-2022-present/Trossulus/californianus_painted.csv")
alisha_norm <- alisha %>% filter(SA.perio.perc > 50) %>% select(shell.W)
alisha <- alisha %>% filter(shell.W %in% alisha_norm$shell.W)
### Error in initial dataset where shell.wt is interchanged with shell.W (weight for width). Correct this.
alisha2 <- alisha %>% select(c(shell.W, shell.wt))
colnames(alisha2) <- c("shell.wt","real.shell.wt")
cali_paint <- cali %>% filter(treatment == "shell.californianus.painted") %>% filter(shell.wt %in% alisha_norm$shell.W) %>% right_join(alisha2, by = "shell.wt")
### G for alisha's data was previously normalized by width (on accident) and not weight. Correct this.
cali_paint$real.G <- (cali_paint$G*cali_paint$shell.wt)/cali_paint$real.shell.wt
cali_paint$shell.wt <- cali_paint$real.shell.wt
cali_paint$G <- cali_paint$real.G
cali_paint <- cali_paint %>% select(-c(real.shell.wt, real.G))
cali <- rbind(cali_unpaint, cali_paint)

### Bind it to trossulus samples
tross <- rbind(trossulus, cali)

### Another data error! Trossulus painted shell length is switched with shell weight
tross$real.shell.wt <- ifelse(tross$treatment2 == "shell.trossulus.painted", tross$shell.length, tross$shell.wt)
tross$real.shell.length <- ifelse(tross$treatment2 == "shell.trossulus.painted", tross$shell.wt, tross$shell.length)
tross$real.G <- (tross$G*tross$shell.wt)/tross$real.shell.wt
tross$G <- tross$real.G
tross$shell.wt <- tross$real.shell.wt
tross$shell.length <- tross$real.shell.length

### Plot all data
ggplot(tross, aes(OmegaAragonite, G, color = as.factor(treatment2))) + geom_point() + 
  geom_smooth(data = subset(tross, treatment2 != "shell.californianus.painted" & treatment2 != "shell.trossulus.painted"), aes(OmegaAragonite, G, color = as.factor(treatment2)), formula = y ~ log(x), se = FALSE) +
  scale_color_manual(values=c("#98BAD9","#59638F","#E8101f","#de9a1b")) +
  ylab("Calcification rate") +
  scale_x_continuous(breaks = seq(0, 8, len = 5)) +
  scale_y_continuous(breaks = seq(-0.45, 0, len = 7)) +
  theme_classic()
  
nrow(tross %>% filter(species == "mytilus_trossulus" & treatment2 == "shell.trossulus")) # 45 unpainted trossulus with standard deviation within 10
nrow(tross %>% filter(species == "mytilus_trossulus" & treatment2 == "shell.trossulus.painted")) # 28 painted trossulus with standard deviation within 10
nrow(tross %>% filter(species == "mytilus_californianus" & treatment2 == "shell.californianus")) # 46 californianus (no standard deviation provided)
nrow(tross %>% filter(species == "mytilus_californianus" & treatment2 == "shell.californianus.painted")) #25 painted californianus with periostracum > 50%

### Bootstrapping: compare tross to californianus at omega < 1 (unpainted)
sub <- tross %>% filter((treatment2 == "shell.trossulus.painted" & OmegaAragonite <1) & !is.na(G))
n = length(sub$G)
B = 10000
result = rep(NA, B)
for (i in 1:B) {
  boot.sample = sample(n, replace = TRUE)
  result[i] = mean(sub$G[boot.sample])
}
with(sub, mean(G) + c(-1, 1) * 2 * sd(result))
mean(sub$G)
nrow(sub)

### Bootstrapping plots
boot <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Trossulus/bootstrap_californianus_LOWOMEGA.csv")
boot$Group <- as.factor(boot$Group)
boot$X_axis[1:4] <- c("Californianus unsealed", "Trossulus unsealed", "Californianus sealed", "Trossulus sealed")
ggplot(boot[1:4,], aes(x = Group, y = Y_axis, fill = X_axis)) +
  geom_bar(position = position_dodge(), stat = "identity",
           colour = "black",
           size = 0.3) +
  geom_errorbar(aes(ymin = Y_axis-sd, ymax = Y_axis+sd), width =0.3, position = position_dodge(.9)) +
  ylim(0,0.3) +
  theme_classic() +
  scale_fill_manual(values=c("#98BAD9","#59638F","#E8101f","#de9a1b")) +
  labs(fill = "Treatment group")

### One-way ANOVA
trossSub3 <- tross %>% filter(treatment2 == "shell.dissolution.californianus" | treatment2 == "shell.dissolution.trossulus")
aov_onew <- aov(G ~ species, data = trossSub3)
plot(aov, which = 2)

### Two-way ANOVA
# Prepare data
tross$treatment3 <- ifelse(tross$treatment2 =="shell.trossulus" | tross$treatment2 == "shell.californianus", "unpainted", "painted")
trossSub <- tross %>% filter(OmegaAragonite < 1 & OmegaAragonite > 0.4) # Filter for < omega and within range of californianus painted
table(trossSub$species, trossSub$treatment3)

# Check that data meets the normality assumption of ANOVA
aov <- aov(G ~ species * treatment3,
           data = trossSub)
plot(aov, which = 2)
hist(aov$residuals) # residuals exhibit normality
# Check that data meets the homogeneity of variance assumption of ANOVA
plot(aov, which = 3)
leveneTest(aov) # there is no homogeneity of variance; this assumption is violated. Here, for one-way, would use Kruskal-Wallis - but this doesn't work for two-way.
# aov2 <- scheirerRayHare(G ~ species + treatment3, data=trossSub) # The two-way analogue to Kruskal-Wallis. This non-parametric test is a good step when doing two-way ANOVA and doesn't meet assumptions.

# Log-transform dependent variable to see if this helps data meet ANOVA assumptions
trossSub$D <- -trossSub$G # Calculate dissolution rate as the inverse of calcification (so that you can remove negatives and log-transform)
aov <- aov(log(D) ~ species * treatment3,
           data = trossSub)
plot(aov, which = 2) # Points fall primarily along reference line, indicating normality
hist(aov$residuals)
shapiro.test(residuals(aov)) # p-value = 0.01057; however, if sample size > 50, q-q plot is preferred diagnostic
leveneTest(aov) # p-value = 0.7758, i.e., no significant difference between variances across groups (meets homogeneity of variance assumption)

#linear regression (tested without log transformation and did not pass diagnostics) to account for differences in omega (continuous independent variables not appropriate to ANOVA).
tross$D <- -tross$G
trossSub2 <- tross %>% filter(OmegaAragonite < 1)
table(trossSub2$species, trossSub2$treatment3)
mod <- lm(log(D) ~ (species * treatment3) + OmegaAragonite, data = trossSub2)
plot(mod) # Diagnostics look fine
summary(mod)

# Tabulate group means
aggregate(x= tross$shell.length,
          # Specify group indicator
          by = list(tross$treatment2),      
          # Specify function (i.e. mean)
          FUN = max)
