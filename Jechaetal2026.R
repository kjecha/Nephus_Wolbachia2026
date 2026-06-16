# load required packages
library(tidyverse)
library(reshape2)
library(ggpubr)
library(EnvStats)
library(lme4)
library(MASS)



# Figure 1 - Wolbachia titer ####
titeregg.av <- read.csv("Wolbtiter.csv")

titeregg.av <- titeregg.av %>%
  mutate(relab = case_when(
    try((Sex.mode == "Parthenogenesis") ~ "Parthenogenetic Females"),
    try((Sex.mode == "Sex") ~ "Sexual Females")
  ))

titeregg.av$timepoint <- factor(titeregg.av$timepoint,
  levels = c(2, 8, 18), ordered = TRUE
)

fig1 <- ggplot(data = titeregg.av, aes(x = timepoint, y = propwolb, fill = Treatment)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_dodge(width = .75), alpha = 0.5, size = 2) +
  ylim(0, 10) + # 1 outlier excluded
  scale_fill_manual(values = c(
    "white", "gray"
  )) +
  xlab("Time point") +
  ylab("Wolbachia titer ratio") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  facet_grid(. ~ relab) +
  scale_color_gradient(low = "gray", high = "black") +
  guides(fill = guide_legend(order = 1))

# png(filename = "fig1.png", width = 1200, height = 700)
plot(fig1)
# dev.off()

# Figure 2 - egg number and development ####
melt.develop <- read.csv("Nephus_egg_exp.csv")
melt.develop$Group <- factor(melt.develop$Group, levels = c("ACU", "ACM", "ATU", "ATM", "STM", "SCM"))
melt.develop$Treatment <- factor(melt.develop$Treatment, levels = c("Control", "Antibiotic"))
melt.develop$Mate <- factor(melt.develop$Mate, levels = c("Unmated", "Mated"))
melt.develop.asex <- melt.develop |> filter(Sex.mode == "Parthenogenesis")
melt.develop.sex <- melt.develop |> filter(Sex.mode == "Sex")

facet_labels <- c("A.", "B.", "C.", "D.")
names(facet_labels) <- c("1", "2", "3", "4")


fig2 <- ggplot(data = (melt.develop), aes(x = Timepoint, y = Eggs, group = ID)) +
  stat_summary(aes(y = DevRate * 10, group = 1, linetype = str_wrap("Average proportion of eggs that develop/hatch", 20)),
    color = "gray20",
    alpha = 0.4, fun.y = mean, geom = "line", group = 1, size = 1, show.legend = T
  ) +
  stat_summary(aes(y = (Eggs), group = 1, linetype = str_wrap("Average number of eggs laid", 20)),
    color = "gray20",
    fun.y = mean, geom = "line", group = 1, size = 1, show.legend = T
  ) +
  facet_wrap(. ~ Treatment ~ Sex.mode, ncol = 2, nrow = 2, axes = "all") +
  geom_vline(xintercept = 8.5, linetype = "dashed", size = 1) +
  theme_gray(base_size = 20) +
  theme(
    legend.position = "top", legend.key.size = unit(1.5, "cm"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(
      fill = "white",
      colour = "white",
      size = 0.5, linetype = "solid"
    ),
    axis.line = element_line(size = 1, colour = "gray11"),
    strip.background = element_blank(), strip.text = element_blank()
  ) +
  xlab("Time point") +
  stat_n_text(y.pos = 10, angle = 45, size = 5) +
  scale_y_continuous(
    name = "Eggs laid", breaks = c(2, 4, 6, 8, 10), # limits=c(0,20),
    sec.axis = sec_axis(~ . / 10, name = "Proportion of eggs that develop/hatch")
  ) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)) +
  guides(linetype = guide_legend(title = ""))

# png(filename = "fig2.png", width = 1200, height = 700)
plot(fig2)
# dev.off()


# Supp figure 1: titer outlier ####
titeregg.av.out <- titeregg.av
titeregg.av.out$propwolb <- (titeregg.av.out$WSP.Conc...cp.µL...dPCR.reaction. / titeregg.av.out$ACTIN.Conc...cp.µL...dPCR.reaction.) / titeregg.av.out$Original.Sample.Conc..Qubit..ng.uL.

titer.out <- ggplot(data = titeregg.av.out, aes(x = timepoint, y = propwolb, fill = Treatment)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_dodge(width = .75), alpha = 0.5, size = 2) +
  scale_fill_manual(values = c(
    "white", "gray"
  )) +
  xlab("Time point") +
  ylab("Wolbachia titer ratio") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_blank()) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_color_gradient(low = "gray", high = "black") +
  guides(fill = guide_legend(order = 1))
titer.out

# Supp figure 2: titer vs egg number/development ####

eggs.melt.develop <- melt.develop[c("ID", "Timepoint", "Eggs")]
eggs.develop <- dcast(eggs.melt.develop, ID ~ Timepoint)
eggs.develop$AverageEggs <- rowMeans(eggs.develop[2:19], na.rm = TRUE)
# titeregg.av.out <- titeregg.av
# titeregg.av.out[18, 14] <- NA # remove outlier

egnum <- ggplot(data = (titeregg.av |> filter(Treatment == "Control" & Sex.mode == "Parthenogenesis")), aes(x = propwolb, y = AverageEggs)) +
  geom_point() +
  geom_smooth(method = "glm")


dev <- ggplot(data = (titeregg.av |> filter(Treatment == "Control" & Sex.mode == "Parthenogenesis")), aes(x = propwolb, y = AverageDev)) +
  geom_point() +
  geom_smooth(method = "glm")

ggarrange(egnum, dev, nrow = 2)



# Statistics ####
melt.develop <- read.csv("Nephus_egg_exp.csv")
melt.develop$Group <- factor(melt.develop$Group, levels = c("ACU", "ACM", "ATU", "ATM", "STM", "SCM"))
melt.develop$Treatment <- factor(melt.develop$Treatment, levels = c("Control", "Antibiotic"))
melt.develop$Mate <- factor(melt.develop$Mate, levels = c("Unmated", "Mated"))
melt.develop.asex <- melt.develop |> filter(Sex.mode == "Parthenogenesis")
melt.develop.sex <- melt.develop |> filter(Sex.mode == "Sex")

eggs.stm.1 <- filter(melt.develop.sex, between(Timepoint, 3, 8)) |> filter(Group == "STM")
eggs.scm.1 <- filter(melt.develop.sex, between(Timepoint, 3, 8)) |> filter(Group == "SCM")
eggs.ac.1 <- filter(melt.develop.asex, between(Timepoint, 3, 8)) |> filter(grp == "AC")
eggs.at.1 <- filter(melt.develop.asex, between(Timepoint, 3, 8)) |> filter(grp == "AT")

sum.stm <- eggs.stm.1 %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Antibiotic")) |>
  cbind("Sex" = c("Sex"))

sum.scm <- eggs.scm.1 %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Control")) |>
  cbind("Sex" = c("Sex"))

sum.at <- eggs.at.1 %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Antibiotic")) |>
  cbind("Sex" = c("Asex")) |>
  cbind("Mate" = c(rep("Mated", 26), rep("Unmated", 30)))

sum.ac <- eggs.ac.1 %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Control")) |>
  cbind("Sex" = c("Asex")) |>
  cbind("Mate" = c(rep("Mated", 26), rep("Unmated", 30)))

sex.sum <- rbind(sum.scm, sum.stm)
asex.sum <- rbind(sum.ac, sum.at)

# 1. "During the treatment period of the experiment, the treated sexual females did not exhibit any significant change in egg production, egg development rate, or egg hatching rate" ####
sex.sum$Treatment <- factor(sex.sum$Treatment, levels = c("Control", "Antibiotic"))
glm(Eggs ~ Treatment, data = sex.sum, family = poisson) |> summary()
glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ Treatment,
  data = sex.sum, family = quasibinomial
) |> summary()

# 2. a) "We found no significant difference between females that were paired with a male and those without, therefore the two groups were pooled for further analyses." ####
#  b) "the parthenogenetic females given a control diet exhibited a small but significant increase in the number, but no significant change in development rateof their eggs over time"
# *c) "We found that egg production and egg development rate among treated, Wolbachia-infected, parthenogenetic females significantly decreased from the beginning to the end of the antibiotic treatment period"
asex.sum$Treatment <- factor(asex.sum$Treatment, levels = c("Control", "Antibiotic"))
asex.sum$Mate <- factor(asex.sum$Mate, levels = c("Unmated", "Mated"))

glm.nb(Eggs ~ Mate, data = filter(asex.sum, Treatment == "Antibiotic")) |> summary()
glm.nb(Eggs ~ Mate, data = filter(asex.sum, Treatment == "Control")) |> summary()
glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ Mate, data = filter(asex.sum, Treatment == "Antibiotic"), family = quasibinomial) |> summary()
glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ Mate, data = filter(asex.sum, Treatment == "Control"), family = quasibinomial) |> summary()
glm.nb(Eggs ~ Treatment, data = asex.sum) |> summary()
glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ Treatment,
  data = asex.sum, family = quasibinomial
) |> summary()

# 3. no strong effect of age or batch####
asex.sum$age <- melt.develop.asex$age[match(asex.sum$ID, melt.develop.asex$ID)]
asex.sum$batch <- melt.develop.asex$Batch[match(asex.sum$ID, melt.develop.asex$ID)]

asex.sum$batch <- factor(asex.sum$batch)
glm.nb(Eggs ~ Treatment + age + batch, data = asex.sum) |> summary()
glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ Treatment + age + batch,
  data = asex.sum, family = quasibinomial
) |> summary()

# 4. *"the total eggs laid were significantly higher in the second half (time points 14-18) of the observation period than the first half (time points 9-13)" ####
eggs.at.9_13 <- filter(melt.develop.asex, between(Timepoint, 9, 13)) |> filter(grp == "AT")
eggs.at.14_18 <- filter(melt.develop.asex, between(Timepoint, 14, 18)) |> filter(grp == "AT")
sum.at.9_13 <- eggs.at.9_13 %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Antibiotic")) |>
  cbind("Sex" = c("Asex")) |>
  cbind("Mate" = c(rep("Mated", 26), rep("Unmated", 30))) |>
  cbind("Half" = c("First"))
sum.at.14_18 <- eggs.at.14_18 %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Antibiotic")) |>
  cbind("Sex" = c("Asex")) |>
  cbind("Mate" = c(rep("Mated", 26), rep("Unmated", 30))) |>
  cbind("Half" = c("Second"))
all.obs <- rbind(sum.at.9_13, sum.at.14_18)
all.obs$Develop <- all.obs$DevRate * all.obs$Eggs

m1 <- glmer.nb(Eggs ~ Half + (1 | ID), data = all.obs)
summary(m1)

m2 <- glmer(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ Half + (1 | ID), weights = Eggs, family = binomial, data = all.obs)
summary(m2)

# 5. "Wolbachia titer estimates significantly increased with the age since emergence from pupae when the control parthenogenetic females were analyzed" ####
model11td <- glm(propwolb ~ deathage,
  family = gaussian,
  data = (titeregg.av |> filter(Treatment == "Control"))
)

summary(model11td)

# 6. "titer levels in treated, parthenogenetic females were significantly lower than in control females at every time point" ####
titer.R <- titeregg.av |> filter(Sex.mode == "Parthenogenesis")
t_test_result_2 <- t.test(propwolb ~ Treatment, data = (filter(titer.R, timepoint == 2)))
t_test_result_2
t_test_result_8 <- t.test(propwolb ~ Treatment, data = (filter(titer.R, timepoint == 8)))
t_test_result_8
t_test_result_18 <- t.test(propwolb ~ Treatment, data = (filter(titer.R, timepoint == 18)))
t_test_result_18

# 7. a) "In parthenogenetic females, Wolbachia titer estimates demonstrate that the Wolbachia level of treated females decreased over the course of the treatment period (time point 2 and 8)" ####
modelt1 <- glm(propwolb ~ Treatment * timepoint, family = gaussian, data = (filter(titer.R, timepoint == 2 | timepoint == 8)))
summary(modelt1)
# b) "they were still significantly higher than the around zero titer estimates of the uninfected sexual females"
t_test_result_sex <- t.test(propwolb ~ Sex.mode, data = (filter(titeregg.av, timepoint == 8)))
t_test_result_sex
titer.Pc <- titeregg.av |>
  filter(Sex.mode == "Sex") |>
  filter(Treatment == "Control")
mean(titer.Pc$WSP.Conc...cp.µL...dPCR.reaction.)
sd(titer.Pc$WSP.Conc...cp.µL...dPCR.reaction.)
# c) "After the treatment, titer levels rebounded to near control levels by the end of the observation period (time point 18)"
t_test_result_asex <- t.test(propwolb ~ Treatment, data = (filter(titer.R, timepoint == 18)))
t_test_result_asex
# 8. "We saw no effect of Wolbachia titer in the average quantity or average quality of eggs produced by control parthenogenetic female at their last time point"####
eggs.ac.all <- filter(melt.develop.asex, between(Timepoint, 3, 19)) |> filter(grp == "AC")

sum.ac.all <- eggs.ac.all %>%
  group_by(ID) %>%
  summarize(
    DevRate = ((sum(Hatch, na.rm = T) + sum(Develop, na.rm = T)) / sum(Eggs, na.rm = T)),
    Eggs = sum(Eggs, na.rm = T)
  ) |>
  cbind("Treatment" = c("Control")) |>
  cbind("Sex" = c("Asex")) |>
  cbind("Mate" = c(rep("Mated", 26), rep("Unmated", 30)))

titeregg.av.sum <- merge(titer.R, sum.ac.all[c("ID", "Eggs")], by = "ID")
modelave8 <- glm(Eggs ~ propwolb,
  family = poisson,
  data = (titeregg.av.sum |> filter(timepoint == 8))
)
modelave18 <- glm(Eggs ~ propwolb,
  family = poisson,
  data = (titeregg.av.sum |> filter(timepoint == 18))
)
summary(modelave8)
summary(modelave18)

dev.melt.develop <- melt.develop[c("ID", "Timepoint", "DevRate")]
dev.develop <- dcast(dev.melt.develop, ID ~ Timepoint)
dev.develop$AverageDev <- rowMeans(dev.develop[2:19], na.rm = TRUE)
titeregg.avd <- merge(titeregg.av.sum, sum.ac.all[c("ID", "DevRate")], by = "ID")

modelave8 <- glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ propwolb, family = quasibinomial, data = (titeregg.avd |> filter(timepoint == 8)))
modelave18 <- glm(cbind(DevRate * Eggs, Eggs - DevRate * Eggs) ~ propwolb, family = quasibinomial, data = (titeregg.avd |> filter(timepoint == 18)))
summary(modelave8)
summary(modelave18)
