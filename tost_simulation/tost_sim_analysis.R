library(dplyr)
library(ggplot2)

allfiles <- list.files(path = "tost_simulation/output/", pattern = ".(?i)rds")
y <- readRDS(paste0("tost_simulation/output/", allfiles[1]))
for(i in allfiles[-1]){
  yt <- readRDS(paste0("tost_simulation/output/",i))
  y <- rbind(y, yt)
  
}

y <- y %>% filter(n_pts == 20)

y <- y %>% filter(!(b_ac == 0.1 & b_con == 0.1),
                  !(b_ac == 0.5 & b_con == 0.5),
                  !(b_ac == 0.5 & b_con == 2),
                  !(b_ac == 1.5 & b_con == 1.5),
                  !(b_ac == 1 & b_con == 1),
                  !(b_ac == 2 & b_con == 2))

y <- y %>% mutate(id = row_number(),
                  mod_id = rep(1:(nrow(y)/5), each = 5))

# remove where estimation errors, careful here but based on bad D estimates
y_keep_D <- y %>% filter(parm == "D") %>% 
  mutate(drop_sim_D = (is.na(SE.beta) | (beta < -14) | (beta > -8)) | (abs(SE.beta / beta) > 0.5))

y_keep_ac <- y %>% filter(parm == "D.stdGC") %>% 
  mutate(drop_sim_ac = (is.na(SE.beta) | (SE.beta > 4) | (abs(SE.beta / beta) > 4)))

y_keep_con <- y %>% filter(parm == "noneuc.stdGC") %>% 
  mutate(drop_sim_con = (is.na(SE.beta) | (SE.beta > 4) | (abs(SE.beta / beta) > 4)))

y <- y %>% left_join(y_keep_D %>% dplyr::select(id, drop_sim_D), by = "id") %>%
  left_join(y_keep_ac %>% dplyr::select(id, drop_sim_ac), by = "id") %>%
  left_join(y_keep_con %>% dplyr::select(id, drop_sim_con), by = "id") 

y <- y %>% group_by(mod_id) %>% 
  mutate(drop_sim = (any(drop_sim_D == TRUE, na.rm = TRUE) | 
                       any(drop_sim_ac == TRUE, na.rm = TRUE) |
                       any(drop_sim_con == TRUE, na.rm = TRUE))) %>% ungroup()

yf <- y %>% filter(drop_sim == 0)

y_sum <- yf %>% filter(parm %in% c("D.stdGC", "noneuc.stdGC")) %>% droplevels() %>% 
  group_by(mean_r, b_ac, b_con, parm) %>%
  summarize(mean_beta = median(beta, na.rm = TRUE),
            lcl1 = stats::quantile(beta, probs = 0.025, na.rm = TRUE),
            ucl1 = stats::quantile(beta, probs = 0.975, na.rm = TRUE),
            lcl2 = mean_beta - 2 * median(SE.beta, na.rm = TRUE),
            ucl2 = mean_beta + 2 * median(SE.beta, na.rm = TRUE),
            in_cl = mean(((mean_beta > lcl) & (mean_beta < ucl)), na.rm = TRUE),
            cl_cont_zero = mean(((lcl < 0) & (ucl > 0)), na.rm = TRUE),
            n_valid = sum(!(is.na(beta)|is.nan(beta))),
            mean_nuniq = mean(n),
            mean_recap = mean(r)) %>%
  mutate(true_beta = ifelse(parm == "D.stdGC", b_ac, b_con),
         bias = mean_beta - true_beta)

#y_sum$b_ac <- factor(y_sum$b_ac, c(0.3,0.5,1), c("b[D]==0.3","b[D]==0.5", "b[D]==1", "b[D]==2"))
#y_sum$b_con <- factor(y_sum$b_con, c(0.3,0.5), c("b[H]==0.3","b[H]==0.5", "b[H]==2"))

y_sum$b_ac <- factor(y_sum$b_ac, c(0.3,0.5,1,2), c("b[D]==0.3","b[D]==0.5", "b[D]==1", "b[D]==2"))
y_sum$b_con <- factor(y_sum$b_con, c(0.3,0.5), c("b[H]==0.3","b[H]==0.5"))

p1 <- y_sum %>% ggplot(aes(x = mean_r, y = cl_cont_zero, colour = parm)) +
  geom_point() + geom_line() +
  facet_grid(. ~ b_ac + b_con, scales = "fixed", labeller = label_parsed) + 
  xlab("Number of recaptures") + ylab("Proportion of CIs containing zero") +
  scale_colour_discrete(breaks = c("D.stdGC", "noneuc.stdGC"), 
                        labels = c(bquote(hat(b)[D] ~ "(Density coefficient)"), 
                                   bquote(hat(b)[H] ~ "(Habitat use coefficient)"))) +
  theme(legend.title = element_blank(), legend.position = "bottom")
p1
# save
ggsave("tost_simulation/output/tost-sim-cis.png", p1, width = 9, height = 3.5, dpi = 300)

y_sum <- yf %>% filter(parm %in% c("D.stdGC", "noneuc.stdGC")) %>% droplevels() %>% 
  mutate(true_beta = ifelse(parm == "D.stdGC", b_ac, b_con),
         bias = beta - true_beta) %>%
  group_by(mean_r, b_ac, b_con, parm) %>%
  summarize(mean_bias = median(bias, na.rm = TRUE),
            lcl1 = stats::quantile(bias, probs = 0.025, na.rm = TRUE),
            ucl1 = stats::quantile(bias, probs = 0.975, na.rm = TRUE),
            cl_cont_zero = mean(((lcl < 0) & (ucl > 0)), na.rm = TRUE),
            n_valid = sum(!(is.na(beta)|is.nan(beta))),
            mean_nuniq = mean(n),
            mean_recap = mean(r),
            true_beta = first(true_beta),
            perc_bias = mean_bias / true_beta) 

#y_sum$b_ac <- factor(y_sum$b_ac, c(0.3,0.5,1,2), c("b[D]==0.3","b[D]==0.5", "b[D]==1", "b[D]==2"))
#y_sum$b_con <- factor(y_sum$b_con, c(0.3,0.5,2), c("b[H]==0.3","b[H]==0.5", "b[H]==2"))

y_sum$b_ac <- factor(y_sum$b_ac, c(0.3,0.5,1,2), c("b[D]==0.3","b[D]==0.5", "b[D]==1", "b[D]==2"))
y_sum$b_con <- factor(y_sum$b_con, c(0.3,0.5), c("b[H]==0.3","b[H]==0.5"))

p2 <- y_sum %>% ggplot(aes(x = mean_r, y = 100 * mean_bias / true_beta, colour = parm)) +
  geom_point() + geom_line() +
  #  geom_pointrange(aes(ymin = lcl1 / true_beta, ymax = ucl1 / true_beta)) +
  facet_grid(. ~ b_ac + b_con, scales = "fixed", labeller = label_parsed) + 
  xlab("Number of recaptures") + ylab("% Bias") + coord_cartesian(ylim = c(-100,100)) +
  scale_colour_discrete(breaks = c("D.stdGC", "noneuc.stdGC"), 
                        labels = c(bquote(hat(b)[D] ~ "(Density coefficient)"), 
                                   bquote(hat(b)[H] ~ "(Habitat use coefficient)"))) +
  theme(legend.title = element_blank(), legend.position = "bottom")
p2
# save
ggsave("tost_simulation/output/tost-sim-bias.png", p2, width = 9, height = 3.5, dpi = 300)


y_sum2 <- yf %>% filter(parm %in% c("D.stdGC", "noneuc.stdGC")) %>% droplevels() %>% 
  mutate(true_beta = ifelse(parm == "D.stdGC", b_ac, b_con),
         bias = beta - true_beta,
         perc_bias = 100 * bias / true_beta) 

y_sum2$b_ac <- factor(y_sum2$b_ac, c(0.3,0.5,1,2), c("b[D]==0.3","b[D]==0.5", "b[D]==1", "b[D]==2"))
y_sum2$b_con <- factor(y_sum2$b_con, c(0.3,0.5,2), c("b[H]==0.3","b[H]==0.5", "b[H]==2"))
y_sum2$parm <- factor(y_sum2$parm, c("D.stdGC", "noneuc.stdGC"), 
                      c("hat(b)[D] ~ '(Density coefficient)'",
                     "hat(b)[H] ~ '(Habitat use coefficient)'"))
p3 <- y_sum2 %>% 
  filter(!is.na(SE.beta)) %>%
  ggplot(aes(x = factor(mean_r), y = perc_bias)) +
  geom_boxplot(aes(group = factor(mean_r)), na.rm = TRUE, outlier.shape = NA) +
  facet_grid(parm ~ b_ac + b_con, scales = "fixed", labeller = label_parsed) + 
  xlab("Number of recaptures") + ylab("% Error") + coord_cartesian(ylim = c(-200,200)) +
  theme(legend.title = element_blank(), legend.position = "bottom")
p3
# save
ggsave("tost_simulation/output/tost-sim-bias-box.png", p3, width = 9, height = 6, dpi = 300)


x3 <- y_sum2 %>% filter(!is.na(SE.beta), b_con == "b[H]==0.5", b_ac == "b[D]==2", mean_r <= 100) %>%
  select(bias, perc_bias, everything())


y_sum %>% filter(parm == "noneuc.stdGC") %>% ggplot(aes(x = mean_r, y = mean_beta)) +
  geom_point() +
  geom_pointrange(aes(ymin = lcl2, ymax = ucl2)) +
  facet_grid(b_ac ~ b_con, scales = "free")

