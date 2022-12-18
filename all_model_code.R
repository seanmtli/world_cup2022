library(tidyverse)
library(glmnet)
library(extraDistr)
source("biv-pois-functions/lm.bp.R")
source("biv-pois-functions/bivpois.table.R")
source("biv-pois-functions/newnamesbeta.R")
source("biv-pois-functions/pbivpois.R")
source("biv-pois-functions/simple.bp.R")
source("biv-pois-functions/splitbeta.R")
int_matches <- read_csv("data/international_matches.csv")
qatar_teams <- read_delim("data/Qatar2022-teams.csv", delim = ";")
matches_sched <- read_delim("data/matches_schedule.csv", delim = ";")
final_roster_data <- read_csv("data/final_roster_data.csv")
worldcup_2022_matches<-read_csv("data/Fifa_world_cup_matches.csv")
qatar_groupstage_teams <- read.csv("data/FIFA2022_GroupStageTeams.csv")
matches_results_2022 <- read.csv("data/matches_results.csv")
cleaned_matches <- read_csv("data/cleaned_matches.csv")
results_and_odds2 <-read_csv("data/results2.csv")

cups_only <- int_matches %>%
  filter(!str_detect(tournament, "Friendly|qualification"), format(date, "%Y") >= 2000) %>%
  select(-home_team_goalkeeper_score, away_team_goalkeeper_score) %>%
  mutate(rank_differential = away_team_fifa_rank - home_team_fifa_rank, #away - home so positive rank differential means home team is paper
         offense_differential = home_team_mean_offense_score - away_team_mean_offense_score,
         defense_differential = home_team_mean_defense_score - away_team_mean_defense_score,
         midfield_differential = home_team_mean_midfield_score - away_team_mean_midfield_score,
         points_differential = home_team_total_fifa_points - away_team_total_fifa_points,
         year = as.numeric(format(date, "%Y"))) 

cups_only_nona <- na.omit(cups_only)

cups_only_nona <- cups_only_nona%>%
  mutate(goals_differential = home_team_score - away_team_score,
         avg_team_differential = (offense_differential + defense_differential + midfield_differential)/3 )


biv_cups_only_nona <- cups_only_nona %>%
  mutate(neutral_location = ifelse(neutral_location == TRUE, 1, 0))

# HOME FIELD COMPENSATION?
n <- nrow(biv_cups_only_nona)


train_data_biv <- biv_cups_only_nona

#all variables are predictors that have significance based on lasso

#each lambda in fin_mod represents the rate param for home, away, and the 
#correlation between the first two vars

#lambda depends on the beta values (the X coefficient)s 





fin_mod <- lm.bp(l1 = home_team_score~ 
                   offense_differential + defense_differential + midfield_differential,
                 l2 = away_team_score~ offense_differential + defense_differential + midfield_differential,
                 l3 = ~ offense_differential + defense_differential + midfield_differential, 
                 data = train_data_biv)

#qatar_teams
matches_sched[matches_sched == "USA"] <- "United States"


group_matches <- matches_sched %>% 
  filter(phase == "group matches") %>% 
  rename(hometeam = country1, awayteam = country2)

group_matches <- group_matches %>% full_join(qatar_groupstage_teams, by = c('hometeam'='Team')) %>% 
  full_join(qatar_groupstage_teams, by = c('awayteam'='Team')) %>% 
  mutate(
    offense_differential = Offense.x - Offense.y,
    defense_differential = Defense.x - Defense.y,
    midfield_differential = Midfield.x - Midfield.y,
    fifa_rank_differential = Rank.x-Rank.y,
    points_differential = Points.x - Points.y,
    goalkeeper_differential=Goalkeeper.x - Goalkeeper.y,
    possession_differential=avg_possession.x-avg_possession.y,
    shots_on_target_differential=avg_shotsontarget.x-avg_shotsontarget.y,
    shots_differential=avg_shots.x-avg_shots.y,
    yellow_cards_differential=avg_yellowcards.x-avg_yellowcards.y,
    red_card_differential=avg_redcards.x-avg_redcards.y,
    fouls_differential=avg_fouls.x-avg_fouls.y,
    score_differential=avg_goals.x-avg_goals.y) %>% 
  select(hometeam, awayteam, offense_differential, defense_differential, midfield_differential,
         fifa_rank_differential,points_differential,goalkeeper_differential,possession_differential,
         shots_on_target_differential,shots_differential,yellow_cards_differential,red_card_differential,
         fouls_differential,score_differential)

final_model2_data <- matches_results_2022 %>%
  left_join(group_matches, by = c("country1" = "hometeam", "country2" = "awayteam")) %>%
  mutate(neutral_location = ifelse(country1 == "Qatar" | country2 == "Qatar", FALSE, TRUE)) %>%
  .[1:48, ]

final_model1_data <- final_model2_data %>%
  select(date,country1, country2, country1_score, country2_score, 
         offense_differential, defense_differential, midfield_differential) 


test_home <- final_model2_data[, "country1_score"] 
test_away <- final_model2_data[, "country2_score"] 

home_teams <- final_model2_data[, "country1"] 
away_teams <- final_model2_data[, "country2"] 



td1 <- final_model1_data %>%
  select(offense_differential, 
         defense_differential,midfield_differential) %>% 
  mutate(int_col = rep(1,nrow(.))) %>% 
  relocate(int_col)
td2 <- final_model1_data %>%
  select(offense_differential, 
         defense_differential,midfield_differential) %>% 
  mutate(int_col = rep(1,nrow(.))) %>% 
  relocate(int_col)
td3 <- final_model1_data %>%
  select(offense_differential, 
         defense_differential,midfield_differential) %>% 
  mutate(int_col = rep(1,nrow(.))) %>% 
  relocate(int_col)

#get final values for lambdas
final_lambdas <- tibble(lambda_1 = exp(as.matrix(td1) %*% as.matrix(fin_mod$beta1)),
                        lambda_2 = exp(as.matrix(td2) %*% as.matrix(fin_mod$beta2)),
                        lambda_3 = exp(as.matrix(td3) %*% as.matrix(fin_mod$beta3)))

n_test_obs <- nrow(final_lambdas)
n_sims <- 100000
set.seed(1)
all_monte_carlos <- as.data.frame(matrix(NA, nrow = n_sims, ncol = 2*n_test_obs))
#draw randomly from the distirbution of the poisson distributions 
for(i in 1:n_test_obs) {
  # [home-away] score prediction so ex: 2-0 or 3-3
  monte_carlo <- rbvpois(n_sims, final_lambdas$lambda_1[i], 
                         final_lambdas$lambda_2[i], 
                         final_lambdas$lambda_3[i])
  all_monte_carlos[, (2*i-1)] <- monte_carlo[, 1]
  all_monte_carlos[, (2*i)] <- monte_carlo[, 2]
}

test_actuals <- tibble(match_num = seq(1, length(test_home), 1), 
                       home_team = home_teams,
                       away_team = away_teams,
                       home = test_home, away = test_away) %>%
  mutate(outcome = case_when(
    test_home - test_away == 0 ~ "Draw",
    test_home - test_away > 0 ~ "Win",
    test_home - test_away < 0 ~ "Loss"
  ),
  pred_win_prob = NA,
  pred_loss_prob = NA,
  pred_draw_prob = NA,
  )

for(i in 1:n_test_obs) {
  test_actuals[i, "pred_win_prob"] <- 
    (sum(all_monte_carlos[, (2*i - 1)] - all_monte_carlos[, (2*i)] > 0))/n_sims
  test_actuals[i, "pred_loss_prob"] <- 
    (sum(all_monte_carlos[, (2*i - 1)] - all_monte_carlos[, (2*i)] < 0)) /n_sims
  test_actuals[i, "pred_draw_prob"] <- 
    (sum(all_monte_carlos[, (2*i - 1)] - all_monte_carlos[, (2*i)] == 0)) /n_sims  
}

goal_diffs <- tibble(mean_diffs = rep(NA, n_test_obs))
for(i in 1:n_test_obs) {
  goal_diffs[i, "mean_diffs"] <- 
    mean(all_monte_carlos[, (2*i - 1)] - all_monte_carlos[, (2*i)])
}

test_actuals <- test_actuals %>%
  mutate(pred_outcome = case_when(
    pmax(pred_win_prob, pred_loss_prob, pred_draw_prob) == pred_draw_prob ~ "Likely Draw",
    pmax(pred_win_prob, pred_loss_prob, pred_draw_prob) == pred_win_prob ~ "Likely Win",
    pmax(pred_win_prob, pred_loss_prob, pred_draw_prob) == pred_loss_prob ~ "Likely Loss"
    
  ))

# Confusion matrix assuming max_prob is the predicted outcome - bit misleading as performance metric because sometimes probabilities are only diff by 0.05 or so but this will still pick a single outcome. Should find a way to adjust for close wins/losses/draws
confusion_mat <- table(test_actuals$outcome, test_actuals$pred_outcome)

confusion_mat

# Draw accuracy still zero bruh

# win accuracy:
confusion_mat[6]/(confusion_mat[3] + confusion_mat[6])

# Loss accuracy
confusion_mat[2]/(confusion_mat[2] + confusion_mat[5])

# Overall accuracy
(confusion_mat[2] + confusion_mat[6]) /n_test_obs

books_accuracy <- results_and_odds2 %>%
  mutate(pred_outcome = case_when(
    pmax(winProb, loseProb, drawProb) == drawProb ~ "Likely Draw",
    pmax(winProb, loseProb, drawProb) == winProb ~ "Likely Win",
    pmax(winProb, loseProb, drawProb) == loseProb ~ "Likely Loss"
    
  ))

# Confusion matrix assuming max_prob is the predicted outcome - bit misleading as performance metric because sometimes probabilities are only diff by 0.05 or so but this will still pick a single outcome. Should find a way to adjust for close wins/losses/draws
confusion_mat_books <- table(books_accuracy$outcome, books_accuracy$pred_outcome)

confusion_mat_books

# win accuracy:
confusion_mat_books[6]/(confusion_mat_books[3] + confusion_mat_books[6])

# Loss accuracy
confusion_mat_books[2]/(confusion_mat_books[2] + confusion_mat_books[5])

# Overall accuracy
(confusion_mat_books[2] + confusion_mat_books[6]) /n_test_obs