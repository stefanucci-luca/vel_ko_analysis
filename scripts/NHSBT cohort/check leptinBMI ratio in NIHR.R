
# NIHR cohort
wdf

```{r }
# Add new ratio
wdf %>% 
  mutate(lep_over_bmi = LEPT / BMI) -> wdf
```

```{r }
# cehck distribution new vraiable
plot(density(wdf$lep_over_bmi, na.rm = T))
```

```{r }
# correct distribution
wdf %>% 
  # correct for skewness
  mutate(lep_over_bmi = log(lep_over_bmi) ) %>% 
  # standardise
  mutate(lep_over_bmi = lep_over_bmi - mean(lep_over_bmi, na.rm =T) / sd(lep_over_bmi, na.rm = T)) -> wdf
```

```{r }
model = lm(lep_over_bmi ~ COHORT + GENDER + AGE, data = wdf)
parameters::parameters(model)
```