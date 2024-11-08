## EDA considering clinical data
source("00_Support_Functions.R")
source("01_Load_Data.R")

cat("Starting EDA\n")
## Gender Distribution
p1 <- ggplot(df_clinical_gbm, aes(x = gender, fill = gender)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("female" = "pink",
                               "male" = "steelblue"),
                    name = "Gender",
                    labels = c("female" = "F", "male" = "M")) +
  labs(title = "Gender Distribution", x = "", y = "Count") +
  theme(axis.text.x = element_blank())

## Race Distribution
p2 <- ggplot(df_clinical_gbm, aes(x = race, fill = race)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("asian" = steel_red,
                               "black or african american" = steel_blue,
                               "not reported" = gray,
                               "white" = pink),
                    name = "Race") +
  labs(title = "Race Distribution", x = "", y = "Value") +
  coord_flip()

## Vital Status
p3 <- ggplot(df_clinical_gbm, aes(x = vital_status, fill = vital_status)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("Alive" = steel_red,
                               "Dead" = steel_blue),
                    name = "Vital Status") +
  labs(title = "Vital Status Distribution", x = "", y = "Value")

## Age Distribution
p4 <- ggplot(df_clinical_gbm, aes(x = gender, y = age_at_index,
                            fill = gender)) +
  geom_boxplot() +
  scale_fill_manual(values = c("female" = pink,
                               "male" = steel_blue),
                    name = "Gender") +
  geom_jitter(size = 1, alpha = 5) +
  labs(title = "Age Distribution", subtitle = "Median Age: 56",
       x = "", y = "Age")

## Pharmaceutical Treatment
p5 <- ggplot(df_clinical_gbm,
       aes(x = treatments_pharmaceutical_treatment_or_therapy,
           fill = treatments_pharmaceutical_treatment_or_therapy)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("no" = steel_blue,
                               "not reported" = gray,
                               "yes" = steel_red),
                    name = "Treatment") +
  labs(title = "Pharmaceutical Treatment Distribution",
       x = "", y = "Value") +
  coord_flip()

## Radiation Treatment
p6 <- ggplot(df_clinical_gbm,
       aes(x = treatments_radiation_treatment_or_therapy,
           fill = treatments_radiation_treatment_or_therapy)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("no" = steel_blue,
                               "not reported" = gray,
                               "yes" = steel_red),
                    name = "Treatment") +
  labs(title = "Radiation Treatment Distribution", x = "", y = "Value") +
  coord_flip()

## Median Age
median_age <- group_by(df_clinical_gbm, gender) %>% 
  summarise(count = n(),
            mean = mean(age_at_index, na.rm = T),
            sd = sd(age_at_index, na.rm = T),
            median = median(age_at_index, na.rm = T),
            IQR = IQR(age_at_index, na.rm = T))

## Print
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(median_age)

cat("EDA Done\n")