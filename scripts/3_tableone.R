### tableone ###

if (!require("tableone")) install.packages("tableone")
library(tableone)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

date <- '20251011'
time_window_width <- 6
data_dir <- "./data/"

load(paste0(data_dir, "df_", date , "_", time_window_width, "hr.RData"))
df <- df %>% filter(time_window_index == 0)

myVars <- c('age', 'female','hospital_id',
            'infection','trauma','circulatory','digestive','neoplasms',
            'neurological','external_causes','respiratory','others',
            'acute_coronary_syndrome','septic_shock', 
            'bt', 'hr', 'rr', 'mbp', 'spo2', 'hemoglobin', 'ph', 'lactate',
            'label_mechanical_ventilation','noradrenaline','adrenaline', 'dopamine','vasopressin')

catVars <- c('female','hospital_id',
             'infection','trauma','circulatory','digestive','neoplasms',
             'neurological','external_causes','respiratory','others',
             'acute_coronary_syndrome','septic_shock',
             'label_mechanical_ventilation','noradrenaline','adrenaline', 'dopamine','vasopressin')

tab1 <- CreateTableOne(vars = myVars,
                       data = df, 
                       factorVars = catVars) 
tab1

