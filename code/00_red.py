# %% Imports
import json
import numpy as np
import pandas as pd

recode_dict = {
    'unknown': np.nan,
    '': np.nan,
    'hospitalized': 1,
    'live': "excluded",
    'released': "excluded",
    'symptomatic': "excluded",
    'deceased': 1,
    'outpatient': 0,
    '-': np.nan,
    'alive': "excluded",
    'recovered': "excluded",
    'asymptomatic': 0,
    'hospitalized (mild)': "excluded",
    'unknow': np.nan,
    'mild': 0,
    'mild clinical signs without hospitalization': 0,
    'home': 0,
    'cured': "excluded",
    'hospitalized (severe)': 1,
    'hospitalized (moderate)': 1,
    'released, live': "excluded",
    'not hospitalized': 0,
    'hospitalized; stable': 1,
    'hospitalized or to be hospitalized': 1,
    'live, released': "excluded",
    'suspected corona': "excluded",
    'n/a': np.nan,
    'naso-pharyngeal swab': "excluded",
    'discharged': "excluded",
    'ehpad': "excluded",
    'physician network': "excluded",
    'in-hospital': 1,
    'acute upper respiratory infection, unspecified': "excluded",
    'quarantine': 1,
    'ward': 1,
    'unkown': np.nan,
    'hospitalized, live': 1,
    '\ufeffunknown': np.nan,
    'overseas inflow': "excluded",
    'suspected coronavirus infection': "excluded",
    'intensive care unit': 1,
    'recovering': 'excluded',
    'hospitalized (critical)': 1,
    'death': 1,
    'still hospitalized': 1,
    'severe / icu': 1,
    'no clinical signs without hospitalization': 0,
    'no clinical signs': 0,
    'contact with and exposure to other communicable diseases': 0,
    'domestic infection': "excluded",
    'icu': 1,
    'moderate / outpatient': 0,
    'mild / contact exposure / asymptomatic': 0,
    'coronavirus infection': "excluded",
    'inpatient': 1,
    'icu; serious': 1,
    'discharged after recovery': "excluded",
    'hospitalized/released': 1,
    'icd-10 disease: j06.9 acute upper respiratory infection, unspecified': "excluded",
    'hospitalized, deceased': 1,
    'pneumonia (chest x-ray)': 1,
    'quarantined': 0,
    'benigne': 0,
    'hospitalised': 1,
    'hospitalized, released': 1,
    'ehpad_ira': 'excluded',
    'icd-10 disease: z03.8 observation for other suspected diseases and conditions': "excluded",
    'physician': 'excluded',
    'encounter for general adult medical examination': 0,
    'hospitalized/deceased': 1,
    'dama': "excluded",
    'facility quarantine': 0,
    'oro-pharyngeal swab': "excluded",
    'fever': "excluded",
    'acute bronchitis': "excluded",
    'encounter for observation for other suspected diseases and conditions ruled out': 0,
    'na': np.nan,
    'hospitalized (intensive care unit)': 1,
    'symptoms indicative of upper respiratory infection': "excluded",
    'initially hospitalized, but now improved and discharged': "excluded",
    'asymptomatic, identified as positive during preoperation investigation': 0,
    'hospitalization': 1,
    'pneumonia, unspecified organism': 1,
    'mild covid-19': 0,
    'live, acute respiratory infection': "excluded",
    'bronchitis': "excluded",
    'isolation': 0,
    'mild, at home.': 0,
    'recovered and released': "excluded",
    'unknwon': np.nan,
    'stable in quarantine': 0,
    'healthy': 0,
    'mild case': 0,
    'screening': 0,
    'rebased': "excluded",
    'asympomatic': 0,
    'other acute upper respiratory infections of multiple sites': "excluded",
    'interned': "excluded",
    'pneumonia (chest x-ray), not critical': 1,
    'hospsitalized, icu, fully recovered': 1,
    'live, mild symptoms, at home': 0,
    'epidemiology study': 0,
    'moderarte covid-19': "excluded",
    'live, physical examination': "excluded",
    'hospitalized, oxygenotherapy, diarrhea': 1,
    'mild symptoms inpatient for observation': "excluded",
    'asymptomatic/released': 0,
    'mild symptoms (fever, cardiovascular disorders)': 0,
    'hospitaized': 1,
    'uncknown': np.nan
}

pd.DataFrame(recode_dict, index=[0]).T.to_csv("data/2020-10-21_recode.csv")

with open('data/2020-10-21_gisaid-patient.json') as f:
    df = pd.DataFrame(json.loads(line) for line in f)

df["covv_patient_status"] = df["covv_patient_status"].str.strip().str.lower()
df["covv_patient_status"].value_counts().to_csv("data/2020-10-21_status-value-counts.csv")

df = df.assign(
    is_red=df["covv_patient_status"].map(recode_dict),
    pid=df["covv_accession_id"].str.extract(
        r"(\d+)",
        expand=False
    ).astype(int)
).set_index("pid")

# json file information (first section of the Results)
(
    df["is_red"].eq("excluded").sum()  # 4200
    + df["is_red"].isna().sum()  # 148121
    + (df["is_red"].ne("excluded") & df["is_red"].notna()).sum()  # 3637
)  # total: 155958

df = df[df["is_red"].ne("excluded") & df["is_red"].notna()]

cols = ["continent", "country", "region", "city"]
df[cols] = df["covv_location"].str.split("/", expand=True).iloc[:, :-1]

pd.Series(df.index).to_csv("data/2020-10-21_trimmed-pids.txt")
df.to_parquet("data/2020-10-21_green-red.parquet")
df.to_csv("data/2020-10-21_green-red.csv")
