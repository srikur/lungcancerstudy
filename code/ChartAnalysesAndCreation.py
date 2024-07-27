import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from plotnine import annotate, ggplot, aes, geom_line, labs, scale_color_manual, theme_minimal, theme, geom_hline, geom_smooth
from plotnine import scale_x_continuous, scale_y_continuous, stat_smooth
from plotnine import ggplot, aes, geom_bar, facet_wrap, labs, theme, element_text
pd.options.mode.chained_assignment = None
base_folder = "/Users/srikur/Documents/GitHub/lungcancerstudy"

# ----------------- Data Preparation -----------------
df = pd.read_excel(f"{base_folder}/data/lung_cancer_studies_final.xlsx", sheet_name="Sheet1")
df = df[(df["Start Year"] >= 2002) & (df["Start Year"] <= 2021)]
df = df[df["NCT Number"] != "NCT02981108"] # Outlier
df.dropna(subset=["NCT Number"], inplace=True)
df["Other"] = df["Mixed"] + df["Unknown Race"]
df["minority_number"] = df["Total Number"] - df["White"]

# For rows where "White" is not NaN, assign "has_race" to 1, else 0
df["has_race"] = df["White"].notna().astype(int)

df["has_ethnicity"] = df["Hispanic"].notna().astype(int)
df.loc[df["Funder Type"].isin(["NETWORK", "FED"]), "Funder Type"] = "OTHER"
df["Start Year"] = df["Start Year"].astype(int)
lung_cancer_prevalence = pd.read_csv(f"{base_folder}/data/lung_cancer_prevalence.csv")

# ----------------- Prevalence-to-Participation Ratio Denominator Calculations -----------------
total_males = lung_cancer_prevalence[lung_cancer_prevalence["Group"] == "Total"]["Males"].values[0]
total_females = lung_cancer_prevalence[lung_cancer_prevalence["Group"] == "Total"]["Females"].values[0]
total_total = lung_cancer_prevalence[lung_cancer_prevalence["Group"] == "Total"]["Total"].values[0]

p2p_ratio_female_ethnicities = lung_cancer_prevalence["Females"] / total_females * lung_cancer_prevalence["PopPercentage"]
p2p_ratio_male_ethnicities = lung_cancer_prevalence["Males"] / total_males * lung_cancer_prevalence["PopPercentage"]
p2p_ratio_total_ethnicities = lung_cancer_prevalence["Total"] / total_total * lung_cancer_prevalence["PopPercentage"]

# concat the two series into a dataframe with Ethnicity as a column
p2p_ratios = pd.DataFrame({
    'Group': lung_cancer_prevalence['Group'],
    'Female_Ratio': p2p_ratio_female_ethnicities,
    'Male_Ratio': p2p_ratio_male_ethnicities,
    'Total_Ratio': p2p_ratio_total_ethnicities
})

# Calculate overall ppr point statistic for Female, Minority, Asian, Black, White, Native American, Hispanic
group_columns = ['Female', 'minority_number', 'Asian', 'Black', 'White', 'Native American', 'Hispanic']
overall_df = df[df["White"].notna()]
# first calculate participation rate for each group in total across all years
participation_rate = overall_df[group_columns].sum() / overall_df["Total Number"].sum()
# then calculate the ppr for each group individually
pprs = {}
for column in group_columns:
    if column == "minority_number":
        pprs[column] = participation_rate[column] / p2p_ratios["Total_Ratio"][0:4].sum()
    else:
        pprs[column] = participation_rate[column] / p2p_ratios[p2p_ratios["Group"] == column]["Total_Ratio"].values[0]

gender_columns = ["Female", "Male"]
gender_participation_agg = df.groupby("Start Year").agg({k: "sum" for k in gender_columns}).reset_index()
gender_participation_agg = gender_participation_agg.astype(int)
gender_participation_agg["total_participants"] = gender_participation_agg.apply("sum", axis=1) - gender_participation_agg["Start Year"]
gender_participation_agg["female_participation_rate"] = gender_participation_agg["Female"] / gender_participation_agg["total_participants"]
gender_participation_agg["female_ppr"] = gender_participation_agg["female_participation_rate"] / p2p_ratios[p2p_ratios["Group"] == "Female"]["Total_Ratio"].values[0]
gender_participation_agg = gender_participation_agg[(gender_participation_agg["Start Year"] >= 2002) & (gender_participation_agg["Start Year"] <= 2021)]

hispanic_agg = df.groupby("Start Year").agg({"Hispanic": "sum", "Non-His": "sum", "Unknown Ethnicity": "sum"}).reset_index()
hispanic_agg = hispanic_agg.astype(int)
hispanic_agg["total_participants"] = hispanic_agg.apply("sum", axis=1) - hispanic_agg["Start Year"]
hispanic_agg["hispanic_participation_rate"] = hispanic_agg["Hispanic"] / hispanic_agg["total_participants"]
hispanic_agg["hispanic_ppr"] = hispanic_agg["hispanic_participation_rate"] / p2p_ratios[p2p_ratios["Group"] == "Hispanic"]["Total_Ratio"].values[0]
hispanic_agg = hispanic_agg[(hispanic_agg["Start Year"] >= 2002) & (hispanic_agg["Start Year"] <= 2021)]

# ----------------- Combined Female, Minority, Hispanic PPR by year -----------------
combined_ppr = pd.merge(gender_participation_agg, minority_participation_agg, on="Start Year")
combined_ppr = pd.merge(combined_ppr, hispanic_agg, on="Start Year")
combined_ppr = combined_ppr[['Start Year', 'minority_ppr', 'hispanic_ppr', 'female_ppr']]
combined_ppr.columns = ['Start Year', 'Minority', 'Hispanic', 'Female']

# Linear regression
x = combined_ppr["Start Year"]
y = combined_ppr["Minority"]
minority_pval = stats.linregress(x,y)[3]

y = combined_ppr["Hispanic"]
hispanic_pval = stats.linregress(x,y)[3]

y = combined_ppr["Female"]
female_pval = stats.linregress(x,y)[3]

# Plotting
base_plot = (
    ggplot(combined_ppr, aes(x="Start Year")) +
    geom_line(aes(y="Minority", color="'Minority'"), size=1, na_rm=True) +
    geom_line(aes(y="Hispanic", color="'Hispanic'"), size=1, na_rm=True) +
    geom_line(aes(y="Female", color="'Female'"), size=1, na_rm=True) + 
    geom_smooth(aes(y="Minority", color="'Minority'"), method='lm', linetype='dashed', se=False) +
    geom_smooth(aes(y="Hispanic", color="'Hispanic'"), method='lm', linetype='dashed', se=False) +
    geom_smooth(aes(y="Female", color="'Female'"), method='lm', linetype='dashed', se=False) +
    labs(x="Year", y="PPR", title="Female, Hispanic, and Minority PPR by Year") + 
    scale_color_manual(values={"Minority": "mediumturquoise", "Hispanic": "tomato", "Female": "purple"}) +
    theme(legend_title=element_text(text="Group"))
)

# Add p-values to the plot
base_plot += annotate('text', x=2010, y=0.21, label=f"Female p-value: {female_pval.round(3)}", color="purple", ha='left') # type: ignore
base_plot += annotate('text', x=2010, y=0.13, label=f"Hispanic p-value: {minority_pval.round(3)}", color="mediumturquoise", ha='left') # type: ignore
base_plot += annotate('text', x=2010, y=0.05, label=f"Ethnicity p-value: {hispanic_pval.round(3)}", color="tomato", ha='left') # type: ignore
base_plot += geom_hline(yintercept=1, color='green', linetype='dashed')

base_plot.show()

base_plot.save(f"{base_folder}/plots/combined_hispanic_female_minority.png", dpi=1200)
base_plot.save(f"{base_folder}/plots/combined_hispanic_female_minority.svg", dpi=1200)

# ----------------- Surgery Type and Demographic Chart -----------------
minority_df = df[df["White"].notna()]
minority_surgery_participation = minority_df.groupby("Surgery?").agg({"minority_number": "sum", "Total Number": "sum"}).reset_index()
minority_surgery_participation["minority_surgery_participation_rate"] = minority_surgery_participation["minority_number"] / minority_surgery_participation["Total Number"]
minority_surgery_participation["minority_surgery_ppr"] = minority_surgery_participation["minority_surgery_participation_rate"] / p2p_ratios["Total_Ratio"][0:4].sum()

women_surgery_participation = df.groupby("Surgery?").agg({"Female": "sum", "Total Number": "sum"}).reset_index()
women_surgery_participation["female_participation_rate"] = women_surgery_participation["Female"] / women_surgery_participation["Total Number"]
women_surgery_participation["female_ppr"] = women_surgery_participation["female_participation_rate"] / p2p_ratios[p2p_ratios["Group"] == "Female"]["Total_Ratio"].values[0]

df["female_participation_rate"] = df["Female"] / df["Total Number"]
women_surgery_participation_yes = df[df["Surgery?"] == 1.0]
# fillna(0) because there are some NaN values in the minority participation rate column
women_surgery_participation_yes.fillna({"female_participation_rate": 0}, inplace=True)
women_surgery_participation_no = df[df["Surgery?"] == 0.0]
women_surgery_participation_no.fillna({"female_participation_rate": 0}, inplace=True)

hispanic_df = df[df["Hispanic"].notna()]
hispanic_surgery_participation = hispanic_df.groupby("Surgery?").agg({"Hispanic": "sum", "Total Number": "sum"}).reset_index()
hispanic_surgery_participation["hispanic_participation_rate"] = hispanic_surgery_participation["Hispanic"] / hispanic_surgery_participation["Total Number"]
hispanic_surgery_participation["hispanic_ppr"] = hispanic_surgery_participation["hispanic_participation_rate"] / p2p_ratios[p2p_ratios["Group"] == "Hispanic"]["Total_Ratio"].values[0]

dfs_to_merge = [minority_surgery_participation, women_surgery_participation, hispanic_surgery_participation]
merged_df = dfs_to_merge[0]
for surgery_df in dfs_to_merge[1:]:
    merged_df = pd.merge(merged_df, surgery_df, on="Surgery?")

merged_df = merged_df[["Surgery?", "minority_surgery_participation_rate", 
                        "female_participation_rate", "female_ppr", 
                        "hispanic_participation_rate", "hispanic_ppr", 
                        "minority_surgery_ppr"]]

# Reshape the data
melted_df = merged_df.melt(id_vars='Surgery?', value_vars=['minority_surgery_ppr', 'female_ppr', 'hispanic_ppr'])
melted_df['Demographic'] = melted_df['variable'].map({'minority_surgery_ppr': 'Minority', 'female_ppr': 'Female', 'hispanic_ppr': 'Hispanic'})
melted_df['Surgery?'] = melted_df['Surgery?'].map({0: 'No', 1: 'Yes'})

# Create grouped bar chart on a single axis
base_plot = (
    ggplot(melted_df, aes(x='Surgery?', y='value', fill='Demographic')) +
    geom_bar(stat='identity', position='dodge') +
    labs(x='Surgery Offered?', y='Prevalence-to-Participation Ratio (PPR)', title='PPR by Surgery Offered and Demographic') +
    theme(axis_text_x=element_text(rotation=0))
)

base_plot.show()

base_plot.save(f"{base_folder}/plots/ppr_by_surgery_and_demographic.png", dpi=1200)
base_plot.save(f"{base_folder}/plots/ppr_by_surgery_and_demographic.svg", dpi=1200)

# ----------------- Funding Type and Demographic Chart -----------------
minority_df = df[df["White"].notna()]
minority_df["minority_participation_rate"] = minority_df["minority_number"] / minority_df["Total Number"]
minority_funder_participation = minority_df.groupby("Funder Type").agg({"minority_number": "sum", "Total Number": "sum"}).reset_index()
minority_funder_participation["minority_participation_rate"] = minority_funder_participation["minority_number"] / minority_funder_participation["Total Number"]
minority_funder_participation["minority_ppr"] = minority_funder_participation["minority_participation_rate"] / p2p_ratios["Total_Ratio"][0:4].sum()

women_funder_participation = df.groupby("Funder Type").agg({"Female": "sum", "Total Number": "sum"}).reset_index()
women_funder_participation["female_participation_rate"] = women_funder_participation["Female"] / women_funder_participation["Total Number"]
women_funder_participation["female_ppr"] = women_funder_participation["female_participation_rate"] / p2p_ratios[p2p_ratios["Group"] == "Female"]["Total_Ratio"].values[0]

hispanic_df = df[df["Hispanic"].notna()]
hispanic_funder_participation = hispanic_df.groupby("Funder Type").agg({"Hispanic": "sum", "Total Number": "sum"}).reset_index()
hispanic_funder_participation["hispanic_participation_rate"] = hispanic_funder_participation["Hispanic"] / hispanic_funder_participation["Total Number"]
hispanic_funder_participation["hispanic_ppr"] = hispanic_funder_participation["hispanic_participation_rate"] / p2p_ratios[p2p_ratios["Group"] == "Hispanic"]["Total_Ratio"].values[0]

dfs_to_merge = [minority_funder_participation, women_funder_participation, hispanic_funder_participation]
merged_df = dfs_to_merge[0]
for funder_df in dfs_to_merge[1:]:
    merged_df = pd.merge(merged_df, funder_df, on="Funder Type")
merged_df = merged_df[["Funder Type", "minority_ppr", "female_ppr", "hispanic_ppr"]]

# Bar Chart
melted_df = merged_df.melt(id_vars='Funder Type', value_vars = ["minority_ppr", "female_ppr", "hispanic_ppr"])
melted_df['Demographic'] = melted_df['variable'].map({'minority_ppr': 'Minority', 'female_ppr': 'Female', 'hispanic_ppr': 'Hispanic'})
melted_df['Funder Type'] = melted_df['Funder Type'].map({'OTHER': 'Other', 'INDUSTRY': 'Industry', 'NIH': 'NIH', 'OTHER': 'Other'})

base_plot = (
    ggplot(melted_df, aes(x='Demographic', y='value', fill='Funder Type')) +
    geom_bar(stat='identity', position='dodge') +
    labs(x='Funder Type', y='Prevalence-to-Participation Ratio (PPR)', title='PPR by Funder Type and Demographic') +
    theme(axis_text_x=element_text(rotation=0))
)

base_plot.show()

base_plot.save(f"{base_folder}/plots/ppr_by_funder_and_demographic.png", dpi=1200)
base_plot.save(f"{base_folder}/plots/ppr_by_funder_and_demographic.svg", dpi=1200)

# ----------------- Race and Ethnicity Reporting Rate -----------------
reporting_df = df.groupby("Start Year").agg({"has_ethnicity": ["mean"], "has_race": ["mean"]})
reporting_df.columns = ["ethnicity_reporting_rate", "race_reporting_rate"]
reporting_df = reporting_df[reporting_df.index >= 2002].reset_index(drop=False)

# Linear regression for both
x = reporting_df["Start Year"]
y = reporting_df["ethnicity_reporting_rate"]
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
eth_pval = p_value

x = reporting_df["Start Year"]
y = reporting_df["race_reporting_rate"]
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
race_pval = p_value

# ggplot of race and ethnicity reporting rate over time
base_plot = (
    ggplot(reporting_df, aes(x="Start Year")) +
    geom_line(aes(y="ethnicity_reporting_rate", color='"Ethnicity"'), size=1, na_rm=True) +
    geom_line(aes(y="race_reporting_rate", color='"Race"'), size=1, na_rm=True) +
    geom_smooth(aes(y="ethnicity_reporting_rate", color='"Ethnicity"'), method='lm', linetype='dashed', se=False) +
    geom_smooth(aes(y="race_reporting_rate", color='"Race"'), method='lm', linetype='dashed', se=False) +
    labs(x="Year", y="Reporting Rate", title="Race and Ethnicity Reporting Rate by Year") +
    scale_color_manual(values={"Ethnicity": "mediumturquoise", "Race": "tomato"}) +
    theme(legend_title=element_text(text="Reporting Type"))
)

# Add p-values to the plot
base_plot += annotate('text', x=2012.5, y=0.05, label=f"Ethnicity p-value: {'{:.3e}'.format(eth_pval) if eth_pval >= 0.001 else '< 0.001'}", color="mediumturquoise", ha='left') # type: ignore
base_plot += annotate('text', x=2012.5, y=0.10, label=f"Race p-value: {'{:.3e}'.format(race_pval) if race_pval >= 0.001 else '< 0.001'}", color="tomato", ha='left') # type: ignore

base_plot.show()

base_plot.save(f"{base_folder}/plots/race_and_ethnicity_reporting_rate_over_time.png", dpi=1200)
base_plot.save(f"{base_folder}/plots/race_and_ethnicity_reporting_rate_over_time.svg", dpi=1200)