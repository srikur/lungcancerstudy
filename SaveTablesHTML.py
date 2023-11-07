import pandas as pd
from jinja2 import Template

# Sample DataFrame
data = pd.read_stata('esophageal_scraped_tables.dta')
ctg_studies = pd.read_excel("Esophageal Cancer Studies 2.xlsx", sheet_name=0, skiprows=1)

ctg_studies['NCT Number'] = ctg_studies['NCT Number'].astype(str)
# make original_studies['nct'] strings
data['nct'] = data['nct'].astype(str)

# list of ctg-studies NCT Number
ncts = ctg_studies['NCT Number'].tolist()

# filter data to only include ncts in ncts list
data = data[data['nct'].isin(ncts)]

# sort by nct, ascending
data = data.sort_values(by=['nct'])

# Load the HTML template
with open('template.html', 'r') as file:
    template_content = file.read()
template = Template(template_content)

# Create a list of tables as (identifier, HTML_table) pairs
tables = [(row['nct'], row['table']) for _, row in data.iterrows()]

# Render the HTML
html_output = template.render(tables=tables)

# Save the HTML to a file
with open('esophageal_tables.html', 'w') as output_file:
    output_file.write(html_output)