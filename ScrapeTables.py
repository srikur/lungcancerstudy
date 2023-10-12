import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
# import selenium exceptions
from selenium.common.exceptions import *
from dataclasses import dataclass
import os

from bs4 import BeautifulSoup

options = Options()
# options.add_argument('--headless')
options.add_argument('--window-size=1920x1080')
driver = webdriver.Chrome(options=options)

unscraped_links = []
df = pd.DataFrame(columns=['NCT Number', 'Study URL', 'Table'])

@dataclass
class Study:
    nct: str
    link: str
    table: str

def remove_text_between_tags(table_str: str):
    """
    Removes all of the text in between a < and > for each tag, except for the tag type
    """

    # Remove all of the text in between a < and > for each tag, except for the tag type
    new_table = ''
    index = 0
    while index < len(table_str):
        if table_str[index] == '<':
            new_table += table_str[index]
            index += 1
            # add characters until the next space is found
            while (table_str[index] != ' ') and (table_str[index] != '>'):
                new_table += table_str[index]
                index += 1
            while table_str[index] != '>':
                index += 1
        else:
            new_table += table_str[index]
        index += 1
    
    # Remove all \n and spaces from string
    new_table = new_table.replace('\n', '')
    new_table = new_table.replace(' ', '')
    new_table = new_table.replace('<!--', '')
    return new_table

ctg_studies = pd.read_excel('lung_cancer_studies_zipcodes_v2.xlsx')
original_studies = pd.read_stata('studies_scraped_tables_original.dta')

# Remove studies that have already been scraped
ctg_studies = ctg_studies[~ctg_studies['NCT Number'].isin(original_studies['nct'])]
print(f"NCTs to scrape: {len(ctg_studies)}")

# Zip together 'NCT Number' and "Study URL"
studies = list(ctg_studies['NCT Number'])
links = [f"https://clinicaltrials.gov/study/{nct}" for nct in studies]

studies_dc = []
for nct, link in zip(studies, links):
    # 1. Use Selenium to open the link
    # 2. Click on the "Results Posted" tab
    # 3. Click on the "Expand all" button with the attribute data-ga-category="Baseline Characteristics"
    # 4. Extract the first instance of a <table> tag that is a child of a <ctg-sticky-container> tag
    try:
        driver.get(link)
        # Wait until the "Results Posted" tab is clickable
        WebDriverWait(driver, 5).until(
            lambda driver: driver.find_element(By.XPATH, "//*[contains(text(), 'Results Posted')]").is_displayed()
        )
        driver.find_element(By.XPATH, "//*[contains(text(), 'Results Posted')]").click()
        # Wait until the "Expand all" button is clickable. Use XPATH to find the button by its data-ga-category attribute
        WebDriverWait(driver, 5).until(
            lambda driver: driver.find_element(By.XPATH, "//button[@data-ga-action='Baseline Characteristics']")
        )
        driver.find_element(By.XPATH, "//button[@data-ga-action='Baseline Characteristics']").click()
        soup = BeautifulSoup(driver.page_source, 'html.parser')
        table = soup.select_one('ctg-baseline-characteristics').select_one('table').prettify()
        if table:
            # Store in a dataframe with columns "NCT Number", "Study URL", and "Table"
            # First store in dataclass
            # table = remove_text_between_tags(table)
            study = Study(nct, link, table)
            studies_dc.append(study)
            # write to tsv
            # with open('studies_scraped_tables.tsv', 'a') as f:
            #     f.write(f'{study.nct}\t{study.link}\t{study.table}\n')
        else:
            print('No table found for', link)
            unscraped_links.append(link)
            with open('unscraped_ncts.txt', 'a') as f:
                f.write(nct + '\n')
    except Exception as e:
        print(e)
        df = pd.DataFrame(studies_dc)
        df.to_stata('studies_scraped_tables_error.dta', write_index=False, version=118)
        driver.close()

# Save to stata .dta file
df = pd.DataFrame(studies_dc)
df.to_stata('studies_scraped_tables.dta', write_index=False, version=118)