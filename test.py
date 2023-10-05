import table1
import pandas as pd

# read data
df = pd.read_csv('https://api.vitaldb.net/cases')

# add columns
df['opdur'] = df['opend'] - df['opstart']
df['anedur'] = df['aneend'] - df['anestart']
df['hospdur'] = df['dis'] - df['adm']

# remove columns
df.drop(columns=['opstart', 'opend', 'anestart', 'aneend', 'dis', 'adm', 'caseid'], inplace=True)
df = df.loc[:, ~df.columns.str.endswith('id')]

# create table one
df = table1.table1(df, 'department')  # 'death_inhosp'
df.to_csv('table1.csv', index=False, header=False)
print(df)
