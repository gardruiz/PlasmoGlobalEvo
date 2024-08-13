import seaborn as sns
import matplotlib.pyplot as plt
import sys
import pandas as pd

data1= pd.read_csv(sys.argv[1])
data2= pd.read_csv(sys.argv[2])
data1['country']=sys.argv[3]
data2['country']=sys.argv[4]
data = pd.concat([data1, data2])
# Assuming `data` is your DataFrame
print(data.head())
print(data.columns)


# Filter the data for one country, e.g., 'Nigeria'
country = 'Nigeria'
data_country = data[data['country'] == country]


# Create the count plot
plt.figure(figsize=(14, 8))
sns.histplot(data=data2, x='snp_distance', bins=100, hue='country', palette='Set2')



plt.xlabel('SNP Distance Bins')
plt.ylabel('Frequency')
plt.title('Frequency of SNP Distance Values by Country')
plt.legend(title='Country')
plt.tight_layout()
plt.show()

