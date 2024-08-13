import plotly.express as px
import pandas as pd

# List of countries to highlight in different colors
highlight_countries_1 = ['Ghana', 'Mali', 'Gambia', 'Kenya', 'Tanzania','Vietnam', 'Cambodia', 'Myanmar', 'Laos', 'Thailand']
highlight_countries_2 = ['Cameroon', 'Malawi', 'Sudan', 'Papua New Guinea', 'India', 'Bangladesh']
# Create a DataFrame for highlighted countries with corresponding colors
df = pd.DataFrame({
    'country': highlight_countries_1 + highlight_countries_2,
    'Sample': ['Dxy and dN/dS continent'] * len(highlight_countries_1) + ['dN/dS regional'] * len(highlight_countries_2)
})

# Load the map
fig = px.choropleth(df, locations="country",
                    locationmode='country names',
                    color='Sample',
                    hover_name="country",
                    color_discrete_map={'Dxy and dN/dS continent': "lightgreen" , 'dN/dS regional': "lightblue"})


fig.update_layout(
    geo=dict(showframe=False, showcoastlines=False),
    legend=dict(
        font=dict(
            size=34  # Change this value to make the legend text larger or smaller
        )
    )
)

fig.show()

