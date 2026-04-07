import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# New Colormap
newcolors = ['#FF0500','#FF1E00','#FF3800','#FF5100','#FF6B00','#FF8400','#FF9E00','#FFB700','#FFD100','#FFEA00','#F4F40A','#C1C13D','#8E8E70','#5B5BA3','#2828D6','#0000F9','#0000E0','#0000C6','#0000AD','#000093']
newcolors = newcolors[::-1]
cmap = ListedColormap(newcolors)

# Load CSV file
df = pd.read_csv("output.csv")

print("NuType ", df.NuType.min(), df.NuType.max())
print("InitialFlavour ", df.InitialFlavour.min(), df.InitialFlavour.max())
print("FinalFlavour ", df.FinalFlavour.min(), df.FinalFlavour.max())
print("Energy ", df.Energy.min(), df.Energy.max())
print("CosineZ ", df.CosineZ.min(), df.CosineZ.max())
print("Probability ", df.Probability.min(), df.Probability.max())

# Replace with your actual column names
x_col = "CosineZ"
y_col = "Energy"
z_col = "Probability"

# Choose subset
df_sub=df[(df.NuType==1) & (df.InitialFlavour==2) & (df.FinalFlavour==2)]
print("InitialFlavour ", df_sub.InitialFlavour.min(), df_sub.InitialFlavour.max())

# Create scatter plot
plt.figure()
scatter = plt.scatter(
    df_sub[y_col],
    df_sub[x_col],
    c=df_sub[z_col],
    cmap=cmap # "seismic"  # you can change this (e.g., plasma, inferno, coolwarm)
)

# Add colorbar
plt.colorbar(scatter, label=z_col)
plt.xscale('log')

# Labels
plt.xlabel("Energy (GeV)")
plt.ylabel("Cosine Zenith Angle")
# plt.title(f"{y_col} over {x_col} (colored by {z_col})")

# Show plot
plt.savefig('NSI_plot_mu2mu.png')
plt.show()
