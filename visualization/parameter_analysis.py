#!/usr/bin/env python3
"""
Simple parameter space analysis and visualization for Notch EMT model
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

# Set up the plotting style
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# Load the data
data_path = "../Data/parameter_databases/Notch_core_bistability_db.csv"
if not os.path.exists(data_path):
    print(f"Error: Could not find {data_path}")
    exit(1)

df = pd.read_csv(data_path)
print(f"Loaded data with shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")

# Create output directory for figures
output_dir = "../figures/parameter_analysis"
os.makedirs(output_dir, exist_ok=True)

# 1. Distribution of stability types
plt.figure(figsize=(8, 6))
stability_counts = df['stability'].value_counts()
plt.pie(stability_counts.values, labels=stability_counts.index, autopct='%1.1f%%')
plt.title("Distribution of Stability Types")
plt.savefig(f"{output_dir}/stability_distribution.png", dpi=150, bbox_inches='tight')
plt.close()

# 2. Parameter distributions by stability type
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()
params = ['k', 'p', 'pp', 'kk', 'm', 'd']

for i, param in enumerate(params):
    ax = axes[i]
    for stability_type in df['stability'].unique():
        if pd.notna(stability_type):
            data = df[df['stability'] == stability_type][param]
            ax.hist(data, alpha=0.5, label=stability_type, bins=20)
    ax.set_xlabel(param)
    ax.set_ylabel('Count')
    ax.legend()
    ax.set_title(f'Distribution of {param}')

plt.tight_layout()
plt.savefig(f"{output_dir}/parameter_distributions.png", dpi=150, bbox_inches='tight')
plt.close()

# 3. Violin plots for key parameters
fig, axes = plt.subplots(1, 4, figsize=(16, 5))
key_params = ['k', 'p', 'pp', 'kk']

for i, param in enumerate(key_params):
    ax = axes[i]
    # Filter out NaN stability values
    df_clean = df.dropna(subset=['stability'])
    sns.violinplot(data=df_clean, y=param, x='stability', ax=ax)
    ax.set_title(f'{param} by Stability Type')
    ax.set_xlabel('')

plt.tight_layout()
plt.savefig(f"{output_dir}/parameter_violins.png", dpi=150, bbox_inches='tight')
plt.close()

# 4. Correlation heatmap
plt.figure(figsize=(10, 8))
# Select numeric columns only
numeric_cols = ['d', 'm', 'p', 'k', 'pp', 'kk', 'δ', 'α1']
corr_matrix = df[numeric_cols].corr()
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, 
            square=True, linewidths=1)
plt.title('Parameter Correlation Matrix')
plt.tight_layout()
plt.savefig(f"{output_dir}/correlation_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()

# 5. Scatter plot matrix for key parameters
from pandas.plotting import scatter_matrix

key_params_scatter = ['k', 'p', 'pp', 'kk']
# Sample data if too large
if len(df) > 10000:
    df_sample = df.sample(n=10000, random_state=42)
else:
    df_sample = df

scatter_matrix(df_sample[key_params_scatter], alpha=0.2, figsize=(12, 12), diagonal='hist')
plt.suptitle('Scatter Matrix of Key Parameters', y=0.995)
plt.tight_layout()
plt.savefig(f"{output_dir}/scatter_matrix.png", dpi=150, bbox_inches='tight')
plt.close()

print(f"\nGenerated {len(os.listdir(output_dir))} plots in {output_dir}/")
print("Analysis complete!")