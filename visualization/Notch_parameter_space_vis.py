import numpy as np
import matplotlib.pyplot as plt
import altair as alt
import pandas as pd
import seaborn as sns
df = pd.read_csv("../Data/parameter_databases/Notch_parameter_bistability_r|nr_statistics.csv")
df = pd.read_csv("../Data/parameter_databases/Notch_core_bistability_db.csv")


df_test = df.loc[(df['pp'] ==100) & (df['p'] >=90),:]

df_test.boxplot(by =['stability', 'pp'])











sns.set_theme(style="ticks")
# Initialize the figure with a logarithmic x axis
f, ax = plt.subplots(figsize=(12, 8))
sns.boxplot(x="k", y="pp", data=df, whis=[0, 20], width=.6, palette="vlag")
sns.stripplot(x="k", y="p", data=df, size=4, color=".3", linewidth=0)

sns.violinplot(
    df.dropna(),
    x="k", hue="stability",
    multiple="stack",
    palette="light:m_r",
    edgecolor=".3",
    linewidth=.5)


sns.displot(
    data=df.dropna(), x="k", y="kk", col="type",
    col_wrap=4, height=4, aspect=.7,
)


sns.swarmplot(x="p", y="pp", hue="stability", data=df)




alt.data_transformers.disable_max_rows()
alt.Chart(df.dropna()).transform_density(
    'k',
    as_=['k', 'density'],
    extent=[0, 20],
    groupby=['type']
).mark_area(orient='horizontal').encode(
    y='k:Q',
    color='type:N',
    x=alt.X(
        'density:Q',
        stack='center',
        impute=None,
        title=None,
        axis=alt.Axis(labels=False, values=[0], grid=False, ticks=True),
    ),
    column=alt.Column(
        'type:N',
        header=alt.Header(
            titleOrient='bottom',
            labelOrient='bottom',
            labelPadding=0,
        ),
    )
).properties(
    width=100
).configure_facet(
    spacing=0
).configure_view(
    stroke=None
)


alt.Chart(df.dropna()).transform_density(
    'p',
    as_=['p', 'density'],
    extent=[0, 20],
    groupby=['type']
).mark_area(orient='horizontal').encode(
    y='p:Q',
    color='type:N',
    x=alt.X(
        'density:Q',
        stack='center',
        impute=None,
        title=None,
        axis=alt.Axis(labels=False, values=[0], grid=False, ticks=True),
    ),
    column=alt.Column(
        'type:N',
        header=alt.Header(
            titleOrient='bottom',
            labelOrient='bottom',
            labelPadding=0,
        ),
    )
).properties(
    width=100
).configure_facet(
    spacing=0
).configure_view(
    stroke=None
)


alt.Chart(df.dropna()).transform_density(
    'pp',
    as_=['pp', 'density'],
    extent=[0, 20],
    groupby=['type']
).mark_area(orient='horizontal').encode(
    y='pp:Q',
    color='type:N',
    x=alt.X(
        'density:Q',
        stack='center',
        impute=None,
        title=None,
        axis=alt.Axis(labels=False, values=[0], grid=False, ticks=True),
    ),
    column=alt.Column(
        'type:N',
        header=alt.Header(
            titleOrient='bottom',
            labelOrient='bottom',
            labelPadding=0,
        ),
    )
).properties(
    width=100
).configure_facet(
    spacing=0
).configure_view(
    stroke=None
)

alt.Chart(df.dropna()).transform_density(
    'kk',
    as_=['kk', 'density'],
    extent=[0, 20],
    groupby=['type']
).mark_area(orient='horizontal').encode(
    y='kk:Q',
    color='type:N',
    x=alt.X(
        'density:Q',
        stack='center',
        impute=None,
        title=None,
        axis=alt.Axis(labels=False, values=[0], grid=False, ticks=True),
    ),
    column=alt.Column(
        'type:N',
        header=alt.Header(
            titleOrient='bottom',
            labelOrient='bottom',
            labelPadding=0,
        ),
    )
).properties(
    width=100
).configure_facet(
    spacing=0
).configure_view(
    stroke=None
)
