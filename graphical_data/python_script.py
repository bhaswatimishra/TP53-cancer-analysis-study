import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# to read the original mutation data file
data = pd.read_csv("mutations.txt")

# save the data in csv format
data.to_csv("TP53_cBIO_PORTAL_data.csv", index=False)


# Graph 1: WT vs Mutant distribution across cancer studies
def stacked_bar_chart(csv_file):

    # read the csv file
    data = pd.read_csv(csv_file, sep="\t")

    # condition to find if the TP53 column contained WT.It was neccessary for filtration
    condition = (data['TP53'] == "WT")

    # it form a new column and assign the value as WT if condition is true or mutant if condition is false
    data['STATUS'] = np.where(condition, "WT", "mutant") # syntax - np.where(condition, if_true , if_false)

    # find the value for each study id.
  # pd.crosstab(x,y) -> here x refers to raw and y for column
  # .reset_index was used to assign a separate index named 'STATUS' in the dataframe
    report_data = pd.crosstab(data['STUDY_ID'], data['STATUS']).reset_index()

    #save the report to csv
    report_data.to_csv("mutation_report.csv", index=False)

    plot = report_data

    # convert long study ids into short cancer names which will be impactful for graph representation and cleanliness
    name_list = []
    for study in plot['STUDY_ID']:
        name = study.split("_")[0].upper()
        name_list.append(name)

    # set figure size
    plt.figure(figsize=(25, 20))

    x = name_list
    wt_height = plot['WT']
    # bar chart matplotlib syntax -> plt.bar(x,height,label)
    # x represents the data to be shown in the x axis
    # height represent the ylabel of the graph
    # label represent the name you want to give to that bars.since it is a stacked bar chart hence putting a label will be helpful
    plt.bar(x, wt_height, label="WT")

    mutant_height = plot['mutant']
    # bottom represent to stack the mutant data on the WT bars so that stacked bars can be seen
    plt.bar(x, mutant_height, bottom=wt_height, label="mutant")

    #rotation -> rotates the x axis labels to 45 degree to avoid overlapping between labels
    plt.xticks(rotation=45, fontsize=9)

    # add graph title and axis labels
    plt.title("Association plot of TP53")
    plt.xlabel("study ids")
    plt.ylabel("number of variants type")

    # to show the labels on bars for WT and mutant
    plt.legend()

    # Graph display 
    plt.show()


# Graph 2: TP53 mutation frequency across cancer types
def mutation_frequency(csv_file):

    data = pd.read_csv(csv_file)

    # calculate the total sample = WT + mutant
    data['TOTAL'] = data['WT'] + data['mutant']

    #calculate mutation freq percentage
    data['Freq'] = (data['mutant'] / data['TOTAL']) * 100

    # convert study id into short cancer names
    name_list = []
    for study in data['STUDY_ID']:
        name = study.split("_")[0].upper()
        name_list.append(name)

    plt.figure(figsize=(30, 25))

    x = name_list
    height = data["Freq"]

    # plot mutation freq in bar chart
    plt.bar(x, height)

    plt.xticks(rotation=45)
    plt.title("TP53 Mutation Frequency Across Cancer Types")
    plt.xlabel("Study ids")
    plt.ylabel("Mutation Frequency (%)")

    plt.show()


# Graph 3: TP53 hotspot mutation analysis
def hotspot_analysis(csv_file):

  
    data = pd.read_csv(csv_file, sep="\t")

    # reomove the WTraws as we only need to anlyse the muated part
    modified_csv = data[data['TP53'] != 'WT'].copy()

    # remove the SAMPLE_ID column as the column as not needed further
    mod_csv = modified_csv.drop(columns=['SAMPLE_ID'])

    # empty list to store the mutation
    positions = []

    # loop to generate the numeric postion for each mutation
    for mutation in mod_csv['TP53']:

        # used ragex to search for number in mutation string
        match = re.search(r"\d+", mutation) # syntax -> \d -> search for number
                                            # \d+ search for more than 1 numbers which are present in continuous manner

        if match:
            # if number was found store it in the pos
            pos = match.group()
            positions.append(pos)
        else:
            # if no number was found store none
            positions.append(None)

    # add the extracted positions into a new column
    mod_csv['positions'] = positions

    # count the combination of each study id with each position
    # reset index(name = "counts") -> make a separate column named as 'counts' to show howmany times the study ids appeared with same position
    hotspot_data = mod_csv.groupby(['STUDY_ID', 'positions']).size().reset_index(name="counts")
    hotspot_data.to_csv("hotspot_data.csv", index=False)

    # mod_csv['positions'].value_counts() -> count how many times each position was repeated overall
    hotspot_pos = mod_csv['positions'].value_counts().reset_index(name="counts")

    # filter out all the positions which have count value less than 20
    new_data = hotspot_pos[hotspot_pos["counts"] >= 20].copy()

    plt.figure(figsize=(30, 25))

    x = new_data['positions']
    height = new_data['counts']
    plt.bar(x, height)
    plt.xticks(rotation=90)
    plt.title("Hotspot analysis report")
    plt.xlabel("mutation position")
    plt.ylabel("freq")
    plt.show()


# Graph 4: Heatmap of top TP53 hotspots across cancer types
def heatmap(csv_file):

    df = pd.read_csv(csv_file)

    # Top 3 hotspot positions were selected from analysing the bar graph of hotspot mutation analysis
    value_list = [273, 248, 175]

    #filter out only the desired hotspot positions
    top3_data = df[df['positions'].isin(value_list)].copy()
    top3_data.to_csv("new.csv", index=False)

    # convert the data to pivot table format to generate heatmap
    heatmap_data = top3_data.pivot(
        index='STUDY_ID',
        columns='positions',
        values='counts'
    ).fillna(0)


    name_list = []

    for study in heatmap_data.index:
        name = study.split("_")[0].upper()
        name_list.append(name)

    #plot the seaborn heatmap data
    sns.heatmap(
        heatmap_data,
        annot=True,          # show the value of count in each cell
        cmap='Purples',      # choose theme as purple
        linewidths=0.5       # add thin lines between cells
    )

    
    plt.title("Cancer-specific Association of TP53 Hotspots")

    # as the study ids were too large to plot on graph so yticks was used to replace the study id to short cancer name
    plt.yticks(
        ticks=range(len(name_list)),
        labels=name_list
    )

    # adjust the layout for heatmap
    plt.tight_layout()

    
    plt.show()



stacked_bar_chart("TP53_cBIO_PORTAL_data.csv")
mutation_frequency("mutation_report.csv")
hotspot_analysis("TP53_cBIO_PORTAL_data.csv")
heatmap("hotspot_data.csv")
