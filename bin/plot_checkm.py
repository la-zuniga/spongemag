import re
import argparse
import os
import pandas as pd
import plotly.express as px


# read in the files

def read_file(data):
    with open(data, 'r', encoding='utf-8') as file:
        file_lines = file.readlines()
    return file_lines

# paring the file
def parse_report(report):
    header = report[0].strip().split("\t")
    data = [line.strip().split("\t") for line in report[1:]]
    df = pd.DataFrame(data, columns = header)
    num_cols = ["Completeness", "Contamination", "GC_Content", "Genome_Size", "Total_Coding_Sequences"]
    df[num_cols] = df[num_cols].apply(pd.to_numeric)
    return df


# Generate a scatter plot of Completeness vs. Contamination
def plot_quality(df, output_file):
    fig = px.scatter(
        df,
        x="Contamination",
        y="Completeness",
        text="Name",
        title="CheckM2 Quality Metrics",
        labels={"Contamination": "Contamination (%)", "Completeness": "Completeness (%)"},
        hover_data=["GC_Content", "Genome_Size", "Total_Coding_Sequences"],
    )
    fig.update_traces(marker=dict(size=10, opacity=0.7, color="blue"), textposition="top center")

# Save the plot based on file extension
    if output_file.endswith(".html"):
        fig.write_html(output_file)
    elif output_file.endswith((".png", ".jpg", ".jpeg", ".svg", ".pdf")):
        fig.write_image(output_file)
    else:
        raise ValueError("Unsupported output format. Use .html, .png, .jpg, .svg, or .pdf")

# Set up argument parsing
def main():
    parser = argparse.ArgumentParser(description="Make some plots from the checkm2 output.")
    parser.add_argument('-i', '--report', required=True, help="Path to the qauality report.")
    parser.add_argument('-o', '--output', required=True, help="Path to save the plots.")
    # Parse the arguments
    args = parser.parse_args()
    # Call the function with the parsed arguments
    concoct_report = read_file(args.report)
    df = parse_report(concoct_report)
    plot_quality(df, args.output)


# Run the script if it is executed directly
if __name__ == "__main__":
    main()

