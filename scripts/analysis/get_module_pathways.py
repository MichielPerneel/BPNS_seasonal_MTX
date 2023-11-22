## Run this program if you want more information for the obtained modules
## It needs the module content files from the WGCNA analysis

import requests
import os

# Define input and output directories
input_dir = "data/analysis/WGCNA"
output_dir = "data/analysis/WGCNA/information"

# Create the output directory if it doesn't already exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterate over all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith("_content.txt"):
        # Define input and output files for this file
        input_file = os.path.join(input_dir, filename)
        output_file1 = os.path.join(output_dir, filename.replace("_content.txt", "_pathways.txt"))
        output_file2 = os.path.join(output_dir, filename.replace("_content.txt", "_pathway_counts.txt"))
        
        # Create a dictionary to store the KO IDs and their corresponding pathway
        ko_pathway_dict = {}

        # Read the input file
        with open(input_file, "r") as f:
            for line in f:
                # Parse the KO ID from the line
                ko_id = line.split("\t")[0]
                
                # Send a GET request to the KEGG API to get the pathway information
                response = requests.get(f"https://rest.kegg.jp/link/pathway/{ko_id}")
                
                # Parse the pathway information from the response
                pathways = [pathway.split(":")[2] for pathway in response.text.split("\n") if pathway]
                
                # Add all pathways to the dictionary if they contain the letters 'map'
                for pathway in pathways:
                    if 'map' in pathway:
                        ko_pathway_dict[ko_id] = pathway

        # Create a dictionary to store the pathway counts
        pathway_count_dict = {}

        # Read the values from the ko_pathway_dict and count the occurrences of each pathway
        for pathway in ko_pathway_dict.values():
            if pathway not in pathway_count_dict:
                pathway_count_dict[pathway] = 1
            else:
                pathway_count_dict[pathway] += 1

        # Write the output files
        with open(output_file1, "w") as f1, open(output_file2, "w") as f2:
            # Write the header for the first output file
            f1.write("KO ID\tPathway\tName\tClass\n")
            
            # Iterate over the KO IDs in the dictionary and write the corresponding pathway information to the first output file
            for ko_id, pathway in ko_pathway_dict.items():
                # Send a GET request to the KEGG API to get the pathway information
                response = requests.get(f"https://rest.kegg.jp/get/path:{pathway}")
                
                # Parse the pathway name and class from the response
                pathway_name = ""
                pathway_class = ""
                for line in response.text.split("\n"):
                    if line.startswith("NAME"):
                        pathway_name = line.split("        ")[1]
                    elif line.startswith("CLASS"):
                        pathway_class = line.split("    ")[1]
                
                # Write the KO ID, pathway, pathway name, and pathway class to the first output file
                f1.write(f"{ko_id}\t{pathway}\t{pathway_name}\t{pathway_class}\n")
            
            # Write the header for the second output file
            f2.write("Pathway\tCount\n")
            
            # Iterate over the pathway counts and write them to the second output file
            for pathway, count in pathway_count_dict.items():
                f2.write(f"{pathway}\t{count}\n")
