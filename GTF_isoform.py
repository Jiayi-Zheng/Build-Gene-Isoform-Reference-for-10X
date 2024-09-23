#this script functions to set all transcripts isoform of a certain gene (RPGR here) in gtf file as individual genes, thereby enabling isoform detection using 10X pipeline. However, it should be noted that given the limitation of short reads 10X seq, result requires further verification. 



# Define things
directory = "/Users/jiajia/Desktop/" #define directory of things
gtf_file = directory + "genes.gtf" #gtf file from 10 official genome reference folder
rpgr_output = directory + 'temp.txt'
gene = 'RPGR'
output_file = directory + gene + '.gtf'


import subprocess
subprocess.run(f'grep \'"{gene}"\' "{gtf_file}" > "{rpgr_output}"', shell=True)

# Global variables
n_transcript = 0

# Function to modify IDs
def modify_ids(n_transcript, gene_id, transcript_name, hgnc_id, havana_gene):
    # Create new gene_id by modifying ENSG
    new_gene_id = gene_id.replace("ENSG00", f"ENSG{'{:02}'.format(n_transcript)}")
    
    # Use transcript_name for the new gene_name
    new_gene_name = transcript_name
    
    # Modify hgnc_id
    hgnc_prefix = "HGNC:"
    new_hgnc_id = hgnc_id.replace("HGNC:", "HGNC:")
    hgnc_digits = int(new_hgnc_id.split(":")[1]) - n_transcript + 80000  # Subtract from 10000 and change to start from 9
    new_hgnc_id = f"{hgnc_prefix}{hgnc_digits}"
    
    # Modify havana_gene
    new_havana_gene = havana_gene.replace("OTTHUMG00", f"OTTHUMG{'{:02}'.format(n_transcript)}")
    
    return new_gene_id, new_gene_name, new_hgnc_id, new_havana_gene

# Open and read input file
with open(rpgr_output, 'r') as infile:
    lines = infile.readlines()

print('Total ' + str(len(lines)) + ' lines containing target gene.')

#find which lines are those transcripts on
transcript_indices = []

# Iterate over the list with index
for i, line in enumerate(lines):
    # Check if 'HAVANA\ttranscript\t' is in the string
    if 'HAVANA\ttranscript\t' in line:
        transcript_indices.append(i)

print('Total ' + str(len(transcript_indices)) + ' transcript isoforms.')

# Store the modified lines here
modified_lines = []

# Iterate over each line
for i, line in enumerate(lines):
    # Delete any line that contains 'HAVANA gene'
    if i == 0:
        continue  # Skip this line
    
    # Check if it's a 'HAVANA transcript' line
    elif i in transcript_indices:
        # Duplicate the line
        transcript_line = line

        # Extract everything starting from 'gene_id' and save it
        gene_info_index = line.find('gene_id')
        #which character does gene_id starts
        gene_info = line[gene_info_index:].strip()
        
        # Remove gene info from the line
        line = line[:gene_info_index].strip()

        #here we split the transcript line into basic info + attributes
        
        # Split gene_info into its components
        gene_id = gene_info.split('gene_id ')[1].split(';')[0].replace('"', '')
        hgnc_id = gene_info.split('hgnc_id ')[1].split(';')[0].replace('"', '')
        havana_gene = gene_info.split('havana_gene ')[1].split(';')[0].replace('"', '')

        # Extract transcript_name from the gene_info
        transcript_name = gene_info.split('transcript_name ')[1].split(';')[0].replace('"', '')
        
        # Modify the values according to the n_transcript
        new_gene_id, new_gene_name, new_hgnc_id, new_havana_gene = modify_ids(n_transcript, gene_id, transcript_name, hgnc_id, havana_gene)
        
        # Create the new 'gene' line based on the modifications and change 'transcript' to 'gene'
        new_gene_line = line.replace('transcript', 'gene', 1)  # Change 'transcript' to 'gene'
        new_gene_line = f'{new_gene_line}\tgene_id "{new_gene_id}"; gene_version "14"; gene_type "protein_coding"; gene_name "{new_gene_name}"; level 2; hgnc_id "{new_hgnc_id}"; tag "ncRNA_host"; havana_gene "{new_havana_gene}";\n'
        
        # **Insert the new 'gene' line first, then update and append the original transcript line**
        modified_lines.append(new_gene_line)  # New 'gene' line first

        # Update the original transcript line with the new gene_id, gene_name, hgnc_id, and havana_gene
        updated_transcript_line = transcript_line.replace(gene_id, new_gene_id)
        updated_transcript_line = updated_transcript_line.replace('gene_name "RPGR"', f'gene_name "{new_gene_name}"')
        updated_transcript_line = updated_transcript_line.replace(hgnc_id, new_hgnc_id)
        updated_transcript_line = updated_transcript_line.replace(havana_gene, new_havana_gene)
        
        # Append the updated transcript line
        modified_lines.append(updated_transcript_line)
        
        # Now modify all lines in this transcript
        start_index = transcript_indices[n_transcript]+1
        if i != transcript_indices[-1]:
            end_index = transcript_indices[n_transcript+1]-1
        elif i == transcript_indices[-1]:
            end_index = len(lines)-1

        n_transcript += 1

    elif i in range(start_index,end_index+1):
        modified_transcript_line = line
        modified_transcript_line = modified_transcript_line.replace(gene_id, new_gene_id)
        modified_transcript_line = modified_transcript_line.replace("gene_name \"RPGR\"", f"gene_name \"{new_gene_name}\"")
        modified_transcript_line = modified_transcript_line.replace(hgnc_id, new_hgnc_id)
        modified_transcript_line = modified_transcript_line.replace(havana_gene, new_havana_gene)
        modified_lines.append(modified_transcript_line)

# Write the modified content to the new output file
with open(rpgr_output, 'w') as outfile:
    outfile.writelines(modified_lines)

print('Total ' + str(len(modified_lines)) + ' lines printed.')
expected_n = len(lines) + len(transcript_indices) -1
print('Total ' + str(expected_n) + ' lines expected.')

# Step 2: Remove the block containing "RPGR" and find its location
result = subprocess.run(['grep', '-n', f'"{gene}"', gtf_file], stdout=subprocess.PIPE, text=True)
lines = result.stdout.strip().split('\n')
line_numbers = [int(line.split(':')[0]) for line in lines]
start_index = line_numbers[0]-1
end_index = line_numbers[-1]+1

subprocess.run(['sed', '-n', f'1,{start_index}p', gtf_file], stdout=open(output_file, 'w'))
subprocess.run(['cat', rpgr_output], stdout=open(output_file, 'a'))
subprocess.run(['sed', '-n', f'{end_index},$p', gtf_file], stdout=open(output_file, 'a'))

#Check point
result = subprocess.run(['wc', '-l', gtf_file], stdout=subprocess.PIPE, text=True)
nline_gtf = int(result.stdout.split()[0])
result = subprocess.run(['wc', '-l', gtf_file], stdout=subprocess.PIPE, text=True)
nline_new_gtf = int(result.stdout.split()[0])
print(str(nline_gtf) + ' lines in original gtf file.')
print(str(nline_new_gtf) + ' lines in new gtf file.')
expected_n = nline_gtf + len(transcript_indices) -1
print('Total ' + str(expected_n) + ' lines expected in new gtf file.')



#now run cellranger pipeline for build ref and alignment
#cellranger mkref --genome=RPGR --fasta=/usersdata/joyzheng/yard/refdata-gex-GRCh38-2020-A/fasta/genome.fa --genes=/usersdata/joyzheng/GW21/RPGR.gtf
#cellranger count --id GW21 --transcriptome /usersdata/joyzheng/GW21/RPGR --fastqs /usersdata/joyzheng/GW21 --sample GW21 --localcores 24 --localmem 128
