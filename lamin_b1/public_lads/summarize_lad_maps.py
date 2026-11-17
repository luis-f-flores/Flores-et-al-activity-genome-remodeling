
from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)


lad_files = {
    "dorsal_npc": "Dorsal_NPCs_merged_lads.bed",
    "hippocampus": "LaminB1_Hippocampus_LAD_peaks_30wk_mm10.bed",
    "liver": "Liver_merged_lads.bed",
    "mef": "MEF_merged_lads.bed",
    "ventral_npc": "Ventral_NPCs_merged_lads.bed",
    "lung": "combined_E13_Lung_LaminB.bed"
}


total_genome_bases = 2725521370


output_file = "LAD_summary.txt"


with open(output_file, "w") as out:
    
    out.write("Tissue\tTotal_LAD_bases\tTotal_Genome_bases\tGenome_fraction_in_LADs\n")

    
    for tissue, file in lad_files.items():
        total_lad_bases = 0  

        try:
            with open(file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")

                    
                    if len(parts) < 3:
                        continue  
                    
                    try:
                        start, end = int(parts[1]), int(parts[2])
                        total_lad_bases += (end - start)
                    except ValueError:
                        continue  

        except FileNotFoundError:
            print(f"WARNING: File {file} not found. Skipping.")
            continue

        
        genome_fraction = total_lad_bases / total_genome_bases

        
        out.write(f"{tissue}\t{total_lad_bases}\t{total_genome_bases}\t{genome_fraction:.6f}\n")

print(f"LAD summary saved to {output_file}")
