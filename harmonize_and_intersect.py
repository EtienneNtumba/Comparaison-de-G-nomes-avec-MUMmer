import pandas as pd
import sys
import subprocess

# === Vérification des arguments ===
if len(sys.argv) != 5:
    print("Usage: python harmonize_and_intersect.py <chrom_mapping.txt> <bed1> <bed2> <output_bed>")
    sys.exit(1)

# === Récupération des fichiers depuis les arguments ===
mapping_file = sys.argv[1]   # Fichier de correspondance des chromosomes
bed1_file = sys.argv[2]      # Premier fichier BED (H37Rv)
bed2_file = sys.argv[3]      # Deuxième fichier BED (N1274)
output_bed = sys.argv[4]     # Fichier de sortie pour l'intersection

# === Charger le mapping des chromosomes ===
def load_chrom_mapping(mapping_file):
    """Charge le fichier de correspondance des chromosomes sous forme de dictionnaire."""
    mapping = pd.read_csv(mapping_file, sep="\t", header=0)
    return dict(zip(mapping["Chrom_H37Rv"], mapping["Chrom_N1274"]))

# === Harmoniser un fichier BED ===
def harmonize_bed_file(input_bed, output_bed, chrom_mapping):
    """Charge un fichier BED et remplace les noms de chromosomes selon le mapping."""
    bed_df = pd.read_csv(input_bed, sep="\t", header=None, names=["chrom", "start", "end"])

    # Remplacement des noms de chromosomes selon le mapping
    bed_df["chrom"] = bed_df["chrom"].map(chrom_mapping).fillna(bed_df["chrom"])

    # Sauvegarde du fichier BED harmonisé
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
    print(f"[✔] Fichier BED harmonisé généré -> {output_bed}")

# === Étape 1 : Charger le mapping ===
chrom_mapping = load_chrom_mapping(mapping_file)

# === Étape 2 : Harmoniser les fichiers BED ===
harmonized_bed1 = bed1_file.replace(".bed", "_harmonized.bed")
harmonized_bed2 = bed2_file.replace(".bed", "_harmonized.bed")

harmonize_bed_file(bed1_file, harmonized_bed1, chrom_mapping)
harmonize_bed_file(bed2_file, harmonized_bed2, chrom_mapping)

# === Étape 3 : Exécuter l'intersection avec bedtools ===
intersect_command = f"bedtools intersect -a {harmonized_bed1} -b {harmonized_bed2} > {output_bed}"
subprocess.run(intersect_command, shell=True, check=True)
print(f"[✔] Intersection des BED terminée -> {output_bed}")

# === Étape 4 : Vérifications ===
num_lines = sum(1 for _ in open(output_bed))
print(f"[✔] Nombre de lignes dans l'intersection : {num_lines}")
