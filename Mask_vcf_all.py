import vcfpy
import sys
# === Normaliser les noms de chromosomes ===
def normalize_chrom(chrom):
    """Normalise les noms des chromosomes pour qu'ils correspondent entre BED et VCF."""
    chrom = chrom.replace("chr", "").replace("_", "").upper()
    return chrom

# === Charger les régions BED ===
def load_bed_file(bed_file):
    """Charge les régions d'intérêt depuis un fichier BED en tant que dictionnaire {chrom: set(positions)}"""
    bed_regions = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                chrom, start, end = normalize_chrom(fields[0]), int(fields[1]), int(fields[2])
                if chrom not in bed_regions:
                    bed_regions[chrom] = set()
                bed_regions[chrom].update(range(start, end + 1))
    return bed_regions

# === Vérifier si une position est dans le masque BED ===
def is_in_bed(chrom, pos, bed_regions):
    """Vérifie si une position donnée se trouve dans le dictionnaire de régions BED."""
    chrom = normalize_chrom(chrom)  # Normalisation du nom du chromosome
    return chrom in bed_regions and pos in bed_regions[chrom]

# === Modifier REF et ALT en fonction des masques BED ===
def mask_ref_and_alt_in_vcf(vcf_file, mask_h37rv_bed, mask_l3_bed, output_vcf):
    """Masque REF et ALT dans un fichier VCF en fonction des fichiers BED fournis."""
    
    # Charger les fichiers BED
    mask_h37rv_regions = load_bed_file(mask_h37rv_bed)  # Masque pour REF
    mask_l3_regions = load_bed_file(mask_l3_bed)  # Masque pour ALT

    # Ouvrir le fichier VCF en lecture
    reader = vcfpy.Reader.from_path(vcf_file)
    
    # Ouvrir le fichier VCF en écriture compressée
    writer = vcfpy.Writer.from_path(output_vcf, reader.header)

    # Parcourir chaque variant et appliquer les masques
    for record in reader:
        chrom, pos = record.CHROM, record.POS

        # Appliquer le masque REF
        if is_in_bed(chrom, pos, mask_h37rv_regions):
            record.REF = "N"
        
        # Appliquer le masque ALT
        if is_in_bed(chrom, pos, mask_l3_regions):
            record.ALT = [vcfpy.Substitution(type_="SNV", value="N")]

        # Écrire l'entrée modifiée
        writer.write_record(record)

    # Fermeture des fichiers
    reader.close()
    writer.close()

    print(f"[✔] Masquage terminé : fichier généré -> {output_vcf}")

# === Exemple d'utilisation ===
vcf_file = sys.argv[1]  # Fichier VCF d'entrée (compressé)
mask_h37rv_bed = sys.argv[2]  # Masque pour REF
mask_l3_bed = sys.argv[3]  # Masque pour ALT
output_vcf = sys.argv[4]  # Fichier de sortie (compressé)

mask_ref_and_alt_in_vcf(vcf_file, mask_h37rv_bed, mask_l3_bed, output_vcf)
