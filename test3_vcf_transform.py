import vcfpy

# Charger le fichier BED et extraire les régions d'intérêt
def load_bed_file(bed_file):
    bed_regions = []
    with open(bed_file, 'r') as bed:
        for line in bed:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                # Ajouter les régions (chromosome, start, end)
                bed_regions.append((fields[0], int(fields[1]), int(fields[2])))
    return bed_regions

# Vérifier si une position se trouve dans les régions du fichier BED
def is_in_bed(chrom, pos, bed_regions):
    for region in bed_regions:
        if region[0] == chrom and region[1] <= pos < region[2]:
            return True
    return False

# Modifier REF par 'N' et ALT par 'N' dans les régions spécifiées dans le fichier BED
def mask_ref_and_alt_in_vcf(vcf_file, intersect_bed, mask_h37rv_bed, mask_l3_bed, output_vcf):
    intersect_regions = load_bed_file(intersect_bed)
    mask_h37rv_regions = load_bed_file(mask_h37rv_bed)
    mask_l3_regions = load_bed_file(mask_l3_bed)
    
    # Ouvrir le fichier VCF en lecture
    reader = vcfpy.Reader.from_path(vcf_file)
    
    # Créer un objet Alternate avec la valeur 'N'
    n_alt = vcfpy.record.Substitution(type_='SNV', value='N')
    
    # Ouvrir le fichier de sortie en mode écriture
    with open(output_vcf, 'w') as out_file:
        writer = vcfpy.Writer(out_file, reader.header)
        # Parcourir chaque record et modifier REF et ALT si nécessaire
        for record in reader:
            # Appliquer le masque pour intersect.bed : remplacer REF et ALT par 'N'
            if is_in_bed(record.CHROM, record.POS, intersect_regions):
                record.REF = 'N'
                record.ALT = [n_alt]
            # Appliquer le masque pour mask_H37RV.bed : remplacer REF par 'N'
            elif is_in_bed(record.CHROM, record.POS, mask_h37rv_regions):
                record.REF = 'N'
            # Appliquer le masque pour mask_L3.bed : remplacer ALT par 'N'
            elif is_in_bed(record.CHROM, record.POS, mask_l3_regions):
                record.ALT = [n_alt]
            # Écrire le record modifié dans le fichier de sortie
            writer.write_record(record)

# Exemple d'utilisation
vcf_file = sys.argv[1]  # Fichier VCF d'entrée
intersect_bed = sys.argv[2]  # Fichier BED pour intersection
mask_h37rv_bed = sys.argv[3]  # Fichier BED pour H37Rv
mask_l3_bed = sys.argv[4]  # Fichier BED pour L3
output_vcf = sys.argv[4]  # Fichier VCF de sortie

mask_ref_and_alt_in_vcf(vcf_file, intersect_bed, mask_h37rv_bed, mask_l3_bed, output_vcf)
