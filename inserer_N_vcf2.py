import vcfpy
import sys

# === CONFIGURATION ===
VCF_INPUT = sys.argv[1]  # Fichier VCF compressé en entrée
VCF_OUTPUT = sys.argv[2]  # Fichier VCF compressé en sortie
MASK_H37RV = sys.argv[3]  # Masque pour H37Rv (REF)
MASK_LX = sys.argv[4]  # Masque pour Lx (ALT)

# === CHARGER LES MASQUES BED ===
def load_bed_mask(bed_file):
    """Charge les positions du fichier BED en un set de positions."""
    mask_positions = set()
    with open(bed_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            for pos in range(start, end + 1):
                mask_positions.add((chrom, pos))
    return mask_positions

# Charger les masques
mask_h37rv = load_bed_mask(MASK_H37RV)  # Masque REF (H37Rv)
mask_lx = load_bed_mask(MASK_LX)  # Masque ALT (Lx)

# === OUVERTURE DU FICHIER VCF ===
reader = vcfpy.Reader.from_path(VCF_INPUT)
writer = vcfpy.Writer.from_path(VCF_OUTPUT, reader.header)

# === TRAITEMENT DES VARIANTS ===
for record in reader:
    chrom = record.CHROM
    pos = record.POS
    
    # Appliquer le masque H37Rv (modifier REF)
    if (chrom, pos) in mask_h37rv:
        record.REF = "N"

    # Appliquer le masque Lx (modifier ALT)
    if (chrom, pos) in mask_lx:
        if record.ALT and record.ALT[0] != ".":
            record.ALT = [vcfpy.SymbolicAllele("N")]

    writer.write_record(record)

# Fermeture des fichiers
reader.close()
writer.close()

print(f"[✔] Masquage terminé : fichier généré -> {VCF_OUTPUT}")
