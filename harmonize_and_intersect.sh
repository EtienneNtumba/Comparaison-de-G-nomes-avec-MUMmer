#!/bin/bash

# === VÉRIFICATION DES ARGUMENTS ===
if [ "$#" -ne 4 ]; then
    echo "Usage: bash harmonize_and_intersect.sh <reference.vcf> <bed1> <bed2> <output_bed>"
    exit 1
fi

# === RÉCUPÉRATION DES PARAMÈTRES ===
VCF_FILE="$1"       # Fichier VCF de référence
BED1_FILE="$2"      # Premier fichier BED
BED2_FILE="$3"      # Deuxième fichier BED
OUTPUT_BED="$4"     # Fichier de sortie pour l'intersection

# === EXTRAIRE LES CHROMOSOMES DU VCF ===
echo "[1] Extraction des noms de chromosomes du VCF..."
VCF_CHROMS=$(bcftools query -f '%CHROM\n' "$VCF_FILE" | sort -u)

# === AFFICHER LES CHROMOSOMES TROUVÉS ===
echo "[INFO] Chromosomes extraits du VCF :"
echo "$VCF_CHROMS"

# === CRÉER UN FICHIER TEMPORAIRE DE MAPPING ===
MAPPING_FILE="chrom_mapping.txt"
rm -f "$MAPPING_FILE"  # Supprimer l'ancien fichier s'il existe

# === ASSOCIER DYNAMIQUEMENT LES CHROMOSOMES DES BED AVEC LE VCF ===
function generate_chrom_mapping {
    BED_FILE="$1"
    echo "[2] Génération du mapping pour $BED_FILE..."
    
    cut -f1 "$BED_FILE" | sort -u | while read BED_CHROM; do
        BEST_MATCH=$(echo "$VCF_CHROMS" | grep -i "$BED_CHROM" | head -n 1)
        if [ -z "$BEST_MATCH" ]; then
            BEST_MATCH="$BED_CHROM"  # Garder le même nom si aucun match trouvé
        fi
        echo -e "$BED_CHROM\t$BEST_MATCH" >> "$MAPPING_FILE"
    done
}

# Générer le mapping pour les deux fichiers BED
generate_chrom_mapping "$BED1_FILE"
generate_chrom_mapping "$BED2_FILE"

# === FONCTION POUR HARMONISER UN FICHIER BED ===
function harmonize_bed {
    INPUT_BED="$1"
    OUTPUT_BED="$2"
    echo "[3] Harmonisation de $INPUT_BED..."
    
    awk 'NR==FNR{a[$1]=$2; next} {if ($1 in a) $1=a[$1]; print}' "$MAPPING_FILE" "$INPUT_BED" > "$OUTPUT_BED"
    
    # Vérifier et corriger la tabulation
    awk -F' ' 'BEGIN{OFS="\t"} {print $1, $2, $3}' "$OUTPUT_BED" > "${OUTPUT_BED}.fixed"
    mv "${OUTPUT_BED}.fixed" "$OUTPUT_BED"

    # Vérifier si les coordonnées sont bien des nombres
    awk '$2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/' "$OUTPUT_BED" > "${OUTPUT_BED}.filtered"
    mv "${OUTPUT_BED}.filtered" "$OUTPUT_BED"

    echo "[✔] Fichier BED harmonisé généré -> $OUTPUT_BED"
}

# === HARMONISER LES FICHIERS BED ===
HARMONIZED_BED1="${BED1_FILE%.bed}_harmonized.bed"
HARMONIZED_BED2="${BED2_FILE%.bed}_harmonized.bed"

harmonize_bed "$BED1_FILE" "$HARMONIZED_BED1"
harmonize_bed "$BED2_FILE" "$HARMONIZED_BED2"

# === VÉRIFIER SI LES FICHIERS BED HARMONISÉS NE SONT PAS VIDES ===
if [ ! -s "$HARMONIZED_BED1" ] || [ ! -s "$HARMONIZED_BED2" ]; then
    echo "[ERROR] Un des fichiers BED harmonisés est vide ! Vérifiez vos fichiers d'entrée."
    exit 1
fi

# === INTERSECTION AVEC BEDTOOLS ===
echo "[4] Intersection des BED avec bedtools..."
bedtools intersect -a "$HARMONIZED_BED1" -b "$HARMONIZED_BED2" > "$OUTPUT_BED"

# === STATISTIQUES ===
NUM_LINES=$(wc -l < "$OUTPUT_BED")
echo "[✔] Intersection terminée : $OUTPUT_BED avec $NUM_LINES lignes."

# === SUPPRESSION DES FICHIERS TEMPORAIRES ===
rm -f "$MAPPING_FILE"

