#!/bin/bash

# Vérification des arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.vcf.gz> <H37Rv.bed> <Lx.bed>"
    exit 1
fi

# Variables d'entrée
VCF_GZ=$1
H37RV_BED=$2
LX_BED=$3
TEMP_VCF="temp_uncompressed.vcf"
OUTPUT_VCF="masked_final.vcf.gz"

echo "🔍 Étape 1: Vérification des noms des contigs dans VCF et BED..."

# Extraire les noms uniques du VCF et des BED files
bcftools view -H "$VCF_GZ" | cut -f1 | sort | uniq > vcf_contigs.txt
cut -f1 "$H37RV_BED" | sort | uniq > bed_H37Rv_contigs.txt
cut -f1 "$LX_BED" | sort | uniq > bed_Lx_contigs.txt

# Comparer les différences
diff vcf_contigs.txt bed_H37Rv_contigs.txt > rename_map.txt
diff vcf_contigs.txt bed_Lx_contigs.txt >> rename_map.txt

# Vérifier si des différences existent
if [ -s rename_map.txt ]; then
    echo "⚠️ Des différences de noms ont été trouvées entre le VCF et les BED files."
    echo "🛠️ Création d'un fichier de correspondance..."
    awk '{print $2, $1}' rename_map.txt > rename_map_fixed.txt

    echo "🔄 Correction des noms dans les fichiers BED..."
    awk 'NR==FNR {map[$1] = $2; next} {if ($1 in map) $1 = map[$1]; print}' rename_map_fixed.txt "$H37RV_BED" > H37Rv_fixed.bed
    awk 'NR==FNR {map[$1] = $2; next} {if ($1 in map) $1 = map[$1]; print}' rename_map_fixed.txt "$LX_BED" > Lx_fixed.bed
else
    echo "✅ Aucune différence détectée, utilisation des fichiers d'origine."
    cp "$H37RV_BED" H37Rv_fixed.bed
    cp "$LX_BED" Lx_fixed.bed
fi

echo "🔓 Étape 2: Décompression du fichier VCF.gz..."
bcftools view -O v "$VCF_GZ" > "$TEMP_VCF"

echo "📍 Étape 3: Extraction des positions des variantes..."
bcftools view -H "$TEMP_VCF" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' > vcf_positions.bed

echo "🛑 Étape 4: Application du masque H37Rv sur REF..."
bedtools intersect -a vcf_positions.bed -b H37Rv_fixed.bed | awk 'BEGIN{OFS="\t"} {print $1, $3}' > mask_H37Rv_positions.txt

echo "🛑 Étape 5: Application du masque Lx sur ALT..."
bedtools intersect -a vcf_positions.bed -b Lx_fixed.bed | awk 'BEGIN{OFS="\t"} {print $1, $3}' > mask_Lx_positions.txt

echo "✏️ Étape 6: Modification du fichier VCF..."
awk '
    BEGIN {
        while ((getline < "mask_H37Rv_positions.txt") > 0) {
            mask_H37Rv[$1 ":" $2] = 1
        }
        close("mask_H37Rv_positions.txt")
        
        while ((getline < "mask_Lx_positions.txt") > 0) {
            mask_Lx[$1 ":" $2] = 1
        }
        close("mask_Lx_positions.txt")
    }
    /^#/ { print; next }
    {
        key = $1 ":" $2
        if (key in mask_H37Rv) $4 = "N"
        if (key in mask_Lx) $5 = "N"
        print
    }
' "$TEMP_VCF" > masked_final.vcf

echo "📦 Étape 7: Recompression du fichier VCF en VCF.gz..."
bgzip -c masked_final.vcf > "$OUTPUT_VCF"
tabix -p vcf "$OUTPUT_VCF"

# Nettoyage des fichiers temporaires
rm -f "$TEMP_VCF" vcf_positions.bed mask_H37Rv_positions.txt mask_Lx_positions.txt masked_final.vcf rename_map.txt rename_map_fixed.txt

echo "🎉 Processus terminé. Fichier final généré: $OUTPUT_VCF"
