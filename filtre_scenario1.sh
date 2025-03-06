#!/bin/bash
# Étape 5: Résumer les résultats et vérifier le nombre de variants détectés
echo "Résumé des variants..."
S=$1  # The input gzipped VCF file

# Get the base filename without .vcf.gz extension
filename=$(basename "${S}" .vcf.gz)

# Decompress the file and save it with the correct name
gunzip -c "${S}" > "${filename}.vcf"
# Print the filename
echo "Decompressed file: ${filename}.vcf"





gunzip -c ${filename}.vcf.gz > ${filename}.vcf

# Initialiser un tableau associatif pour stocker les résultats
declare -A scenario_counts

# Scenario 1: L'allèle est présent dans le génome de référence et dans le génome query
scenario_counts[Scénario_1]=$(awk '$4 ~ /^N+$/ && $5 ~ /^N+$/' ${filename}.vcf | wc -l)

# Scenario 2: L'allèle est présent dans le génome de référence mais masqué dans le génome query
scenario_counts[Scénario_2]=$(awk '$4 ~ /^[ACTG]+$/ && $5 ~/^\.$/' ${filename}.vcf | wc -l)

# Scenario 3: L'allèle est absent dans le génome de référence mais présent dans le génome query
scenario_counts[Scénario_3]=$(awk '$4 ~ /^[ACTG]+$/ && $5 ~ /^[ACTG]+$/' ${filename}.vcf | wc -l)

# Scenario 4: L'allèle est absent dans les deux génomes
scenario_counts[Scénario_4]=$(awk '$4 ~ /^\.$/ && $5 ~ /^N+$/' ${filename}.vcf | wc -l)

# Scenario 5: L'allèle est présent dans le génome de référence mais absent dans le génome query
scenario_counts[Scénario_5]=$(awk '$4 ~ /^[ACTG]+$/ && $5 ~ /^\.$/' ${filename}.vcf | wc -l)

# Scenario 6: L'allèle est présent dans le génome query mais absent dans la référence
scenario_counts[Scénario_6]=$(awk '$4 ~ /^N+$/ && $5 ~ /^[ACTG]+$/' ${filename}.vcf | wc -l)

# Scenario 7: L'allèle est présent dans les deux génomes mais avec une différence d'intensité ou de présence
scenario_counts[Scénario_7]=$(awk '$4 ~ /^[ACTG]+$/ && $5 ~ /^N+$/' ${filename}.vcf | wc -l)

# Scenario 8: L'allèle est dans une région masquée dans les deux génomes
scenario_counts[Scénario_8]=$(awk '$4 ~ /^N+$/ && $5 ~ /^\.$/' ${filename}.vcf | wc -l)

# Nom du fichier de sortie CSV
csv_file=$2

# Écriture des résultats dans le fichier CSV
echo "Scénario,Nombre" > "$csv_file"
for key in "${!scenario_counts[@]}"; do
    echo "$key,${scenario_counts[$key]}" >> "$csv_file"
done

echo "Analyse terminée. Les résultats sont enregistrés dans $csv_file."
