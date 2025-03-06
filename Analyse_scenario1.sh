#!/bin/bash


# Étape 5: Résumer les résultats et vérifier le nombre de variants détectés
echo "Résumé des variants..."

gunzip -c masked_variants.vcf.gz > masked_variants.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' masked_variants.vcf > variants_summary.txt


# Étape 6: Analyser les 8 scénarios en fonction des résultats obtenus

echo "Analyse des scénarios :"

# Scenario 1: L'allèle est présent dans le génome de référence et dans le génome query
echo "Scénario 1: Présence dans les deux génomes (référence et query)"
grep -P '\t[ACTG]+[\t][ACTG]+' masked_variants.vcf >> scenario_1.txt

# Scenario 2: L'allèle est présent dans le génome de référence mais masqué dans le génome query
echo "Scénario 2: Présence dans référence, masqué dans query"
grep -P '\t[ACTG]+[\t]\.+' masked_variants.vcf >> scenario_2.txt

# Scenario 3: L'allèle est absent dans le génome de référence mais présent dans le génome query
echo "Scénario 3: Absence dans référence, présent dans query"
grep -P '\t\.+[\t][ACTG]+' masked_variants.vcf >> scenario_3.txt

# Scenario 4: L'allèle est absent dans les deux génomes
echo "Scénario 4: Absence dans les deux génomes"
grep -P '\t\.+[\t]\.+' masked_variants.vcf >> scenario_4.txt

# Scenario 5: L'allèle est présent dans le génome de référence mais absent dans le génome query
echo "Scénario 5: Présence dans référence, absence dans query"
grep -P '\t[ACTG]+[\t]\.+' masked_variants.vcf >> scenario_5.txt

# Scenario 6: L'allèle est présent dans le génome query mais absent dans la référence
echo "Scénario 6: Présence dans query, absence dans référence"
grep -P '\t\.+[\t][ACTG]+' masked_variants.vcf >> scenario_6.txt

# Scenario 7: L'allèle est présent dans les deux génomes mais avec une différence d'intensité ou de présence
echo "Scénario 7: Différence d'intensité entre référence et query"
grep -P '\t[ACTG]+[\t][ACTG]+' masked_variants.vcf | grep -v 'N' >> scenario_7.txt

# Scenario 8: L'allèle est dans une région masquée dans les deux génomes
echo "Scénario 8: Allèle dans une région masquée"
grep -P '\t[ACTG]+[\t]\.+' masked_variants.vcf >> scenario_8.txt

echo "Analyse terminée. Les résultats pour chaque scénario sont enregistrés dans les fichiers scenario_1.txt, scenario_2.txt, ..."

# Étape 7: Nettoyage
echo "Nettoyage des fichiers temporaires..."
rm query2ref.sam query_variants.vcf

echo "Tout est prêt. Analyse terminée."

