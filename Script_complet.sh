#!/bin/bash

# === CONFIGURATION ===
REFERENCE="CP010329_H37Rv.fasta"  # Génome de référence
QUERY="N1274.LR.Asm.fasta"        # Génome à comparer
OUTDIR="genome_comparison"         # Dossier de sortie
THREADS=8                          # Nombre de threads utilisés

mkdir -p $OUTDIR
cd $OUTDIR

# === ÉTAPE 1 : Préparation des fichiers ===
echo "[1] Indexation du génome de référence..."
samtools faidx $REFERENCE

# === ÉTAPE 2 : Comparaison des génomes avec Mummer ===
echo "[2] Détection des SNPs et INDELs avec Mummer..."
nucmer --prefix=genome_comparison $REFERENCE $QUERY
show-snps -ClrH genome_comparison.delta > genome_comparison.snps

# === ÉTAPE 3 : Conversion des SNPs en format VCF ===
echo "[3] Conversion en VCF..."
echo "##fileformat=VCFv4.2" > genome_comparison.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> genome_comparison.vcf
awk '{print $NF "\t" $1 "\t.\t" $2 "\t" $3 "\t.\tPASS\t."}' genome_comparison.snps >> genome_comparison.vcf

# === ÉTAPE 4 : Compression et Indexation du VCF ===
echo "[4] Compression et indexation du VCF..."
bgzip genome_comparison.vcf
tabix -p vcf genome_comparison.vcf.gz

# === ÉTAPE 5 : Vérification des résultats ===
echo "[5] Vérification du fichier VCF final..."
bcftools stats genome_comparison.vcf.gz

# === FIN ===
echo "[✔] Pipeline terminé ! Résultats disponibles dans $OUTDIR/genome_comparison.vcf.gz"
