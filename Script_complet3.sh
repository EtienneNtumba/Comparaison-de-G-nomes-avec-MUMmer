#!/bin/bash

# === CONFIGURATION ===
#REFERENCE="/home/p0129674/Data_Analysis/Outbreat1/Reference_genome_fastq/COMPARAISON/CP010329_H37Rv.fasta"  # Génome de référence
REFERENCE=$1
QUERY_DIR="/home/p0129674/Data_Analysis/Outbreat1/Reference_genome_fastq/COMPARAISON/FASTA"         # Dossier contenant les génomes à comparer
OUTDIR="genome_comparison2"         # Dossier de sortie
THREADS=8                           # Nombre de threads utilisés

# Vérifier si le fichier de référence existe
if [ ! -f "$REFERENCE" ]; then
    echo "[ERROR] Fichier de référence $REFERENCE introuvable ! Vérifiez le chemin."
    echo "[INFO] Contenu du répertoire actuel :"
    ls -lh
    exit 1
fi

mkdir -p $OUTDIR
cd $OUTDIR

# === ÉTAPE 1 : Préparation des fichiers ===
echo "[1] Indexation du génome de référence..."
samtools faidx "$REFERENCE" || { echo "[ERROR] Échec de l'indexation du génome ! Vérifiez que $REFERENCE est valide."; exit 1; }

# Vérifier si des fichiers query existent
QUERY_FILES=($QUERY_DIR/*.fasta)
if [ ${#QUERY_FILES[@]} -eq 0 ]; then
    echo "[ERROR] Aucun fichier .fasta trouvé dans $QUERY_DIR ! Vérifiez le chemin."
    echo "[INFO] Contenu de $QUERY_DIR :"
    ls -lh "$QUERY_DIR"
    exit 1
fi

# === ÉTAPE 2 : Comparaison des génomes avec Mummer ===
echo "[2] Détection des SNPs et INDELs avec Mummer pour chaque fichier de requête..."
for QUERY in "${QUERY_FILES[@]}"; do
    BASENAME=$(basename "$QUERY" .fasta)
    PREFIX="${BASENAME}_comparison"
    
    echo "[2.1] Traitement de $QUERY..."
    nucmer --prefix=$PREFIX "$REFERENCE" "$QUERY" || { echo "[ERROR] Échec de nucmer pour $QUERY"; continue; }
    show-snps -ClrH $PREFIX.delta > $PREFIX.snps || { echo "[ERROR] Échec de show-snps pour $QUERY"; continue; }

    # Vérification si le fichier SNPs est vide
    if [ ! -s "$PREFIX.snps" ]; then
        echo "[ERROR] Fichier $PREFIX.snps est vide, passage au suivant."
        continue
    fi

    # === ÉTAPE 3 : Conversion des SNPs en format VCF ===
    echo "[3] Conversion en VCF pour $QUERY..."
    echo "##fileformat=VCFv4.2" > $PREFIX.vcf
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> $PREFIX.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> $PREFIX.vcf
    awk '{print $NF "\t" $1 "\t.\t" $2 "\t" $3 "\t.\tPASS\t.\tGT\t1/1"}' $PREFIX.snps >> $PREFIX.vcf || { echo "[ERROR] Échec de la conversion VCF pour $QUERY"; continue; }

    # === ÉTAPE 4 : Compression et Indexation du VCF ===
    echo "[4] Compression et indexation du VCF pour $QUERY..."
    bgzip $PREFIX.vcf && tabix -p vcf $PREFIX.vcf.gz || { echo "[ERROR] Échec de la compression/indexation pour $QUERY"; continue; }

    # === ÉTAPE 5 : Vérification des résultats ===
    echo "[5] Vérification du fichier VCF final pour $QUERY..."
    bcftools stats $PREFIX.vcf.gz || { echo "[ERROR] Problème avec bcftools stats pour $QUERY"; continue; }

done

# === FIN ===
echo "[✔] Pipeline terminé ! Tous les résultats sont disponibles dans $OUTDIR."
