#!/bin/bash

# === CONFIGURATION ===
REFERENCE="CP010329_H37Rv.fasta"  # Génome de référence
QUERY_DIR=$1         # Dossier contenant les génomes à comparer
OUTDIR="genome_comparison"         # Dossier de sortie
THREADS=8                           # Nombre de threads utilisés

mkdir -p $OUTDIR
cd $OUTDIR

# === ÉTAPE 1 : Préparation des fichiers ===
echo "[1] Indexation du génome de référence..."
samtools faidx $REFERENCE

# === ÉTAPE 2 : Comparaison des génomes avec Mummer ===
echo "[2] Détection des SNPs et INDELs avec Mummer pour chaque fichier de requête..."
for QUERY in $QUERY_DIR/*.fasta; do
    BASENAME=$(basename $QUERY .fasta)
    PREFIX="${BASENAME}_comparison"
    
    echo "[2.1] Traitement de $QUERY..."
    nucmer --prefix=$PREFIX $REFERENCE $QUERY
    show-snps -ClrH $PREFIX.delta > $PREFIX.snps

    # === ÉTAPE 3 : Conversion des SNPs en format VCF ===
    echo "[3] Conversion en VCF pour $QUERY..."
    echo "##fileformat=VCFv4.2" > $PREFIX.vcf
    echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> $PREFIX.vcf
    awk '{print $NF "\t" $1 "\t.\t" $2 "\t" $3 "\t.\tPASS\t."}' $PREFIX.snps >> $PREFIX.vcf

    # === ÉTAPE 4 : Compression et Indexation du VCF ===
    echo "[4] Compression et indexation du VCF pour $QUERY..."
    bgzip $PREFIX.vcf
    tabix -p vcf $PREFIX.vcf.gz

    # === ÉTAPE 5 : Vérification des résultats ===
    echo "[5] Vérification du fichier VCF final pour $QUERY..."
    bcftools stats $PREFIX.vcf.gz

done

# === FIN ===
echo "[✔] Pipeline terminé ! Tous les résultats sont disponibles dans $OUTDIR."
