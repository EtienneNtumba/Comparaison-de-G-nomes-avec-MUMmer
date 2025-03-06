

## 📌 Description

Ce pipeline permet de comparer un génome de référence (H37Rv) avec plusieurs génomes de requête (par exemple, Lx) en utilisant **MUMmer**, **Samtools**, et **BCFtools**. Le pipeline détecte les **SNPs** et **INDELs** entre le génome de référence et les génomes de requête et génère des fichiers **VCF** compressés et indexés.

## 🛠️ Prérequis

Avant d'exécuter ce pipeline, assurez-vous d'avoir installé les outils suivants :

- **MUMmer** (composants `nucmer` et `show-snps`) pour l'alignement et la comparaison des génomes.
- **Samtools** (composant `faidx`) pour la gestion des fichiers fasta.
- **BCFtools** (composants `bcftools`, `bgzip`, `tabix`) pour la gestion des fichiers VCF.
- **Awk** pour la conversion des SNPs en format VCF.
- **Bash** (les scripts sont testés sous Unix/Linux).

## 🔧 Installation des outils

Installez les outils nécessaires en exécutant la commande suivante :

```bash
sudo apt update && sudo apt install -y mummer samtools bcftools tabix


## 📁 Structure des fichiers

Voici l'organisation des fichiers dans le projet :

project/
│── Script_complet3.sh                   # Script principal
│── Reference_genome.fasta               # Génome de référence
│── QUERY_DIR/                           # Dossier des génomes spécifiques
│   ├── L1.fasta
│   ├── L2.fasta
│   ├── ...
│── genome_comparison2/                  # Dossier des résultats générés



🎬 Exécution

Pour exécuter le pipeline, suivez les étapes suivantes :

    Rendez le script principal exécutable :

chmod +x Script_complet3.sh

Lancez le script avec le chemin vers le génome de référence (H37Rv) :

    ./Script_complet3.sh $path_of_ref_genome

📝 Résultats attendus

Après exécution, dans le dossier genome_comparison2/, vous obtiendrez :

    L1_comparison.snps : Liste des SNPs détectés.
    L1_comparison.vcf.gz : Fichier VCF compressé avec les variantes.
    L1_comparison.vcf.gz.tbi : Index du fichier VCF.

Le même résultat sera généré pour d'autres génomes comme L2, L3, etc.
🔒 Étape 2 : Script de Masquage des Variants dans un Fichier VCF

Ce script Python permet de modifier les bases de référence (REF) et alternatives (ALT) dans un fichier VCF en fonction de régions spécifiques définies dans des fichiers BED. Il remplace les bases par "N" si elles se trouvent dans ces régions.
🛠️ Prérequis

Avant d'exécuter ce script, vous devez installer Python 3.x et le module vcfpy.
Installation de vcfpy :

pip install vcfpy

🎬 Exécution

Pour lancer le script, utilisez la commande suivante :

python Mask_vcf_all.py Lx_comparison.vcf.gz Mask_H37Rv.bed Mask_Lx.bed Lx_mask_N.vcf.gz

📁 Structure des fichiers

project/
│── Mask_vcf_all.py                # Script Python pour masquer les variants
│── Lx_comparison.vcf.gz           # Fichier VCF compressé
│── Mask_H37Rv.bed                 # Fichier BED pour le masque REF
│── Mask_Lx.bed                    # Fichier BED pour le masque ALT
│── Lx_mask_N.vcf.gz               # Fichier VCF de sortie avec régions masquées

📊 Étape 3 : Analyse des Variants dans un Fichier VCF Comprimé

Ce script Bash permet de :

    Décompresser un fichier VCF.gz (ex. : Lx_mask_N.vcf.gz).
    Compter les occurrences des variantes selon plusieurs scénarios.
    Générer un fichier CSV résumant les résultats.

🛠️ Prérequis

Avant d'exécuter ce script, assurez-vous d'avoir installé :

    gunzip (inclus dans gzip).
    awk (pour l'extraction des données).
    wc (pour le comptage des lignes).

Installation sous Ubuntu (si nécessaire) :

sudo apt update && sudo apt install gzip gawk

🎬 Exécution

Pour exécuter le script, utilisez la commande suivante :

chmod +x filtre_scenario1.sh
./filtre_scenario1.sh Lx_mask_N.vcf.gz Tableau_scenario.csv

📁 Structure des fichiers

project/
│── filtre_scenario1.sh            # Script Bash pour l'analyse des variants
│── Lx_mask_N.vcf.gz               # Fichier VCF compressé après masquage
│── Tableau_scenario.csv          # Fichier CSV généré avec les résultats des scénarios

📄 Résumé

Ce projet permet de :

    Comparer un génome de référence (H37Rv) avec plusieurs génomes de requête à l'aide de MUMmer, Samtools et BCFtools.
    Masquer certaines régions génomiques dans les fichiers VCF.
    Analyser les variants détectés selon plusieurs scénarios et résumer les résultats dans un fichier CSV.



















