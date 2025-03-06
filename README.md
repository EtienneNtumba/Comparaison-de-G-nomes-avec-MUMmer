

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

