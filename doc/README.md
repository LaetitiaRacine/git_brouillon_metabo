# README : Git_Metabo_Analysis
# Analyse bioinformatique du projet METABO multiomic

********************************************************************************
********************************************************************************

## I. Architecture du dossier 

L'arborescence du dossier de travail est la suivante :   
- dossier bin : contient les scripts   
- dossier data : contient les données brutes  
- dossier html : contient les rapports html issus de chaque script  
- dossier exp : contient les fichiers de sortie des scripts (graphiques, tableaux)  
- dossier src : contient les codes sources permettant de lancer tous les scripts  
  
Le dossier data contient des sous-dossiers selon le type de données (scATACseq, scRNAseq ou bulk ATACseq).   
Les fichiers dans le dossier data ont été renommés pour plus de clarté (ajout du nom de la drogue dans chaque fichier).   
  
Le dossier exp contient des sous-dossiers du nom du script auquel ils sont associés.  
Dans chaque sous-dossier, on trouve d'autres dossiers nommés avec la date du jour de leur création.
  
Les fichiers html sont datés du jour de leur création. 


## II. Description des donnees brutes

Ce dossier data contient des données multi-omiques dans différentes conditions de stress métabolique.   
Les données seront analysées séparement dans un premier temps puis confrontées lorsque cela est possible.  
  
### scATAC-seq 10X

Manipulation du 25 octobre 2021  
**4 conditions** : CTRL, DON, 2DG, AOA  
**Timing** : 24h  
Fichiers au format .fastq pré-traités par la plateforme de bioinformique Imagine.   
On démarre l'analyse à partir des fichiers de sortie de Cell Ranger.     
L'analyse s'effectue sous R.  
Initialement, les données ont été récupérées dans des dossiers issus de Cell Ranger et nommés full_cellranger_outs/ATAC_NomDrogue. Pour plus de simplicité, on a copié dans le dossier data seulement les fichiers utiles pour les futures analyses. De plus, le nom de chaque drogue a été ajouté au début du nom du fichier. Enfin, les fichiers barcodes.tsv, matrix.mtx et peaks.bed qui étaient rangés dans un sous-dossier filtered_peak_bc_matrix ont été sortis du sous-dossier et on a ajouté le sous-dossier dans le nom de chaque fichier individuel.   
  
  
### scRNA-seq 10X couplé à CITE-seq

Manipulations du 27 mai 2021 et du 04 juillet 2022  
**10 conditions** : CTRL x2, DON, 2DG, AOA, CTRLaK, DONaK, 2DGaK, AOAaK, VPA    
**Timing** : 96h    
Fichiers au format .fastq pré-traités par la plateforme de bioinformique de Cochin.     
On démarre l'analyse à partir des fichiers de sortie de Cell Ranger.    
L'analyse s'effectue sous R.    
  
  
### bulk ATAC-seq méthode FAST-ATAC

Manipulations du 23 juillet 2019 et du 29 octobre 2019    
**7 conditions** : Xvivo, CTRL, DON, 2DG, AOA, CTRLaK, VPA  
**Timings**: 00h, 03h, 06h, 12h, 24h    
Fichiers au format .fastq pré-traités par la plateforme du CEA.    
On démarre l'analyse à partir des fichiers .bam.   
Une partie de l'analyse s'effectue en utilisant Snakefile et la suite sous R.      
  
  
## III. Description des analyses 




