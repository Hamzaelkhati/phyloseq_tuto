DADA2 Tuto
================

  - [Partie1: importation des données
    MiSeq\_SOP](#partie1-importation-des-données-miseq_sop)
  - [Partie2: Filtrer et découper](#partie2-filtrer-et-découper)
      - [Obtention des listes des fichiers fastq sens et anti-sens (fnFs
        et
        FnRs)](#obtention-des-listes-des-fichiers-fastq-sens-et-anti-sens-fnfs-et-fnrs)
      - [plots de qualité pour fnFs et
        fnRs](#plots-de-qualité-pour-fnfs-et-fnrs)
      - [modèle d’erreur](#modèle-derreur)
  - [partie3: Construire une table de séquence et supprimer des
    chimères](#partie3-construire-une-table-de-séquence-et-supprimer-des-chimères)
      - [Construire une table de
        séquence](#construire-une-table-de-séquence)
      - [supprimer les chimères](#supprimer-les-chimères)
      - [Suivre les lectures dans le
        pipeline](#suivre-les-lectures-dans-le-pipeline)
  - [Partie4: Attribuer une taxonomie](#partie4-attribuer-une-taxonomie)
  - [Partie5: Évaluer la précision](#partie5-évaluer-la-précision)

``` r
library(rmarkdown)
library(knitr)
```

\#ce tutoriel nécessite le packages dada2 et ggplot2 déja télécharger
dans le fichier 00\_install\_packages

``` r
library(dada2); packageVersion("dada2")
```

    ## Loading required package: Rcpp

    ## [1] '1.18.0'

``` r
library(ggplot2)
```

\#Les données utilsés ici sont des séquences d’amplicon Illumina Miseq
2x250 très chevauchantes de la région V4 du gène 16S, les 360
échantillons de matiere fécale sont prélevé de 12 souris, Ces données
sont téléchargées à partir de l’emplacement de données suivant et
décompressées( 01\_data\_import
(<https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip>))
sous le nom de MiSeq\_SOP

# Partie1: importation des données MiSeq\_SOP

\#l’emplacement de mes fichiés correspond au
\~/phyloseq\_tuto/MiSeq\_SOP

``` r
path <- "~/phyloseq_tuto/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

# Partie2: Filtrer et découper

## Obtention des listes des fichiers fastq sens et anti-sens (fnFs et FnRs)

\#La premier et 2ème ligne de code sert a ordonner les fichiés fastq
sens et anti-sens dans le meme ordre \_R1\_001.fastq et \_R2\_001.fastq
avec path comme chemin de spécifique, et la 2ème c’est pour faire sortir
des exemples de noms de fichies fastq en supposant que le nom est
SAMPLENAME\_XXX.fastq

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

## plots de qualité pour fnFs et fnRs

\#le plot de qualité montre les qualités de séquence sens est antisens
dont la moyenne est en vert, la médiane la ligne orange pleine et les
quartiles sont les lignes orange pointillées, ce code est pour affiché
le score de qualité pour les 2 premieres échantillons

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\#\# Filtration des données

\#noms de fichiers pour les fichiers fastq.gz filtrés,
\_F\_filt.fastq.gz pour les forwards, \_R\_filt.fastq.gz pour les revers

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

\#ici on a tronquer les lectures sens à la position 245 car les lectures
sens(forwards) maintiennent une qualité élevée tout au long avec une
dimunition de score de qualités , et les lectures anti-sens (Revers)sont
tronqués à la position 160 car le score de qualité diminue aprés cette
position (TruncLen=c()). Nous choisissons également de couper les 10
premiers nucléotides de chaque lecture car ces positions de base sont
particulièrement susceptibles de contenir des erreurs pathologiques.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

## modèle d’erreur

\#l’estimation des taux d’erreur des échantillon ce fait par un modèle
d’erreur paramétrique ( err), tout d’abord l’algorithme doit commencer
par une estimation initiale pour laquelle les taux d’erreur maximums
possibles dans ces données sont utilisés, par exemple pour errF 33514080
de bases dans 139642 reads de 20 échantillons vont etre utilisé pour
apprenez le taux d’erreur, par la suite il est utile de visualiser les
taux d’erreur estimés par la fonction plotErrors pour errF et errR.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

\#Les erreurs de substitutions, quand les clusters sont trop proches,
peuvent être dues des erreurs de lecture. Les points gris sont les
erreurs de substitutions. Pour les substitutions (A donne C, deuxième
carré de graphique), dada2 pense que ce sont des erreurs car il a
comparé à une autre séquence proche. Chaque point donne la probabilité
d’avoir une substitutions de A donne C en fonction du Q score de la
position. La courbe en noir représente la régression. On voit que le
modèle suit les points. Plus le Q score est élevé plus le taux de
substition est faible. On va appliquer ce modèle de substitution à tout
notre jeu de données afin de le corriger ainsi que de corriger les
erreurs. On peut ensuite considérer chaque séquence.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
\# partie3: Inférence d’échantillon

\#l’algorithme d’inférence de l’échantillon de base aux données de
séquence filtrées et découpées. dadaFs et dadaRs sont issus de la
fonction dada. Ils prennent les reads Foward et Reverse filtrés sur
lesquelles on a mis des filtres mais aussi les erreurs Foward ou
Reverse. On obtnient normalement des reads corrigé.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

\#Inspection de l’ dada-class qui porte des diagnostics et des
informations sur la qualité de chaque variante de séquence débruitée,
exemple de 3ème echantillon. \#’algorithme DADA2 a déduit 128 variantes
de séquence à partir des séquences uniques de 1979.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# partie3: Construire une table de séquence et supprimer des chimères

\#ici on va crée la «table OTU» commune, c’est-à-dire une table de
caractéristiques échantillon par séquence valorisée par le nombre de
fois que chaque séquence a été observée dans chaque échantillon.
\#Merger permet de refusionner les séquences Foward et Reverse pour
reconstituer le V4.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6540 paired-reads (in 107 unique pairings) successfully merged out of 6891 (in 197 pairings) input.

    ## 5028 paired-reads (in 101 unique pairings) successfully merged out of 5190 (in 157 pairings) input.

    ## 4986 paired-reads (in 81 unique pairings) successfully merged out of 5267 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2754 (in 108 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3646 paired-reads (in 55 unique pairings) successfully merged out of 4109 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6514 (in 198 pairings) input.

    ## 3968 paired-reads (in 91 unique pairings) successfully merged out of 4388 (in 187 pairings) input.

    ## 14233 paired-reads (in 143 unique pairings) successfully merged out of 15355 (in 352 pairings) input.

    ## 10528 paired-reads (in 120 unique pairings) successfully merged out of 11165 (in 278 pairings) input.

    ## 11154 paired-reads (in 137 unique pairings) successfully merged out of 11797 (in 298 pairings) input.

    ## 4349 paired-reads (in 85 unique pairings) successfully merged out of 4802 (in 179 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7191 (in 187 pairings) input.

    ## 4426 paired-reads (in 67 unique pairings) successfully merged out of 4603 (in 127 pairings) input.

    ## 4576 paired-reads (in 101 unique pairings) successfully merged out of 4739 (in 174 pairings) input.

    ## 6092 paired-reads (in 109 unique pairings) successfully merged out of 6315 (in 173 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

\#Le tableau séquences montre les séquences combinées (Foward et
Reverse). Dans le second tableau, la colonne nmatch montre le nombre de
nucléotide qui match dans la régions de l’overlaps, et nindel vaut 0 car
il n’y a pas d’indels dans Illumina. Prefer montre quelle séquence est
prise entre R1 et R2 s’il y a des mismatsh au niveau des chevauchements.
Pour choisir, l’ordinateur va regarder le meilleur Q score entre R1 et
R2.

## Construire une table de séquence

\#A partir de la fonction mergers, on peut construire une table de
séquence. La dimension de la table c’est 20 colonnes et 293 lignes. la
fonction MakeSequenceTable construit une table analogue à une table
d’OTU. C’est une matrice d’observation échantillon/séquence. \#létape
2 correspond a l’alignement du R1 avec le R2, on peut peut faire une
vérification pour voir si on obtient bien des séquences aux alentours
de 250. On a les nombres de séquences que l’on retrouve dans chaque
classes (1 dans 251, 88 dans 252, etc).

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

## supprimer les chimères

\#Ici on va détecter si le début d’une séquence appartient à une
bactérie et la fin à une autre car ce qui est sûr c’est que le début de
la séquence est dans notre jeu de données et la fin aussi. Le ratio
entre le nombre de séquences sans les chimères divisés par le nombre
avec les chimères est de 96%. Donc ce sont des occurences rares (les
chimères).

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

``` r
1-sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.03596257

## Suivre les lectures dans le pipeline

\#ici on va suivre l’évolution du nombre de read à chaque étape. Et cela
permet de vérifier qu’il n’y a pas eu trop de perte. Cette dernière peut
être du à la filtration, le dénoising, quand on a mergé et on a aussi
enlevé des chimères.

\#les resultats pour l’éch.1: On passe de 7793 à 6528, pur l’éch.2: de
5299 à 5017

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6976      6979   6540    6528
    ## F3D1    5869     5299      5227      5239   5028    5017
    ## F3D141  5958     5463      5331      5357   4986    4863
    ## F3D142  3183     2914      2799      2830   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4151      4228   3646    3507

# Partie4: Attribuer une taxonomie

\#Pour faire l’assignation taxonomique, il faut une base de données de
référence. Pour ce faire on utilise un algorithme d’assignation
taxonomique, on va utiliser Silva 132 (01\_data\_import),la méthode de
Blast qui va etre utiliser et le résultat d’une assignation taxonomique
c’est une classification hétérogène avec des séquences qui vont
s’arrêter au genre, à la classe, à la famille ou même au phylum.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/phyloseq_tuto/tax/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

\#Avec un autre jeu de données pour l’assignation taxonomique, on peut
essayer de voir s’il y a 100% de similarité avec d’autres séquences pour
avoir des résultats jusqu’à l’espèce. Avec la commande précédente, on ne
va que jusqu’au genre. On va donc ajouter les espèces dans le tableau
taxa (01\_data\_import).

``` r
taxa <- addSpecies(taxa, "~/phyloseq_tuto/tax/silva_species_assignment_v138.fa.gz")
```

\#taxa.print prend la table taxa et enlève les noms de lignes de la
table taxa.print par ce que ce sont les séquences.

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

# Partie5: Évaluer la précision

\#Dans le jeu de données, il y avait une Mock community (souches connues
dans des proportions connues). C’est le seul éléments de données qui
permet d’avoir un recul sur le pipeline.

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

\#Ici On va chercher toutes les séquences qui vont être notées Mock. On
va ensuite les trier par abondances et regarder celles qui ne sont pas
dans les ASV. Les unqs.mock qui sont supérieurs à 0 vont etre enlever

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

\#A travers ce code on va chercher les vraies séquences qui sont dans le
fichier HMP\_MOCK.v35.fasta et on va les mettre dans l’objet Mock.ref.
le résultat montre qu’on a 20 ASV différents dans la communauté mock et
20 sont des matchs exacts des séquences de références. On peut donc dire
que DADA2 a bien fonctionné pour ces séquences là. On a pu tester le
programme de dada2 sur les corrections.

``` r
save.image(file = "02_data-analysis-with-DADA2_FinalENV")
```
