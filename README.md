# Mapping_project

1. API de manipulation de fichiers FASTA et FASTQ

OPTIONS:
--help Print usage and exit.
--fastq, -q Si il s’agit d’un FastQ.
--num, -n Numéro des séquences à utiliser avec bord supérieur non inclus.
Exemple usage : -n1,2
--seq, -s Sous-sequence avec extrémités incluses.
Exemple usage : -s1,2
--rev, -r Inverser la/les séquences.
--cmp, -c Sequence complementaire.
--count, -w Seulement compter les séquences.
--entete, -e Analyse les entetes.
--concat, -t Concatenation des séquences.
--info, -i Recupère toutes les informations sur les séquences.
--pref, -p Supprime préfixe d’une longueur donnée.
--suff, -f Supprime suffixe d’une longueur donnée.
--polya, -f Supprime queue polyA d’une longueur minimale donnée. (Defaut: 4)
--qual, -u Supprime les parties de mauvaise qualité selon une qualité
minimale et une longueur de sequence minimale.
(Defaut: Lg = 15 et Qu = 20)
Exemple usage : -u10,20 avec 10 la longueur et 20 la qualité

2. API de création et d'utilisation d'une table des suffixes 

OPTIONS:
--help Print usage and exit.
--num, -n Numéro des séquences à utiliser avec bord supérieur non inclus.
Défaut : seulement la première.
--concat, -t Concatenation des séquences.
--print, -p Imprime le tableau ordonneé de suffixes.
--prints, -i Imprime la table des suffixes.
--search, -s Cherche une ou plusieurs séquences.
--enc Si la séquence doit etre encodee.
--occ, -o Recuperer le facteur de longueur k qui a le plus d’occurrences
dans la sequence.
Recuperation du 1er facteur seulement si plusieurs.
--dist, -d Recuperer les facteurs distincts de longueur k
--fact, -f Recuperer le ieme facteur de longueur k de la sequence
Exemple : -fi,k
(Defaut : 1,1)
--mot, -m Recuperer le mot de longueur K à la position i
Exemple : -mi,k
(Defaut : 1,1)

3. Algorithme de Mapping 

- Indexation du génome de référence
- Découpage des lectures en k-mers, profils de localisation et support (Recherche exacte de k-mers, Recherche partielle de k-mers, 
Localisation des lectures et détection d’une suite)
