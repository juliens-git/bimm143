class13\_markdown
================

R Markdown
----------

Load bio3d library

``` r
library(bio3d)
```

Store 1hsg as file.name

``` r
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

Read pdb file and store as hiv

``` r
hiv <- read.pdb(file.name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Question 1: The two non-protein resid values are HOH (127), MK1 (1). Resid corresponds to the amino acids and other things (water and ligand molecules) in the protein structure. The full listing in the protein sequence below.

Isolate protein and ligand and create files

``` r
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

Question 2: In ADT (Auto Dock Tools) a little hard to see docking site so add H atoms, H is usually omitted because

Question 3: The charges make sense based upon biochemical amino acid interactions

Processing the all.pdbqt to a PDB format file that can be loaded into VMD

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

Question 4: The docks look in the right region and seem the fit well, however there are many different configurations. The first file configuration looks very clean.

Calculate the RMSD between docking results and crystal strucutre

``` r
# res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

Question 5: They look good except only the first configuration of crystal binding is within 1A RMSD for all atoms.

Question 6: To determine the RMSD for heavy atoms only, I would create new variables that are trimmed to only contain non H atoms. Then I would run the RMSD.

``` r
res.noh <- trim.pdb(res, "noh")
ori.noh <- trim.pdb(ori, "noh")
rmsd(ori.noh,res.noh)
```

    ##  [1]  0.458 11.021 10.374  4.301 10.891  3.717  5.764  3.791  5.498 10.759
    ## [11]  4.224  6.308 10.889  8.776
