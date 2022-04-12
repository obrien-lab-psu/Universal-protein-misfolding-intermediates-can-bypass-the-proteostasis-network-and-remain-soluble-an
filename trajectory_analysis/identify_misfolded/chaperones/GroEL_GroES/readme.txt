DATA SOURCE
------------
Kerner is the original report, which identified 252 substrate proteins that interact with GroEL.
Fujiwara is a follow-up that found 4 additional substrate proteins (though their main goal was classifying the proteins identified by Kerner).
Niwa identified an additional 20 substrate proteins based on a more refined assay able to detect smaller concentrations of proteins.
Together, these 252+4+20 = 276 proteins represent the current knowledge concerning protein interactions with GroEl in E. coli

NAME SCHEME
------------
The Kerner data set is in terms of the uniprot entry name already, so no modification is needed. 
The 4 additional proteins from Fujiwara were in terms of the gene name, so I converted them to uniprot entry name as necessary.
The Niwa data set was in terms of the gene name so these were also converted to the uniprot entry name as needed.
Proteins from my data set (myproteins.dat) are listed in terms of their uniprot entry name.
One protein, uniprot name yjjv, is in both the Niwa and Kerner data sets despute Niwa claiming their 20 proteins are novel substrates.

The file "all_substrate_proteins.dat" contains the uniprot entry names of the combined lists from Fujiwara, Niwa, Kerner.

OUTPUT
------------
The file myproteins_that_are_substrates.dat are the interesction of the lists myprotein.dat and all_substrate_proteins.dat 
