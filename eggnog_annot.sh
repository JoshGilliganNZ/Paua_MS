trinity_v3_nop1.Trinity.fasta.transdecoder_mod.pep##make new env
##activate env
##install
#mamba install bioconda::eggnog-mapper
##
#modify the line below, then run
emapper.py --cpu 64 --sensmode more-sensitive -i CHANGETOMOD.PEP  -o SMART_DESCRIPTIVE_NAME_HERE
