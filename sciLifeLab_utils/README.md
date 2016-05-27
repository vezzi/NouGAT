# NGI specific scripts
`scilifelab_denovo`, automated convenience for NGI users. 

### Requirements
* Fastq sequence data in the folder structure as delivered by [NGI Stockholm](https://portal.scilifelab.se/genomics/)
* Needs to be run on [UPPMAX](http://www.uppmax.uu.se/). 


The commandline tool `scilifelab_denovo` supports customizable parameter defaults for all it's sub-commands. It attempts to read a [configuration file](../config_files/scilifelab.conf) found in the file `$HOME/.nougat/scilifelab.conf`.
```bash
mkdir $HOME/.nougat
cp ../config_files/scilifelab.conf $HOME/.nougat/

# Now edit the template
vim $HOME/.nougat/scilifelab.conf
```
