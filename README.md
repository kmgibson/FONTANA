# FONtANA: FOrensic Next-generation ANAlysis
## Forensics Microhaplotype Pipeline
### For Daniele Podini's Forensic Lab
---
---

<br />

## Beginning - Prior to running FONTANA
Please see [Colonial One's help page](https://colonialone.gwu.edu/getting-access/). For working through the terminal, please complete the [Codecademy's Command Line Tutorial](https://www.codecademy.com/learn/learn-the-command-line). For help related to GW's HPC (Colonial One), please see [our help page](https://gwcbi.github.io/HPC/).


<br />

### Beginning on Colonial One
See Colonial One's [wiki page](https://colonialone.gwu.edu/available-filesystems/). You will run *all* FONTANA related stuff on `lustre` (**NOT** your home directory). You can access `lustre` by such: `cd /lustre/podinigrp`.

You will then want to create a directory for each person individually. Do so by doing `mkdir -p $name`, so for example I would do `mkdir -p kgibson` and now I would have a directory named kgibson in this directory.

#### Starting:

See this [help page](https://colonialone.gwu.edu/quick-start/). Where it shows `$name` insert your colonialone id. Mine is `kmgibson`.
```
ssh $name@login.colonialone.gwu.edu
$enter_password
# here we are logging into Colonial one

cd /lustre/groups/podinigrp # here we are moving into the correct directory
mkdir -p $name/samples # here we are creating a directory for you and creating a samples directory
cd $name # here we are moving into your directory
```
---
## Downloading/Uploading data to Colonial One

1. Download the `BAM` files from the server.
2. Upload them to Colonial One like so:
This is done on your own computer. Not within Colonial One.
```
rsync -P -av "MP" kmgibson@login3.colonialone.gwu.edu:/lustre/groups/podinigrp/kgibson/samples

```
This is saying upload directory "MP" to Colonial One at this location. You will need to change the directory name to whatever yours is and change kgibson to whatever your Colonial One ID is. So fill in $name and $directory in the one below:
```
rsync -P -av "$directory" $name@login3.colonialone.gwu.edu:/lustre/groups/podinigrp/$name/samples
```

If you have done this successfully, your files on your computer in `$directory` will be uploaded to Colonial One under the directory `samples`.






3. Getting a list of all samples, renaming the files and creating a directory for each sample.
```
cd samples &&

for f in *.bam; do 
    name=$(echo $f |awk -F[_.] '{print $2}')
    samp=$(samtools view -H $f | grep "SM:" | cut -f10 | sort | uniq | awk -F[:.] '{print $2}')
    echo -e "${name}\t${samp}" >> dirsamp.txt # this gets a list of all samples
    echo $f
done

cat dirsamp.txt # view this list

# this loop creates directories and renames the file
for f in *.bam; do
    name=$(echo $f |awk -F[_.] '{print $2}')
    samp=$(samtools view -H $f | grep "SM:" | cut -f10 | sort | uniq | awk -F[:.] '{print $2}' | cut -d " " -f1)
    echo -e "$f\t$name\t$samp"
    mkdir -p Ion${name}_${samp}/00_raw
    mv $f Ion${name}_${samp}/00_raw
    mv Ion${name}_${samp}/00_raw/$f Ion${name}_${samp}/00_raw/original.bam
done

ls -d Ion* > list.txt # creates a list of just the directory names. This file is used in FONtANA.
```

---

## Starting Fontana
You need to allocate a node for yourself, because I have not made this automatic (into a job) yet. See this [help page](https://github.com/gwcbi/HPC/blob/master/interactive_jobs.md) for getting an interactive node.

Basic commands are below. You can replace the `defq` option with either `short` or `debug`.
```
sinfo    # look for idle nodes.
salloc -p defq -N 1 -t 500    # request a node
squeue   # find the node you got
ssh node$num     # insert number where $num is, the number is on the allocated node you obtained.

# You need to move back to the correct directory
cd /lustre/podinigrp/$name/
```

Loading an environment for FONtANA.
```
module use /groups/podinigrp/shared/modulefiles
module unload python
module load miniconda3/4.3.27.1
source activate fontana
```

A dry run of FONtANA, to see what the program is going to do.
```
snakemake Snakemake
```

Running FONtANA:
```
snakemake Snakemake
```



---

/groups/podinigrp/shared/apps/miniconda3/4.3.27.1
<br />

---
---
---

#### For Future:
Create documentation with https://readthedocs.org/
