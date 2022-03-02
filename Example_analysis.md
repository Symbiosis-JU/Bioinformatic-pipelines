# Introduction to Symbiosis Evolution Group bioinformatics pipeline for amplicon sequencing data analysis
We hope that this script will help you navigate through the analyzes of an example amplicon dataset - a dozen of libraries from our Greenland project.

## Contents ---
1. Some basic Linux commands.
2. Accessing example data
3. Workflow overview



##1. Before we start let's get familiar with some Linux commands!
- `pwd` --- where are you? (prints the PATH to your current position).
- `ls` --- listing directory contents.
  - `ls -l` --- lists directory contents while displaying their characteristics  
- `cd` --- changing directories:
  - `cd Workshop` --- change working directory to the directory "Workshop" that is in the current directory
  - `cd ..` --- will move you one level up in your directories tree 
  - `cd` --- by default, typing just 'cd' will take you to your home directory
- `cp` --- copying file, need to be followed by item you want to copy and a path to the directory:
  - `cp /path/to/file.txt .` --- copies a remote "file.txt" to your current working directory, symbolized by a dot (".")
  - `cp Sequences.fasta Workshop/` --- copies a file "Sequences.fasta" from the current working directory to subdirectory "Workshop"
  - `cp -r /path/to/directory ~` --- copy recursively the whole directory/file structure from a remote location to your home directory ("~")
- `mkdir` --- make new directory
  - `mkdir MyNewDirectory` --- will create a directory with the requested name
- `mv` --- move, need to be followed by the name of the item you want to move and a path to the destination directory. Can also be used to rename items.
  - `mv file.txt Workshop/` --- will move a file to the "Workshop" subdirectory
  - `mv OldName.txt NewName.txt` --- will rename your text document
- `rm` --- delete item
  - `rm FileNotNeededAnymore.txt`
  - `rm *.fasta` --- will remove all files with extension "fasta". Asterisk ("*") 

### ...and some additional Linux tools! 
- `screen` --- very usefull tool that 'can be used to multiplexes a physical console between several processes' by creating virtual sessions that you can connect to or disconnect from, as desired. Generally, you want to run your proccesses within a screen to ensure processes are not disrupted by network connection issues! ):
  - `screen -S MyNewSession` --- will create a new session with the requested name. Make it informative!
  - `screen -r MyRunningSession` --- will re-attach you to your session.
  - `screen -ls` --- will list all the session.
  - `ctr + a + d` --- will detach you from a session, without killing it
- `htop` --- starts a 
- `gunzip` --- will uncompress your gzip-compressed files (.fastq.gz ---> .fastq)

**Lets use some of those beautiful commends in action!**



##2. Copying example data to your folder.
- First, log in to your account on *azor* cluster.
- Then, copy the prepared sample data to the directory of your choosing (we recommend using your home directory):
```
cp -r /mnt/matrix/symbio/workshop_march_2022 ~/
```
- Now you have folder "workshop_march_2022" containing R1.fastq and R2.fastq files for each sample.
**Those are uncompressed files, remember that when you obtaine your files from sequencing facility they will be with compressed using gzip (extension: ".fastq.gz") and will need to be decompressed prior to use!
- Go to the copied directory, display the contents:
```
cd ~/workshop_march_2022
ls
ls -l
```



## Running splitting (MultiPISS) script.

First we want to split our libraries into bins representing our target genes and from each library delete sequences with tags uncharacteristic to its well.
To do that we are going to use our [MultiPISS.py](https://github.com/Symbiosis-JU/Bioinformatic-pipelines/blob/main/multiPISS.py) script:
1. Click into link above and copy it (by clicking the icon 'copy raw contents' in the right upper corner of the box).
2. use command:
``` nano MultiPISS.py```.
This will create an empty file named MultiPISS.py. Paste the script and exit file by ctr + x with saving the changes.
3. Make this script executable by using command: ```chmod +x MultiPISS.py```.

**Now you can run script!**
Type ```./MultiPISS.py``` (```./``` --- indicates that this file is in the current directory).

Oh no! You got an ERROR!:


```
ERROR! CHECK YOUR INPUT PARAMETERS!
Please provide:
1) sample list with information about your libraries created in following manner (tab-separated):
Sample_name Sample_name_R1.fastq	Sample_name_R2.fastq
Please remember to first un-gzip your .gz files!!!
2) FULL path to the directory with R1 and R2 fiels for all the amplicon libraries that you want to analyse e.g.:
/home/Data/For/Nature/Publication/)
shortcuts such as "./" are unlikely to work
3) output directory path.
4) Information whether the last two characters of your sample name indicate well number
1=True, 0=False
If you claim 1 but the last two characters are not numbers, it may create an error!
5) Number of cores to use
```

As you can see, we need to provide some input parameters for our script.

#### Sample_list:
To create a sample_list in a manner that script requires, just type following command in the directory with your R1 and R2 files:
```
for file in *_R1.fastq; do
    SampleName=`basename $file _R1.fastq `
    SampleNameMod=$(echo "$SampleName" | sed 's/-/_/g' | sed 's/_S[0-9]\+$//g')
    echo $SampleNameMod "$SampleName"_R1.fastq "$SampleName"_R2.fastq >> sample_list.txt
done
```
**Now we are ready to run this script for real!**
Use following command:
```
./MultiPISS.py sample_list.txt ~/workshop_march_2022 ~/workshop_march_2022/splitted 0 20
```
This command will run our script using created sample_list, output it without visualization in 'splitted' subdirectory using 20 cores.

Type: 
```
cd ~/workshop_march_2022/splitted
```

Now you can see that you have four new subdirectories: **COI_trimmed  incorrect_untrimmed  V12_trimmed  V4_trimmed**
- COI_trimmed --- containes mitochondrial COI reads with trimmed primers sequences,
- V12_trimmed --- containes bacterial 16S V1-V2 reads with trimmed primers sequences,
- V4_trimmed --- containes bacterial 16S V4 reads with trimmed primers sequences,
- incorrect_untrimmed --- with sequences that were unrecognize by the script with primers still attached.

**Lets start with COI analysis!**

## LSD script
This script joins F and R (R1 and R2) reads, passes only those of high quality. 
Next it converts fastq to fasta file, dereplicate and denoise sequences in each library seperately.
Joins all the libraries into one table and assigns all the sequences to taxonomy.
**This is a first step of analysis of bacterial 16S and COI data!**

Again, go to our repository and copy [LSD.py](https://github.com/Symbiosis-JU/Bioinformatic-pipelines/blob/main/LSD.py).
Copy it.
Use ```nano LSD.py``` to create an empty file in the directory tou have your marker gene reads (for COI: ~/workshop_march_2022/splitted/COI_trimmed).
Close it with saving changes and make this file executable with ```chmod +x LSD.py```.

Again if you will type just ```./LSD.py``` you will get an error about putting some parameters. In this case:
```
ERROR! CHECK YOUR INPUT PARAMETERS!
Please provide:
1) sample list with information about your libraries created in following manner:
Sample_name Sample_name_R2.fastq Sample_name_R2.fastq
Please remember to first un-gzip your .gz files!!!
2) path to the directory with R1 and R2 fiels for all the amplicon libraries that you want to analyse e.g.:
/home/Data/For/Nature/Publication/)
3) type of data, please indicate if you are going to analyse 16SV4, 16SV1-V2 (Bacterial) or COI.""")
```
So, again we need a sample list. But this time as an input we are using an output of MultiPISS.
Now our imput files for COI analysis looks like this:
```
GRE2290_F_COI.fastq
GRE2290_R_COI.fastq
```
We will use simmilar (**but not the same**) commend as before:
```
for file in *_F_COI.fastq; do
    SampleName=`basename $file _F_COI.fastq `
    SampleNameMod=$(echo "$SampleName" | sed 's/-/_/g' | sed 's/_S[0-9]\+$//g')
    echo $SampleNameMod "$SampleName"_F_COI.fastq "$SampleName"_R_COI.fastq >> sample_list_COI.txt
done
```
**Running LSD.py:**
As this part may take a while **remember to use screen session!**:
```
screen -S LSD
```
and then:
```
./LSD.py sample_list_COI.txt ~/workshop_march_2022/splitted/COI_trimmed COI
```
Now you have **A LOT** of files. But for now the most important for you is:
```
zotu_table_expanded.txt
```
This is a table with all the COI sequences in all your libraries assigned to a various taxonomic levels.
You can use it for your analysis
OR
You can use it as an imput for [MAO.py](https://github.com/Symbiosis-JU/Bioinformatic-pipelines/blob/main/MAO.py) script!

## MAO script
This script is simple, but briliant at the same time.
It uses ```zotu_table_expanded.txt``` of COI data as an in input and produces:
- barcode.txt that contains info about most abundand COI barcode, taxonomy and bacteria presence
- barcode.fasta, containing most abundand Eucaryotic COI barcode per each library/sample
- euc5.txt, containing information about Eucaryotic COI sequences that represents at least 5% of total Eucaryotic reads per sample
- bac.txt, containing information about bacterial COI sequences

Just do our trick with creating an ampty file with ```nano MAO.py```, paste the script and close the file with saving.
Make script executable with ```chmod +x MAO.py``` and run it with ```./MAO.py zotu_table_expanded.txt```.

