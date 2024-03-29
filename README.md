# CHAPAO
**CHAPAO** (Compressing  Alignments  using  Hierarchical  and  Probabilistic  Approach) is a novel reference-based technique for compressing MSA files. The name of this method is inspired from a bengali word **“চাপাও”**, which means “compress” or “squeeze”, and is pronounced the same way as CHAPAO (/chɑ:paʊ/). This is to our knowledge the first application of the reference-based technique for compressing MSAs. Unlike conventional reference-based methods where an “extra” sequence (not included in the input sequences to compress) is used as the reference, we used a novel hierarchical referencing technique where a suitable subset of the input sequences in the MSA file is used as reference sequences. CHAPAO offers a substantial improvement in compression gain over the existing best alternate methods for both general purpose compression algorithms (zip, Bzip2, gzip) and special purpose compression algorithm (MFCompress).

Details in [CHAPAO: Likelihood and hierarchical reference-based representation of biomolecular sequences and applications to compressing multiple sequence alignments](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0265360)

# Dependencies 
Python 3.0 or later

# File Descripttion
| First Name    | Uses          |
| ------------- | ------------- |
| compress.py   | Python script for compression  |
| decompress.py |  Python script for decompression  |
| compress.exe  | Windows executable for compression  |
| decompress.exe  | Windows executable for decompression  |

<sup>Windows sequrity system might consider those executables as malwares. If there is any problem in downloading the files, please follow the these alternative links for [python script](https://drive.google.com/drive/folders/1P4cwftp9Pr7oZ87wGLmu_cxTeTDQDVgV?usp=sharing) and [executables](https://drive.google.com/drive/folders/1A10qPlQ4lXjHEr-zKdzGCZiuguB8mUDh?usp=sharing) </sup>
# Usage 
### Input File Format
Input files should be in FASTA or PHYLIP format.

### Compressing Single File

**Using Python Script**
```bash
python3 compress.py PATH_OF_THE_FILE WINDOW_SIZE OVERLAP_AMOUNT
```
Example:
```bash
python3 compress.py /home/Desktop/DATA/avian/chr1_96_s.fasta 30 28
```

**Using EXE**
```bash
compress.exe PATH_OF_THE_FILE WINDOW_SIZE OVERLAP_AMOUNT
```
Example:
```bash
compress.exe /home/Desktop/DATA/avian/chr1_96_s.fasta 30 28
```

This will create a folder with extention **.chapao** which is the output of the compression algorithm. There will be two separate files in the folder, namely metadata.txt and ref.txt.

**WINDOW_SIZE** and **OVERLAP_AMOUNT** are two important hyperparameters in our algorithm.  With the increase of **WINDOW_SIZE**, the compression gain generally increases at the cost of more compression time. **OVERLAP_AMOUNT** should be less than **WINDOW_SIZE** and its effect is similar to  **WINDOW_SIZE**.

### Guidelines for hyperparameter selection
For smaller files (<100MB) **WINDOW_SIZE** of 40-50 and **OVERLAP_AMOUNT** of 35-48 will give a high compression ratio within reasonable amounts of time.
For larger files (>100MB) **WINDOW_SIZE** of 5-10 and **OVERLAP_AMOUNT** of 3-8 should be used to compress files within reasonable amounts of time.



### Decompressing:

**Using Python Script**

```bash
python3 decompress.py PATH_OF_COMPRESSED_FOLDER OUTPUT_FORMAT(f/p)
```

OUTPUT_FORMAT = 'f' - to produce decompressed files in **FASTA** format

OUTPUT_FORMAT = 'p' - to produce decompressed files in **PHYLIP** format

Example:
```bash
python3 decompress.py /home/Desktop/DATA/avian/chr1_96_s.fasta.chapao/ f
```

This will create the decompressed file named "Decompressed.txt" inside the same folder.

**Using EXE**

```bash
decompress.exe PATH_OF_COMPRESSED_FOLDER OUTPUT_FORMAT(f/p)
```

Here 'f' and 'p' have the same meaning mentioned above.

Example:
```bash
decompress.exe /home/Desktop/DATA/avian/chr1_96_s.fasta.chapao/ f
```
This will also create the decompressed file named "Decompressed.txt" inside the same folder.

# Sample Dataset
Links of several MSA datasets, which are used in the development and testing of this technique, can be found [here](https://drive.google.com/file/d/1pIBNOJEmtN-sjwPjj_bLs8byx3dcEkek/view?usp=sharing).
### Contributor
**Abdullah Aman Tutul**

**Sifat Muhammad Abdullah**
## License
[MIT](https://choosealicense.com/licenses/mit/)

## Bug reports: ashiqbuet14@gmail.com
