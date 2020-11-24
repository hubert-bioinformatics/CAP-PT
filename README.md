## CAP-PT
Scripts to help you enter the answer for CAP-PT (Proficiency Testing)

## Table of Contents
* [Intallation](#installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [Credits](#credits)
* [License](#license)


## <a name="installation">Installation</a>
### Requirements
* The server state that Xvfb(X Virtual Frame Buffer) has already been installed is needed to run the script


## <a name="usage">Usage</a>

* Basic
  * command
    * python cap_vcf_parsing_v2.py \
      -v [snpeff annotated vcf] \
      -c [CAP-PT excel file] \
      -r [bam file] \
      -b [target region bed file]
  * input
    * snpeff annotated vcf: a vcf file annotated with snpeff is required for matching a transcript id.
    * CAP-PT excel file: an excel file including transcript id, gene, region, and so on from CAP.
    * bam file: a bam file for making a snapshot for the problem.
    * target region bed file: a bed file for bioinformatics analysis to know whether a region of problem is in the analysis target region or not.
  * output
    * the answer file (tsv format)
    
    <br>
    
    [![the answer file](https://github.com/hubert-bioinformatics/CAP-PT/blob/main/README_images/CAP-PT_parse2.png)](https://github.com/hubert-bioinformatics/CAP-PT/blob/main/README_images/CAP-PT_parse2.png)
    <br>
    
    * IGV snapshot files: 
    
    <br>
    
    [![snapshot image](https://github.com/hubert-bioinformatics/CAP-PT/blob/main/README_images/9_2_29416572_29416573_h2000.png)](https://github.com/hubert-bioinformatics/CAP-PT/blob/main/README_images/9_2_29416572_29416573_h2000.png)
    <br>
    
  * demo play
  
  <br>
  
  [![domo play](https://github.com/hubert-bioinformatics/CAP-PT/blob/main/README_images/cap-pt_parsing.gif)](https://github.com/hubert-bioinformatics/CAP-PT/blob/main/README_images/cap-pt_parsing.gif)
  <br>

## <a name="contributing">Contributing</a>


Welcome all contributions that can be a issue report or a pull request to the repository.


## <a name="credits">Credits</a>


hubert (Jong-Hyuk Kim)


## <a name="license">License</a>

Licensed under the MIT License.

