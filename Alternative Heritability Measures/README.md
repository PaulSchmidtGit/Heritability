# Alternative Estimation Methods for H² on Entry-Mean Basis
There are several alternative estimation methods for (broad-sense) heritability on an entry-mean basis. Please see our paper for more information:

P. Schmidt, J. Möhring, J. Rath and H.-P. Piepho. 2018. Estimating broad-sense heritability with unbalanced data from agricultural cultivar trials. Crop Science **forthcoming**

### Based on 3 different mixed model functions
* `asreml()` of the R-package [ASReml-R](https://www.vsni.co.uk/software/asreml-r/) Version 3.0
* `mmer2()`  of the R-package [sommer](https://cran.r-project.org/web/packages/sommer/index.html)
* `PROC MIXED` in SAS

### Application
**For both R packages** you will find example analyses that you can copy-paste into R and run immediately (given all required packages are installed).

**For SAS** it works similarly, yet you will see that the example analyses make use of SAS `%MACROs`. These macros are also provided in the SAS folders. You do not need to copy-paste or download the macros in order to run the example analyses, since they are included automatically via their URL and a `proc http` command at the top of each code.

### Work in progress
Note that at the moment we do not provide code for all methods described in our paper and/or for all of the three mixed model functions. You can find an overview here:

H² Method | `asreml()` | `mmer2()` | `PROC MIXED` | 
:--- | :---: | :---: | :---: |
Standard |  |  |  |
Holland |  |  |  |
Piepho | X | °° | X |
Cullis | X | X | X |
Oakey | X | X | X |
BLUP~BLUE | X | °° | X |
sumdiv | X | °° | X |
Simulated | ° | X | X |

° *in asreml-R version 3.0 it is not possible to extract the required estimated matrices*  <br />
°° *not yet possible to obtain BLUEs in sommer package*
