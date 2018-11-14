## Alternative Estimation Methods for H² on Entry-Mean Basis in SAS
There are several alternative estimation methods for (broad-sense) heritability on an entry-mean basis. Please see our paper for more information:

> P. Schmidt, J. Hartung, J. Rath and H.-P. Piepho. 2018. Estimating broad-sense heritability with unbalanced data from agricultural cultivar trials. Crop Science **forthcoming**

### General Application
For each method, you will find
* a `%MACRO` and
* an example with where this `%MACRO` is applied to a simple dataset. 

You can simply copy-paste the example code into SAS and run it (given you are connected to the internet). This works, because the respective `%MACRO` is run directly from this github page via a `proc http` command at the beginning of each example. Everything else (i.e. dataset, modelling procedure etc.) is provided in the example code.

#### Note on additional, supporting `%MACRO`s
Due to the limitations of SAS' Output Delivery System, some of the estimated tables/matrices you can obtain from `PROC MIXED` need further processing in order to obtain the desired estimate. Specifically the estimated variance-covariance matrices of the random (genotype) effects obtained via `G=` and the estimated variance-covariance matrix of the genotype BLUPs contained in the `MMEqSol=` need to be extracted before being used in `PROC IML`. Since these are used in more than one of the alternative H2 methods, we decided to create separate `%MACRO`s whose only purpose is the preprocessing/formatting of the `G=` and `MMEqSol=` datasets. They are called `%getC22g`, and `%getGFD`. Finally, there is a third `%MACRO` named `%getGamma`, which calculates the Gamma matrix [see Eq. 13 in Piepho & Möhring (2007)] required for H2 Sim. 
These additional `%MACRO`s are all contained in the SAS-file "MACROS getC22g getGFD getGamma.sas". You will find that within the `%MACRO`s `%H2_cullis`, `%H2_oakey` and `%H2RSim`, the respective additional `%MACRO`s are run from github and executed without further need for action by the user. Thus, the addtional macros are included inside some heritability macros.

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/10/SAS_logo_horiz.svg/1200px-SAS_logo_horiz.svg.png" data-canonical-src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/10/SAS_logo_horiz.svg/1200px-SAS_logo_horiz.svg.png" height="150" />
