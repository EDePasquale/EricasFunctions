#################
#               #
# Package Setup #
#               #
#################

#Modified from: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/


#Step 0: Packages you will need
#The packages you will need to create a package are devtools and roxygen2. 
#I am having you download the development version of the roxygen2 package.

install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)



#Step 1: Create your package directory
#You are going to create a directory with the bare minimum folders of R 
#packages. I am going to make a cat-themed package as an illustration.

#Because I made an R Studio project for this package, it advised that I shouldn't 
#run the following code:

#setwd("parent_directory")
#create("cats")

#Instead, I created 2 empty folders in the project called "man" and "R" and an
#empty file called "DESCRIPTION" with no extension. I filled this DESCRIPTION file
#with the package name, title description, version, etc.



#Step 2: Add functions
#Literally just save them in the R folder you created



#Step 3: Add documentation
#The way it works is that you add special comments to the beginning of each 
#function, that will later be compiled into the correct format for package 
#documentation. The details can be found in the roxygen2 documentation: 
#https://github.com/klutometis/roxygen#roxygen2



#Step 4: Process your documentation
#Now you need to create the documentation from your annotations earlier. 
#You’ve already done the “hard” work in Step 3. Step 4 is as easy doing this:
document()



#Step 5: Install!
setwd("..")
install("Ericas_Functions")



#Bonus) Step 6: Make the package a GitHub repo
#First put it on GitHub then...
install_github('Ericas_Functions','EDePasquale')
