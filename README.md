# MixtureERGMs

This C++ projects implements variational algorithms for mixture models of exponential random graphs models.

The corresponding research paper was published in Annals of Applied Statistics, 7:2. The draft version can be found in the directory **papers**

## Compiling the code:

1. Unzip the file

2. Go to the unzipped directory

3. Download boost library from http://www.boost.org/users/download/ and unzip the file. No compilation for boost library is needed.

4. Specify directories in the file "makefile":
	a. HBOOST : path to where you unzip the boost library
	b. HRLIB : path to R header files, e.g. /usr/share/R/include
	c. HMYDIR : path to the source folder, i.e to the upzipped directory + /src
	d. RLIB = path to the directory containing compiled libraries of R, e.g. /usr/lib/R/lib

5. Run "make"

## Run the code:

1. Examples are in the file exampleCommands.txt

2. Sample data sets are in directory data:

...a. **Epinions** from Paolo Massa and Kasper Souren of trustlet.org.

Massa, P., and Avesani, P. (2007), “Trust metrics on controversial users: balancing be-
tween tyranny of the majority and echo chambers,” International Journal on Semantic
Web and Information Systems.

...b. **PolBlogs**:

Adamic, L. and Glance, N. (2005). The political blogosphere and the 2004 US election:
Divided they blog. In Proceedings of the 3rd International Workshop on Link Discovery
36–43. ACM Press, New York.

