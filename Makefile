# Compilateur utilisé
CC=mpicxx

CCseq=g++
#g++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -std=c++11
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -std=c++11

# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler

SRC = main1.cpp GradConj.cpp Problem.cpp BC.cpp Output.cpp Readfile.cpp

SRC = main.cpp GradConj.cpp Problem.cpp BC.cpp Output.cpp Readfile.cpp

SRC1 = main1.cpp GradConj.cpp Problem.cpp BC.cpp Output.cpp Readfile.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

seq: $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)
# à changer le compilateur pour le seq

par: $(SRC1)
	$(CC) $(SRC1) $(CXX_FLAGS) -o $(PROG)
# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
