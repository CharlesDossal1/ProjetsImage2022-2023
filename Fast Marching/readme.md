# Qu'est-ce que l'algorithme de Dijkstra et le Fast Marching ?

L'algorithme de Dijkstra est un algorithme de recherche de chemin dans
une image qui permet de résoudre le problème du plus court chemin entre
deux points sur une image. La valeur prise par chacun des pixels
représente une topographie de terrain, ou une viscosité.

L'algorithme de Dijkstra est un exemple d'algorithme glouton, car il prend la décision la plus avantageuse localement à chaque étape dans l'espoir de trouver une solution globale optimale. Il est utilisé pour trouver les plus courts chemins entre les pixels d'une image, qui peuvent par exemple représenter des réseaux de routes.

Le Fast Marching est une amélioration de l'algorithme de Dijsktra, permettant de calculer des distances *L²* dans une image. Le Fast Marching permet une représentation plus réaliste des distances.

# Contenu du repo

Dans ce repo, vous trouverez un notebook python expliquant point par
point comment utiliser Dijkstra et Fast Marching. Un exemple jouet est
disponible dans le notebook et dans le fichier python
*exemple_jouet_franc.py*.