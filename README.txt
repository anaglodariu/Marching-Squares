TEMA1 APD
- alocarile de memorie pentru image, scaled_image si grid sunt facute inaintea crearii threadurilor
- pentru a nu se folosi variabile globale, se creeaza o structura "my_arg"
= se creeaza un vector de structuri "my_arg" pentru a nu avea race conditions, in care doua sau 
mai multe threaduri ar modifica simultan acelasi loc din memorie (cum ar fi modificarea valorii id-ului)
= structura "my_arg" are ca campuri pointeri la bariera, imaginea initiala, contour_map, imaginea 
rescalata si grid, pentru ca acestea vor fi share-uite de threaduri
- la creearea threadurilor este data ca argument functia "thread_function", in care sunt paralelizate 
functiile de rescale, sample_grid si march
- pentu sincronizarea threadurilor sunt folosite 2 bariere intre functiile de rescale, sample_grid
si march, deoarece este necesara terminarea tuturor calculelor unei functii inaintea inceperii celeilalte
= fiind date ca pointeri in structura "my_arg", contour_map, scaled_image si grid, schimbarile ce au loc asupra
lor vor fi vazute de toate threadurile
= campurile din structura care isi vor schimba valoarea mai tarziu in cod sunt date ca pointeri,
pentru ca schimbarile sa fie vazute de toate threadurile
- paralelizarea functiilor rescale_image, sample_grid si march se face prin impartirea loop-ului exterior in segmente
pe baza numarului de threaduri, folosind variabilele start si end
= atunci cand e necesara scalarea imaginii initiale, functia rescale_image este cea mai costisitoare din punct de vedere al timpului 
- abia dupa ce se da join la threaduri, se afiseaza imaginea finala si se elibereaza toata memoria folosita
- de notat ca atunci cand e nevoie de rescale la image, trebuie dezalocata si imaginea initiala, nu doar cea scalata