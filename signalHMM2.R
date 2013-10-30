#model hmm do sekwencji sygnalowej
#mamy cztery typy sygnałów
#1) aminokwasy z dodatnim ładunkiem K,R,H
#2) aminokwasy hydrofobowe V,I,L,M,F,W,C
#3) aminokwasy polarne bez ładunku S,T,N,Q
#4) reszta D,E,Y,G

#mamy cztery stany ukryte
#1) n-region
#2) h-region
#3) c-region
#4) region odcięcia

#stany ukryte produkują sygnały zgodnie z rozkładem wielomianowym
#zmiana stanów zgodnie z zadaną macierzą przejść (dofitowane parametry rozkładu wykładniczego)

#w takim modelu pakiet depmixS4 powinien sobie sam poradzić




