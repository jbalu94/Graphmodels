# Véletlen gráfmodellek

E mappa, a 2019 tavaszán írt szakdolgozatom keretében készített szimulációk és megjelenítők forráskódját tartalmazza.

A szakdolgozat megtalálható itt: 

A dolgozatban olyan véletlen növekvő gráfmodelleket vizsgáltam, melyben valamilyen téreblei struktúra jelen van, elméleti és gyakorlati szempontból egyaránt.

A forráskódok szabadon felhasználhatók saját felelősségre.


## Modellek
Mind a 4 féle modellből lehetőségünk van generálni gráfokat és elmenteni a gráfok tulajdonságait. 
A main() függvény megváltoztatásával állíthatjuk be a kívánt paramétereket, és választhatjuk ki a modellt 
Megadhatunk egy paraméterteret, melyet úgy járunk, be hogy minden lehetőséget kipróbálunk, és mindegyik beállítás mellett T db gráfot generálunk. Az eredményeket .txt fájlokba mentem ki, amit aztán ipython notebook segítségével dolgozunk fel.

### Preferential Attachment modell
Használat:
A main() függvényben a megfelelő paraméterek mellett a következő függvényt kell meghívjuk:
simulatePAandSave(path,T,N,M,deltas,gen);

### Geometric Preferential Attachment modell
#### F = F0
Ez az eset megegyezik a PA modellel.

#### F = F1
Használat:
A main() függvényben a megfelelő paraméterek mellett a következő függvényt kell meghívjuk:
simulateGPA1andSave(path,T,N,M,R,alphas,deltas,gen);

#### F = F2
Használat:
A main() függvényben a megfelelő paraméterek mellett a következő függvényt kell meghívjuk:
simulateGPA2andSave(path,T,N,M,betas,alphas,deltas,gen);

#### F = F3
Használat:
A main() függvényben a megfelelő paraméterek mellett a következő függvényt kell meghívjuk:
simulateGPA3andSave(path,T,N,M,betas,alphas,deltas,gen);

### Spatial Preferential Attachment modell
Használat:
A main() függvényben a megfelelő paraméterek mellett a következő függvényt kell meghívjuk:
simulateSPAandSave(path,T,N,P,A1,A2,gen);

### Hyperbolic Preferential Attachment modell
Használat:
A main() függvényben a megfelelő paraméterek mellett a következő függvényt kell meghívjuk:
simulateHPAandSave(path,T,N,M,betas,deltas,gen);


## Eredmények megjelenítése
A processszakdogafiles.ipynb  fájl segítségével könnyedén megjeleníthetjük a generált eredményeket, néhány példa található, hogy hogyan.


