# ProteinPow(d)er

Een eiwit heeft een belangrijke regelfunctie in het menselijk lichaam en bestaat uit een keten van verschillende aminozuren: H-amino's, C-amino's en P-amino's. Deze ketens moeten op een specifieke manier uitgevouwen om goed te functioneren; verkeerde uitvouwingen staan aan de basis van ernstige ziektes. De verschillende amino's maken verbindingen met hun indirecte, naastliggende 'buren'. H-amino's liggen graag naast andere H-amino's of naast C-amino's en andersom. H-H verbindingen en H-C/ C-H verbindingen zijn gelijkwaardig sterk. De sterkst mogelijke verbinding bestaat tussen twee C amino's. Om de medische wetenschap te ondersteunen, wordt er gezocht naar een eiwit met meeste verbindingen tussen de aminozuren en daarmee de beste uitvouwing(sscore).

## Aan de slag

### Vereisten

Deze codebase is volledig geschreven in Python 3.9.6. In requirements.txt staan alle benodigde packages om de code succesvol te draaien. Deze zijn te installeren via pip dmv. de volgende instructie:

```
pip install -r requirements.txt
```

Of via conda:

```
conda install --file requirements.txt
```

### Gebruik
#### Stap 1
In main.py kies eerst een eiwit bestaande uit een combinatie van C, H, P. 
#### Stap 2
Kies vervolgens een algoritme naar keuze:
 - RandomSolution: Dit algoritme genereert een willekeurige uitgevouwen keten en kiest steeds een willekeurige richting voor de aminozuren om uit te vouwen. Alleen de uitgevouwen ketens die zichzelf niet kruisen (dwz. valide uitvouwingen) worden behouden.
 - HillClimber: Dit algoritme begint met een willekeurige uitgevouwen keten (een valide uitgevouwen keten gegenereerd door RandomSolution). Vervolgens worden er kleine aanpassingen aan de keten gemaakt. Het algoritme kiest altijd de beste wijziging (de verbetering van de score) en blijft dit doen totdat er geen betere oplossing meer wordt gevonden, wat resulteert in een lokaal optima. HillClimber kan niet uit een lokaal optima ontsnappen.
 - SimulatedAnnealing: Dit algoritme begint ook met een willekeurige uitgevouwen keten, maar in plaats van altijd de beste wijziging te kiezen, wordt er af en toe een slechtere uitvouwing geaccepteerd. Dit gebeurt volgens een afkoelingsschema, waarbij de kans om een slechtere uitvouwing te accepteren afneemt naarmate het algoritme vordert. Dit stelt het algoritme in staat om uit lokale optima te ontsnappen en uiteindelijk een betere algehele oplossing elders te vinden.
 - DepthFirst: Dit algoritme gaat de aminozuren in de eiwitketen één voor één af, waarbij het telkens de meest belovende uitvouwing kiest op basis van de hoogste score. Wanneer er geen geldige uitvouwing meer mogelijk is zonder dat de keten zelf kruist, maakt het algoritme een stap terug (LIFO: last in, first out) naar het vorige aminozuur en probeert daar de op één na beste uitvouwing. Het herhaalt dit proces totdat de volledige keten succesvol is uitgevouwen. Om sneller tot oplossingen te komen, wordt er een heuristiek toegepast die het eiwit in kleine stukjes knipt (chunks) en per chunk scores berekent ipv. per aminozuur. 
 - BreadthFirst: Dit algoritme werkt in de basis vergelijkbaar met DepthFirst, alleen werkt dit algoritme volgens een FIFO (first in, first out) principe. Dit betekent dat het algoritme eerst de huidige mogelijke uitvouwingsstaten volledig onderzoekt voordat het verdergaat naar de volgende. Dit zorgt voor een grondigere maar computationeel intensievere zoektocht naar de optimale oplossing. Om de zoekruimte beheersbaar te houden, wordt er hier ook gebruikgemaakt van het heuristiek die het eiwit in kleinere stukjes (chunks) verdeelt en scoort.

#### Stap 3
Bepaal de runtijd van het experiment in seconden.
#### Stap 4 
Run vervolgens het experiment door het aanroepen van:
```
PYTHONPATH=code python main.py
```

Zie de progressie van het experiment in de terminal en vervolgens de beste uitvouwing in een 3d omgeving, gevolgd door een histogram van alle gevonden uitvouwingsscores. Verder worden de resultaten en andere gegevens uit het experiment opgeslagen in een folder genaamd: 'experiment_results', met voor elke run van het experiment een aparte file. Deze folder wordt automatisch aangemaakt tijdens het experiment. 

### Structuur

De hierop volgende lijst beschrijft de belangrijkste mappen en files in het project, en waar je ze kan vinden:

- **/code**: bevat alle code van dit project
  - **/code/algorithms**: bevat de code voor algoritmes
  - **/code/classes**: bevat de drie benodigde classes voor deze case
  - **/code/visualisation**: bevat de code voor de visualisatie

## Auteurs
- Robin Wierda
- Sem Loogman

