# ProteinPow(d)er

Een eiwit heeft een belangrijke regelfunctie in het menselijk lichaam en bestaat uit een keten van verschillende aminozuren: H-amino's, C-amino's en P-amino's. Deze ketens moeten op een specifieke manier ontvouwen worden om goed te functioneren; verkeerde ontvouwingen staan aan de basis van ernstige ziektes. De verschillende amino's maken verbindingen met hun indirecte, naastliggende 'buren'. H-amino's liggen graag naast andere H-amino's of naast C-amino's en andersom. H-H verbindingen en H-C/ C-H verbindingen zijn gelijkwaardig sterk. De sterkst mogelijke verbinding bestaat tussen twee C amino's. Om de medische wetenschap te ondersteunen, wordt er gezocht naar een eiwit met meeste verbindingen tussen de aminozuren en daarmee de beste ontvouwing(sscore).

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
 - RandomSolution: Dit algoritme genereert een willekeurige ontvouwde keten en kiest steeds een willekeurige richting voor de aminozuren om te ontvouwen. Alleen de ontvouwde ketens die zichzelf niet kruisen (dwz. valide ontvouwingen) worden behouden.
 - HillClimber: Dit algoritme begint met een willekeurige ontvouwde keten (een valide ontvouwde keten gegenereerd door RandomSolution). Vervolgens worden er kleine aanpassingen aan de keten gemaakt. Het algoritme kiest altijd de beste wijziging (de verbetering van de score) en blijft dit doen totdat er geen betere oplossing meer wordt gevonden, wat resulteert in een lokaal optima. HillClimber kan niet uit een lokaal optima ontsnappen.
 - SimulatedAnnealing: Dit algoritme begint ook met een willekeurige ontvouwde keten, maar in plaats van altijd de beste wijziging te kiezen, wordt er af en toe een slechtere ontvouwing geaccepteerd. Dit gebeurt volgens een afkoelingsschema, waarbij de kans om een slechtere ontvouwing te accepteren afneemt naarmate het algoritme vordert. Dit stelt het algoritme in staat om uit lokale optima te ontsnappen en uiteindelijk een betere algehele oplossing elders te vinden.
 - DepthFirst: Dit algoritme gaat de aminozuren in de eiwitketen één voor één af, waarbij het telkens de meest belovende ontvouwing kiest op basis van de hoogste score. Wanneer er geen geldige ontvouwing meer mogelijk is zonder dat de keten zelf kruist, maakt het algoritme een stap terug (LIFO: last in, first out) naar het vorige aminozuur en probeert daar de op één na beste ontvouwing. Het herhaalt dit proces totdat de volledige keten succesvol is ontvouwen. Om sneller tot oplossingen te komen, wordt er een heuristiek toegepast die het eiwit in kleine stukjes knipt (chunks) en per chunk scores berekent ipv. per aminozuur. 
 - BreadthFirst: Dit algoritme werkt in de basis vergelijkbaar met DepthFirst, alleen werkt dit algoritme volgens een FIFO (first in, first out) principe. Dit betekent dat het algoritme eerst de huidige mogelijke ontvouwingsstaten volledig onderzoekt voordat het verdergaat naar de volgende. Dit zorgt voor een grondigere maar computationeel intensievere zoektocht naar de optimale oplossing. Om de zoekruimte beheersbaar te houden, wordt er hier ook gebruikgemaakt van het heuristiek die het eiwit in kleinere stukjes (chunks) verdeelt en scoort.

#### Stap 3
Bepaal de runtijd van het experiment in seconden.
#### Stap 4 
Run vervolgens het experiment door het aanroepen van:
```
PYTHONPATH=code python main.py
```

Zie de progressie van het experiment in de terminal en vervolgens de beste ontvouwing in een 3d omgeving, gevolgd door een histogram van alle gevonden ontvouwingsscores. Verder worden de resultaten en andere gegevens uit het experiment opgeslagen in een folder genaamd: 'experiment_results', met voor elke run van het experiment een aparte file. Deze folder wordt automatisch aangemaakt tijdens het experiment. 

### Structuur

De hierop volgende lijst beschrijft de belangrijkste mappen en files in het project, en waar je ze kan vinden:

- **/code**: bevat alle code van dit project
  - **/code/algorithms**: bevat de code voor algoritmes
  - **/code/classes**: bevat de drie benodigde classes voor deze case
  - **/code/visualisation**: bevat de code voor de visualisatie

## Auteurs
- Sem Loogman

