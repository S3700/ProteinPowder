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

In main.py kies een eiwit, algoritme en runtime (in seconden). Run vervolgens het experiment door het aanroepen van:
```
PYTHONPATH=code python main.py
```

Zie de progressie van het experiment in de terminal en vervolgens de beste uitvouwing in een 3d omgeving, gevolgd door een histogram van alle gevonden uitvouwingsscores.

### Structuur

De hierop volgende lijst beschrijft de belangrijkste mappen en files in het project, en waar je ze kan vinden:

- **/code**: bevat alle code van dit project
  - **/code/algorithms**: bevat de code voor algoritmes
  - **/code/classes**: bevat de drie benodigde classes voor deze case
  - **/code/visualisation**: bevat de code voor de visualisatie

## Auteurs
- Robin
- Sem Loogman

