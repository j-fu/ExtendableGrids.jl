Hallo,


> .... bei mir reift die Idee im Kopf...
> 
> Grundsätzlich haben Timo und ich folgendes abgesprochen:
> 
> - Es gibt eine Struktur, die  nennt sich
>    SuperGrid, GridCollection oder so, die enthält die Liste der Knotenkoordinaten
>    und mehrere SubGrids, die sparse Adjazenzen gespeichert sind.
> 
> - Für die einzelnen Grids gibt es jetzt verschiedene Möglichkeiten. Die Extreme sind m.E.
>    a) volle Flexibilität:  Elementtyp, Dimension (im Typ codiert) und Region können
>       innerhalb eines SubGrids variieren
>    b) Typ, Dimension und Region sind Attribut des SubGrids und alle Elemente eines gegebenen
>       SubGrids sind gleich (man bräuchte dann auch keine Sparse-Strukturen)
>    In gewisser Weise entsprechen die SubGrids den partitions in pdelib, die Struktur
>    würde auch für die Parallelisierung nutzbar sein.
> 
>    Voraussetzungslos würde ich  eigentlich zu b) tendieren, möglicherweise mit der Ausnahme, dass alle Elemente
>    eines SubGrids dieselbe Dimension haben.
>       Am meisten rückwärtskompatibel wäre eine Variante für die
> Dimension und Typ für
>    alle Elemente eines SubGrids gleich sind, nur Regionen können verschieden sein. Da könnte man sich
>    die Sparsity sparen, und die jetzigen Gitter würden schon mal in dieses Schema passen.
> 
>    Was denkst Du ? - Nach dem, was Du gebaut hast, denke ich dass Deine Bedürfnisse hier
>    mit entscheidend sind.


Ja, ich denke b) würde mir auch etwas besser gefallen. Bei mir sieht
es im Moment so aus, dass mein Mesh-Objekt ein Attribut
elementtypes::Array{AbstractElementType,1} hat, wo die Element-Typen
für jedes subgrid drin stehen könnten. Im Moment habe ich jedenoch
keine echten subgrids und daher steht da in dem Array nur ein Eintrag,
der für das ganze Mesh gilt. Aber man könnte das jetzt
natürlicherweise dahin wachsen lassen, dass ich aus meinem Mesh ein
Array{Mesh,1} mache und dann das j-te SubMesh (oder SubGrid ...) zu
dem j-ten ElemenType gehört und dann alles zusammen mit den
coords4nodes in die OberKlasse (SuperGrid/CompoundGrid) wandert. (Die
Finiten-Elemente muss ich dann auch so umbauen, dass es ein
FinitesElement pro SubGrid gibt. )

Andererseits gibt es zum Beispiel die virtual element methods, die ja
auf beliebigen Ansammlungen von Polygonen arbeiten, d.h. da kann im
Prinzip jede Zelle des Gitters unterschiedlich viele Ecken habe. Das
wäre dann eher mit a) und einer Sparse-Variante für die
Knoten-Zell-Beziehung zu vereinbaren, denke ich. (Allerdings habe ich
momentan nicht vor, virtual element methods zu implementieren, wäre
aber andererseits auch schön, sich diese Option nicht unnötig zu
verbauen). Man kann natürlich auch innerhalb von b) die Polygone aus
dem Gitter aufteilen in SubGrids mit Dreiecken, Vierecken, Fünfecken
etc. (ist es dann wahrscheinlich doch eher unwahrscheinlich, dass da
jemand Zwanzig-oder-mehr-Ecke machen will.)

Also ich votiere daher mal für b).

Und dann wäre da noch die Frage, wie wir die ganzen Adjazenzen
speichern. Ich brauche jetzt zum Beispiel alles Mögliche wie

nodes4faces

faces4cells

cells4faces

normals4faces

area4faces

volumes4cells

sign4cells

bfaces

bregions (4bfaces)

Das möche ich alles irgendwo im Grid speichern, wenn möglich. Im
Moment habe ich daher eine mutable struct, wo es diese ganzen Felder
gibt, die aber erstmal leer initialisert werden und nur bei Bedarf
ausgerechnet werden (nicht jedes Finite-Element benötigt jedes dieser
Felder). Das ist wahrscheinlich nicht die feine englische Art, aber
ich wußte mir erstmal nicht anders zu helfen. Wie löst man das
elegant? Timo hatte mir mal von einer Art Rucksack erzählt, wo man das
alles reintun könnte. Wie würde soetwas aussehen und wie stelle ich
sicher, dass dann der Finite-Element-Code noch weiß, wo er nachschauen
muss? (Und mutable wäre das dann trotzdem? Es sei denn man rechnet
einfach alles im Konstruktor des Meshes aus.)

Man müsste auch nochmal überlegen, wie man die Knoten usw. clever
abspeichert. Ich nehme jetzt immer Array{Float64,...}
bzw. Array{Int64,...}, aber da kennt ihr eventuell bessere Typen? Und
dann war ja da noch die Sache mit dem transponieren, damit die
Nachbarn auch im Speicher Nachbarn sind, d.h. man müßte dann im
SuperGrid die Knoten als 2xAnzahlKnoten-Array speichern und im SubGrid
die DreiecksKnoten als 3xAnzahlDreiecke, richtig? Und das Float64 soll
vermutlich auch variabel gehalten werden?

