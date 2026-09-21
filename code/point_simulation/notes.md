# Note tecniche — scena sintetica

Documentazione di `synthetic_point_cloud.py`, `change_region.py`,
`camera_locus.py`, `visibility_filter.py`, `detection_sweep.py`,
`knowledge.py`, `measurement.py`, `planner.py`, `trajectory.py`,
`plot_scene.py` e `scene_viewer.py`. Il codice
porta solo commenti da una riga; il ragionamento, le formule e la taratura
stanno qui.

## Perché una superficie e non un volume

Un sensore di profondità misura superfici: non restituisce mai punti
all'interno di un solido. La scelta non è estetica, decide se il resto della
pipeline ha senso. Occlusione, angolo di incidenza e confine fra noto e
ignoto sono proprietà che esistono solo se la scena ha superfici e c'è un
punto di vista. Con punti sparsi in un volume nessuna di queste è definita.

## Perché un'isosuperficie di un campo casuale

L'alternativa era un catalogo di forme parametriche (sfera, toro, C, L) con
parametri casuali. Scartata: le forme restano quelle che hai nominato, e la
varietà è apparente. L'isosuperficie di un campo casuale dà invece
gratuitamente, dalla casualità e non dalla costruzione:

- **topologia arbitraria** — manici, cavità, componenti staccate;
- **concavità e auto-occlusione**, che è il motivo per cui la selezione di
  viste non è banale: su un corpo convesso ogni vista vede la silhouette e il
  problema si risolve da solo;
- **multi-scala con una manopola sola**, il decadimento spettrale;
- **riproducibilità** dal seed.

Un tentativo intermedio, scartato, era un tubo generato spazzando una sezione
circolare lungo una curva casuale di Fourier. Funzionava ma produceva sempre
un filo ingrossato: topologicamente un tubo, mai un oggetto con massa.

## Il campo

$$f(x) = \sum_k a_k \sin(w_k \cdot x + p_k) - c(x), \qquad a_k = |w_k|^{-\text{decay}}$$

con $w_k$ il vettore d'onda del modo $k$ e $p_k$ la sua fase. L'oggetto è la
regione $f(x) > \text{level}$, quindi la superficie campionata è
$f(x) = \text{level}$.

Gradiente, analitico e calcolato dalla stessa proiezione del valore:

$$\nabla f(x) = \sum_k a_k \cos(w_k \cdot x + p_k)\, w_k$$

**Direzioni** uniformi sulla sfera, ottenute normalizzando un vettore
gaussiano. Campionare gli angoli uniformemente non va: addensa ai poli.

**Frequenze** log-uniformi sulla banda, così ogni ottava è equamente
rappresentata. Uniformi in frequenza affollerebbero l'estremo alto e
lascerebbero sotto-campionate le feature grandi, che sono quelle che
definiscono la forma complessiva.

**Ampiezze** che scalano come la frequenza elevata a `-decay`. È tutto il
significato della manopola: `decay` è la pendenza dello spettro di potenza, e
scambia la dimensione delle feature con la loro rugosità.

## Normalizzazione a deviazione unitaria

Il campo viene riscalato perché abbia deviazione standard unitaria. Senza,
la sua ampiezza dipende da `n_modes` e da `decay`, quindi `STRENGTH`
significherebbe una cosa diversa per ogni configurazione e il termine di
chiusura andrebbe ritarato ogni volta. Con la normalizzazione `STRENGTH` è
misurato in deviazioni standard del campo, e `decay` resta una manopola
indipendente.

## Termine radiale di chiusura

$$c(x) = \text{STRENGTH} \cdot s(x)^2, \qquad s(x) = \max\left(0, \frac{|x| - r_{\text{in}}}{R - r_{\text{in}}}\right)$$

con $r_{\text{in}} = \text{CLOSURE} \cdot R$. Gradiente:

$$\nabla c(x) = \frac{2\,\text{STRENGTH}\, s(x)}{(R - r_{\text{in}})\,|x|}\, x$$

Senza questo termine l'isosuperficie viene tagliata aperta dove attraversa la
sfera della scena, e l'oggetto resta senza un interno ben definito — cosa che
serve a valle per i test dentro/fuori nell'occlusione.

Il `max(0, ...)` è ciò che lo rende **a un lato**: $s$ è identicamente zero
dentro $r_{\text{in}}$, quindi il termine e il suo gradiente svaniscono lì e
la forma non viene disturbata. Una prima versione usava una potenza quarta su
tutto il raggio: dominava il campo e schiacciava ogni oggetto a palla. È uno
dei due errori che hanno prodotto i blob convessi del primo tentativo.

## Isolivello dal quantile

L'oggetto è la regione $f(x) > \text{level}$, quindi la frazione di volume
della scena occupata è la frazione di campioni sopra il livello. Prendere il
quantile $1 - \text{fill}$ di $f$ su un campione uniforme nel volume dà
allora esattamente il livello che lascia dentro la frazione richiesta,
qualunque sia la distribuzione di $f$ per quel campo particolare.

## Campionamento: il guscio a spessore costante

La distanza al primo ordine da $x$ all'isosuperficie è

$$d(x) = \frac{|f(x) - \text{level}|}{|\nabla f(x)|}$$

Dividere la differenza di campo per la pendenza trasforma quindi il test in
una distanza vera. Conta per la densità: un guscio a differenza di campo
costante sarebbe spesso dove il campo è piatto e sottile dove è ripido, e i
punti si ammasserebbero nelle zone piatte. Con un guscio a spessore costante
escono distribuiti uniformemente sulla superficie — verificato sotto.

Il campionamento è per rigetto, a lotti. Il tasso di accettazione non è noto
a priori perché dipende da quanta superficie ha quel campo dentro la scena,
quindi si chiede quello che manca e si ripete. Sovra-campionamento ×8: il
guscio è sottile, la maggior parte dei candidati cade lontano.

## Proiezione di Newton

Newton sull'equazione scalare $f(x) = \text{level}$. Il gradiente è una riga
sola e non un sistema quadrato, quindi si usa la sua pseudo-inversa:

$$x \leftarrow x - \frac{f(x) - \text{level}}{|\nabla f(x)|^2}\, \nabla f(x)$$

Il passo si muove lungo la normale all'insieme di livello, che è anche la via
più breve per raggiungerlo. Convergenza quadratica, da cui il numero piccolo
di passi.

Il termine radiale è ripido vicino al bordo, quindi lì un passo di Newton può
spingere un punto fuori dalla sfera della scena. Quei punti vengono scartati
e non riportati dentro: clampare li lascerebbe fuori dalla superficie.

## Taratura delle costanti

| costante | valore | ruolo |
|---|---|---|
| `FREQ_MIN` | 3.0 | Estremo basso della banda, in unità di 1/R. Fissa la dimensione della feature più grande. |
| `FREQ_MAX` | 16.0 | Estremo alto. |
| `CLOSURE` | 0.60 | Frazione del raggio entro cui il termine radiale è nullo. |
| `STRENGTH` | 3.0 | Deviazioni standard del campo raggiunte dal termine radiale al bordo. |
| `SHELL` | 0.02 | Semispessore del guscio di campionamento, in frazioni di R. |
| `STEPS` | 4 | Passi di Newton. |
| `LEVEL_PROBE` | 20000 | Campioni per stimare l'isolivello. |
| `NORM_PROBE` | 4000 | Campioni per stimare la deviazione del campo. |
| `MAX_BATCHES` | 200 | Lotti prima di rinunciare. |

Note sulla taratura, cioè in che direzione muoverle:

- **`FREQ_MIN` è critica.** Con la lunghezza d'onda più lunga molto maggiore
  del raggio il campo è quasi costante sulla scena e ogni oggetto esce come
  un blob convesso. Il primo tentativo aveva `FREQ_MIN = 1.5`, cioè
  lunghezza d'onda ~4R: tutti blob. Deve essere comparabile alla scena.
- **`SHELL`** sottile significa pochi candidati sopravvissuti, spesso
  significa che Newton deve viaggiare più a lungo e può attraversare in un
  ramo vicino della superficie.
- **`STEPS = 4`** è generoso: dal guscio sopra il residuo arriva alla
  precisione macchina in due o tre passi.
- **`LEVEL_PROBE` e `NORM_PROBE`** stimano un quantile e un momento di una
  funzione liscia sulla palla: qualche migliaio di campioni è già molto più
  della precisione necessaria.
- **`MAX_BATCHES`** lo raggiunge solo un campo patologico, o un `fill` così
  estremo che l'isosuperficie è quasi vuota.

## Verifiche

**I punti stanno sulla superficie.** Residuo $|f - \text{level}| /
|\nabla f|$ sui punti restituiti: mediana 5.7e-17, massimo 2.2e-13. Newton
converge alla precisione macchina.

**La densità superficiale è uniforme.** Coefficiente di variazione della
distanza al vicino più prossimo: 0.521, contro 0.523 di un processo di
Poisson uniforme sul piano. I punti sono un campione uniforme sulla
superficie, che è ciò che il guscio a spessore costante doveva garantire.

**L'auto-occlusione è reale.** Misurata con l'operatore di hidden point
removal di Katz et al. (2007), calibrato su una sfera: misura 0.314 contro il
valore analitico esatto $(1 - 1/d)/2 = 0.300$ da distanza $d = 2.5$, quindi
lo strumento è affidabile. Frazione di punti visibile da una singola vista,
media su 30 pose × 5 seed:

| scena | media | minimo |
|---|---|---|
| sfera (convessa, riferimento) | 0.314 | 0.300 |
| isosuperficie `fill=0.30` | 0.181 | 0.079 |
| isosuperficie `fill=0.15` | 0.153 | 0.065 |
| isosuperficie `decay=0.8` | 0.114 | 0.053 |

Una vista singola vede fra un terzo e un ottavo di quello che vedrebbe su un
corpo convesso, e ci sono pose che ne vedono il 5%. È la condizione che rende
la selezione di viste non banale.

Risultato non ovvio: **`decay` è la manopola di occlusione più forte**, molto
più di quanto suggerisca l'impressione visiva. Abbassarlo aggiunge dettaglio
fine che si auto-occlude.

## Scena prima e dopo il cambiamento

`generate_scene_pair` produce la stessa scena **due volte**: come la conosce
il modello obsoleto e com'è adesso, con oggetti aggiunti, rimossi o spostati.
È ciò che rende il simulatore capace di mettere alla prova l'argomento per cui
si lavora a punti.

### Perché due scene e non un'etichetta

Una prima versione (`change_region.py`, sotto) marcava come cambiati punti
che esistevano già, cioè simulava una superficie *nota* che cambia. È
esattamente il caso che i metodi basati su Fisher, FisherRF in testa, già
risolvono: i Gaussiani ci sono e cambiano valore. Il caso che motiva la tesi è
un oggetto **aggiunto**, che non ha Gaussiani e quindi ha informazione di
Fisher identicamente nulla. Con l'etichetta il pianificatore conosce le
posizioni dei punti cambiati prima di averli visti, cioè esattamente
l'informazione che nel caso reale manca. In più l'etichetta ignora
l'occlusione — un oggetto nuovo nasconde lo sfondo, uno tolto lo rivela — e
rende la rimozione rilevabile come un'aggiunta, mentre per accorgersi che
qualcosa è sparito bisogna vedere *attraverso* il posto dove stava.

L'etichetta resta il modello giusto per un solo caso: un cambiamento di
**aspetto** con la stessa geometria, come una superficie ridipinta.

Misurato, la differenza è netta. Con due camere di sweep l'etichettatura
dava il 46% del cambiamento visto, gli oggetti veri il **25%**; con quattro,
96% contro **71%**. Un oggetto vero ha facce rivolte verso la superficie su
cui poggia e un retro, che una toppa etichettata non ha.

### Composizione con i campi impliciti

Ogni parte della scena — lo sfondo e ogni oggetto — è un campo con
l'isolivello già sottratto, quindi solida dove è positiva. Allora:

- **unione**: $g_{\text{scena}}(x) = \max(g_{\text{sfondo}}(x), g_1(x), g_2(x), \dots)$;
- **rimozione**: si toglie un termine dal massimo;
- **spostamento**: si valuta il campo dell'oggetto in coordinate
  trasformate, $g(R^\top (x - t))$, con gradiente $R \nabla g$.

La superficie dell'unione è corretta da sola anche dove un oggetto tocca o
compenetra lo sfondo. Un oggetto è un campo casuale come lo sfondo, nella
propria palla di raggio `object_radius` e con `OBJECT_FILL` più alto, così
esce corposo.

Uno **spostamento è una rimozione in una posa più un'aggiunta in un'altra**,
con la stessa forma. Si perde il fatto che sia lo stesso oggetto: per
pianificare non conta, perché bisogna guardare sia dov'era sia dov'è, ma per
aggiornare il GS sì, perché un oggetto spostato ha già i suoi Gaussiani e
basterebbe trasformarli rigidamente. È una semplificazione dichiarata.

### Il cambiamento si misura, non si etichetta

Il cambiamento è la **differenza geometrica** fra le due scene:

- **aggiunto** = punto della scena nuova lontano dalla superficie vecchia;
- **rimosso** = punto della scena vecchia lontano dalla superficie nuova;

con la distanza al primo ordine $|g| / |\nabla g|$ e soglia `TOLERANCE`
mezza spaziatura. Due conseguenze escono giuste senza programmarle:

- uno **spostamento breve** segna solo la parte dell'oggetto che non si
  sovrappone alla posa vecchia;
- un oggetto appoggiato **rimuove la superficie che copre**, e toglierlo la
  **rivela**, che risulta aggiunta: è informazione che il modello obsoleto
  non ha mai avuto.

La differenza è solo geometrica: un cambio d'aspetto su superficie che resta
al suo posto non viene visto.

### Appoggio e densità

Un oggetto poggia sullo sfondo: si sceglie un punto di superficie con normale
rivolta verso l'alto (componente verticale oltre `UPWARD`) e lo si fa
scendere lungo la normale, provando le distanze in `CONTACT`, finché almeno
`BURIED` dei suoi punti affonda nella superficie. Gli oggetti restano
separati di almeno due raggi, tutte le pose comprese. Si richiede solo che il
**centro** stia nella sfera della scena: lo sfondo la riempie fino al bordo,
le superfici rivolte in alto stanno lì vicino, e imporre che l'oggetto intero
ci stia dentro lasciava 18 ancore utili su 720 e il posizionamento falliva.
Un oggetto può quindi sporgere dalla sfera della scena al massimo del proprio
raggio.

Lo sfondo si campiona una sola volta ed è condiviso fra le due scene, così i
punti che non cambiano coincidono esattamente. Gli oggetti hanno la stessa
**densità** di punti dello sfondo — altrimenti la densità diventerebbe un
peso implicito nell'obiettivo. L'area di una superficie si misura durante il
campionamento stesso: il guscio di semispessore $h$ contiene un volume
$2hA$, quindi la frazione di candidati che ci cade dà $A$.

### Costanti

| costante | valore | ruolo |
|---|---|---|
| `OBJECT_FILL` | 0.45 | Frazione della propria palla occupata da un oggetto. |
| `OBJECT_PROBE` | 256 | Punti campionati prima, per misurare l'area di un oggetto. |
| `MIN_OBJECT` | 50 | Punti minimi di un oggetto, anche se piccolo. |
| `UPWARD` | 0.3 | Componente verticale minima della normale su cui poggiare. |
| `CONTACT` | 1.0 … 0.2 | Distanze di discesa provate, in raggi dell'oggetto. |
| `BURIED` | 0.02 | Frazione dell'oggetto affondata che conta come contatto. |
| `ATTEMPTS` | 200 | Pose provate prima di rinunciare a un oggetto. |
| `TOLERANCE` | 0.5 | Distanza che conta come cambiamento, in spaziature. |

### Verifiche

- **Nessun cambiamento** → le due scene sono identiche punto per punto.
- **Aggiunta e rimozione sono inverse esatte.** Con lo stesso seme lo stesso
  oggetto nella stessa posa, una volta messo e una volta tolto:

| caso | old | new | aggiunti | rimossi |
|---|---|---|---|---|
| nessuno | 20000 | 20000 | 0 | 0 |
| 1 aggiunto | 20000 | 20702 | 877 | 179 |
| 1 rimosso | 20702 | 20000 | 179 | 877 |
| 1 spostato | 20702 | 21018 | 1236 | 918 |
| 2 aggiunti, 1 rimosso, 1 spostato | 21481 | 22538 | 3003 | 1938 |

  I 179 sono lo sfondo coperto dall'oggetto: aggiungerlo li rimuove,
  toglierlo li rivela.
- **Densità**: spaziatura mediana 0.0096 sugli oggetti contro 0.0103 su
  tutta la scena.
- `generate_point_cloud` è **identica al bit** a prima della
  ristrutturazione, quindi tutti i numeri delle sezioni precedenti restano
  validi.

## Regione cambiata

> **Superata** dalla scena prima e dopo il cambiamento. `change_region.py`
> non è più usato da nessun modulo; resta descritto qui finché il file
> esiste.

Implementata in `change_region.py`, separata dal generatore. Lo stato viene
assegnato **dopo** la generazione e **a caso**, come nel MATLAB dove
`create_change_cluster` è un file distinto da `create_splats`. Il vantaggio
della separazione è che la stessa scena si riusa con cambiamenti diversi
muovendo solo il seed della regione, e i due seed sono indipendenti.

La regione è una **palla geodetica** su un grafo dei k vicini pesato con la
distanza euclidea: si estrae un punto seme a caso e si marcano tutti i punti
entro `reach` da esso, misurando la distanza lungo la superficie e non
attraverso lo spazio vuoto.

$$\text{reach} = \text{extent} \cdot \max_i |x_i - \bar{x}|$$

Il raggio della nuvola come riferimento serve perché la dimensione di una
regione non dipenda dalla scala a cui la scena è stata generata.

Con più regioni, un seme candidato viene rifiutato se le regioni già piazzate
lo raggiungono entro `SEPARATION` volte il reach, altrimenti due regioni si
fondono e ne ottieni meno di quante chieste — lo stesso vincolo di distanza
minima fra i centri che c'è in `create_change_cluster`.

### Perché la geodetica, e perché la motivazione ovvia è sbagliata

La motivazione naturale sarebbe: su una superficie con una fessura le due
pareti sono spazialmente vicine ma lontane lungo la superficie, quindi una
palla euclidea prende due lembi staccati e produce una Ω con una struttura di
occlusione assurda. **Misurato, non succede.** Confrontando le due metriche a
parità di seme e di raggio su 24 configurazioni (`fill` da 0.10 a 0.30,
`decay` da 1.0 a 1.6, 12 semi ciascuna):

| | geodetica | euclidea |
|---|---|---|
| punti inclusi in più dall'euclidea | — | 12–20% |
| di cui oltre 2× reach lungo la superficie | — | 12 su ~9000 |
| componenti connesse della regione | 1 sempre (20 semi) | max 2, una volta |

I punti che solo l'euclidea include stanno a distanza geodetica 1.0–1.3 volte
il reach: sono un bordo, non un lembo lontano. La differenza fra le due
metriche è che la palla euclidea sborda un po' oltre il confine dove la
superficie curva, niente di topologico.

La ragione per cui la geodetica resta la scelta è più modesta: **garantisce**
una regione connessa sulla superficie, mentre l'euclidea si è spezzata in due
componenti in un caso su 24. Per una tesi in cui Ω deve essere una regione
cambiata ben definita, una garanzia vale le trenta righe in più. In aggiunta,
il grafo kNN non è infrastruttura sprecata: serve comunque per stimare le
normali con una PCA locale quando arriverà il modello di sensore.

Se in futuro il costo del grafo diventasse un problema, passare all'euclidea
è un cambio di due righe e sposta il 15% dei punti al bordo della regione.

### Costanti

| costante | valore | ruolo |
|---|---|---|
| `NEIGHBOURS` | 10 | Vicini per punto nel grafo su cui corre la geodetica. |
| `SEPARATION` | 2.0 | Distanza minima fra semi, in unità di reach. |
| `ATTEMPTS` | 200 | Candidati estratti prima di rinunciare a una regione. |

Nota d'uso: `extent` è un raggio, non una frazione di punti. Su una nuvola di
20000 punti a `fill=0.30`, `extent=0.25` marca circa il 2% dei punti e
`extent=0.15` circa lo 0.7%, perché l'area di una calotta di raggio
geodetico r cresce come r² mentre l'area totale della superficie è molto
maggiore. Il MATLAB parametrizzava invece per numero di punti (`n_min`,
`n_max`): se in fase di esperimenti ti serve controllare la frazione di
scena cambiata piuttosto che la dimensione fisica del cambiamento, conviene
aggiungere quella modalità.

## Locus delle camere

Implementato in `camera_locus.py`. Le camere hanno il centro su una sfera
attorno all'origine e guardano tutte verso di essa. Il raggio del locus è
**indipendente** da quello che contiene i punti della scena: il primo è una
proprietà dello spazio di lavoro del braccio, il secondo dell'oggetto.

### Convenzione delle pose

La stessa del MATLAB, per restare compatibile con il modo in cui `render.m`
usa le pose. Rotazione **world-to-camera** in convenzione OpenCV, con l'asse
z della camera in avanti e l'asse y verso il basso, e traslazione come
**centro della camera in coordinate mondo**:

$$x_{\text{cam}} = R \, (x_{\text{world}} - c)$$

Le righe di $R$ sono, in ordine, `right`, `down`, `forward`. La direzione
`forward` è quella che va dal centro della camera verso l'origine. La `right`
si ottiene come prodotto vettoriale fra `forward` e la verticale del mondo,
il che fissa il rollio lasciando l'orizzonte orizzontale; sopra i poli
`forward` è parallela alla verticale, il prodotto vettoriale si annulla e si
ripiega su un asse alternativo.

La matrice di calibrazione non entra nel posizionamento: viene trasportata
nel risultato perché il gruppo di camere sia autosufficiente per proiettare,
non perché serva a collocarle.

### Perché un reticolo di Fibonacci e non una griglia

La richiesta era una "griglia", ma una griglia vera in latitudine e
longitudine addensa le camere ai poli e lascia rado l'equatore, e il problema
peggiora al crescere del numero di camere. Il reticolo di Fibonacci —
passi uguali in quota, angolo d'oro in azimut — distribuisce quasi
uniformemente. Spaziatura angolare al vicino più prossimo, in gradi:

| camere | metodo | media | min | max | max/min |
|---|---|---|---|---|---|
| 100 | Fibonacci | 19.28 | 17.78 | 19.73 | **1.1** |
| 98 | lat-long | 16.46 | 5.68 | 25.71 | 4.5 |
| 200 | Fibonacci | 13.80 | 12.55 | 14.31 | **1.1** |
| 200 | lat-long | 11.49 | 2.80 | 17.78 | 6.3 |
| 400 | Fibonacci | 9.64 | 8.87 | 9.89 | **1.1** |
| 406 | lat-long | 7.91 | 1.39 | 12.34 | 8.9 |

Il rapporto fra spaziatura massima e minima resta 1.1 per Fibonacci a
qualunque numero di camere, mentre per la griglia passa da 4.5 a 8.9. Conta
perché il locus è lo spazio su cui la selezione ottimizza: una
discretizzazione non uniforme introduce un bias, e gruppi di camere quasi
coincidenti ai poli portano informazione ridondante e peggiorano il
condizionamento della matrice di Fisher.

### Verifiche

Tutte esatte alla precisione macchina, su 200 camere:

- centri sulla sfera, raggio costante a 1e-16;
- rotazioni ortonormali, $R R^\top = I$ a 4e-16, e destrorse, $\det R = +1$;
- **l'origine si proietta nel punto principale di ogni camera** con scarto
  massimo 9e-14 pixel. È la verifica che conta, perché passa attraverso la
  calibrazione e valida insieme orientamento, convenzione e intrinseci;
- origine sempre davanti alla camera, $z_{\text{cam}} > 0$, cioè la
  convenzione OpenCV è quella giusta e non è specchiata.

### Costanti

| costante | valore | ruolo |
|---|---|---|
| `GOLDEN` | $\pi(3 - \sqrt{5})$ | Angolo d'oro, passo in azimut del reticolo. |
| `UP` | $(0,0,1)$ | Verticale del mondo, fissa il rollio. |
| `ASIDE` | $(1,0,0)$ | Verticale di ripiego per una camera sopra il polo. |

### Banda di elevazione

`min_elevation` e `max_elevation` limitano il locus a una fascia, per
escludere le pose che il braccio non può raggiungere — tipicamente il
sottosuolo, se l'oggetto sta su un tavolo. Sono in **gradi** sopra
l'equatore e non in radianti come nel MATLAB: sono parametri che si scrivono
a mano in uno script, e `-30` si legge mentre `-0.52` no. I default
`(-90, 90)` sono la sfera intera, cioè nessun vincolo.

Il reticolo campiona uniformemente in **quota** e non in elevazione, che è
ciò che lo rende quasi uniforme: su una sfera l'area di una fascia è
proporzionale alla sua altezza, non all'angolo che sottende. Per questo la
banda si traduce in $z \in [\sin(\text{min}), \sin(\text{max})]$, e
l'uniformità si conserva dentro la fascia:

| banda | elevazione ottenuta | spaziatura max/min |
|---|---|---|
| [-90, 90] | [-84.3, 84.3] | 1.14 |
| [-30, 90] | [-29.8, 85.0] | 1.12 |
| [0, 90] | [0.1, 85.9] | 1.12 |
| [-20, 70] | [-19.8, 69.5] | 1.04 |
| [20, 60] | [20.1, 59.9] | 1.04 |

## Passata di detection

`detection_sweep.py` simula la passata economica che precede la
pianificazione: poche pose, distribuite il più lontano possibile, e la parte
del cambiamento che vedono davvero. Quella è la regione cambiata **nota** al
pianificatore; la regione generata da `change_region.py` diventa il **ground
truth**. La differenza fra le due è il cambiamento che lo sweep non ha visto,
cioè ciò che l'esplorazione deve ancora trovare. È il primo pezzo del modello
di acquisizione: rende Ω stimata invece che assunta.

### Dispersione, non copertura

"Più lontane possibile" e "copertura massima" sono obiettivi diversi. La
copertura vera è un problema di max-coverage sulla matrice di visibilità,
risolvibile con il greedy a garanzia $(1 - 1/e)$, ma va calcolata su una
geometria: sulla nuvola vera sarebbe barare, perché userebbe informazione
che lo sweep non ha. La versione legittima la calcolerebbe sul modello
obsoleto, che è noto — ma nel simulatore attuale vecchio e nuovo hanno la
stessa geometria, e il cambiamento è solo un'etichetta sui punti. Finché è
così, la dispersione geometrica è la scelta onesta: non guarda la scena.

### Farthest point sampling sul locus esistente

Le pose si scelgono con un farthest point sampling greedy sulle candidate
del reticolo di Fibonacci, invece che da un reticolo nuovo con poche camere.
Così lo sweep è un **sottoinsieme delle candidate**, e il pianificatore sa
quali viste sono già state prese — da escludere o da trattare come
informazione a priori $F_0$. Rispetta anche da solo la banda di elevazione.
La prima posa è la più alta, il che rende la scelta deterministica senza
seed. La distanza è quella euclidea fra i centri, monotona con quella
angolare sulla sfera.

Distanza minima fra le camere scelte, contro un reticolo nuovo di $k$
camere:

| $k$ | FPS sul locus | Fibonacci nuovo |
|---|---|---|
| 4 | 4.33 | 3.56 |
| 6 | 3.14 | 2.89 |
| 8 | 2.55 | 2.50 |
| 12 | 1.67 | 2.05 |

Per pochi punti il FPS disperde meglio; da una dozzina in su peggio, perché è
greedy e pesca da pose predefinite. Per uno sweep, che per definizione ha
poche camere, va bene.

### Detection per punto, con due misure diverse

Ciò che lo sweep osserva lo decide il proxy della detection in
`knowledge.py`, sotto. Lo sweep guarda la **scena nuova**; la vecchia è solo
ciò che il modello ricorda, confrontato con ciò che le camere vedono.

- Un'**aggiunta** si rileva vedendo un punto nuovo: `visible_points` sulla
  scena nuova.
- Una **rimozione** si rileva vedendo *oltre* il posto dove stava un punto
  vecchio: `vacated_points`, vedi il filtro di visibilità.

Nessuna delle due si estende a tutta la regione: vederne un lembo dice che lì
c'è cambiamento, non fin dove arriva, e il resto rimane ignoto. Quel bordo fra
osservato e non osservato è esattamente dove andranno i punti fantasma
dell'esplorazione.

Riserva dichiarata: nella realtà un punto visto non basta per una detection.
Il render-and-compare deve produrre un residuo sopra il rumore, quindi serve
vedere abbastanza della regione. Per ora è ignorato.

### Le pose ammissibili, superate dal pianificatore

`feasible_cameras` teneva le pose che vedono almeno `min_points` aggiunte
note, sugli occlusori creduti. Era un filtro al posto di un pianificatore, e
ora che il pianificatore c'è non è più usata: il pianificatore guarda tutte le
pose non ancora prese e le ordina per guadagno di informazione. Quando la
soglia filtrava ancora, sulla scena etichettata con sweep rado toglieva 10
pose su 94 che guardavano solo cambiamento mai rilevato.

### Le pose già prese: $F_0$ ed esclusione esplicita

Le pose già acquisite contano due volte. Primo, ciò che hanno misurato entra
come informazione a priori $F_0$: è la `information` accumulata in
`Knowledge`, e fa cercare al pianificatore viste *diverse* da quelle già
prese. Secondo, vanno **escluse esplicitamente** dalle candidate.

Qui c'era un errore nelle note precedenti, che sostenevano che con $F_0$ una
posa già presa rendesse zero al margine e non servisse escluderla. È falso:
nel modello di Fisher una misura ripetuta è indipendente e aggiunge
informazione, $\log\det(A + F) > \log\det(A)$ finché $F \neq 0$. Il
rendimento cala, non si annulla. Misurato: rifare una posa dello sweep rende
**39.4**, non zero. Senza l'esclusione il pianificatore rischierebbe di
rifare viste che ha già.

### Cosa vede lo sweep

Cambiamento misto — 2 aggiunti, 1 rimosso, 1 spostato — con 3003 punti
aggiunti e 1938 rimossi:

| $k$ | aggiunti visti | rimossi visti oltre |
|---|---|---|
| 1 | 25.3% | 66.6% |
| 2 | 25.4% | 75.7% |
| 4 | 71.2% | 80.0% |
| 6 | 81.1% | 85.1% |
| 10 | 93.0% | 85.9% |

**Le rimozioni si confermano con poche viste, le aggiunte no.** Tolto un
oggetto, *tutti* i suoi punti lungo il raggio — davanti e dietro — sono
spazio vuoto, e una vista sola li attraversa tutti. Un'aggiunta invece la
si ricostruisce una faccia alla volta. È un'asimmetria che il pianificatore
dovrà rispettare. La quota di rimossi che resta mai confermata, circa il
14%, è in buona parte sfondo coperto dagli oggetti nuovi, che non si può
vedere attraverso un oggetto.

La seconda camera del farthest point sampling non aggiunge quasi nulla
(25.3% → 25.4%). Il motivo è una proprietà del farthest point sampling su una
banda: partendo dalla camera più alta, le successive sono le più lontane, e
su una calotta i punti più lontani fra loro stanno sul **cerchio di bordo**,
che è il più largo. Con la banda da −10° a 85° le prime quattro camere hanno
elevazione 82.4°, −9.7°, −5.8°, −9.2°: tre su quattro guardano la scena di
taglio dal bordo basso, da dove le superfici rivolte in alto su cui poggiano
gli oggetti si vedono male. Lo sweep è quindi disperso sul locus ma non nelle
direzioni di vista che contano; è un argomento per la copertura sul modello
obsoleto invece della sola dispersione, quando ci sarà.

## Conoscenza del sistema e nuvola unita

### Lo scope: la detection è data per fatta

Nella pipeline reale l'ingresso è il GS vecchio con la sua nuvola; lo sweep
produce viste da confrontare con i render del GS (canale fotometrico) e una
nuvola grezza da confrontare con quella del GS (canale geometrico). In questa
fase la detection **non si simula**: si assume fatta, per isolare il
pianificatore dai suoi errori e misurarli separatamente più avanti.

"Fatta" vuol dire **perfetta su ciò che le viste osservano**, e ignara del
resto. L'alternativa — perfetta e completa, cioè nota anche sulle facce mai
viste — è stata scartata: nel reale i punti aggiunti vengono dalla nuvola
grezza dello sweep e non esistono dove nessuno ha guardato, e con una
detection completa il pianificatore tornerebbe a conoscere posizioni che non
può avere.

### Il proxy della detection

`observe` in `knowledge.py`. Una vista rivela i punti aggiunti che vede
(`visible_points` sulla scena nuova) e conferma i rimossi oltre cui vede
(`vacated_points`). Lo stesso proxy serve lo sweep e poi ogni vista che il
pianificatore acquisirà: la conoscenza cresce con le viste prese. Tre
proprietà verificate:

- **monotona**: si aggiunge soltanto, niente viene dimenticato;
- **indipendente dall'ordine**: osservare le viste A e poi B dà esattamente
  ciò che dà osservarle insieme, che è ciò che serve in un ciclo di
  acquisizione;
- **equivalente allo sweep precedente**, al bit, sulle stesse pose.

`Knowledge` tiene anche le pose già prese (`acquired`): è da lì che verrà la
$F_0$ del pianificatore.

### La nuvola unita

`merge` produce la nuvola su cui lavora il pianificatore: la nuvola vecchia —
ciò che il modello obsoleto conosce — più i punti aggiunti osservati finora,
ciascuno con uno stato:

- **invariato**: punto vecchio che il sistema crede ancora lì, compresi i
  rimossi non ancora confermati;
- **rimosso**: punto vecchio oltre cui una vista ha visto;
- **aggiunto**: punto nuovo che una vista ha visto.

I punti nuovi non si possono segnare sulla nuvola del GS perché lì non
esistono: per questo l'uscita è un'unione e non una nuvola etichettata. Da
ricampionare sono rimossi e aggiunti; **occludono** gli invariati e gli
aggiunti (`occluding`), non i rimossi confermati.

### Il pianificatore vede attraverso ciò che sa

Le occlusioni si calcolano sugli **occlusori creduti**, non sulla scena nuova
vera, che il sistema non ha. Dove il cambiamento non è osservato il modello
creduto sbaglia in due direzioni opposte: un oggetto tolto e non confermato
continua a occludere, un pezzo aggiunto e non visto non occlude.

Sulla scena mista (seed 3) dopo uno sweep di 4 camere restano **388** punti
rimossi creduti presenti e **864** aggiunti ignoti. L'effetto sulla
pianificazione:

| scena | ammissibili, occlusori creduti | con geometria vera | errore sui punti noti visti per camera |
|---|---|---|---|
| seed 3, misto, 20000 punti | 119 | 120 | medio 33.6, massimo 189 |
| seed 0, `run_scene`, 2000 punti | **96** | 103 | medio 8.6, massimo 33 |

L'errore va soprattutto in una direzione: il pianificatore **sovrastima** ciò
che una camera vede, fino a 189 punti in più, perché le parti aggiunte che non
ha ancora visto mancano come occlusori e crede di vedere attraverso di esse.

## Modello di misura RGB

`measurement.py`. Il sensore è una camera RGB: una vista misura *dove* un
punto cade nell'immagine, non a che distanza sta. L'informazione di Fisher di
un'immagine sulla posizione del punto è

$$F = \frac{1}{\sigma^2} J^\top J, \qquad J = \frac{\partial u}{\partial x} = \begin{pmatrix} f_x / z & 0 & -f_x x / z^2 \\ 0 & f_y / z & -f_y y / z^2 \end{pmatrix} R$$

con $(x, y, z)$ il punto in coordinate camera, $R$ la rotazione
world-to-camera e $\sigma$ = `PIXEL_NOISE` in pixel. È lo Jacobiano esatto
della proiezione, non l'approssimazione $(f / \sigma z)^2 (I - d d^\top)$:
fuori asse i due autovalori non nulli non sono uguali.

$F$ ha **rango 2** e il suo nucleo è esattamente il raggio di vista: una vista
sola non colloca un punto in profondità. Servono due viste da direzioni
diverse, e l'autovalore minimo della somma cresce come il quadrato del seno
dell'angolo fra i raggi. È la triangolazione, ed è ciò che fa cercare al
pianificatore **baseline** invece di viste vicine e frontali — l'opposto di
quello che farebbe con un sensore di profondità.

Verifiche su 500 punti: autovalore minimo 3e-11, $F$ applicata al raggio
1.5e-11, gli altri due autovalori con mediana 71780 e 73036 contro
$(f / z)^2 = 71111$ atteso, e Jacobiano uguale alle differenze finite a 7e-11.
Due viste bastano per il rango pieno.

## Pianificatore

`planner.py`. Sceglie le prossime viste massimizzando l'obiettivo D-ottimo
sui punti aggiunti noti,

$$J(S) = \sum_j \log\det\Big(A_j + \sum_{i \in S} F_{ij}\Big), \qquad A_j = F_{0j} + \frac{I}{(\text{PRIOR} \cdot r)^2}$$

con $F_{0j}$ l'informazione già raccolta sul punto e un prior debole, di
deviazione pari al raggio della nuvola, che rende invertibile un punto visto
da un lato solo. Il log det è submodulare, quindi il **greedy** — una vista
alla volta, quella con il guadagno marginale più alto — è entro il fattore
$1 - 1/e$ dall'ottimo. È l'approccio di FisherRF; il rilassamento convesso
del progetto MATLAB è stato abbandonato.

Il pianificatore vede la scena **solo attraverso la nuvola unita**: quali
punti vedrebbe una candidata lo predice sugli occlusori creduti, mai sulla
scena vera. L'informazione passata invece viene da ciò che le viste hanno
davvero misurato, cioè la `information` che `observe` accumula: il passato è
osservato, il futuro è predetto.

Affina **solo il cambiamento già noto**. Le parti mai viste non sono
bersagli; però una vista pianificata che le inquadra le rivela, e nel round
successivo diventano bersagli. È un'esplorazione per effetto collaterale,
non cercata. Quella esplicita, con i punti fantasma sulla frontiera, non c'è
ancora.

### Dal conteggio di viste al percorso

La prima versione sceglieva `k` viste alla volta con costo pari al loro
numero, a orizzonte recedente. È stata sostituita dal pianificatore a
**tratti**, sotto; i risultati delle due tabelle seguenti sono di quella
prima versione.

### Cosa fa, sulla scena di `run_scene`

| round | viste | aggiunti noti | rimossi confermati | deviazione peggiore mediana |
|---|---|---|---|---|
| sweep | 0, 1, 7, 119 | 244 / 258 | 103 / 207 | 1.0000 R |
| 1 | 88, 25 | 256 / 258 | 122 / 207 | 0.0025 R |
| 2 | 104, 106 | 257 / 258 | 135 / 207 | 0.0019 R |
| 3 | 96, 93 | 257 / 258 | 138 / 207 | 0.0017 R |
| 4 | 112, 38 | 257 / 258 | 142 / 207 | 0.0016 R |

Dopo lo sweep la deviazione mediana lungo la direzione peggiore è il prior:
la maggior parte dei punti è stata vista da una vista sola e non ha
profondità, perché le camere del farthest point sampling sono troppo lontane
fra loro per sovrapporsi. Il primo round la porta a 0.0025 R: il
pianificatore sceglie subito viste che triangolano. È il comportamento
previsto per una camera RGB, e conferma che uno sweep RGB andrebbe fatto a
coppie di viste vicine.

### Contro le baseline

Scena mista (seed 3, 8000 punti di sfondo, 1133 aggiunti veri), sweep di 4,
poi 8 viste a orizzonte 2, misurato sulla verità:

| politica | cambiamento noto | deviazione peggiore mediana | ben triangolati (< 1% R) |
|---|---|---|---|
| solo sweep | 68.7% | 1.0000 R | 16.4% |
| **greedy** | **98.4%** | **0.0021 R** | **88.1%** |
| dispersione (FPS che continua) | 94.9% | 0.0023 R | 80.4% |
| casuale, media di 3 | 94.1% | 0.0026 R | 78.3% |

Il greedy vince su tutte e tre le misure, anche sulla scoperta che non
cerca. Ma i margini sono modesti — otto punti di ben triangolati sulla
dispersione — ed è **una scena sola**: è un primo segnale, non un risultato.
Per scriverlo servono molte scene, più budget, e intervalli di confidenza.

### Costo

Questi tempi sono di prima dell'ottimizzazione del depth buffer: il greedy
impiegava 96 s contro i 2 s della dispersione, quasi tutti nella visibilità.
Restano due ottimizzazioni possibili se servirà: la visibilità cambia poco
fra un round e l'altro, quindi si potrebbe riusare, e il greedy pigro
sfrutta la submodularità per non ricalcolare tutti i guadagni.

## Traiettoria

### Modello di costo

Il costo non è più il numero di viste ma il **tempo**: il viaggio lungo il
percorso, più un tempo fisso per ogni immagine, perché la camera rallenta o
si ferma per scattare. I budget sono due, **tempo** e **numero di
immagini**. Le pose attraversate fra una sosta e l'altra sono quindi
economiche: costano solo il tempo della foto, nessuna deviazione.

### Il moto, con il braccio astratto

`trajectory.py`. La camera resta sul locus, rivolta al centro, e si muove
**interpolando azimut ed elevazione** invece che lungo un grande cerchio:
somiglia a una base che ruota e a una spalla che si alza, e non esce mai
dalla banda di elevazione, cosa che un grande cerchio fra due pose basse su
lati opposti potrebbe fare. L'azimut gira dalla parte corta. Una posa del
reticolo è **attraversata** se il percorso le passa a meno di mezza
spaziatura del reticolo (`PASSING`).

Verifiche su 300 tratti a caso: estremi esatti, raggio costante, **nessuna
uscita dalla banda**; lunghezza in mediana 1.03 volte il grande cerchio,
fino a 1.61 per i tratti che passano vicino alla sommità, dove interpolare
l'azimut allunga la strada. Un tratto attraversa in mediana 7 pose del
reticolo, estremi compresi.

### Il pianificatore a tratti

`plan_leg`. Dalla posa corrente sceglie la prossima **sosta** e quali pose
attraversate **fotografare**, massimizzando l'efficienza del tratto,
guadagno di informazione diviso costo, dentro entrambi i budget residui. Il
costo è **normalizzato sui budget residui**: la frazione di tempo che il
tratto consuma più la frazione di immagini, ciascuna rispetto a quanto ne
resta,

$$c = \frac{t_{\text{viaggio}} + n \cdot t_{\text{foto}}}{T_{\text{residuo}}} + \frac{n}{I_{\text{residue}}}$$

così il budget più scarso pesa di più: con tempo da parte il pianificatore
massimizza il guadagno per immagine e il viaggio diventa quasi gratis, con
immagini da parte il guadagno per secondo. La
sosta si fotografa sempre; una posa attraversata solo se alza l'efficienza
del tratto. Le pose già prese non possono essere soste — tornarci apposta è
uno spreco — ma si possono rifotografare passandoci sopra.

Con un costo che dipende dall'ordine delle soste il greedy sull'efficienza
**perde la garanzia** $1 - 1/e$ del caso a cardinalità: qui è un'euristica.
Il braccio parte dall'ultima posa dello sweep; lo sweep non consuma budget.

### Cosa fa

**Con il costo in secondi** (prima versione), sulla scena di `run_scene`:
dodici tratti da una foto sola, ciascuno al vicino più prossimo sul reticolo
— circa 14°, meno di un secondo di viaggio, zero pose intermedie.

**Con il costo normalizzato**, tre regimi di budget, velocità 1 e 2 s per
foto:

| scena | regime | tratti | foto per tratto | viaggio medio | ben triangolati |
|---|---|---|---|---|---|
| `run_scene`, 2000 punti | equilibrato, 40 s e 12 immagini | 11 | 1.00 | 1.41 | 99.2% |
| | tempo abbondante, 200 s e 12 | 12 | 1.00 | 2.20 | 99.2% |
| | immagini abbondanti, 40 s e 50 | 13 | 1.00 | 1.07 | 99.2% |
| mista, 8000 punti | equilibrato | 10 | 1.00 | 1.77 | 85.2% |
| | **tempo abbondante** | 12 | 1.00 | **3.33** | **90.8%** |
| | immagini abbondanti | 11 | 1.00 | 1.44 | 84.0% |

Sulle soste la normalizzazione fa ciò che deve: col tempo abbondante il
pianificatore viaggia il doppio, perché il viaggio costa poco, e triangola
meglio. Ma in tutti i casi fa **una foto per tratto**.

### Perché le foto di passaggio non si usano

Non è un bug. Valutando ogni sosta con e senza le foto di passaggio, sulla
scena mista: la sosta scelta rende 4060, una foto di passaggio in media 1250,
un terzo. La sosta è la posa **migliore** fra 115, quelle di passaggio sono
pose qualsiasi che capitano sulla strada. Finché una foto al volo costa
quanto una sosta, spendere un'immagine lì abbassa l'efficienza, e la regola
la scarta correttamente. La prima ipotesi — che il pianificatore strisciasse
perché il vicino è il modo più economico di fare una foto — valeva solo per
la scena di `run_scene`: su quella mista le soste migliori sono tratti lunghi
4–5 unità.

### Sosta e foto al volo

Il pezzo che manca è nel modello: la camera **si ferma o rallenta**, e
rallentare costa meno che fermarsi. Prova con un pianificatore a parte,
identico tranne che per due tempi distinti — 2 s la sosta, variabile la foto
al volo. A 2 s riproduce esattamente i numeri del codice vero, quindi misura
la stessa cosa.

| regime | foto al volo | foto per tratto | immagini | deviazione mediana | ben triangolati |
|---|---|---|---|---|---|
| equilibrato | 2.0 s | 1.00 | 10 | 0.0019 R | 85.2% |
| | 0.5 s | 1.33 | 12 | 0.0017 R | 84.7% |
| | 0.2 s | 1.22 | 11 | 0.0020 R | 92.6% |
| immagini abbondanti | 2.0 s | 1.00 | 11 | 0.0017 R | 84.0% |
| | 0.5 s | 2.00 | 16 | 0.0016 R | 88.9% |
| | **0.2 s** | **3.00** | **21** | **0.0013 R** | **91.4%** |

Quando il limite è il **tempo**, le foto al volo pagano chiaramente: tre per
tratto, il doppio delle immagini nello stesso tempo, dall'84% al 91% di ben
triangolati. Quando il limite sono le **immagini** cambia poco e in modo non
netto, com'è ragionevole: se le immagini sono contate conviene spenderle
sulle pose migliori. Una scena sola: è un'indicazione. I due tempi distinti
non sono ancora nel codice.

### Costo di calcolo

Un tratto costa 0.6 s sulla scena di `run_scene` (2000 punti) e 1.3 s su una
da 20000; una missione di dodici tratti 8 s. Prima dell'ottimizzazione del
depth buffer, descritta nel filtro di visibilità, la stessa missione
impiegava tre minuti, con esattamente gli stessi tratti. Il grosso resta la
visibilità: a ogni tratto si ricalcola quella di tutte le 120 pose sugli
occlusori creduti.

## Visualizzazione

`plot_scene.py` disegna scena, regione cambiata e locus. Restituisce la
figura senza mostrarla, così il chiamante decide: `run_scene.py` fa
`plt.show()`, e un domani `make_figures.py` può passarla a `viz.save_fig()`
per farla finire in `5_report/figures/auto/` e renderla rigenerabile da
`rebuild.sh`. Se il disegno stesse dentro lo script sarebbe da
copiaincollare, e da quel momento avresti due versioni della stessa figura
che divergono.

### Perché i frustum e non punti con una freccia

Un frustum mostra tre cose in una: dove sta la camera, dove guarda, e quanto
campo vede. I punti con una freccia mostrano le prime due, e il campo di
vista — che è ciò che determina quali parti della scena una posa può
osservare — resta invisibile.

I corner si ottengono retroproiettando gli angoli dell'immagine alla
profondità del frustum,

$$x_{\text{cam}} = \left(\frac{u - c_x}{f_x} d, \; \frac{v - c_y}{f_y} d, \; d\right)$$

e portandoli nel mondo con $x_{\text{world}} = R^\top x_{\text{cam}} + c$.
La dimensione dell'immagine si assume pari al doppio del punto principale in
entrambe le direzioni, cosa esatta per una calibrazione centrata.

Verifica: i quattro corner di ogni frustum, riproiettati attraverso `K`,
cadono esattamente sugli angoli dell'immagine — (0,0), (500,0), (500,500),
(0,500) con `resolution=500`. Valida insieme la geometria del frustum e la
trasformazione camera-mondo.

Nelle tre proiezioni ortogonali le camere sono disegnate come soli centri:
i frustum coprirebbero proprio la struttura che le proiezioni servono a
mostrare. In compenso è lì che si vede la banda di elevazione, come una
calotta vuota sotto la scena.

### Limite noto

Con il locus nella stessa figura, i limiti degli assi si allargano al raggio
del locus e l'oggetto diventa piccolo nella vista 3D. È il prezzo di
mostrare entrambi nello stesso riquadro: le proiezioni restano leggibili,
ed è lì che si giudica la forma.

## Filtro di visibilità

`visibility_filter.py` dice quali punti una camera vede e quali le nasconde
una superficie più vicina. Richiede **solo i punti**: niente mesh, niente
ricostruzione, niente normali fornite da fuori, perché lo stesso codice deve
girare su una scena sintetica e su una scansione reale.

### Il metodo

Depth buffer con impronte. I punti si proiettano attraverso `K` e ciascuno
copre un disco invece di un pixel solo: una nuvola è un insieme di campioni,
e senza impronta lo sfondo trapela dai buchi fra l'uno e l'altro. Un punto è
visto quando niente di più vicino della sua profondità, oltre lo spessore
della sua toppa, copre il suo pixel.

La scelta rispetto ad HPR (Katz, Tal e Basri, SIGGRAPH 2007) è semantica: la
domanda "quali punti sono coperti rispetto alla camera" *è* la domanda di un
depth buffer, quindi questo è il modello di sensore e non un'approssimazione
geometrica. In più costa O(N) per camera invece di un inviluppo convesso, e
il suo parametro si stima dal dato invece di essere tarato a mano.

### Ogni lunghezza viene dalla spaziatura locale

Niente parametri di scena. La spaziatura $s_j$ è la distanza mediana di ogni
punto ai suoi $k$ vicini — mediana e non media, perché ignora l'outlier
isolato. Da lì:

- **Impronta**: un campione rappresenta una toppa larga quanto la sua
  distanza dai vicini, quindi il raggio in pixel è
  $\text{FOOTPRINT} \cdot f s_j / z_j$.
- **Tolleranza di profondità**: una toppa inclinata copre un intervallo di
  profondità, e un punto non deve essere nascosto dalla propria superficie.
  Lo spessore è $s_j / \cos\theta$, con $\theta$ l'incidenza ricavata dalla
  normale, e `GRAZING` mette un pavimento sul coseno perché una toppa vista
  di taglio non pretenda spessore infinito.
- **Normali**: stimate dalla nuvola stessa con una PCA locale sugli stessi
  $k$ vicini — autovettore dell'autovalore minore. Non servono in input.

Resta $k$, che è adimensionale e non dipende dalla scala: è una scelta di
robustezza statistica, non una taratura.

### Le due costanti adimensionali, e come sono state fissate

`FOOTPRINT` e `MARGIN` non si possono dedurre: vanno calibrate. Ma una volta
sola, contro la verità, non scena per scena. Sulle sintetiche la verità
esiste ed è esatta — l'oggetto è l'insieme di livello di un campo analitico,
quindi si marcia lungo il raggio dalla camera a ogni punto e si guarda se
entra nell'oggetto.

| footprint | margin | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|---|
| 0.35 | 0.5 | 0.893 | 9.96% | 0.78% |
| 0.50 | 0.5 | 0.952 | 3.60% | 1.22% |
| **0.70** | **1.0** | **0.965** | **2.42%** | **1.08%** |
| 1.00 | 1.0 | 0.965 | 1.75% | 1.75% |
| 1.40 | 0.5 | 0.921 | 0.38% | 7.54% |

Il compromesso si legge bene: impronta troppo piccola e lo sfondo trapela
(falsi visibili al 10%), troppo grande e i punti si occludono a vicenda
(falsi occlusi al 7.5%). L'ottimo è piatto fra 0.5 e 1.0, quindi la scelta è
robusta e non su misura.

### La taratura trasferisce?

Congelati i due valori, riverificati su scene **mai usate per tararli** e su
camere diverse:

| punti | fill | decay | seed | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|---|---|---|
| 20000 | 0.30 | 1.6 | 3 | 0.983 | 0.47% | 1.22% |
| 20000 | 0.15 | 1.3 | 5 | 0.988 | 0.67% | 0.55% |
| 20000 | 0.45 | 2.0 | 11 | 0.989 | 0.55% | 0.58% |
| **6000** | 0.30 | 1.6 | 7 | **0.960** | **1.73%** | **2.29%** |
| 50000 | 0.30 | 1.6 | 2 | 0.987 | 0.86% | 0.42% |

Due letture. La prima: `fill` da 0.15 a 0.45 e `decay` da 1.3 a 2.0 non
spostano nulla, quindi la taratura non è su misura per una forma.

La seconda conta di più: **l'unica variabile che muove l'errore è la
densità**. La nuvola rada a 6000 punti è il caso peggiore, con entrambi i
tipi di errore circa raddoppiati, mentre da 20000 a 50000 non cambia niente.
È il comportamento atteso — meno punti, spaziatura maggiore, impronta più
grossolana — ed è anche la ragione per cui ci si può fidare sul reale: una
scansione vera è rada proprio dove il range è lungo e l'incidenza radente, ed
è lì che il filtro degrada per primo. Su un fattore 8 di densità l'accordo
scende solo dal 98.7% al 96.0%.

Attenzione a non confrontare questi numeri con il 96.5% della tabella di
taratura: sono camere diverse, e la varianza fra pose è dello stesso ordine
di quella fra scene.

### Un errore che il metodo ha smascherato da solo

La prima calibrazione dava un tasso di falsi occlusi **bloccato al 15%
qualunque parametro**. Un parametro che non muove l'errore vuol dire che
l'errore non viene da lì: erano punti non occlusi ma **fuori
dall'inquadratura**. Il ray-marching misura solo l'occlusione, il filtro
risponde a "questa camera lo vede", che include il campo di vista. Due
domande diverse confrontate fra loro.

Vale la pena registrarlo perché è il modo tipico in cui una validazione
mente: non sbagliando il numero, ma misurando una cosa diversa da quella che
si crede. L'invarianza rispetto ai parametri è stata il campanello.

Ne è uscita anche una diagnosi sui parametri della scena: con focale 800 e
locus a raggio 2, **l'oggetto sottende 792 pixel su un'immagine di 500**. La
camera non lo inquadra intero. È una condizione legittima per la
pianificazione di viste, ma deve essere una scelta e non una svista.

### Quanto regge quando la nuvola si sporca

Le due perturbazioni che separano una scansione vera da una sintetica pulita,
applicate alla nuvola e rimisurate contro la verità esatta. La verità viene
ricalcolata sulle posizioni **perturbate**: si marcia dalla camera al punto
dov'è finito, non a dov'era.

Rumore gaussiano isotropo, con σ in unità di spaziatura locale:

| σ | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|
| 0.00 | 0.973 | 0.63% | 2.05% |
| 0.25 | 0.964 | 0.76% | 2.82% |
| 0.50 | 0.938 | 1.09% | 5.12% |
| 1.00 | 0.908 | 1.92% | 7.25% |

Outlier sparsi uniformemente nel volume, come frazione dei punti aggiunta
alla nuvola (quindi il totale cresce, e gli spuri sono valutati anch'essi):

| quota | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|
| 0.0% | 0.973 | 0.63% | 2.05% |
| 0.2% | 0.972 | 0.64% | 2.17% |
| 1.0% | 0.935 | 0.62% | 5.86% |
| 3.0% | 0.938 | 0.58% | 5.60% |

**Degrada con grazia.** Nessun crollo: rumore pari a un quarto della
spaziatura costa un punto di accuratezza, e serve rumore pari all'intera
distanza fra punti vicini per scendere al 91%. Gli outlier sotto lo 0.2%
sono innocui; a partire dall'1% costano quattro punti, poi saturano.

**Sbaglia sempre dalla parte giusta.** In entrambi i casi il danno è quasi
tutto in falsi *occlusi*, mentre i falsi visibili restano sotto il 2%. Il
filtro diventa conservativo, cioè dichiara non visto ciò che forse vedrebbe.
Per la pianificazione di viste è il bias corretto: non si sostiene mai di
aver osservato qualcosa che non è stato osservato, e l'informazione stimata
è un limite inferiore.

**Il meccanismo è lo stesso in tutti e tre i casi** — rada, rumorosa,
sporca. Il rumore gonfia le distanze ai k vicini, quindi la spaziatura
stimata cresce, quindi le impronte crescono e i punti si auto-occludono di
più. Un outlier davanti a una superficie occlude ciò che le sta dietro. Una
nuvola rada ha spaziatura maggiore per costruzione. Tutte e tre le strade
portano a impronte troppo grosse, il che spiega perché l'errore vada sempre
nella stessa direzione.

**Ne segue che una funzione di rumore a sé non serve**: genererebbe una
difficoltà che il filtro assorbe già, e di cui ora si conosce il costo.

**Quello che resta non verificato è la parzialità**, e non è verificabile
perturbando. Rumore e outlier peggiorano un numero; una nuvola acquisita da
poche pose cambia la *semantica*, perché "occluso" e "mai catturato"
diventano indistinguibili e il filtro risponde a una domanda diversa da
quella che si crede di porgli. Serve il modello di acquisizione.

### Vedere attraverso: il test di rimozione

`vacated_points` rileva una rimozione proiettando il punto vecchio nella
camera e confrontandolo con il depth buffer della scena **nuova** in quel
pixel. Se la superficie osservata sta più lontano del punto, oltre lo
spessore della sua toppa, la camera vede oltre: lì adesso c'è spazio vuoto.
Un pixel vuoto vale infinito, quindi anche non vedere niente conta come
oltre. Sul reale è il confronto fra la profondità renderizzata dal modello
vecchio e quella osservata.

Usa lo stesso rendering di `visible_points`, spezzato in `_to_image` e
`_render`. La ristrutturazione è **identica al bit** sulla visibilità, quindi
le tarature sopra restano valide.

Validato contro la verità esatta — marcia lungo il raggio nel campo della
scena nuova, e controlla che fra camera e punto ci sia solo spazio libero —
su 1757 punti rimossi e 6 camere: accordo **91.5%**, **falsi svuotati
0.0%**, mancati 8.5%.

Lo stesso bias conservativo della visibilità: non dice mai che qualcosa è
sparito quando c'è ancora, e l'errore sta tutto nei mancati. Solo il 23% dei
mancati sta a meno di due spaziature dalla superficie nuova, quindi non sono
soprattutto punti che toccavano ciò che c'era dietro. La causa è
l'**impronta**: i dischi della superficie in primo piano sbordano oltre la
silhouette e coprono pixel da cui in realtà si vede oltre.

| `FOOTPRINT` | accordo | falsi svuotati | mancati |
|---|---|---|---|
| 0.35 | 0.909 | 3.2% | 5.9% |
| 0.50 | 0.913 | 0.1% | 8.6% |
| 0.70 | 0.882 | 0.0% | 11.8% |
| 1.00 | 0.845 | 0.0% | 15.5% |

Per la rimozione l'ottimo sarebbe 0.5. Resta 0.7, lo stesso della
visibilità: tre punti di accordo non valgono una seconda costante tarata su
una scena sola, e a 0.7 i falsi svuotati sono zero, che è la proprietà che
conta — un pianificatore che crede sparito qualcosa che c'è non ci
guarderebbe più. (Questa tabella è della prima versione del posizionamento;
l'accordo a 0.7 con quella attuale è il 91.5% sopra.)

### Il costo del depth buffer

Il profilo del pianificatore dava il 100% del tempo in `_depth_buffer`, con
837 mila chiamate a `np.minimum.at` per 122 render: quasi settemila passate
a render. La prima versione dipingeva le impronte con una passata per ogni
posizione del disco, dimensionando il ciclo sull'impronta **più grande della
nuvola**. In una nuvola rada qualche punto isolato ha un'impronta di ~40
pixel, e allora ogni render scorreva 83×83 posizioni per tutti.

Due passi, entrambi **identici al bit** sulla visibilità e sul test di
rimozione, verificati sulla scena di `run_scene` e sul riferimento da 20000
punti da cui vengono le tarature sopra — il minimo per pixel non dipende
dall'ordine in cui lo si calcola:

| versione | 120 visibilità + 12 test di rimozione |
|---|---|
| passate dimensionate sull'impronta massima | 12.9 s |
| punti raggruppati per estensione dell'impronta, una sola riduzione | 4.2 s |
| ciclo sui punti e sul loro disco, compilato con numba | **0.56 s** |

Il limite rimasto è fisico: con 2000 punti l'impronta mediana è di 15 pixel
e un render scrive 1.7 milioni di pixel per 2120 punti, sette volte
l'immagine, perché serve a sigillare la superficie. Il ciclo compilato non
costruisce array di indici per quei milioni di scritture, ed è lì che stava
il costo della versione numpy.

`numba` è quindi una dipendenza. È già installato nel Python di sistema, e il
venv lo eredita con `--system-site-packages`.

### Cosa rompe sul reale e non sul sintetico

- **Outlier.** Un punto isolato nel vuoto occlude tutto ciò che sta dietro.
  Va rimosso a monte, con la stessa statistica kNN, ma è un passo dichiarato
  perché **modifica il dato**, non qualcosa da nascondere nel filtro.
- **Cuciture di registrazione.** Scansioni allineate imperfettamente lasciano
  due copie della stessa superficie a qualche millimetro: l'una occlude
  l'altra e produce auto-occlusione fantasma.
- **"Occluso" è relativo a ciò che è stato acquisito.** Un punto mai
  catturato non è occluso, è assente. La funzione risponde a una domanda più
  stretta di quanto il nome suggerisca.
- **La calibrazione trasferisce solo se le sintetiche somigliano alle reali**
  nelle proprietà che contano: densità variabile, rumore, parzialità. Oggi le
  sintetiche sono complete, prive di rumore e uniformi, quindi la taratura
  sopra è ottimistica. Renderla onesta richiede il modello di acquisizione.

### Costanti

| costante | valore | ruolo |
|---|---|---|
| `NEIGHBOURS` | 12 | Vicini su cui si stimano spaziatura e normale. |
| `FOOTPRINT` | 0.7 | Raggio dell'impronta, in spaziature locali. |
| `GRAZING` | 0.1 | Pavimento sul coseno di incidenza. |
| `MARGIN` | 1.0 | Tolleranza, in spessori di toppa. |
| `NEAR` | 1e-6 | Profondità sotto cui il punto è dietro la camera. |

### Visualizzatore interattivo

`scene_viewer.py` logga su **rerun**, che mostra tutto in una finestra sola.
Gli stadi sono passi di una **timeline** che si scorre ruotando liberamente:

1. `scene/old` — la scena come la conosce il modello obsoleto;
2. `scene/new` — la scena com'è adesso, nello stesso colore, così dove non
   cambia nulla i due strati coincidono;
3. `scene/added` in rosso sopra la nuova, `scene/removed` in viola sopra la
   vecchia;
4. `cameras/candidate` — tutto il locus, in grigio;
5. un passo per round: `cameras/sweep` in blu al primo, poi un passo per
   tratto con `cameras/planned` in verde e il percorso della camera in
   `trajectory/leg_NN`, e a ogni round `known/added` in giallo e
   `known/removed` in azzurro sopra il cambiamento vero, con la conoscenza
   di quel momento. Scorrendo la timeline si vede la conoscenza crescere; il
   rosso e il viola che restano visibili alla fine sono cambiamento che
   nessuno ha osservato.

Spegnendo `scene/old` resta il mondo com'è adesso.

L'ordine è quello causale: le ammissibili vengono dopo lo sweep perché sono
scelte a partire da esso.

Niente viene mai cancellato: ogni passo aggiunge, e ciascuna entità si spegne
a mano dal pannello laterale. Una prima versione cancellava le candidate al
passo delle ammissibili, ma il viewer apre sull'ultimo istante della
timeline, dove a quel punto non esistevano più.

Le pose dello sweep e quelle pianificate sono sottoinsiemi del locus, quindi i loro frustum
coincidono con quelli grigi nella stessa posa. Due geometrie identiche
lasciano al renderer la scelta, e non sceglie in modo affidabile l'ultima:
per questo sono disegnati con una linea più spessa e un piano immagine più
lungo del 2%, quanto basta a vincere il pareggio senza vedersi.

### Tre trappole di rerun

Nessuna delle tre stava nel codice della scena: erano convenzioni di rerun
date per scontate invece che verificate.

- **`image_plane_distance` vale 1.0 unità di mondo per default.** Con il
  locus a raggio 3 e la scena a raggio 1, 120 frustum lunghi un'unità
  convergono sull'oggetto e lo seppelliscono. Ora è il 10% del raggio del
  locus.
- **Il layout è memorizzato per `application_id`.** Cambiando la struttura
  delle entità con lo stesso nome, un'entità nuova può non comparire nel
  pannello laterale finché il layout non viene reimpostato a mano. Per
  questo `view_scene` manda un blueprint esplicito a ogni registrazione.
- **Il segno di un raggio ne cambia l'unità.** Positivo è in unità di
  scena, negativo in punti dell'interfaccia. Uno spessore di linea `2.5`
  significava linee più grosse della scena intera; quello voluto è `-2.5`,
  cioè 2.5 pixel.

Il vantaggio concreto sul matplotlib: `Pinhole` è una primitiva nativa, quindi
il frustum lo disegna il viewer partendo da `K`, senza ricalcolare gli otto
spigoli a mano. E mplot3d non fa depth buffering per frammento ma ordina per
artista, quindi i punti dietro trapelano davanti e la nuvola sembra una
foschia; rerun ha un depth buffer vero.

**Convenzioni, verificate e non assunte.** Rerun vuole la trasformata
dell'entità nello spazio del padre, cioè **world-from-camera**, che è la
trasposta della rotazione world-to-camera prodotta dal locus. Il `Pinhole`
assume RDF — x a destra, y in basso, z in avanti — la stessa OpenCV usata
ovunque qui. Controllo: `rotation.T @ [0,0,1]` coincide con la direzione
verso l'origine a errore **0.00e+00**, quindi la trasposta è quella giusta e
non serve alcuna conversione.

**Versioni.** L'SDK Python e il viewer devono combaciare o il viewer rifiuta
la connessione. Qui il viewer è lo snap 0.33.1, quindi serve
`rerun-sdk==0.33.1` e non l'ultima. Sta in `.venv/` alla radice del repo,
creato con `--system-site-packages` per ereditare numpy, scipy e matplotlib
già presenti nel Python di sistema: `rerun-sdk` è l'unica cosa aggiunta.
Da cui il modo di lanciare, `../../.venv/bin/python run_scene.py`.

Nota d'API: `set_time_sequence` non esiste più in 0.33, si usa
`set_time(timeline, sequence=n)`.

**`save_path`** scrive un `.rrd` invece di aprire la finestra. Serve a
guardare una scena più tardi, a passarla a qualcun altro, e a far girare il
codice dove non c'è un display.

## Limiti noti

- **Manici veri, cioè buchi passanti, sono rari** ai parametri di default. Le
  forme sono blob lobati con fessure profonde e a volte componenti staccate,
  non oggetti topologicamente complessi. Si ottengono abbassando insieme
  `fill` e `decay`, ma non sono garantiti a ogni seed.
- **La detection è un proxy perfetto su ciò che si vede**: niente nuvola
  grezza, niente rumore di ricostruzione, niente canale fotometrico. È
  un'ipotesi di scope, da togliere quando si simulerà la detection vera.
- **Il modello obsoleto è una scansione completa.** La scena vecchia contiene
  tutta la superficie, anche quella che in realtà nessuna vista aveva mai
  raggiunto. Serve il modello di acquisizione per renderla parziale.
- **Il cambiamento è solo geometrico**: un cambio d'aspetto su superficie
  che resta al suo posto non viene rilevato.
- **Uno spostamento perde l'identità dell'oggetto**, vedi sopra.
- **Gli oggetti possono sporgere dalla sfera della scena** fino al proprio
  raggio.
- **Il generatore non produce l'ignoto.** Una nuvola ignota non è una cosa
  che si genera, è il risultato di un'acquisizione: serve un modello di
  sensore che da questa superficie e da alcune pose produca la nuvola
  parziale, rumorosa e bucata che il pianificatore vede davvero. Questo
  modulo fornisce solo il ground truth.
