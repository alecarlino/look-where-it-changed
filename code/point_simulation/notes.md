# Note tecniche — scena sintetica

Documentazione di `synthetic_point_cloud.py`, `camera_locus.py`,
`visibility_filter.py`, `detection_sweep.py`, `knowledge.py`,
`measurement.py`, `planner.py`, `trajectory.py` e `scene_viewer.py`. Il codice
porta solo commenti da una riga; il ragionamento, le formule e la taratura
stanno qui.

## Perché una superficie e non un volume

Un sensore di profondità misura superfici: non restituisce mai punti
all'interno di un solido. La scelta non è estetica, decide se il resto della
pipeline ha senso. Occlusione, angolo di incidenza e confine fra noto e
ignoto sono proprietà che esistono solo se la scena ha superfici e c'è un
punto di vista. Con punti sparsi in un volume nessuna di queste è definita.

## Perché l'insieme di escursione di un processo gaussiano

L'alternativa era un catalogo di forme parametriche (sfera, toro, C, L) con
parametri casuali. Scartata: le forme restano quelle che hai nominato, e la
varietà è apparente. Un tentativo intermedio, anch'esso scartato, spazzava
una sezione circolare lungo una curva casuale di Fourier: produceva sempre
un filo ingrossato, topologicamente un tubo, mai un oggetto con massa.

L'insieme di escursione di un campo casuale — la regione dove supera una
soglia — dà invece gratuitamente, dalla casualità e non dalla costruzione:

- **topologia arbitraria** — manici, cavità, componenti staccate;
- **concavità e auto-occlusione**, che è il motivo per cui la selezione di
  viste non è banale: su un corpo convesso ogni vista vede la silhouette e il
  problema si risolve da solo;
- **due parametri con un nome standard**, la lunghezza di scala e la
  rugosità;
- **riproducibilità** dal seed.

Ogni pezzo del generatore è una **costruzione pubblicata**, per poterlo
citare invece di doverlo giustificare: random Fourier features (Rahimi e
Recht, 2007), kernel Matérn, sua densità spettrale e funzione media non nulla
(Rasmussen e Williams, *Gaussian Processes for Machine Learning*, 2006, §2.7
e §4.2), insiemi di escursione (Adler e Taylor, *Random Fields and
Geometry*, 2007), intersezione come minimo (Ricci, 1973), proiezione sulla
superficie lungo il gradiente (Witkin e Heckbert, 1994), distanza al primo
ordine (Taubin, 1991). Un generatore precedente, costruito ad hoc, è
descritto in fondo a questa parte.

## Il campo: random Fourier features con kernel Matérn

$$f(x) = \sqrt{\tfrac{2}{K}} \sum_{k=1}^{K} \cos(\omega_k \cdot x + b_k), \qquad b_k \sim U(0, 2\pi), \qquad \omega_k \sim p(\omega)$$

con gradiente analitico, dalla stessa proiezione del valore,

$$\nabla f(x) = -\sqrt{\tfrac{2}{K}} \sum_{k=1}^{K} \sin(\omega_k \cdot x + b_k)\, \omega_k .$$

Per il teorema di Bochner, se le frequenze $\omega_k$ si estraggono dalla
densità spettrale $p(\omega)$ di un kernel stazionario, $f$ è un campione
approssimato di un processo gaussiano con quel kernel, a media nulla e
**varianza unitaria per costruzione** (Rahimi e Recht, 2007). L'ampiezza è la
stessa per tutti i termini: lo spettro lo decide la distribuzione delle
frequenze, non il loro peso.

Il kernel è il **Matérn**, con due parametri:

- la lunghezza di scala $\ell$ = `length` × R, la dimensione delle feature;
- la rugosità $\nu$ = `smoothness`: più è bassa, più la superficie è
  irregolare.

La sua densità spettrale in 3D è proporzionale a
$(2\nu/\ell^2 + |\omega|^2)^{-(\nu + 3/2)}$, cioè una **t di Student
multivariata** con $2\nu$ gradi di libertà e scala $1/\ell$, che si estrae
in una riga:

$$\omega_k = \frac{z_k}{\ell} \sqrt{\frac{2\nu}{g_k}}, \qquad z_k \sim \mathcal{N}(0, I), \quad g_k \sim \chi^2_{2\nu} .$$

$K$ = `n_modes` = 256: più termini rendono la distribuzione marginale più
vicina alla gaussiana, a costo lineare.

## La media che scende verso il bordo

$$g(x) = f(x) + m(x) - u, \qquad m(x) = -\text{MEAN\_DROP}\,\frac{|x|^2}{R^2}$$

Un processo gaussiano può avere una funzione media non nulla (Rasmussen e
Williams, §2.7). Con una media che scende verso il bordo il campo sta sopra
soglia soprattutto al centro, e l'oggetto si concentra dentro la scena.

Senza media la costruzione è corretta ma **non realistica**. Un campo
stazionario tocca il bordo con probabilità pari a `fill` ovunque sulla
sfera, quindi l'intersezione con la sfera lo taglia estesamente e la scena
diventa una **sfera intagliata**: con `length` 0.3, dal 25 al 34% della
superficie era calotta sferica, cioè una superficie grande, convessa e
visibile da tutto il locus. Al crescere della caduta della media (prova a
parte, con $\nu$ = 1.5, 5 semi):

| `MEAN_DROP` | calotta | componenti oltre il 2% dei punti | `fill` misurato |
|---|---|---|---|
| 0 | 27.1% | 2.0 | 0.328 ± 0.078 |
| 2 | 9.0% | 1.6 | 0.318 ± 0.071 |
| **4** | **0.8%** | **1.0** | **0.309 ± 0.055** |
| 6 | 0.0% | 1.0 | 0.308 ± 0.040 |

A 4 le calotte sono ridotte a pochi punti e resta un oggetto solo; oltre, il
campo si comprime verso una palla. L'intersezione con la sfera resta come
rete di sicurezza per la rara parte che ci arriva.

## Il livello dalla frazione attesa

Il campo è gaussiano a varianza unitaria, quindi un punto a frazione $r$ del
raggio sta dentro l'oggetto con probabilità $1 - \Phi(u - m(r))$. Il livello
$u$ è quello per cui la media di questa probabilità sul volume della palla,
pesato da $3r^2\,dr$, vale `fill`:

$$\text{fill} = 3 \int_0^1 \big(1 - \Phi(u - m(r))\big)\, r^2\, dr$$

risolto numericamente per $u$ (`quad` e `brentq`). È la parametrizzazione
degli insiemi di escursione per soglia (Adler e Taylor, 2007).

`fill` è quindi la frazione **attesa**, non quella di ogni scena. La palla
contiene solo circa tre lunghezze di correlazione, quindi una singola
realizzazione oscilla molto: chiedendo 0.30 si misura **0.318 ± 0.083** su 20
campi, 0.170 ± 0.067 chiedendo 0.15 e 0.465 ± 0.088 chiedendo 0.45 — un
leggero eccesso sistematico di circa 0.02. Il generatore precedente, che
fissava il livello col quantile empirico, dava il `fill` esatto per ogni
scena; qui lo si è scambiato con una formula citabile e con più varietà fra
scene.

## La chiusura: intersezione con una sfera

$$g_{\text{scena}}(x) = \min\big(g(x),\; R - |x|\big)$$

L'intersezione di due solidi è il minimo dei loro campi (Ricci, 1973), la
stessa operazione, duale, dell'unione come massimo usata per gli oggetti.
Dove l'insieme di escursione raggiunge la sfera, è la sfera a chiuderlo, e
l'oggetto ha sempre un interno ben definito — che serve a valle per i test
dentro/fuori nell'occlusione. Con la media che scende la calotta è dello
0–3% della superficie.

## Campionamento: il guscio a spessore costante

La distanza al primo ordine da $x$ alla superficie è (Taubin, 1991)

$$d(x) = \frac{|g(x)|}{|\nabla g(x)|}$$

Dividere il valore del campo per la pendenza trasforma quindi il test in una
distanza vera. Conta per la densità: un guscio a valore di campo costante
sarebbe spesso dove il campo è piatto e sottile dove è ripido, e i punti si
ammasserebbero nelle zone piatte. Con un guscio a spessore costante escono
distribuiti uniformemente sulla superficie — verificato sotto.

I candidati si estraggono in una palla **più larga di un guscio** di quella
della scena. La calotta sta esattamente sulla sfera: estraendo solo dentro,
il guscio attorno a lei esisterebbe da un lato solo, e le calotte
riceverebbero metà dei candidati delle altre superfici, con densità e area
sottostimate della metà.

Il campionamento è per rigetto, a lotti. Il tasso di accettazione non è noto
a priori perché dipende da quanta superficie ha quel campo dentro la scena,
quindi si chiede quello che manca e si ripete. Sovra-campionamento ×8: il
guscio è sottile, la maggior parte dei candidati cade lontano.

## Proiezione di Newton

Newton sull'equazione scalare $g(x) = 0$ (Witkin e Heckbert, 1994). Il
gradiente è una riga sola e non un sistema quadrato, quindi si usa la sua
pseudo-inversa:

$$x \leftarrow x - \frac{g(x)}{|\nabla g(x)|^2}\, \nabla g(x)$$

Il passo si muove lungo la normale all'insieme di livello, che è anche la via
più breve per raggiungerlo. Convergenza quadratica, da cui il numero piccolo
di passi.

Vicino allo **spigolo** dove la superficie del campo incontra la calotta, il
minimo cambia componente da un passo all'altro e Newton può rimbalzare senza
atterrare: in una prima versione un punto finiva a 0.47 R dalla superficie.
Dopo i passi di Newton si tengono quindi solo i punti atterrati, a distanza
al primo ordine sotto `LANDED` × R. Il campionamento è già per rigetto, e lo
spigolo è una linea: scartarli non cambia la densità altrove.

## Parametri e costanti

| parametro | default | ruolo |
|---|---|---|
| `n_modes` | 256 | Termini delle random Fourier features. |
| `length` | 0.30 | Lunghezza di scala del Matérn, in raggi della scena. |
| `smoothness` | 2.5 | Rugosità $\nu$ del Matérn. |
| `fill` | 0.30 | Frazione attesa della scena dentro l'oggetto. |

| costante | valore | ruolo |
|---|---|---|
| `MEAN_DROP` | 4.0 | Caduta della media dal centro al bordo, in deviazioni standard. |
| `SHELL` | 0.02 | Semispessore del guscio di campionamento, in raggi. |
| `STEPS` | 4 | Passi di Newton. |
| `LANDED` | 1e-9 | Distanza dalla superficie che conta come atterrato, in raggi. |
| `MAX_BATCHES` | 200 | Lotti prima di rinunciare. |

- **La lunghezza di scala va confrontata con la scena.** Con feature molto
  più grandi del raggio il campo è quasi costante e ogni oggetto esce come
  un blob convesso: è l'errore che ha prodotto i blob del generatore
  precedente, e vale identico qui per `length`.
- **`SHELL`** sottile significa pochi candidati sopravvissuti, spesso
  significa che Newton deve viaggiare più a lungo e può attraversare in un
  ramo vicino della superficie.
- **`STEPS = 4`** è generoso: dal guscio il residuo arriva alla precisione
  macchina in due o tre passi, spigolo a parte.

### Perché $\nu$ = 2.5

La rugosità scambia difficoltà della scena contro affidabilità del filtro di
visibilità. Filtro a impronta 0.7 e margine 1.0, 3 scene × 4 camere, contro
la visibilità esatta:

| $\nu$ | accordo del filtro | falsi visibili | falsi occlusi | visibile da una vista |
|---|---|---|---|---|
| 0.5 | 84.2% | 4.7% | 11.1% | 27% |
| 1.5 | 90.9% | 3.8% | 5.2% | 35% |
| **2.5** | **93.9%** | **2.8%** | **3.3%** | **37%** |
| 5.0 | 94.3% | 2.6% | 3.1% | 43% |

Più la superficie è rugosa, più la scena è difficile e meno il filtro è
affidabile: pieghe più piccole dell'impronta, e normali stimate peggio dal
piano locale. A 2.5 il filtro torna vicino al 95% che aveva sul generatore
precedente e la difficoltà resta la stessa; oltre, si guadagna quasi niente
in accuratezza e le scene diventano più facili. Anche fisicamente oggetti
reali sono lisci alla scala della distanza fra i punti, e una rugosità come
$\nu$ = 0.5 non rappresenta nulla di plausibile.

## Verifiche

**Il campo ha la covarianza dichiarata.** 4000 campi indipendenti, coppie di
punti a distanza $r$, contro il Matérn $\nu$ = 5/2 analitico
$k(r) = (1 + \sqrt5 r/\ell + 5r^2/3\ell^2)\,e^{-\sqrt5 r/\ell}$ con $\ell$ = 0.3:

| $r$ | 0.0 | 0.1 | 0.2 | 0.3 | 0.5 | 0.8 |
|---|---|---|---|---|---|---|
| misurata | 1.000 | 0.914 | 0.742 | 0.519 | 0.235 | 0.050 |
| analitica | 1.000 | 0.916 | 0.728 | 0.524 | 0.225 | 0.048 |

Varianza 1.037, media 0.004. È la verifica che rende la formula citabile: il
campo è davvero ciò che il riferimento dice.

**I punti stanno sulla superficie.** Residuo massimo sotto 1e-9 su 5 scene da
20000 punti, cioè `LANDED`.

**La densità superficiale è uniforme.** Coefficiente di variazione della
distanza al vicino più prossimo fra 0.523 e 0.531, contro 0.523 di un
processo di Poisson uniforme sul piano. Anche sulle calotte, grazie alla
palla allargata: la spaziatura lì è il 94–99% di quella altrove.

**Un oggetto solo, con poca calotta.** Calotta dallo 0 al 3.1% della
superficie; una sola componente sopra il 2% dei punti, due in un seme su
cinque.

**Nessuna superficie irraggiungibile.** Nessun punto è invisibile da tutte le
120 pose di una sfera intera di camere: niente bolle chiuse dentro l'oggetto,
che sarebbero bersagli che nessuna vista può raggiungere.

**L'auto-occlusione è reale.** Frazione di superficie visibile da una vista,
con la **visibilità esatta** — marcia lungo il raggio nel campo — su 120 pose
a raggio 3: dal 36 al 40% su tre scene. Il generatore precedente, sulle stesse
camere, dava il 41.6%: le scene nuove sono un po' più difficili.

**Una lezione di metodo.** In una prima misura l'auto-occlusione delle scene
nuove era stata stimata con l'operatore HPR (Katz, Tal e Basri, 2007), che
dava il 4–9% di superficie visibile da una vista. Era sbagliato: su queste
superfici l'HPR concorda con la visibilità esatta solo al 60–70%, contro il
92% del filtro di visibilità. L'HPR era stato calibrato su una sfera, dove è
esatto, ma sulle superfici rugose non regge. Da allora il riferimento è solo
la visibilità esatta, e i numeri dell'HPR di quella prima misura sono stati
scartati.

### Il generatore precedente

Il primo generatore era una costruzione **ad hoc**: somma di sinusoidi con
direzioni isotrope, frequenze log-uniformi fra 3/R e 16/R e ampiezze che
decadevano come la frequenza alla potenza `decay`, normalizzata a varianza
unitaria; chiusura con un termine radiale quadratico nullo dentro il 60% del
raggio; livello dal quantile empirico. Funzionava — superfici a residuo
1e-12, oggetti lobati e concavi — ma nessun pezzo aveva un riferimento da
citare, e il termine di chiusura era del tutto inventato. È stato sostituito
dalla costruzione sopra, che a parità di camere dà scene equivalenti e un po'
più auto-occludenti. Tutti i numeri a valle sono stati rifatti sul generatore
nuovo.

## Scena prima e dopo il cambiamento

`generate_scene_pair` produce la stessa scena **due volte**: come la conosce
il modello obsoleto e com'è adesso, con oggetti aggiunti, rimossi o spostati.
È ciò che rende il simulatore capace di mettere alla prova l'argomento per cui
si lavora a punti.

### Perché due scene e non un'etichetta

Una prima versione marcava come cambiati punti
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

La differenza è strutturale e non dipende da un numero: un oggetto vero ha
facce rivolte verso la superficie su cui poggia e un retro, occlude ciò che
gli sta dietro e rivela ciò che copriva, e la sua geometria non è nota prima
di essere vista. Una toppa etichettata non ha nessuna di queste proprietà.
Sul generatore precedente la differenza si vedeva anche nei numeri — con
quattro camere di sweep l'etichettatura dava il 96% del cambiamento visto,
gli oggetti veri il 71% — ma quel confronto non è stato rifatto sul
generatore attuale, perché la versione etichettata non è più usata.

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
| 1 aggiunto | 20000 | 20869 | 889 | 52 |
| 1 rimosso | 20869 | 20000 | 52 | 889 |
| 1 spostato | 20869 | 20881 | 961 | 927 |
| 2 aggiunti, 1 rimosso, 1 spostato | 21755 | 22768 | 3020 | 2060 |

  I 52 sono lo sfondo coperto dall'oggetto: aggiungerlo li rimuove,
  toglierlo li rivela.
- **Densità**: spaziatura mediana 0.0117 sugli oggetti contro 0.0101 su
  tutta la scena. Gli oggetti escono circa il 15% più radi dello sfondo,
  mentre col generatore precedente erano il 7% più fitti. La causa non è
  stata indagata: il numero di punti di un oggetto si fissa da un'area
  stimata su un primo campione, ed è il primo candidato da controllare.

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
pianificatore; il cambiamento vero della coppia di scene è il **ground
truth**. La differenza fra le due è il cambiamento che lo sweep non ha visto,
cioè ciò che l'esplorazione deve ancora trovare. È il primo pezzo del modello
di acquisizione: rende Ω stimata invece che assunta.

### Dispersione, non copertura

"Più lontane possibile" e "copertura massima" sono obiettivi diversi. La
copertura vera è un problema di max-coverage sulla matrice di visibilità,
risolvibile con il greedy a garanzia $(1 - 1/e)$, ma va calcolata su una
geometria: sulla nuvola vera sarebbe barare, perché userebbe informazione
che lo sweep non ha. La versione legittima la calcolerebbe sul modello
obsoleto, che è noto. Non è implementata: la dispersione geometrica è la
scelta più semplice che non guarda la scena, e resta un confronto da fare.

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

Cambiamento misto — 2 aggiunti, 1 rimosso, 1 spostato, sfondo di 20000 punti
— con 3020 punti aggiunti e 2060 rimossi:

| $k$ | aggiunti visti | rimossi visti oltre |
|---|---|---|
| 1 | 50.9% | 70.3% |
| 2 | 67.4% | 74.3% |
| 4 | 83.0% | 80.7% |
| 6 | 90.1% | 84.3% |
| 10 | 97.6% | 85.0% |

**Con poche viste le rimozioni si confermano prima delle aggiunte.** Tolto un
oggetto, *tutti* i suoi punti lungo il raggio — davanti e dietro — sono
spazio vuoto, e una vista sola li attraversa tutti; un'aggiunta invece la si
ricostruisce una faccia alla volta. Con una camera sola si conferma il 70%
delle rimozioni e si vede il 51% delle aggiunte. Con più viste le aggiunte
raggiungono e superano le rimozioni, che si fermano attorno all'85%: la
parte mai confermata è in buona parte sfondo coperto dagli oggetti nuovi, e
attraverso un oggetto non si può vedere.

Una proprietà del farthest point sampling su una banda da tenere presente:
partendo dalla camera più alta, le successive sono le più lontane, e su una
calotta i punti più lontani fra loro stanno sul **cerchio di bordo**, che è
il più largo. Con la banda da −10° a 85° le prime quattro camere hanno
elevazione 82.4°, −9.7°, −5.8°, −9.2°: tre su quattro guardano la scena di
taglio dal bordo basso. Lo sweep è disperso sul locus ma non nelle direzioni
di vista; è un argomento per la copertura sul modello obsoleto invece della
sola dispersione, quando ci sarà.

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
creduto sbaglia in due direzioni opposte, e per costruzione solo in queste
due:

- un oggetto tolto e non confermato **continua a occludere**, quindi il
  pianificatore crede di vedere **meno** di quanto vedrebbe;
- un pezzo aggiunto e non ancora visto **non occlude**, quindi crede di
  vedere **di più**.

Sulla scena mista (seed 3, 20000 punti) dopo uno sweep di 4 camere restano
**397** punti rimossi creduti presenti e **514** aggiunti ignoti. Per ogni
camera del locus, il numero di punti aggiunti noti che il pianificatore crede
di vedere contro quelli che vedrebbe davvero: errore medio 27.1 punti, fino a
**53 in eccesso** e fino a **106 in difetto**. Su questa scena prevale la
sottostima, cioè l'effetto dei rimossi non confermati.

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
**tratti** della sezione sulla traiettoria, e i suoi risultati non sono più
riportati: erano misurati su codice e generatore che non esistono più.

### Cosa fa, sulla scena di `run_scene`

Sfondo di 2000 punti, un oggetto aggiunto, uno rimosso, uno spostato; sweep di
4 camere, poi 40 s e 12 immagini, velocità 1 e 2 s per foto:

| tratto | immagini | aggiunti noti | rimossi confermati | deviazione peggiore mediana |
|---|---|---|---|---|
| sweep | 119, 0, 7, 1 | 170 / 179 | 191 / 255 | 0.0034 R |
| 1 | 56 | 179 / 179 | 199 / 255 | 0.0031 R |
| 2 | 98 | 179 / 179 | 199 / 255 | 0.0023 R |
| 6 | 69 | 179 / 179 | 203 / 255 | 0.0019 R |
| 10 | 103, 111 | 179 / 179 | 207 / 255 | 0.0014 R |
| 11 | 114 | 179 / 179 | 207 / 255 | 0.0013 R |

Undici tratti, 36.8 s e tutte le 12 immagini; uno solo usa una foto di
passaggio. Il primo tratto completa la scoperta delle aggiunte, e da lì il
pianificatore affina la triangolazione, dimezzando la deviazione. Restano 48
punti rimossi creduti presenti: il pianificatore non ha le rimozioni come
bersaglio, e le conferma solo quando capita.

**Quanto triangola lo sweep dipende dalla scena.** Qui dopo le 4 viste dello
sweep la deviazione mediana è già 0.0034 R, cioè la maggior parte dei punti è
vista da due direzioni. Sulla scena mista sotto invece resta il prior,
1.0 R: la maggior parte dei punti è vista da una vista sola e non ha
profondità, perché le camere del farthest point sampling sono troppo lontane
fra loro per sovrapporsi. Uno sweep a coppie di viste vicine renderebbe il
risultato meno dipendente dalla scena.

### Contro le baseline

Scena mista (seed 3, sfondo di 8000 punti, 1188 aggiunti veri), sweep di 4,
poi 40 s e 12 immagini, velocità 1 e 2 s per foto. Le baseline si muovono con
lo stesso modello di moto e fotografano solo la sosta: la **dispersione**
sceglie la posa più lontana da quelle già prese, la **casuale** una posa
qualsiasi. Tutto misurato sulla verità:

| politica | cambiamento noto | deviazione peggiore mediana | ben triangolati (< 1% R) | tempo | immagini |
|---|---|---|---|---|---|
| solo sweep | 82.2% | 1.0000 R | 29.8% | — | — |
| **greedy** | **98.8%** | **0.0019 R** | **93.0%** | 37.6 s | 9 |
| casuale, media di 3 | 95.3% | 0.0027 R | 84.0% | | |
| dispersione | 95.1% | 0.0028 R | 66.0% | 36.9 s | 5 |

Il greedy vince su tutte e tre le misure, anche sulla scoperta che non
cerca: nove punti di ben triangolati sul casuale, ventisette sulla
dispersione. La dispersione paga il modello di costo: le pose più lontane
costano molto viaggio, e nello stesso tempo scatta solo 5 immagini. Il greedy
ne usa 9 su 12 perché finisce prima il tempo. È ancora **una scena sola**:
per scriverlo come risultato servono molte scene e intervalli di confidenza.

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

Tre regimi di budget, velocità 1 e 2 s per foto:

| scena | regime | tratti | foto per tratto | viaggio medio | tempo | immagini | deviazione mediana | ben triangolati |
|---|---|---|---|---|---|---|---|---|
| `run_scene`, 2000 punti | equilibrato, 40 s e 12 immagini | 11 | 1.09 | 1.16 | 36.8 s | 12 | 0.0013 R | 100% |
| | tempo abbondante, 200 s e 12 | 12 | 1.00 | 2.09 | 49.1 s | 12 | 0.0012 R | 100% |
| | immagini abbondanti, 40 s e 50 | 12 | 1.08 | 0.97 | 37.7 s | 13 | 0.0013 R | 99.4% |
| mista, 8000 punti | equilibrato | 8 | 1.12 | 2.46 | 37.6 s | 9 | 0.0019 R | 93.0% |
| | **tempo abbondante** | 11 | 1.09 | **3.63** | 63.9 s | 12 | 0.0016 R | **97.3%** |
| | immagini abbondanti | 8 | 1.25 | 2.22 | 37.7 s | 10 | 0.0018 R | 94.6% |

La normalizzazione fa ciò che deve: col tempo abbondante il pianificatore
viaggia di più, perché il viaggio costa poco, e sulla scena mista triangola
meglio. Con 200 s non li usa tutti: finiscono prima le immagini. Sulla scena
di `run_scene` i regimi quasi non si distinguono, perché è già satura: tutto
il cambiamento è ben triangolato in ogni caso.

### Perché le foto di passaggio si usano poco

Fra 1.00 e 1.25 foto per tratto: quasi sempre solo la sosta. Non è un bug.
La sosta è la posa **migliore** fra tutte quelle raggiungibili, mentre quelle
di passaggio sono pose qualsiasi che capitano sulla strada. Finché una foto
al volo costa quanto una sosta, spendere un'immagine lì di solito abbassa
l'efficienza del tratto, e la regola la scarta correttamente. La prima
ipotesi, che il pianificatore facesse tratti brevi perché il vicino è il modo
più economico di fare una foto, è smentita dai numeri: sulla scena mista il
viaggio medio è 2.2–3.6 unità, cioè soste lontane.

### Sosta e foto al volo

Il pezzo che manca è nel modello: la camera **si ferma o rallenta**, e
rallentare costa meno che fermarsi. Prova con un pianificatore a parte,
identico tranne che per due tempi distinti: 2 s la sosta, variabile la foto
al volo. A 2 s riproduce esattamente i numeri del codice vero, quindi misura
la stessa cosa. Scena mista:

| regime | foto al volo | foto per tratto | immagini | deviazione mediana | ben triangolati |
|---|---|---|---|---|---|
| equilibrato | 2.0 s | 1.12 | 9 | 0.0019 R | 93.0% |
| | 0.5 s | 1.50 | 12 | 0.0017 R | 95.8% |
| | 0.2 s | 1.50 | 12 | 0.0017 R | 95.8% |
| immagini abbondanti | 2.0 s | 1.25 | 10 | 0.0018 R | 94.6% |
| | 0.5 s | 2.57 | 18 | 0.0015 R | 95.3% |
| | **0.2 s** | **3.00** | **24** | **0.0013 R** | **98.1%** |

Quando il limite è il **tempo** le foto al volo pagano chiaramente: tre per
tratto, più del doppio delle immagini nello stesso tempo, dal 94.6% al 98.1%
di ben triangolati. Nel regime equilibrato aiutano finché ci sono immagini:
a 0.5 s il budget di 12 è già esaurito, e scendere a 0.2 s non cambia nulla.
Una scena sola: è un'indicazione. I due tempi distinti non sono ancora nel
codice.

### Costo di calcolo

La missione di `run_scene`, sweep e undici tratti, impiega circa 9 s. Prima
dell'ottimizzazione del depth buffer, descritta nel filtro di visibilità, una
missione simile impiegava tre minuti, con esattamente gli stessi tratti. Il
grosso resta la visibilità: a ogni tratto si ricalcola quella di tutte le 120
pose sugli occlusori creduti. Se servirà, la si può riusare fra un tratto e
l'altro, perché cambia poco, e il greedy pigro può sfruttare la
submodularità per non ricalcolare tutti i guadagni.

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

Griglia su 3 scene (seed 3, 5, 11, 8000 punti) e 4 camere della banda,
accordo medio:

| footprint | margin 0.5 | margin 1.0 | margin 2.0 |
|---|---|---|---|
| 0.35 | 0.847 | 0.833 | 0.786 |
| 0.50 | 0.925 | 0.923 | 0.877 |
| **0.70** | 0.914 | **0.939** | 0.902 |
| 1.00 | 0.834 | 0.919 | 0.900 |
| 1.40 | 0.748 | 0.875 | 0.889 |

All'ottimo, 0.7 e 1.0, i falsi visibili sono il 2.77% e i falsi occlusi il
3.28%. Il compromesso si legge bene: impronta troppo piccola e lo sfondo
trapela (a 0.35 e margin 2.0 i falsi visibili sono il 20.9%), troppo grande
e i punti si occludono a vicenda (a 1.4 e margin 0.5 i falsi occlusi sono il
25.1%). Il margin sposta l'errore da un tipo all'altro nello stesso modo.
L'ottimo è piatto fra 0.5 e 1.0 di impronta, quindi la scelta è robusta e
non su misura. È anche lo stesso ottimo trovato con il generatore precedente,
su scene di forma diversa.

**Quanto dipende dalla forma.** Con impronta e margin fissi, al variare della
regolarità $\nu$ del campo: 84.2% a $\nu$ = 0.5, 91.0% a 1.5, 93.9% a 2.5,
94.3% a 5.0. Una superficie più ruvida ha più bordi e più pieghe, dove il
disco di un campione sborda oltre la silhouette.

### La taratura trasferisce?

Congelati i due valori, riverificati su una scena **mai usata per tararli**
(seed 7) e su 3 camere, variando un parametro alla volta rispetto a 8000
punti, fill 0.30, lunghezza 0.3:

| variazione | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|
| fill 0.15 | 0.926 | 3.48% | 3.90% |
| fill 0.45 | 0.913 | 4.01% | 4.72% |
| lunghezza 0.2 | 0.890 | 4.20% | 6.77% |
| lunghezza 0.4 | 0.941 | 2.73% | 3.22% |
| **3000 punti** | **0.901** | **4.04%** | **5.83%** |
| 30000 punti | 0.957 | 2.19% | 2.14% |

Due letture. La prima: `fill` quasi non sposta nulla, mentre la lunghezza di
correlazione sì, per la stessa ragione di $\nu$ sopra: una lunghezza minore
vuol dire una superficie più increspata, con più bordi.

La seconda: **la densità muove l'errore**, ed entrambi i tipi insieme. A
3000 punti l'accordo scende al 90.1%, a 30000 sale al 95.7%. È il
comportamento atteso — meno punti, spaziatura maggiore, impronta più
grossolana — ed è anche la ragione per cui ci si può fidare sul reale: una
scansione vera è rada proprio dove il range è lungo e l'incidenza radente, ed
è lì che il filtro degrada per primo, senza crollare.

Attenzione a non confrontare questi numeri al decimale con il 93.9% della
taratura: sono scene e camere diverse, e la varianza fra pose è dello stesso
ordine delle differenze in tabella.

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

Scena di taratura seed 3, 8000 punti, 2 camere. Rumore gaussiano isotropo,
con σ in unità di spaziatura mediana:

| σ | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|
| 0.00 | 0.944 | 2.76% | 2.88% |
| 0.25 | 0.919 | 2.90% | 5.24% |
| 0.50 | 0.849 | 3.24% | 11.84% |
| 1.00 | 0.775 | 2.84% | 19.70% |

Outlier sparsi uniformemente nel volume, come frazione dei punti aggiunta
alla nuvola (quindi il totale cresce, e gli spuri sono valutati anch'essi):

| quota | accordo | falsi visibili | falsi occlusi |
|---|---|---|---|
| 0.0% | 0.944 | 2.76% | 2.88% |
| 0.2% | 0.924 | 2.74% | 4.92% |
| 1.0% | 0.883 | 2.51% | 9.18% |
| 3.0% | 0.826 | 2.40% | 14.99% |

**Degrada senza crollare, ma non poco.** Rumore pari a un quarto della
spaziatura costa due punti e mezzo di accordo; rumore pari all'intera
distanza fra punti vicini lo porta al 77%. Anche lo 0.2% di outlier costa due
punti, e l'1% sei.

**Il danno va tutto dalla parte giusta.** Sulla nuvola pulita i due errori si
equivalgono, circa il 3% ciascuno: il filtro **non** è conservativo di per
sé. Ma tutto l'errore che le perturbazioni aggiungono è in falsi *occlusi*,
mentre i falsi visibili restano fermi fra il 2.4% e il 3.2%. Sporcando la
nuvola il filtro diventa conservativo, cioè dichiara non visto ciò che forse
vedrebbe. Per la pianificazione di viste è il bias corretto: l'informazione
stimata cala invece di gonfiarsi. Resta però un 3% di falsi visibili di
base, che nessuna perturbazione crea e nessuna toglie.

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
su 1280 punti rimossi (due oggetti, 20000 punti) e 6 camere: accordo
**88.7%**, **falsi svuotati 0.0%**, mancati 11.3%.

Qui sì il bias è conservativo anche sulla nuvola pulita: non dice mai che
qualcosa è sparito quando c'è ancora, e l'errore sta tutto nei mancati. La
causa è l'**impronta**: i dischi della superficie in primo piano sbordano
oltre la silhouette e coprono pixel da cui in realtà si vede oltre. Con il
generatore precedente un'impronta di 0.5 dava tre punti di accordo in più,
al prezzo di qualche falso svuotato; resta 0.7, lo stesso della visibilità,
perché una seconda costante tarata su una scena sola non li vale, e perché
falsi svuotati a zero è la proprietà che conta — un pianificatore che crede
sparito qualcosa che c'è non ci guarderebbe più.

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

## Visualizzazione

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

- **Manici veri, cioè buchi passanti, non sono controllati.** Le forme sono
  blob lobati con fessure profonde e a volte due componenti, non oggetti
  topologicamente complessi per costruzione. Una `length` minore e una
  `smoothness` minore rendono la superficie più increspata e più ricca di
  pieghe, ma quanti manici producano non è misurato.
- **Il filtro di visibilità sbaglia di circa il 3% in entrambi i versi**
  anche su una nuvola pulita, e peggiora su superfici più ruvide. HPR, che
  si era provato prima, su queste superfici non regge: non va usato né come
  filtro né come verità di riferimento, che è invece la marcia lungo il
  raggio nel campo implicito.
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
