# Come generiamo le nuvole sintetiche, spiegato semplice

Spiegazione passo passo di `synthetic_point_cloud.py`: cosa fa, con quali
formule, da dove vengono e perché si fa così. Le note tecniche complete, con
le tarature e le verifiche, stanno in `notes.md`. Questo file serve a capire.

---

## 0. L'idea in una frase

> Inventiamo una "nebbia" casuale nello spazio, diciamo che l'oggetto è dove
> la nebbia è più densa di una soglia, e mettiamo i punti **sulla pelle**
> dell'oggetto.

Tutto il resto serve a rendere precisa questa frase:

- come si fa una nebbia casuale **con le proprietà giuste** (§2);
- come si evita che l'oggetto scappi fuori dalla scena (§3–§5);
- come si mettono i punti sulla pelle **in modo uniforme** (§6–§8);
- come si aggiungono, tolgono e spostano oggetti per creare un cambiamento
  (§9–§12).

---

## 1. Perché una pelle e non un volume

Una camera di profondità o un laser non vedono mai dentro un solido:
restituiscono punti sulle **superfici**. Se mettessimo punti sparsi in un
volume, concetti come "questo punto è nascosto da quello davanti" o "da qui
lo vedo di taglio" non avrebbero senso. Ma sono proprio le cose di cui il
pianificatore di viste ha bisogno. Quindi: **solido → pelle → punti sulla
pelle**.

### Il trucco del campo implicito

Un solido si descrive con una funzione $g(x)$ che a ogni punto dello spazio
dà un numero:

- $g(x) < 0$: sono **dentro**;
- $g(x) = 0$: sono **sulla pelle**;
- $g(x) > 0$: sono **fuori**.

Esempio facile, una palla di raggio $R$: $g(x) = |x| - R$. Al centro vale
$-R$ (dentro), lontano è positivo (fuori), e si annulla esattamente sulla
sfera.

È la convenzione delle *signed distance function*. Attenzione: la
letteratura sugli insiemi di escursione usa il segno opposto, cioè
l'oggetto è dove il campo **supera** un livello. Sono la stessa cosa a meno
di un segno.

Tre vantaggi:

1. dire se un punto è dentro o fuori costa una valutazione;
2. **unire** e **intersecare** solidi è banale (§5 e §10);
3. il **gradiente** $\nabla g$ punta verso l'esterno, perpendicolare alla
   pelle: è già la normale uscente.

Questa rappresentazione si chiama *superficie implicita* o *insieme di
livello*.

### Da $g$ a $f$: chi è chi

$g$ è **il solido**, $f$ è **il suo ingrediente casuale**. La nebbia $f(x)$
del §2 da sola non è un oggetto: è un numero che oscilla attorno a 0 in
tutto lo spazio. Il campo del solido si costruisce da lei in tre mosse:

$$g(x) = \underbrace{u}_{\text{soglia, §4}} - \underbrace{m(x)}_{\text{cupola, §3}} - \underbrace{f(x)}_{\text{nebbia, §2}}, \qquad g_{\text{scena}}(x) = \max\big(g(x),\, |x| - R\big) \;\; \text{(§5)}$$

così "$g < 0$" vuol dire esattamente "nebbia più cupola supera la soglia".
Il segno meno davanti a $f$ e $m$ è tutto qui: si guarda di quanto il campo
**manca** alla soglia, e mancare di meno di zero vuol dire superarla.
Esempio con il valore vero di default, $u = -1.64$ (§4), e la stessa
nebbia $f = -0.5$ in tre posti:

- al centro, $m = 0$: $g = -1.64 - 0 + 0.5 = -1.14$, dentro;
- a metà raggio, $m = -1$: $g = -1.64 + 1 + 0.5 = -0.14$, dentro di poco;
- a tre quarti, $m = -2.25$: $g = -1.64 + 2.25 + 0.5 = 1.11$, fuori.

La stessa nebbia è dentro al centro e fuori verso il bordo: è la cupola.

Dal §6 in poi si lavora solo con $g$: guscio, Newton, distanza e unione non
sanno che dentro ci sono onde. Nel codice è la riga di `_component`
`value = level - mean(distance / radius) - value`, dove `value` entra come
$f$ ed esce come $g$.

---

## 2. La nebbia: un campo casuale fatto di onde

### Cosa vogliamo dalla nebbia

Una funzione $f(x)$ casuale che:

- sia **liscia**, non rumore bianco, così la pelle è una superficie e non
  polvere;
- abbia **un'unica taglia tipica** delle gobbe, regolabile, così possiamo
  dire "oggetti grandi un terzo della scena";
- abbia **una rugosità regolabile**;
- sia uguale in tutte le direzioni e in tutti i punti (isotropa e
  stazionaria);
- abbia valori con una distribuzione nota, così si può calcolare la soglia.

L'oggetto matematico con tutte queste proprietà è un **processo gaussiano**
(GP): in ogni punto $f(x)$ è una gaussiana, e la somiglianza fra $f(x)$ e
$f(y)$ dipende solo dalla distanza $|x - y|$ tramite una funzione chiamata
**kernel** o covarianza $k(r)$.

- $k(0) = 1$: varianza 1.
- $k(r)$ scende con la distanza: punti vicini hanno valori simili, punti
  lontani sono indipendenti.

### Il kernel Matérn

Usiamo il kernel **Matérn**, lo standard in statistica spaziale
(Rasmussen e Williams 2006, §4.2). Ha due manopole:

- $\ell$ (`length`): la **taglia** delle gobbe. A distanza $\ell$ la
  correlazione è scesa a circa metà.
- $\nu$ (`smoothness`): la **rugosità**.
  - $\nu = 0.5$: il campo è continuo ma spigoloso ovunque, come una roccia
    frattale.
  - $\nu = 1.5$: derivabile una volta, in media quadratica.
  - $\nu = 2.5$: derivabile due volte, in media quadratica.
  - $\nu \to \infty$: liscissimo, il kernel diventa gaussiano.

  Vale per il processo teorico. Il campione con 256 onde è sempre
  infinitamente derivabile: $\nu$ decide quanta energia sta nelle rughe
  piccole, e il §16 spiega perché questo conta.

Noi usiamo $\nu = 2.5$, e per quel valore il kernel ha una forma chiusa:

$$k(r) = \left(1 + \frac{\sqrt5\, r}{\ell} + \frac{5 r^2}{3 \ell^2}\right) e^{-\sqrt5\, r / \ell}$$

Perché 2.5 e non altro: è il compromesso fra scene difficili (rugose, piene
di pieghe) e un filtro di visibilità affidabile. La tabella è in `notes.md`.
In breve, a 0.5 il filtro sbaglia il 16% dei punti, a 2.5 il 6%, e oltre si
guadagna poco.

### Come si estrae un campione: onde casuali

Estrarre esattamente un GP su milioni di punti costa troppo: servirebbe una
matrice gigante. Il trucco (Rahimi e Recht 2007, **random Fourier
features**) è sommare tante onde piane con direzione, frequenza e fase
casuali:

$$f(x) = \sqrt{\frac{2}{K}} \sum_{k=1}^{K} \cos(\omega_k \cdot x + b_k)$$

- $K$ = `n_modes` = 256 onde;
- $\omega_k$ è il **vettore d'onda**: la direzione dice verso dove oscilla,
  la lunghezza quanto fitta;
- $b_k \sim U(0, 2\pi)$ è la fase, cioè da dove parte l'onda.

**Immagine mentale:** lanci 256 sassi in uno stagno 3D e guardi
l'interferenza di tutte le onde. Dove si sommano in cresta c'è un picco,
dove si cancellano c'è una valle.

**Perché funziona (teorema di Bochner).** Un kernel stazionario è la
trasformata di Fourier di una densità di probabilità $p(\omega)$, detta
*densità spettrale*. Se estrai le frequenze $\omega_k$ da quella densità, la
somma di coseni ha **esattamente** quel kernel come covarianza. Per $K$
grande è anche gaussiana, per il teorema del limite centrale.

**Perché il peso $\sqrt{2/K}$.** Un coseno con fase uniforme ha varianza
$1/2$. Sommando $K$ termini indipendenti con peso $\sqrt{2/K}$ viene
$K \cdot \tfrac{2}{K} \cdot \tfrac12 = 1$: **varianza unitaria per
costruzione**. Serve al §4.

### Da dove si estraggono le frequenze per il Matérn

La densità spettrale del Matérn in 3D è

$$p(\omega) \propto \left(\frac{2\nu}{\ell^2} + |\omega|^2\right)^{-(\nu + 3/2)}$$

che è esattamente una **t di Student multivariata** con $2\nu$ gradi di
libertà e scala $1/\ell$. Una t si estrae dividendo una gaussiana per la
radice di una chi-quadro indipendente:

$$\omega_k = \frac{z_k}{\ell}\sqrt{\frac{2\nu}{g_k}}, \qquad z_k \sim \mathcal N(0, I_3), \qquad g_k \sim \chi^2_{2\nu}$$

Nel codice (`_random_field`):

```python
scale = np.sqrt(2.0 * smoothness / rng.chisquare(2.0 * smoothness, n_modes))
waves = rng.standard_normal((n_modes, 3)) * scale[:, None] / (length * radius)
```

**Intuizione su $\nu$.** La t ha code pesanti: ogni tanto estrae una
frequenza molto alta, cioè un'onda molto fitta, cioè una piccola ruga. Meno
gradi di libertà ($\nu$ piccolo) vuol dire code più pesanti e più rughe. Per
$\nu \to \infty$ la chi-quadro divisa per i gradi di libertà tende a 1, le
frequenze diventano gaussiane e le rughe spariscono.

### Il gradiente, gratis

Serve per Newton (§7) e per le normali. Si deriva il coseno:

$$\nabla f(x) = -\sqrt{\frac{2}{K}} \sum_{k} \sin(\omega_k \cdot x + b_k)\, \omega_k$$

Nel codice valore e gradiente escono dalla stessa proiezione
`x @ waves.T + phase`, quindi il gradiente non costa quasi niente in più.

### Controllo che sia vero

Su 4000 campi indipendenti, la correlazione misurata fra punti a distanza
$r$ coincide con la formula del Matérn 5/2: per esempio 0.519 misurato
contro 0.524 teorico a $r = 0.3$. La varianza misurata è 1.037. Quindi il
campo **è** quello che la letteratura dice, e si può citare.

---

## 3. La media a cupola: tenere l'oggetto al centro

### Il problema

Il campo del §2 è uguale ovunque. Se prendiamo "dove $f$ supera una
soglia", l'oggetto è sparso uniformemente e **tocca il bordo della scena
dappertutto**. Tagliandolo con la sfera (§5) uscirebbe una sfera piena di
buchi, una specie di groviera: il 25–34% della pelle sarebbe calotta
sferica, cioè una superficie grande, liscia, convessa, visibile da tutte le
camere. Poco realistica e troppo facile.

### La soluzione

Un GP può avere una **funzione media** non nulla (Rasmussen e Williams,
§2.7). Aggiungiamo al campo una **cupola**, cioè una parabola rovesciata:
massima (0) al centro, sempre più negativa verso il bordo. Lungo un raggio:

```
m(r)
  0 ┤●●●●●
    │      ●●●●
 -1 ┤          ●●●
    │             ●●●
 -2 ┤                ●●
    │                  ●●
 -3 ┤                    ●●
    │                      ●
 -4 ┤                       ●
    └──────────────────────────
     centro     metà      bordo
```

Immagina le onde della nebbia appoggiate sopra questa collina, e allaga
tutto fino all'altezza $u$: l'isola che emerge è l'oggetto. Le onde in cima
emergono, quelle sui fianchi restano sott'acqua.

$$m(x) = -\text{MEAN\_DROP} \cdot \frac{|x|^2}{R^2}, \qquad \text{MEAN\_DROP} = 4$$

- Al centro $m = 0$ e il campo è quello di prima.
- Al bordo $m = -4$: il campo è spinto giù di **4 deviazioni standard**,
  e superare la soglia lì è quasi impossibile.

**Immagine mentale:** la nebbia è più densa al centro della stanza e si
dirada verso le pareti. L'oggetto nasce al centro e raramente arriva al
muro.

Perché 4: con 0 la calotta è il 27% della pelle, con 2 il 9%, con 4 lo
0.8%, e resta un oggetto solo. Con 6 l'oggetto si schiaccia verso una palla.

**Cosa è citato e cosa è nostro.** Dal libro viene solo l'idea che un GP
possa avere una media non nulla, anche costruita da funzioni di base fisse
come i polinomi; in geostatistica è il *trend*. La **forma** parabolica e il
**valore** 4 sono una scelta nostra: la parabola è la funzione più semplice,
liscia e uguale in ogni direzione che vale 0 al centro e scende al bordo, e
4 è il valore più basso a cui le calotte spariscono nella prova di
`notes.md`. Nel report va scritta come scelta di modellazione giustificata
da quella prova, non come formula presa da altri.

Nel codice, dentro `_component`:

```python
value = value + mean(distance / radius) - level
gradient = gradient - 2.0 * MEAN_DROP * x / radius ** 2
```

La seconda riga è la derivata di $-4|x|^2/R^2$, che vale $-8x/R^2$.

---

## 4. La soglia: quanto grande deve essere l'oggetto

Il campo completo è

$$g(x) = u - m(x) - f(x)$$

e l'oggetto è dove $g < 0$, cioè dove $f + m$ supera $u$. Resta da scegliere
il **livello** $u$. Non lo
scegliamo a mano: chiediamo "voglio che l'oggetto occupi in media il 30%
della scena" (`fill` = 0.30) e ricaviamo $u$.

### Il ragionamento

Prendiamo un punto a distanza $r$ dal centro, con $r$ in frazione di $R$.

- $f$ lì è una gaussiana standard (media 0, varianza 1, grazie al §2).
- Il punto è dentro se $f > u - m(r)$.
- Quindi è dentro con probabilità

$$P_{\text{dentro}}(r) = 1 - \Phi\big(u - m(r)\big)$$

dove $\Phi$ è la cumulativa della gaussiana standard (`norm.cdf`).

La frazione di volume attesa è la media di questa probabilità su tutta la
palla. In una palla il guscio a raggio $r$ ha volume proporzionale a $r^2$:
la frazione di volume fra $r$ e $r + dr$ è $3r^2\,dr$. Quindi

$$\text{fill} = 3 \int_0^1 \Big(1 - \Phi\big(u - m(r)\big)\Big)\, r^2\, dr$$

Si risolve per $u$ numericamente:

- `quad` calcola l'integrale;
- `brentq` trova lo zero di "integrale − fill" fra −10 e 10.

Questo modo di descrivere un oggetto casuale, cioè "la regione dove un campo
gaussiano supera un livello", si chiama **insieme di escursione** (Adler e
Taylor 2007).

### Cosa rappresenta $u$

È il **livello dell'acqua** nell'immagine dell'isola del §3: acqua alta,
oggetto piccolo; acqua bassa, oggetto grande. Valori veri, con la cupola da
4:

| `fill` | $u$ | dentro al centro | a metà raggio | al bordo |
|---|---|---|---|---|
| 0.15 | −0.84 | 80% | 44% | 0.1% |
| **0.30** | **−1.64** | **95%** | **74%** | **0.9%** |
| 0.45 | −2.26 | 99% | 90% | 4% |

Il livello è **negativo**: la cupola abbassa tutto il campo, e l'acqua deve
scendere sotto lo zero per lasciar emergere il 30%. Senza cupola basterebbe
$u \approx +0.5$. L'oggetto viene quindi fatto di un nucleo quasi pieno al
centro, con la pelle frastagliata dalle onde nella fascia intermedia: è lì
che nascono pieghe e cavità.

### Onestà: `fill` è una media, non una promessa

Una scena singola contiene poche gobbe, circa tre lunghezze di correlazione
per raggio, quindi oscilla molto. Chiedendo 0.30 si ottiene
**0.318 ± 0.083**, con un piccolo eccesso sistematico di circa 0.02. Il
vecchio generatore dava il fill esatto a ogni scena ma non era citabile. È
uno scambio voluto.

---

## 5. La chiusura: tagliare con la sfera

L'oggetto deve essere **chiuso**, cioè avere un dentro e un fuori ben
definiti, e stare nella palla di raggio $R$. La media del §3 lo tiene quasi
sempre dentro, ma "quasi" non basta. Allora lo si interseca con la palla:

$$g_{\text{scena}}(x) = \max\big(g(x),\; |x| - R\big)$$

**Perché il massimo è l'intersezione (Ricci 1973).** Un punto è dentro
entrambi i solidi se entrambi i campi sono negativi, cioè se il **più
grande** dei due è negativo. Dove l'oggetto arriva al muro, è il muro
($|x| - R$) a chiuderlo, e lì la pelle diventa un pezzo di sfera, la
*calotta*.

Il gradiente è quello del campo che "vince" il massimo: dove vince il muro è
$x/|x|$, cioè punta fuori.

```python
wall = distance - radius
cut = wall > value
return np.where(cut, wall, value), np.where(cut[:, None], outward, gradient)
```

Con la media del §3 le calotte sono lo 0–3% della pelle.

---

## 6. Mettere i punti sulla pelle: il guscio

Ora abbiamo $g$ e vogliamo $N$ punti **uniformi** sulla superficie $g = 0$.
Non si possono estrarre direttamente, quindi si procede in tre mosse:

1. **tira** punti a caso nella palla;
2. **tieni** quelli vicini alla pelle, dentro un guscio sottile;
3. **spingi** quelli tenuti esattamente sulla pelle (§7).

### Punti a caso in una palla (`_uniform_in_ball`)

- Direzione: una gaussiana 3D normalizzata, uniforme sulla sfera perché la
  gaussiana è isotropa.
- Raggio: $R \cdot U^{1/3}$. Il volume entro $r$ cresce come $r^3$: per
  essere uniformi nel volume la cumulativa del raggio è $(r/R)^3$, e
  invertendola si ottiene la radice cubica.

### Il guscio a spessore costante

Quanto è lontano un punto dalla pelle? Al primo ordine (Taubin 1991):

$$d(x) \approx \frac{|g(x)|}{|\nabla g(x)|}$$

**Perché dividere per la pendenza.** Immagina una collina. Il valore di $g$
è l'altezza, la pelle è la linea dove l'altezza vale zero. Se il pendio è
ripido, anche un'altezza grande corrisponde a pochi passi dalla linea; se è
piatto, anche un'altezza piccola corrisponde a molti passi. Dividendo per la
pendenza si passa dall'altezza ai **passi**, cioè a una distanza vera.

Teniamo i punti con $d(x) \le h$, dove $h$ = `SHELL` · $R$ = 0.02 $R$.

**Perché conta.** Se tenessimo i punti con $|g| \le$ costante, il guscio
sarebbe spesso dove il campo è piatto e sottile dove è ripido, e i punti si
ammucchierebbero nelle zone piatte. Con spessore costante in distanza, ogni
pezzo di pelle riceve punti in proporzione alla sua **area**: densità
uniforme.

### Il dettaglio della palla più larga

I candidati si tirano in una palla di raggio $R(1 + \text{SHELL})$, un filo
più larga della scena. Le calotte stanno **esattamente** sulla sfera $R$:
tirando solo dentro, il guscio attorno a loro esisterebbe da un lato solo, e
le calotte avrebbero metà dei punti.

### Rigetto a lotti

Non sappiamo in anticipo quanti candidati cadranno nel guscio: dipende da
quanta pelle ha quell'oggetto. Quindi si tira un lotto, si prende quello che
si trova, si tira un altro lotto per quelli che mancano, e così via.

- Ogni lotto ha $8\times$ i punti mancanti, e almeno 4096.
- Si rinuncia dopo `MAX_BATCHES` = 200 lotti.

---

## 7. La spinta sulla pelle: Newton

Un punto nel guscio è **vicino** alla pelle, non **sopra**. Lo si spinge con
Newton sull'equazione $g(x) = 0$ (Witkin e Heckbert 1994):

$$x \leftarrow x - \frac{g(x)}{|\nabla g(x)|^2}\, \nabla g(x)$$

**Come funziona.** È Newton in una dimensione, "valore diviso derivata",
fatto lungo la direzione del gradiente, cioè perpendicolare alla pelle, che è
la strada più corta. La correzione è $d(x)$, la stessa distanza del §6,
nella direzione giusta.

Converge in modo quadratico: l'errore si eleva al quadrato a ogni passo.
Bastano `STEPS` = 4 passi per arrivare alla precisione della macchina.

### Il controllo di atterraggio

Sullo **spigolo** dove la pelle del campo incontra la calotta, il massimo del
§5 cambia vincitore da un passo all'altro. Newton rimbalza fra le due
superfici e può non atterrare mai: in una prima versione un punto era finito
a 0.47 $R$ dalla pelle. Quindi dopo Newton si tengono solo i punti con
$d(x) \le$ `LANDED` · $R$ = $10^{-9} R$. Lo spigolo è una linea: scartare
quei pochi punti non cambia la densità altrove.

Si scartano anche i punti finiti fuori dalla sfera oltre l'arrotondamento.

---

## 8. Regalo: l'area della pelle

Mentre campioniamo, misuriamo anche l'**area** $A$ della superficie, che
serve per le densità (§9).

Il guscio ha semispessore $h$, quindi volume $\approx 2hA$. I candidati sono
uniformi nella palla grande di volume $V$, quindi la frazione che cade nel
guscio è $2hA/V$:

$$A \approx \frac{\text{nel guscio}}{\text{tirati}} \cdot \frac{V}{2h}$$

Dall'area si ricava la **spaziatura** media fra i punti,
$s = \sqrt{A/N}$: ogni punto "possiede" un quadratino di area $A/N$. Serve
per la tolleranza del §12.

**Verifica di uniformità.** Il coefficiente di variazione della distanza al
vicino più prossimo è 0.52–0.53, lo stesso di punti messi uniformemente a
caso su un piano. Vale anche sulle calotte.

---

## 9. Gli oggetti che cambiano

Lo **sfondo** è quanto costruito fin qui. Un **oggetto** si costruisce
**con la stessa ricetta**, in piccolo, con `_object_shape`:

- palla di raggio `object_radius` · $R$ = 0.25 $R$;
- `OBJECT_FILL` = 0.45, più pieno dello sfondo, così esce corposo e non
  filiforme;
- stessi `length` e `smoothness`, relativi al suo raggio.

### Stessa densità dello sfondo

Ogni punto deve rappresentare la stessa area ovunque, altrimenti un oggetto
fitto peserebbe di più nell'obiettivo del pianificatore senza motivo. Quindi:

1. si campionano prima `OBJECT_PROBE` = 256 punti, e intanto si misura
   l'area $A_{\text{obj}}$ (§8);
2. il numero giusto è $N \cdot A_{\text{obj}} / A_{\text{sfondo}}$, con un
   minimo di `MIN_OBJECT` = 50;
3. se servono più punti se ne campionano altri, se ne servono meno si
   taglia.

**Da sapere:** misurati, gli oggetti escono circa il 15% più radi dello
sfondo (spaziatura 0.0117 contro 0.0101). La causa non è stata ancora
indagata; il primo sospetto è la stima d'area sul primo campione da 256.

---

## 10. Appoggiare un oggetto (`_place_object`)

Un oggetto vero non galleggia: poggia su qualcosa.

1. **Scegli dove.** Un punto dello sfondo con la normale rivolta verso
   l'alto: componente verticale della normale uscente > `UPWARD` = 0.3,
   cioè una superficie su cui qualcosa potrebbe stare. La normale uscente
   è $\nabla g / |\nabla g|$, perché il gradiente punta fuori.
2. **Gira a caso.** Una rotazione casuale dalla fattorizzazione QR di una
   matrice gaussiana 3×3; se il determinante esce −1, cioè uno specchio, si
   cambia segno.
3. **Fai scendere.** Metti il centro sopra il punto, lungo la normale, a
   distanza 1.0 · raggio, poi 0.8, 0.6, 0.4, 0.2 (`CONTACT`). Ti fermi
   alla prima distanza in cui almeno il 2% dei punti dell'oggetto
   (`BURIED`) è affondato nello sfondo: l'oggetto **tocca**.
4. **Regole:**
   - il **centro** deve stare nella sfera della scena, quindi l'oggetto
     può sporgere al massimo del suo raggio;
   - due oggetti devono avere i centri distanti almeno due raggi.
5. Se dopo `ATTEMPTS` = 200 tentativi non ci si riesce, errore: ci sono
   troppi oggetti o sono troppo grandi.

### Spostare un campo (`_placed`)

L'oggetto è stato costruito centrato nell'origine. Per metterlo in posa
(rotazione $Q$, centro $t$) non si ricostruisce niente: si valuta il campo
originale nelle coordinate dell'oggetto,

$$g_{\text{in posa}}(x) = g\big(Q^\top (x - t)\big), \qquad \nabla g_{\text{in posa}}(x) = Q\, \nabla g\big(Q^\top(x - t)\big)$$

e i punti si spostano con $Qp + t$.

**Nota:** la QR senza la correzione dei segni sulla diagonale di $R$ non dà
una rotazione esattamente uniforme (Mezzadri 2007). Per mettere oggetti in
pose varie va bene, ma non è una distribuzione uniforme sulle rotazioni.

---

## 11. Comporre la scena: unione come minimo

$$g_{\text{scena}}(x) = \min\big(g_{\text{sfondo}}(x),\; g_1(x),\; g_2(x), \dots\big)$$

È il duale del §5: un punto è dentro l'unione se è dentro **almeno uno**,
cioè se il **più piccolo** dei campi è negativo (Ricci 1973).

- **Aggiungere** un oggetto: un termine in più nel minimo.
- **Toglierlo**: un termine in meno.
- **Spostarlo**: toglierlo da una posa e aggiungerlo in un'altra, con la
  **stessa forma**. Si perde l'informazione che è lo stesso oggetto. Per
  pianificare non importa, perché bisogna guardare sia dov'era sia dov'è.

### I punti della scena (`_assemble`)

I punti di ogni parte si tengono solo se stanno **fuori da tutte le altre
parti**. Un punto dello sfondo finito dentro un oggetto appoggiato non è più
pelle: l'oggetto lo copre. Viceversa per la parte di oggetto affondata nello
sfondo.

### Due scene

`generate_scene_pair` costruisce:

- **old**: sfondo + oggetti rimossi + spostati nella posa vecchia;
- **new**: sfondo + oggetti aggiunti + spostati nella posa nuova.

Lo sfondo è **campionato una volta sola**: i punti che non cambiano sono
identici bit per bit nelle due scene.

- `seed` decide lo sfondo.
- `change_seed` decide forme e pose degli oggetti.

Si può quindi tenere la stessa stanza e cambiare solo cosa succede.

---

## 12. Cos'è cambiato: si misura, non si dichiara

Non diciamo noi "questi punti sono cambiati": lo **misuriamo**.

- **aggiunto** = punto di *new* distante più di mezza spaziatura
  (`TOLERANCE` = 0.5 · $s$) dalla pelle di *old*;
- **rimosso** = punto di *old* distante più di mezza spaziatura dalla pelle
  di *new*;

con la distanza al primo ordine del §6 sul campo unione dell'altra scena.

Così escono giusti da soli, senza programmarli, due casi sottili:

- un oggetto appoggiato **copre** un pezzo di sfondo. Aggiungerlo fa
  risultare **rimossi** i punti coperti; toglierlo li **rivela**, e
  risultano **aggiunti**, perché sono superficie che il modello vecchio non
  aveva mai visto. Nell'esempio di `notes.md` sono i 52 punti;
- uno **spostamento corto** segna solo la parte dell'oggetto che non si
  sovrappone alla posa vecchia.

**Limite:** è solo geometria. Una superficie che cambia colore ma non forma
non risulta cambiata.

---

## 13. Tutto il flusso in un colpo d'occhio

```
seed ──► onde casuali (Matérn, §2)
          + cupola (media, §3)
          − livello u da fill (§4)
          ∩ sfera (max, §5)                        = campo dello sfondo g
          │
          ├─► tira nella palla, tieni il guscio (§6)
          │   Newton sulla pelle, controlla atterraggio (§7)
          │                                        = N punti + area A (§8)
          │
change_seed ─► per ogni oggetto: stessa ricetta in piccolo (§9)
               appoggia su superficie in su (§10)
          │
          ├─► old = min(sfondo, rimossi, spostati prima)   (§11)
          ├─► new = min(sfondo, aggiunti, spostati dopo)
          │
          └─► added / removed = distanza dall'altra pelle > s/2   (§12)
```

---

## 14. Manopole

| parametro | default | cosa cambia a girarlo |
|---|---|---|
| `n_points` | — | punti sullo sfondo, cioè la densità |
| `radius` | — | raggio della scena $R$ |
| `n_modes` | 256 | onde sommate: più onde, campo più gaussiano, più lento |
| `length` | 0.30 | taglia delle gobbe, in raggi: più piccolo, più pieghe e più occlusione |
| `smoothness` | 2.5 | rugosità $\nu$: più piccolo, superficie più rugosa |
| `fill` | 0.30 | frazione attesa della scena occupata dall'oggetto |
| `n_added`, `n_removed`, `n_moved` | 1, 0, 0 | quanti oggetti per tipo |
| `object_radius` | 0.25 | raggio di un oggetto, in raggi della scena |
| `seed`, `change_seed` | — | riproducibilità di sfondo e cambiamento |

Le costanti interne (`MEAN_DROP`, `SHELL`, `STEPS`, …) e il perché dei loro
valori sono nelle tabelle di `notes.md`.

---

## 15. Cosa le scene hanno di buono e cosa no

**Buono, verificato**

- La covarianza del campo coincide con il Matérn teorico.
- I punti stanno sulla pelle con errore sotto $10^{-9} R$, e la densità è
  uniforme.
- C'è **auto-occlusione** vera: da una vista si vede solo il 36–40% della
  pelle, quindi scegliere le viste conta.
- Nessuna superficie è irraggiungibile da tutte le camere.

**Limiti dichiarati**

- `fill` vale solo in media.
- Gli oggetti escono un po' più radi dello sfondo, per una causa non
  indagata.
- Il cambiamento è solo geometrico.
- Uno spostamento perde l'identità dell'oggetto.
- Il numero di buchi passanti (manici) non è controllato.
- **La nuvola è perfetta**: completa, senza rumore, senza buchi. Una nuvola
  reale viene da un'acquisizione parziale. Questo generatore dà la
  **verità**, non ciò che un sensore vedrebbe.

---

## 16. Perché proprio questa $f$ e non un'altra

Qui la scelta di $f$ si giustifica come una catena: prima i **requisiti**, che
vengono dal problema e non da $f$; poi ogni scelta come l'unica, o la più
semplice, che li soddisfa, con le alternative e il requisito su cui cadono.
Alla fine resta quello che è una **convenzione** o una **scelta di
modellazione**, da dichiarare come tale.

### 16.1 I requisiti

Il simulatore è un banco di prova per un pianificatore di viste. Da questo
seguono sei requisiti sull'oggetto:

- **R1 — varietà vera.** Concavità, cavità, più componenti, forme che nessuno
  ha disegnato. Un pianificatore valutato su forme scelte a mano è valutato
  sul gusto di chi le ha scelte.
- **R2 — difficoltà regolabile con parametri interpretabili.** Per dire su
  quali scene si è testato e per rendere il test più difficile in modo
  controllato: la taglia delle pieghe, che decide quanta auto-occlusione c'è,
  e la rugosità, che decide quanto è affidabile il filtro di visibilità.
- **R3 — nessun posto e nessuna direzione privilegiati**, se non quelli
  imposti apposta. Altrimenti il risultato dipende da dove capita il
  cambiamento rispetto a un riferimento arbitrario.
- **R4 — una funzione vera, valutabile in qualunque punto, sempre uguale, con
  il gradiente.** La verità di visibilità marcia lungo raggi arbitrari, Newton
  valuta punti arbitrari, il cambiamento valuta la scena vecchia sui punti
  della nuova. Tutti devono vedere **la stessa** superficie, esattamente.
- **R5 — una distribuzione dei valori nota**, così la dimensione
  dell'oggetto si fissa con una formula (§4) e non per tentativi.
- **R6 — una costruzione standard, con proprietà dimostrate e citabili.**

### 16.2 Perché un campo casuale e non un catalogo di forme

| alternativa | dove cade |
|---|---|
| Forme parametriche con parametri casuali (sfera, toro, C, L) | **R1**: la varietà è apparente, le topologie sono quelle elencate |
| Sezione trascinata lungo una curva casuale (provata) | **R1**: sempre un tubo, mai un oggetto con massa |
| Somma di bolle in posizioni casuali (*metaball*, Blinn 1982) | **R1** in parte: forme a grappolo di sfere. **R5**: il valore del campo dipende da quante bolle si sovrappongono, niente formula per `fill` |
| Modelli CAD o scansioni reali (ShapeNet, ecc.) | **R2**: nessuna manopola di difficoltà. **R4**: la verità esatta richiede una mesh chiusa e pulita. Restano un'alternativa **complementare** per il realismo, da lavoro futuro |

L'insieme di escursione di un campo casuale soddisfa R1 senza sforzo: la
topologia non si costruisce, viene dalla casualità. Aumentando la taglia
relativa della scena rispetto alle pieghe compaiono da sole più componenti,
più manici, più cavità.

### 16.3 Perché gaussiano

| alternativa | dove cade |
|---|---|
| Rumore di Perlin o simplex (Perlin 1985), lo standard della grafica | **R5**: la distribuzione dei valori non è gaussiana e non ha forma chiusa. **R3**: è costruito su un reticolo, e resta un'anisotropia legata agli assi. **R6**: è una tecnica, non un modello con proprietà enunciabili |
| Campi non gaussiani (chi-quadro, t, log-normali) | **R5** in parte: la formula per `fill` esiste ma cambia modello per modello. Aggiungono asimmetrie che niente nel problema richiede |
| **Processo gaussiano** | È determinato **solo** da media e covarianza: tutte le manopole stanno lì (**R2**). I valori sono gaussiani, quindi `fill` ha la formula del §4 (**R5**). Esiste una teoria completa dei suoi insiemi di escursione (Adler e Taylor 2007) (**R6**) |

C'è anche un argomento di principio: fra tutte le distribuzioni con una media
e una covarianza date, la gaussiana è quella a **massima entropia** (Cover e
Thomas 2006, cap. 12), cioè quella che aggiunge meno ipotesi oltre alle due
che si sono scelte. Scegliendo la covarianza si dice quanto sono grandi e
quanto sono rugose le pieghe; scegliendo la gaussiana si dice di non
presumere nient'altro.

### 16.4 Perché stazionario e isotropo

Stazionario e isotropo vuol dire

$$\operatorname{Cov}\big(f(x), f(y)\big) = k(\lvert x - y\rvert),$$

cioè la somiglianza fra due punti dipende solo dalla loro distanza. È la
traduzione letterale di **R3**. Ha due conseguenze utili.

- **Tutta la struttura spaziale è esplicita.** L'unica cosa che distingue il
  centro dal bordo è la cupola $m$, che si dichiara e si giustifica a parte.
  Non c'è struttura nascosta nel campo.
- **La legge della scena è invariante per rotazioni attorno al centro**,
  perché $f$ è isotropa e $m$ dipende solo da $\lvert x\rvert$. Quindi le
  camere sul locus sferico non hanno orientazioni privilegiate: nessun
  risultato dipende da come sono orientati gli assi.

È però una **scelta di modellazione**, e va dichiarata come tale: le scene
vere **non** sono isotrope. La gravità crea pavimenti orizzontali, pareti
verticali, oggetti appoggiati. Il generatore produce scene **generiche**, non
stanze. L'unico elemento non isotropo è voluto: gli oggetti si appoggiano su
superfici rivolte verso l'alto.

### 16.5 Perché la covarianza Matérn

Un processo stazionario e isotropo è determinato dal suo kernel $k(r)$, o
equivalentemente dalla sua **densità spettrale** $p(\omega)$, che dice quanta
energia c'è a ogni frequenza. È lo spettro a decidere la forma: quanto sono
grandi le pieghe e quante rughe piccole ci sono sopra.

Il Matérn ha densità spettrale in 3D

$$p(\omega) \propto \left(\frac{2\nu}{\ell^2} + \lvert\omega\rvert^2\right)^{-(\nu + 3/2)}$$

che si legge in due pezzi:

- per $\lvert\omega\rvert \ll 1/\ell$ è **piatta**: tutte le pieghe più
  grandi di $\ell$ hanno la stessa energia, quindi **$\ell$ è la taglia**
  oltre la quale non c'è più struttura;
- per $\lvert\omega\rvert \gg 1/\ell$ scende come una **legge di potenza**,
  $\lvert\omega\rvert^{-(2\nu + 3)}$: **$\nu$ decide quanto velocemente
  spariscono le rughe piccole**.

Le due manopole sono quindi **separate**: una sposta dove inizia la discesa,
l'altra ne decide la pendenza. È esattamente **R2**.

| alternativa | spettro | dove cade |
|---|---|---|
| Gaussiano (*squared exponential*) | scende come $e^{-\ell^2\lvert\omega\rvert^2/2}$, più veloce di ogni potenza | **R2**: un parametro solo, nessuna rugosità. Liscio in modo irrealistico: Stein (1999) lo sconsiglia proprio per questo per i dati spaziali |
| Esponenziale | legge di potenza con $\nu = 1/2$ | rugoso a ogni scala: vedi 16.6 |
| Rational quadratic | miscela di gaussiani a scale diverse | il secondo parametro mescola scale invece di controllare la rugosità; meno standard per i dati spaziali |
| **Matérn** | piatto, poi potenza | due manopole separate; **contiene gli altri come casi limite** (esponenziale a $\nu = 1/2$, gaussiano per $\nu \to \infty$); è il kernel raccomandato per la statistica spaziale (Stein 1999; Rasmussen e Williams §4.2) |

La coda a legge di potenza ha anche un argomento di realismo: le superfici
reali hanno tipicamente spettri di rugosità a legge di potenza su un ampio
intervallo di scale (Persson et al. 2005), mentre la coda esponenziale del
kernel gaussiano non ha un corrispettivo fisico.

### 16.6 Perché $\nu = 2.5$: teoria, forma chiusa, esperimento

Qui la scelta ha **tre** giustificazioni indipendenti che puntano allo stesso
valore.

**La teoria: $\nu > 2$ perché la curvatura sia definita.** Un processo
Matérn è derivabile $n$ volte in media quadratica se e solo se $\nu > n$
(Rasmussen e Williams §4.2.1). Nelle random Fourier features questo si vede
direttamente. Le frequenze sono una t di Student con $2\nu$ gradi di
libertà, che ha momenti finiti solo di ordine minore di $2\nu$, e

- la varianza del **gradiente** di $f$ è proporzionale a
  $\mathbb E\lvert\omega\rvert^2$, finita solo se $\nu > 1$;
- la varianza delle **derivate seconde**, cioè della curvatura della
  superficie, è proporzionale a $\mathbb E\lvert\omega\rvert^4$, finita
  solo se $\nu > 2$.

Verificato estraendo un milione di frequenze e confrontando le due metà del
campione:

| $\nu$ | $\mathbb E\lvert\omega\rvert^2$, due metà | $\mathbb E\lvert\omega\rvert^4$, due metà |
|---|---|---|
| 0.5 | 1.8·10⁵ contro 8.8·10⁹ | ~10¹⁵ contro ~10²⁵ |
| 1.5 | 8.8 contro 9.0 | 3570 contro 4213 |
| **2.5** | **5.0 contro 5.0** | **127 contro 125** |

Un momento infinito si riconosce così: la media non si stabilizza, e la
domina la frequenza più grande che si è estratta. Cosa vuol dire per la
scena:

- con $\nu \le 1$ le **normali** dipendono dalle poche onde più fitte
  estratte per caso;
- con $1 < \nu \le 2$ le normali sono stabili ma la **curvatura** no;
- con $\nu > 2$ entrambe hanno statistiche che convergono al crescere di $K$.

In altre parole: con $\nu \le 2$ il dettaglio fine della superficie è deciso
dal **troncamento** a 256 onde, cioè da un artefatto numerico. Con
$\nu > 2$ è deciso dal **modello**. La curvatura non è un dettaglio: la
usano la stima delle normali con la PCA locale, la distanza al primo ordine
(il cui errore è proporzionale alla curvatura) e il guscio di campionamento
(la cui uniformità è valida fino a termini di ordine spessore × curvatura).

**La forma chiusa: il primo semi-intero oltre 2.** Per $\nu$ semi-intero
il Matérn si scrive come esponenziale per polinomio, senza funzioni di
Bessel. Rasmussen e Williams indicano 3/2 e 5/2 come i valori di uso
comune. 5/2 è il più piccolo semi-intero con $\nu > 2$: la formula del §2
che si confronta con la covarianza misurata è quella.

**L'esperimento: dove il filtro smette di migliorare.** Accordo del filtro
di visibilità con la verità, e frazione visibile da una vista, cioè la
difficoltà:

| $\nu$ | accordo | visibile da una vista |
|---|---|---|
| 0.5 | 84.2% | 27% |
| 1.5 | 90.9% | 35% |
| **2.5** | **93.9%** | **37%** |
| 5.0 | 94.3% | 43% |

Da 1.5 a 2.5 il filtro guadagna tre punti; da 2.5 a 5 ne guadagna mezzo,
mentre la scena diventa più facile di sei punti. Il gomito cade fra 1.5 e
2.5, cioè a cavallo della soglia $\nu = 2$ oltre la quale la teoria dice che
la curvatura è definita: **i tre argomenti sono coerenti**. Con quattro
valori di $\nu$ provati è una coerenza, non una dimostrazione che il gomito
stia esattamente a 2.

### 16.7 Perché $\ell = 0.3\,R$

$\ell$ decide quanta topologia c'è. Per un campo gaussiano stazionario con
varianza 1 in 3D, il numero atteso di "pezzi di topologia" per unità di
volume, cioè la densità della caratteristica di Eulero dell'insieme di
escursione, è (Adler e Taylor 2007)

$$\rho(u) = \frac{\lambda_2^{3/2}}{(2\pi)^2}\,(u^2 - 1)\,e^{-u^2/2}, \qquad \lambda_2 = -k''(0) = \frac{\nu}{(\nu - 1)\,\ell^2}$$

con $\lambda_2$ il **secondo momento spettrale**, che per il Matérn con
$\nu = 5/2$ vale $5/(3\ell^2)$ (verificato numericamente). Quindi a livello
fissato il numero di componenti, manici e cavità scala come
$(R/\ell)^3$:

- $\ell \gg R$: un solo blob convesso, ogni vista vede quasi tutto;
- $\ell \ll R$: una spugna di bolle piccole, più piccole della spaziatura
  dei punti;
- $\ell = 0.3\,R$: circa tre pieghe per raggio, cioè un oggetto con lobi,
  fessure e auto-occlusione vera (36–40% visibile da una vista), ancora ben
  campionato da qualche migliaio di punti.

Anche questo argomento chiede $\nu > 1$: con $\lambda_2$ infinito la
formula esplode, e il numero di pezzi di topologia dipenderebbe dal
troncamento.

Due precisazioni. La formula vale per un campo senza media; con la cupola
vale localmente, con $u$ sostituito da $u - m(x)$. E $\ell = 0.3$ non è
ottimizzato: è una scelta di scala, e la trasferibilità del filtro è stata
verificata anche a 0.2 e a 0.4 (`notes.md`).

### 16.8 Perché media zero e varianza uno

Non è una restrizione, è una scelta di **unità di misura**. Moltiplicare
$f$ per $\sigma$ e spostarla di $\mu$ dà l'insieme

$$\{\sigma f + \mu + m - u \ge 0\} = \Big\{f + \tfrac{m}{\sigma} - \tfrac{u - \mu}{\sigma} \ge 0\Big\},$$

lo stesso che si ottiene con $f$ standard, cupola $m/\sigma$ e livello
$(u - \mu)/\sigma$. Nessuna forma si perde. Il vantaggio è che tutto si
legge in deviazioni standard: `MEAN_DROP` = 4 vuol dire che al bordo il
campo è spinto giù di quattro deviazioni.

### 16.9 Perché le random Fourier features

| alternativa | dove cade |
|---|---|
| GP esatto con Cholesky sui punti | **R4**: il campo esiste solo sui punti fissati prima, non dove serve dopo (Newton, raggi, l'altra scena); inoltre costa $O(N^3)$ |
| Griglia con FFT, poi interpolazione | **R4**: verità non più esatta, gradiente approssimato, periodicità della griglia |
| Approccio SPDE su mesh (Lindgren, Rue e Lindström 2011) | **R4**: il campo vive sulla mesh, e tra i nodi è interpolato |
| **Random Fourier features** | una funzione vera, valutabile ovunque e sempre uguale, infinitamente derivabile, con gradiente analitico, $O(K)$ per punto. La covarianza è esattamente quella del Matérn per il teorema di Bochner (Rahimi e Recht 2007) |

Il prezzo va dichiarato: per $K$ finito i valori sono gaussiani solo
**approssimativamente**, per il teorema del limite centrale, e la covarianza
di un singolo campione oscilla attorno a quella teorica con errore
$O(1/\sqrt K)$. Con $K = 256$ la covarianza misurata coincide con quella
teorica entro 0.014 (§2), mentre `fill` risulta 0.318 invece di 0.30.

### 16.10 Cosa resta una scelta

| scelta | natura | come si giustifica |
|---|---|---|
| Campo gaussiano | principio | massima entropia date media e covarianza |
| Stazionario e isotropo | **modellazione** | neutralità (R3); le scene vere non lo sono, va nei limiti |
| Matérn | standard | due manopole separate, contiene gli altri kernel, raccomandato da Stein |
| $\nu = 2.5$ | teoria + esperimento | curvatura definita per $\nu > 2$, primo semi-intero, gomito del filtro |
| $\ell = 0.3\,R$ | **scala** | topologia $\propto (R/\ell)^3$; trasferibilità verificata a 0.2 e 0.4 |
| Media 0, varianza 1 | convenzione | nessuna perdita di generalità |
| Random Fourier features | tecnica | unico metodo che dà una funzione valutabile ovunque |
| $K = 256$ | numerico | covarianza verificata; `fill` in eccesso di 0.02 |

**In quattro frasi, per il report:**

> Serve un oggetto casuale di topologia arbitraria, quindi un insieme di
> escursione di un campo casuale. Il campo è gaussiano, perché è la scelta a
> massima entropia una volta fissate media e covarianza e perché rende
> analitica la dimensione dell'oggetto; stazionario e isotropo, per non
> favorire alcun punto o direzione. La covarianza è Matérn, l'unica famiglia
> standard che separa la scala delle pieghe dalla loro rugosità, con
> $\nu = 5/2$, il più piccolo valore in forma chiusa per cui la curvatura
> della superficie è definita. Il campo è campionato con random Fourier
> features, che lo rendono una funzione liscia, valutabile ovunque con il
> suo gradiente e con la covarianza esatta.

---

## 17. Da dove viene ogni pezzo

| pezzo | riferimento |
|---|---|
| Somma di coseni con frequenze casuali, varianza unitaria (§2) | A. Rahimi, B. Recht, *Random Features for Large-Scale Kernel Machines*, NeurIPS 2007 |
| Kernel Matérn, densità spettrale, teorema di Bochner, funzione media (§2–§3) | C. E. Rasmussen, C. K. I. Williams, *Gaussian Processes for Machine Learning*, MIT Press 2006, §2.7 e §4.2 |
| Insieme di escursione e soglia (§4) | R. J. Adler, J. E. Taylor, *Random Fields and Geometry*, Springer 2007 |
| Intersezione come massimo, unione come minimo (§5, §11) | A. Ricci, *A Constructive Geometry for Computer Graphics*, The Computer Journal 1973 |
| Distanza al primo ordine $\lvert g\rvert/\lvert\nabla g\rvert$ (§6, §12) | G. Taubin, *Estimation of Planar Curves, Surfaces, and Nonplanar Space Curves Defined by Implicit Equations*, IEEE TPAMI 1991 |
| Proiezione di Newton sulla superficie implicita (§7) | A. Witkin, P. Heckbert, *Using Particles to Sample and Control Implicit Surfaces*, SIGGRAPH 1994 |
| Matérn raccomandato per i dati spaziali, critica del kernel gaussiano (§16) | M. L. Stein, *Interpolation of Spatial Data: Some Theory for Kriging*, Springer 1999 |
| Gaussiana a massima entropia (§16) | T. M. Cover, J. A. Thomas, *Elements of Information Theory*, 2ª ed., Wiley 2006, cap. 12 |
| Rumore di Perlin, alternativa scartata (§16) | K. Perlin, *An Image Synthesizer*, SIGGRAPH 1985 |
| Metaball, alternativa scartata (§16) | J. F. Blinn, *A Generalization of Algebraic Surface Drawing*, ACM TOG 1982 |
| Campi Matérn come SPDE su mesh, alternativa scartata (§16) | F. Lindgren, H. Rue, J. Lindström, *An Explicit Link between Gaussian Fields and Gaussian Markov Random Fields: the SPDE Approach*, JRSS-B 2011 |
| Spettri di rugosità a legge di potenza delle superfici reali (§16) | B. N. J. Persson et al., *On the Nature of Surface Roughness with Application to Contact Mechanics, Sliding Friction, Rubber Friction and Adhesion*, J. Phys.: Condens. Matter 2005 |
| Rotazioni casuali da QR, e la correzione dei segni (§10) | F. Mezzadri, *How to Generate Random Matrices from the Classical Compact Groups*, Notices of the AMS 2007 |
| Cambiamenti simulati per la change detection su nuvole di punti (§12) | I. de Gélis, S. Lefèvre, T. Corpetti, *Change Detection in Urban Point Clouds: An Experimental Comparison with Simulated 3D Datasets* (Urb3DCD), Remote Sensing 2021 |
