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

- $g(x) > 0$: sono **dentro**;
- $g(x) < 0$: sono **fuori**;
- $g(x) = 0$: sono **sulla pelle**.

Esempio facile, una palla di raggio $R$: $g(x) = R - |x|$. Al centro vale
$R$ (dentro), lontano è negativo (fuori), e si annulla esattamente sulla
sfera.

Tre vantaggi:

1. dire se un punto è dentro o fuori costa una valutazione;
2. **unire** e **intersecare** solidi è banale (§5 e §10);
3. il **gradiente** $\nabla g$ punta verso l'interno, perpendicolare alla
   pelle, e ci dà le normali gratis.

Questa rappresentazione si chiama *superficie implicita* o *insieme di
livello*.

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
  - $\nu = 1.5$: derivabile una volta.
  - $\nu = 2.5$: derivabile due volte.
  - $\nu \to \infty$: liscissimo, il kernel diventa gaussiano.

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

## 3. La media che scende: tenere l'oggetto al centro

### Il problema

Il campo del §2 è uguale ovunque. Se prendiamo "dove $f$ supera una
soglia", l'oggetto è sparso uniformemente e **tocca il bordo della scena
dappertutto**. Tagliandolo con la sfera (§5) uscirebbe una sfera piena di
buchi, una specie di groviera: il 25–34% della pelle sarebbe calotta
sferica, cioè una superficie grande, liscia, convessa, visibile da tutte le
camere. Poco realistica e troppo facile.

### La soluzione

Un GP può avere una **funzione media** non nulla (Rasmussen e Williams,
§2.7). Aggiungiamo al campo una conca, cioè una parabola rovesciata:

$$m(x) = -\text{MEAN\_DROP} \cdot \frac{|x|^2}{R^2}, \qquad \text{MEAN\_DROP} = 4$$

- Al centro $m = 0$ e il campo è quello di prima.
- Al bordo $m = -4$: il campo è spinto giù di **4 deviazioni standard**,
  e superare la soglia lì è quasi impossibile.

**Immagine mentale:** la nebbia è più densa al centro della stanza e si
dirada verso le pareti. L'oggetto nasce al centro e raramente arriva al
muro.

Perché 4: con 0 la calotta è il 27% della pelle, con 2 il 9%, con 4 lo
0.8%, e resta un oggetto solo. Con 6 l'oggetto si schiaccia verso una palla.

Nel codice, dentro `_component`:

```python
value = value + mean(distance / radius) - level
gradient = gradient - 2.0 * MEAN_DROP * x / radius ** 2
```

La seconda riga è la derivata di $-4|x|^2/R^2$, che vale $-8x/R^2$.

---

## 4. La soglia: quanto grande deve essere l'oggetto

Il campo completo è

$$g(x) = f(x) + m(x) - u$$

e l'oggetto è dove $g > 0$. Resta da scegliere il **livello** $u$. Non lo
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

$$g_{\text{scena}}(x) = \min\big(g(x),\; R - |x|\big)$$

**Perché il minimo è l'intersezione (Ricci 1973).** Un punto è dentro
entrambi i solidi se entrambi i campi sono positivi, cioè se il **più
piccolo** dei due è positivo. Dove l'oggetto arriva al muro, è il muro
($R - |x|$) a chiuderlo, e lì la pelle diventa un pezzo di sfera, la
*calotta*.

Il gradiente è quello del campo che "vince" il minimo: dove vince il muro è
$-x/|x|$.

```python
wall = radius - distance
cut = wall < value
return np.where(cut, wall, value), np.where(cut[:, None], -outward, gradient)
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

Sullo **spigolo** dove la pelle del campo incontra la calotta, il minimo del
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
   è $-\nabla g / |\nabla g|$, perché il gradiente punta dentro.
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

## 11. Comporre la scena: unione come massimo

$$g_{\text{scena}}(x) = \max\big(g_{\text{sfondo}}(x),\; g_1(x),\; g_2(x), \dots\big)$$

È il duale del §5: un punto è dentro l'unione se è dentro **almeno uno**,
cioè se il **più grande** dei campi è positivo (Ricci 1973).

- **Aggiungere** un oggetto: un termine in più nel massimo.
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
          + conca (media, §3)
          − livello u da fill (§4)
          ∩ sfera (min, §5)                        = campo dello sfondo g
          │
          ├─► tira nella palla, tieni il guscio (§6)
          │   Newton sulla pelle, controlla atterraggio (§7)
          │                                        = N punti + area A (§8)
          │
change_seed ─► per ogni oggetto: stessa ricetta in piccolo (§9)
               appoggia su superficie in su (§10)
          │
          ├─► old = max(sfondo, rimossi, spostati prima)   (§11)
          ├─► new = max(sfondo, aggiunti, spostati dopo)
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

## 16. Da dove viene ogni pezzo

| pezzo | riferimento |
|---|---|
| Somma di coseni con frequenze casuali, varianza unitaria (§2) | A. Rahimi, B. Recht, *Random Features for Large-Scale Kernel Machines*, NeurIPS 2007 |
| Kernel Matérn, densità spettrale, teorema di Bochner, funzione media (§2–§3) | C. E. Rasmussen, C. K. I. Williams, *Gaussian Processes for Machine Learning*, MIT Press 2006, §2.7 e §4.2 |
| Insieme di escursione e soglia (§4) | R. J. Adler, J. E. Taylor, *Random Fields and Geometry*, Springer 2007 |
| Intersezione come minimo, unione come massimo (§5, §11) | A. Ricci, *A Constructive Geometry for Computer Graphics*, The Computer Journal 1973 |
| Distanza al primo ordine $\lvert g\rvert/\lvert\nabla g\rvert$ (§6, §12) | G. Taubin, *Estimation of Planar Curves, Surfaces, and Nonplanar Space Curves Defined by Implicit Equations*, IEEE TPAMI 1991 |
| Proiezione di Newton sulla superficie implicita (§7) | A. Witkin, P. Heckbert, *Using Particles to Sample and Control Implicit Surfaces*, SIGGRAPH 1994 |
| Rotazioni casuali da QR, e la correzione dei segni (§10) | F. Mezzadri, *How to Generate Random Matrices from the Classical Compact Groups*, Notices of the AMS 2007 |
| Cambiamenti simulati per la change detection su nuvole di punti (§12) | I. de Gélis, S. Lefèvre, T. Corpetti, *Change Detection in Urban Point Clouds: An Experimental Comparison with Simulated 3D Datasets* (Urb3DCD), Remote Sensing 2021 |
