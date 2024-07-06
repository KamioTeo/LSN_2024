# LSN_2024
## Esercizi del Laboratorio di Simulazione Numerica
Ogni cartella Lab contiene la risoluzione di ogni relativa unità di laboratorio, con la generazione dei dati per ogni esercizio utilizzando C++ e relativa analisi riportata su file jupyter. <br>
Ogni file C++ fa riferimento al codice più generale presente nella cartella "RandomNumberGenerator", che contiene le funzioni per generare numeri casuali distribuiti in base alla necessità.

## Note sugli esercizi:
Dato che ho seguito il corso nell'anno accademico 2022-2023, le esercitazioni 4, 6, 7 sono eseguite con il vecchio codice.
Nella cartella MDNVE-MCNVT è presente il codice per le esercitazioni 4 e 7.

### Esercizio 4
Nel file cpp MDNVE-MCNVT per cambiare stato di aggregazione va cambiato il valore in ingresso alla funzione Input(), subito dopo il main. 

Ho aggiunto dei parametri booleani in ingresso nel codice MD-MC:
- PrintAcceptance: se 1 stampa runtime a schermo il valore di accettazione
- PrintEqValues: se 1 stampa su file esterno i risultati di temperatura ed energia potenziale durante l'equilibrazione (Equilibration Values). Utile per il plot dell'equilibrazione.
- PrintEnergy: se 1 stampa su file esterno i valori dell'energia potenziale dopo l'equilibrazione

Alla fine del file Jupyter, è presente il codice per l'animazione 3D della disposizione delle posizioni nella struttura 2p. Il video già renderizzato è presente nella cartella stessa: ".\camera_rotation.mp4"

### Esercizio 6
Rispetto al codice originale, sono stati aggiunti dei parametri d'ingresso impostabili nel file input.dat:
- restart: Se 1, viene caricata l'ultima configurazione raggiunta, come da richiesta, da file config.final. Se 0, ricomincia la simulazione da una configurazione di spin casuale.
- Tsequence: Se 1, stampa un file con i valori dei dati termodinamici al variare della temperatura, da 0.5 a 2 (a passi di 0.1) come da richiesta.<br>
  Se 0, li stampa alla temperatura impostata da file di input.<br>
  Se non è né 0 né 1, stampa le informazioni sull'energia misurata (media a blocchi e incertezza statistica) al variare del numero di blocchi utilizzati (da 10 a 100 a passi di 10), fissati 
  10^6 step.<br>
- Equilibration: Se 1, non esegue il codice ma stampa un file con tutti i valori relativi all'energia del sistema per ogni step eseguito. Serve a monitorare la velocità di equilibrazione del     sistema.<br>
  Se 0 non equilibra (in questo esercizio non è necessario) ed esegue il resto del codice come descritto sopra.

è anche presente un file Python in cui ho implementato la simulazione del modello di Ising 2D, con video della magnetizzazione. Il video già renderizzato è presente al percorso ".\Ising_2D\IsingModel.mp4"

### Esercizio 7
- La base del codice è la stessa dell'esercizio 4, impostando su file di input iNVET=1 <br>
- Le tail corrections sono state aggiunte solamente alla parte di algoritmo di Metropolis, includendola per completezza anche al calcolo dell'energia potenziale per il campionamento della distribuzione di Boltzmann. Quest'ultima non è in realtà necessaria perché si lavora con differenze di energia, quindi la correzione si semplifica. <br>
- Per l'analisi della g(r) sono stati impostati neq = 200, nblk = 100 e nstep = 1000

### Esercizio 8
- Nella cartella è presente il codice "delta_finder.cpp" che, per un intervallo fornito di mu e sigma, calcola il valore migliore di delta da utilizzare nell'algoritmo di Metropolis (quello che fornisce un rapporto di accepted/attempted di circa il 50\%), per ogni coppia di valori mu e sigma con cui testare. <br>
- Per ognuna di queste coppie di valori mu e sigma viene calcolato il valor medio dell'energia, specificando dal file "delta_finder_input.in" il numero di step per equilibrare il sistema, il numero di step da effettuare dopo l'equilibrazione e il numero di blocchi in cui dividerle. <br>
- Dallo stesso file si può specificare il range di mu e sigma con cui utilizzare l'algoritmo e il loro step d'incremento. <br>
delta_finder.cpp stampa i valori di mu, sigma, delta, dell'energia calcolata con metodo a blocchi e relativo errore, della probabilità di accettazione e il numero di tentativi fatti durante la ricerca del delta ottimale. <br>
- Nel file "SA.cpp" è implementato il codice che esegue l'algoritmo del Simulated Annealing, i cui parametri sono caricati dal file "SA_input.in". Oltre al numero di step per equilibrare il campionamento della funzione d'onda e il numero di blocchi da utilizzare, è possibile definire il range di temperatura, dalla massima alla minima, e lo step di decrescita. Si inserisce la coppia di parametri mu e sigma da cui far partire l'algoritmo (trovati con il codice precedente) e il valore di delta utile all'algoritmo SA per il campionamento della distribuzione di Boltzmann. <br>
- L'ultimo parametro booleano permette di utilizzare i valori di mu e sigma inseriti per analizzare quella specifica configurazione, senza effettuare l'algoritmo SA. In questo modo è possibile salvare i dati della distribuzione test trovata con gli specifici valori di mu e sigma che minimizzano l'energia.

### Esercizio 9
Da file esterno è possibile caricare in ordine:
- Numero città da visitare (default 34)
- shape, se 1 le città sono casualmente distribuite su una circonferenza, se 0 IN un quadrato; il diametro del cerchio e il lato del quadrato sono lunghi 1
- dimensione della popolazione
- Numero di generazioni

Il codice esegue 4 diverse mutazioni in sequenza, ognuna con probabilità settabile nel codice (p_mut e p_cross rispettivamente). <br>
- La prima generazione è creata permutando casualmente la sequenza iniziale. <br>
- Ogni generazione successiva è costruita a partire da quella precedente aggiungendo 1 o 2 sequenze per volta (in base a se avviene o no il crossover, per il quale tengo i due figli). <br>
- Il crossover tra una sequenza selezionata e un'altra avviene solo sotto certe condizioni, dunque il secondo genitore candidato viene preso a partire dal primo individual della popolazione precedente e poi a seguire provando gli altri (finché non avviene).
- Per fare in modo che il secondo genitore sia una buona configurazione, ogni popolazione viene riordinata in ordine crescente, in modo che i canditati siano sempre buoni. <br>
- Vi è inoltre una probabilità di sopravvivenza del 50\% che tiene le configurazioni figlie anche se sono più lunghe del genitore, sia per le mutazioni che per il crossover.
- Con questa probabilità, nel crossover se il secondo figlio ha lunghezza maggiore del primo genitore aggiungo quest'ultimo e il primo figlio, altrimenti aggiungo entrambi i figli.

### Esercizio 10
ATTENZIONE: per motivi di spazio, i dati relativi al TSP con città Americane non è stato caricato su GitHub. <br>
Gli input modificabili sono identici a quelli dell'esercitazione 9:<br>
Vi sono le seguenti modifiche dopo parallelizzazione dell'algoritmo genetico con prodocollo MPI:
- Il Random Number Generator viene modificato per aggiungere un costruttore MPI, che utilizza un seed diverso per ogni core
- Il codice viene ora compilato con mpicxx ed eseguito con mpiexec
- Per eseguire il codice digitare "make run", il numero di processi può essere cambiato dal Makefile alla prima variabile "N_PROCESSES"
- Impostando nel file Input.it shape=2 carico un file di posizioni arbitrarie da chiamare "Positions.dat"
- Settando shape=3 carico le coordinate delle città americane. (Il calcolo della distanza è svolto in un altro modo perché le coordinate sono in longitudine e latitudine. Inoltre nel file delle città americane i loro nomi sono stati modificati per evitare spaziature indesiderate)
- Settando shape=4 carico le coordinate delle città italiane.
- Il numero di processi massimo utilizzabile è 10 e il minimo 2
- Il rank 0 genera ANCHE le coppie di core che comunicano
(uno per generare le coppie di core che si scambiano gli individual e gli altri due che se li scambiano. Se sono 3 sicuramente gli altri due se li scambiano sempre)
- Modificato il PrintDistance (funzione che stampa tutte le distanze per ogni generazione) in modo di farlo solo alla fine, cosi evito di aprire 1000 stream
- Il file con le sequenze minime e le distanze lo stampo per ogni core (avrò $2 \\times $ rank files)
