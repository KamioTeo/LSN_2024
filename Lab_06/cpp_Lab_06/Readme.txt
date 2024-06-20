Rispetto al codice originale, sono stati aggiunti alcuni parametri d'ingresso impostabili nel file input.dat:

1) restart: 
Se 1, fa partire le simulazioni da una configurazione specifica caricata da file esterno (viene caricata l'ultima configurazione raggiunta, come da richiesta, da file config.final).
Se 0, ricomincia la simulazione da una configurazione di spin casuale.

2) Tsequence:
Se 1, stampa un file con i valori dei dati termodinamici al variare della temperatura, da 0.5 a 2 (a passi di 0.1) come da richiesta.
Se 0, li stampa alla temperatura impostata da file di inpuT.
Se non è né 0 né 1, stampa le informazioni sull'energia misurata (media a blocchi e incertezza statistica) al variare del numero di blocchi utilizzati (da 10 a 100 a passi di 10), fissati 10^6 step.
  
3) Equilibration
Se 1, non esegue il codice ma stampa un file con tutti i valori relativi all'energia del sistema per ogni step eseguito. Serve a monitorare la velocità di equilibrazione del sistema.
Se 0 non equilibra ed esegue il resto del codice come descritto sopra.

