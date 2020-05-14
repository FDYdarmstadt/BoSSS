
## BoSSS Database Configuration for working on Supercomuters 
... and other remote computers.


### Manual Deployment on other computers

Früher gab es nur eine Datenbank-Pfad pro Rechung; da musste man den Dantenbank-Pfad immer manuell ändern, wenn man z.B. vom eigenen PC aud den Lichtenberg oder den HPC cluster gewechselt hat


The  `AppControl.AlternateDbPaths` should help working with multiple machines; 
E.g., let
object `C` be some `AppControl` instance:
```
C.AlternateDbPaths = new[] {
  (@"\\dc1\userspace\kikker\cluster\cluster_db\ConfinedCylinder_Drag", "hpccluster"), // Kikkers's account
  (@"d:\Users\kummer\default_bosss_db", "terminal03"),     // Kummer's account on Terminal03
  (@"c:\Users\florian\default_bosss_db", "rennmaschin"),   // Florian's old Laptop
  (@"c:\Users\flori\default_bosss_db", "stormbreaker")     // Florian's new Laptop
};
```
D.h. damit kannst du für verschiedene Computer verschiedene Datenbanken bestimmen.


### Deployment through the job management in BoSSSpad 

Diese alternativen Pfade sind für sich ganz praktisch, 
aber mit dem Job-Manger war es immer noch Gefrickel.
Wenn du mit dem Job-Manager arbeitest, 
solltest du dich um die Datenbank-Pfade gar nicht kümmern müssen.
Dafür müssen aber einige Sachen richtig eingestellt sein:
1. die Konfiguration der Batch-Queues: `~/.BoSSS/etc/BatchProcessorConfig.json`
   Hier ein Beispiel:
   ```
   {
   "AllQueues": [
    {
      "$type": "BoSSS.Application.BoSSSpad.MiniBatchProcessorClient, BoSSSpad",
      "DeploymentBaseDirectory": "C:\Users\flori\AppData\Local\BoSSS-LocalJobs",
      "DeployRuntime": false
    },
    {
      "$type": "BoSSS.Application.BoSSSpad.SlurmClient, BoSSSpad",
      "Username": "fk69umer",
      "ServerName": "lcluster3.hrz.tu-darmstadt.de",
      "PrivateKeyFilePath": "C:\Users\flori\.ssh\id_rsa",
      "DeploymentBaseDirectory": "X:\bosss_deploy",
      "DeployRuntime": false,
      "AllowedDatabasesPaths": [ "X:\bosss_db_lichtenberg" ],
      "SlurmAccount": "project01217",
      "DeploymentBaseDirectoryAtRemote": "/home/fk69umer/bosss_deploy"
    }
   ]
   }
   ```
   Interresant ist hier die zweite, also der Slurm-Client für den Lichtenberg:
   ```
   "DeploymentBaseDirectory": "X:\bosss_deploy"    // Daten müssen im lokal gemouteten Verzeichniss X:\bosss_deploy liegen
   "AllowedDatabasesPaths": [ "X:\bosss_db_lichtenberg" ],  // Die verwendete Datenbank muss auf X:\bosss_db_lichtenberg liegen
   "DeploymentBaseDirectoryAtRemote": "/home/fk69umer/bosss_deploy"  // Adresse des DeploymentBaseDirectory auf dem Lichtenberg
   ```
2. Die Lichtenberg-Verzeichnisse müssen lokal gemounted sein, 
   z.B. mit ExpanDrive oder mit SSHFS-Win
   D.h. ich habe auf meinem PC ein Laufwerk `X:`, das meinem 
   Lichtenberg-Homeverzeichnis `/home/fk69umer` entspricht
   (Wenn man auch noch mit dem Lichtenberg-Scratch arbeiten möchte, 
   braucht man einen zweites Laufwerk, z.b. `Y:`)
3. Die Datenbank muss ihre eigene Adresse auf verschiedenen Computern kennen. 
   Das erreicht man über die 'AlternatePaths.txt' Datei, die im 
   Datrenbank-Verzeichniss selbst liegt
   ```
   X:\bosss_deploy
   ├───data
   ├───grids
   ├───sessions
   │   ├───1522e363-903b-4106-a34e-1d9521a0ded6
   │   ├───1f36290b-a6c8-4e65-925c-b94a6df52372
   │   ├───231bcdbd-346e-4daf-9760-cf6c2b935434
   │   ├───4db65cfe-d1ee-408e-ac6d-b1674af72bce
   │   └───e26158ef-4c0e-45b7-95af-afc63d6652fa
   └───timesteps
   └───AlternatePaths.txt
   ```
   der Inhalt von `AlternatePaths.txt` ist eine Auflistung von Verteichniss-Pfaden, mit optionalen Filtern für die Maschine
   ```
   X:\bosss_db_lichtenberg,stormbreaker
   /home/fk69umer/bosss_db_lichtenberg,
   ```
   PS: das Komma am Ende ist kein Typo!
