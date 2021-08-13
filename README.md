# EventIdentificationUsingTMVA
The SiFi-CC machine learning model identifies true Compton events. Different models available in TMVA compete and their performances are compared. Moreover, the energy regression is available to correct the total deposited energies of the predicted Compton events.


Prerequisites
------------------------------------------------
* Minimum Required ROOT version: 6.13/01, but better 6.18.04 or later

Sources
-------

Sources repository:
```
https://github.com/SiFi-CC/Event-Identification-Using-TMVA/

```
To get sources run:

```
git clone https://github.com/SiFi-CC/Event-Identification-Using-TMVA/

```
Building and running
-------------------------
```
cd EventIdentificationUsingTMVA
in the EventIdentificationUsingTMVA directory run:
mkdir build
cd build
cmake ../source
make
./CompId ./path_to_EI.txt
```
