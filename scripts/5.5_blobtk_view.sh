Run on my local machine

conda activate btk
scp alicebalard@curta.zedat.fu-berlin.de:/scratch/alicebalard/outData/blobtools/Z1Z12assembly/*json /home/alice/Documents/Erika-chytridProject/Z1Z12/.

To observe:

blobtools host ~/Documents/Erika-chytridProject/

To plot:

blobtools view --plot --param plotColour=kingdom --out blobKingdomZ1Z2 Z1Z12/

