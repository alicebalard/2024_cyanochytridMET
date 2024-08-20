Run on my local machine

conda activate btk
scp alicebalard@curta.zedat.fu-berlin.de:/scratch/alicebalard/outData/blobtools/Z1Z12/*json /home/alice/Documents/Erika-chytridProject/Z1Z12/.

To observe:

blobtools host ~/Documents/Erika-chytridProject/

To plot:

blobtools view --plot --out blobKingdomZ1Z2 Z1Z12/

To output a table:

blobtools filter --table Z1Z12.table.tsv --table-fields gc,length,bestsumorder_kingdom,bestsumorder_family,bestsumorder_genus Z1Z12
blobtools filter --table In1In12.table.tsv --table-fields gc,length,bestsumorder_kingdom,bestsumorder_family,bestsumorder_genus In1In12

