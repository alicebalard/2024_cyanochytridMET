Run on my local machine

conda activate btk
scp alicebalard@curta.zedat.fu-berlin.de:/scratch/alicebalard/outData/blobtools/FINALTran/*json /home/alice/Documents/Erika-chytridProject/GIT/2024_cyanochytridMET/ignoreThinkpad/FINALTran/.

To observe:

blobtools host ~/Documents/Erika-chytridProject/

To plot:

blobtools view --plot --out blobKingdomFINALTran FINALTran/

To output a table:

blobtools filter --table FINALTran.table.tsv --table-fields gc,length,bestsumorder_kingdom,bestsumorder_family,bestsumorder_genus FINALTran
