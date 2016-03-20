mkdir -p ACDY8
mkdir -p AVPR1A
mkdir -p EPHA7
mkdir -p GALM
mkdir -p GATA2
mkdir -p PCDH7
mkdir -p PCDHA1
mkdir -p SLC6A4
mkdir -p TRPA1
mkdir -p UGT8.2
mkdir -p UGT8
mkdir -p UNC5C
mkdir -p ZDHHC11


cd ACDY8
python ../../../to_midi/fasta2midi.py --fasta ACDY8.fasta
cd ..

cd AVPR1A
python ../../../to_midi/fasta2midi.py --fasta AVPR1A.fasta
cd ..

cd EPHA7
python ../../../to_midi/fasta2midi.py --fasta --fasta EPHA7.fasta
cd ..

cd GALM
python ../../../to_midi/fasta2midi.py --fasta --fasta GALM.fasta
cd ..

cd GATA2
python ../../../to_midi/fasta2midi.py --fasta GATA2.fasta
cd ..

cd PCDH7
python ../../../to_midi/fasta2midi.py --fasta PCDH7.fasta
cd ..

cd PCDHA1
python ../../../to_midi/fasta2midi.py --fasta PCDHA1.fasta
cd ..

cd SLC6A4
python ../../../to_midi/fasta2midi.py --fasta SLC6A4.fasta
cd ..

cd TRPA1
python ../../../to_midi/fasta2midi.py --fasta TRPA1.fasta
cd ..

cd UGT8.2
python ../../../to_midi/fasta2midi.py --fasta UGT8.2.fasta
cd ..

cd UGT8
python ../../../to_midi/fasta2midi.py --fasta UGT8.fasta
cd ..

cd UNC5C
python ../../../to_midi/fasta2midi.py --fasta UNC5C.fasta
cd ..

cd ZDHHC11
python ../../../to_midi/fasta2midi.py --fasta ZDHHC11.fasta
cd ..

