#!/bin/bash

NRUNS=$1
SAMPLE_PATH=$2
HEPMC_NAME=$3

DELPHES_PATH=~/softwares/MG5_aMC_v2_6_7/Delphes

echo -e "#!/bin/bash
#cd $SAMPLE_PATH/Events/run_04_\$1/
#pigz -d tag_1_pythia8_events.hepmc.gz

cd $DELPHES_PATH
./DelphesHepMC /usr/local/share/delphes_cards/delphes_card_CMS_skimmed.tcl $SAMPLE_PATH/Events/run_04_\$1/delphes_\$1.root $SAMPLE_PATH/Events/run_04_\$1/$HEPMC_NAME

cd $SAMPLE_PATH/Events/run_04_\$1/
pigz -1 $HEPMC_NAME
" > run_delphes.sh

for i in $(eval echo "{0..$NRUNS}")
do
bash run_delphes.sh $i &
done
wait
rm run_delphes.sh
