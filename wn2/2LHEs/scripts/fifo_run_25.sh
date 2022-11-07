#!/bin/bash

#export PYTHIA8=/home/soumyadip/Package/Pythia/pythia8306/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHIA8/lib
#export RIVET_ANALYSIS_PATH=$PWD
#export RIVET_DATA_PATH=$PWD

#source /home/soumyadip/Package/Root/Root62208build/bin/thisroot.sh   # root6.22
#source /home/soumyadip/Package/Rivet/Rivet315/rivetenv.sh            # rive 3.1.5

#Run_dir="$PWD"
#ANALYSIS="$PWD"
#export RIVET_ANALYSIS_PATH="/home/soumyadip/Package/Rivet/Rivet315/bin/rivet"

<<comment
rm -rf *.fifo  # delete existing pipes
rm -rf *.log   # delet all logs

mkfifo hepmc.fifo  # initiate pipe

./main89 main89unlops_singleLHE.cmnd hepmc.fifo >> pythia.log &
rivet hepmc.fifo -a MC_ZINC -a MC_ZJETS -o out.yoda &>> rivet.log

rm -rf *.fifo
echo DONE
comment

rm -rf hepmc25.fifo  # delete existing pipes
rm -rf pythia25.log
rm -rf rivet25.log

mkfifo hepmc25.fifo  # initiate pipe

./main89 main89unlops25.cmnd hepmc25.fifo >> pythia25.log &

#rivet hepmc.fifo -a MC_JETS -a ATLAS_2020_I1808726 -o jets.yoda &>> rivet.log
#rivet hepmc.fifo -a MC_JETS -a CMS_2018_I1682495 -o jets.yoda &>> rivet.log
#rivet hepmc.fifo --ignore-beams -a MC_ZKTSPLITTINGS -a MC_ZJETS -a MC_ZINC -a MC_XS -a MC_TAUS -a MC_PHOTONS -a MC_MUONS -a MC_JETS -a MC_ELECTRONS -a CMS_2014_I1305624 -a CMS_2016_I1459051 -a CMS_2017_I1605749 -a CMS_2017_I1635889 -a CMS_2018_I1667854 -a CMS_2018_I1682495 - CMS_2019_I1753680 -a CMS_2020_I1837084 -a ATLAS_2012_I1124167 -a ATLAS_2015_I1393758 -a ATLAS_2018_I1634970 -a ATLAS_2019_I1740909 -a ATLAS_2019_I1724098 -o QCD_4j_LO_18062022.yoda &>> rivet.log

rivet hepmc25.fifo --event-timeout 172800 -a MC_ZKTSPLITTINGS -a MC_ZJETS -a MC_ZINC -a MC_XS -a MC_TAUS -a MC_PHOTONS -a MC_JETS -a MC_ELECTRONS -a CMS_2017_I1635889 -a CMS_2018_I1667854 -a CMS_2019_I1753680 -a CMS_2020_I1837084 -o dyee012j_5f_NLO_UNLOPS25.yoda &>> rivet25.log

rm -rf hepmc25.fifo
echo DONE
