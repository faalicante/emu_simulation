for i in $(seq 0 122); do
	mv $i partitions
	#rm -r $i
	#cp -r 0 $i
	#ln -s -f /eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_BRICK11_RUN1/8461734/$i/sndLHC.Ntuple-TGeant4.root $i/
	#sed -i "s/XPART/$i/" $i/SND2FEDRA_masterconv_light.sh
done
