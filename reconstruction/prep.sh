sim_path=/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_BRICK11_RUN1/8461734
for i in $(seq 0 122); do
	cp -r 0 $i
	ln -s -f $sim_path/$i/sndLHC.Ntuple-TGeant4.root $i/
	sed -i "s/XPART/$i/" $i/SND2FEDRA_masterconv_light.sh
done
