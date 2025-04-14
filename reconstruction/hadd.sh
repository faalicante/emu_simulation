b=(11 21 31 41 51) #ehi, I have written an array in BASH!
for i in $(seq 1 4); do
	brick="$(printf "b%0*d" 6 ${b[$i]})"
	for p in $(seq 1 60); do
		plate="$(printf "p%0*d" 3 $p)"
		mkdir -p $brick/$plate
		hadd $brick/$plate/${b[$i]}.$p.0.0.cp.root {0..122}/$brick/$plate/${b[$i]}.$p.0.0.cp.root 
	done
done
