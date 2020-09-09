let "x = 0";

#for filename in dataset/*.mtx; do
for filename in polybench/*.dot; do
        # echo "running "$filename;
        # date;

        echo "./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --refinement 0 --inipart 6 -n 10 &" > detailtimemetis\_$filename\_149
        ./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --refinement 0 --inipart 6 -n 10 >> detailtimemetis\_$filename\_149 &
        echo "./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --inipart 6 -n 10 &" > detailtimemetis\_$filename\_150
        ./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --inipart 6 -n 10 >> detailtimemetis\_$filename\_150 &
        let "x += 1";
        if [[ $(( $x % 70 )) == 0 ]]; then
                echo "wait ${x}";
                wait
        fi
done
echo "waiting last batch";
wait

let "x = 0";

for filename in dataset/*.mtx; do
# for filename in polybench/*.dot; do
        # echo "running "$filename;
        # date;

        echo "./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --refinement 0 --inipart 6 -n 10 &" > detailtimemetis\_$filename\_149
        ./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --refinement 0 --inipart 6 -n 10 >> detailtimemetis\_$filename\_149 &
        echo "./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --inipart 6 -n 10 &" > detailtimemetis\_$filename\_150
        ./exe/rMLGP $filename 2 --conpar 0 --co_stop_level 0 --inipart_nrun 1 --inipart 6 -n 10 >> detailtimemetis\_$filename\_150 &
        let "x += 1";
        if [[ $(( $x % 70 )) == 0 ]]; then
                echo "wait ${x}";
                wait
        fi
done
echo "waiting last batch";
wait