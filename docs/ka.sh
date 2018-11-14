
launch -c python ../bin/tps.py -T 0.6 --field 0.00 --steps 1000 -t 90.0 -i ../data/ka_rho1.2_T0.6.xyz '/tmp/tps/ka/s{field:.3f}/output'
launch -c python ../bin/tps.py -T 0.6 --field 0.02 --steps 2000 -t 90.0 -i ../data/ka_rho1.2_T0.6.xyz '/tmp/tps/ka/s{field:.3f}/output'
launch -c python ../bin/tps.py -T 0.6 --field 0.04 --steps 4000 -t 90.0 -i ../data/ka_rho1.2_T0.6.xyz '/tmp/tps/ka/s{field:.3f}/output'
launch -c python ../bin/tps.py -T 0.6 --field 0.08 --steps 8000 -t 90.0 -i ../data/ka_rho1.2_T0.6.xyz '/tmp/tps/ka/s{field:.3f}/output'

{
for f in /tmp/tps/ka/s* ; do
    echo $(grep -E '^# field' $f/output.log | cut -d : -f 2 | tail -n 1) \
         $(awk '(!/#/&&(NR>200)){print $4}' $f/output.thermo | ave)
done
} > /tmp/tps/ka/K_s.txt
