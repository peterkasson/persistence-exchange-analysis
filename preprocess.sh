#!/bin/sh

a="a CN1 | a CN2 | a CN3 | a N | a CA | a CB | a OA | a P | a OB | a OC | a OD \n
name 12 HEADGROUPS\n
\"HEADGROUPS\" & \"INNER_TOP\"\n
name 13 TOP_HEADGROUPS\n
\"HEADGROUPS\" & \"INNER_BOT\"\n
name 14 BOT_HEADGROUPS\n
\"TOP_HEADGROUPS\" | \"BOT_HEADGROUPS\"\n
name 15 INNER_HEADGROUPS\n
\"INNER_HEADGROUPS\" | \"SOL_INT\"\n
name 16 ANALYZE\n
quit"


for i; do
    echo $i
    echo $a | make_ndx -f conf.gro -n index.ndx -o index-analyze.ndx
    cat index-analyze.ndx solint-remaining.ndx > index-analyze-remaining.ndx
done

