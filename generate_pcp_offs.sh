#!/bin/bash
target="$HOME/self/pi3/dds-web/polyh/pcp/off"
mkdir -p $target
args="-x 70"
if [ -z "$n_val" ]; then
	n_val=5
fi
[ -z $overwrite ] || args="$args --overwrite"
if [ "$n_val" == 3 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
fi
if [ "$n_val" == 4 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
fi
if [ "$n_val" == 5 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -H $n_val 2 ${target}/pcp_${n_val}_2_h.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
fi

if [ "$n_val" == 6 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args    $n_val 2 ${target}/pcp_${n_val}_2_alt.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
fi

if [ "$n_val" == 7 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 6 ${target}/pcp_${n_val}_6.off
fi

if [ "$n_val" == 8 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 6 ${target}/pcp_${n_val}_6.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 7 ${target}/pcp_${n_val}_7.off
fi

if [ "$n_val" == 9 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 6 ${target}/pcp_${n_val}_6.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 7 ${target}/pcp_${n_val}_7.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 8 ${target}/pcp_${n_val}_8.off
fi

if [ "$n_val" == 10 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 6 ${target}/pcp_${n_val}_6.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 7 ${target}/pcp_${n_val}_7.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 8 ${target}/pcp_${n_val}_8.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 9 ${target}/pcp_${n_val}_9.off
fi

if [ "$n_val" == 11 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 6 ${target}/pcp_${n_val}_6.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 7 ${target}/pcp_${n_val}_7.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 8 ${target}/pcp_${n_val}_8.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 9 ${target}/pcp_${n_val}_9.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 10 ${target}/pcp_${n_val}_10.off
fi

if [ "$n_val" == 12 ]; then
	python n_m_cupolaic_prismatoids.py $args -s $n_val 1 ${target}/pcp_${n_val}_1.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 2 ${target}/pcp_${n_val}_2.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 3 ${target}/pcp_${n_val}_3.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 4 ${target}/pcp_${n_val}_4.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 5 ${target}/pcp_${n_val}_5.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 6 ${target}/pcp_${n_val}_6.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 7 ${target}/pcp_${n_val}_7.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 8 ${target}/pcp_${n_val}_8.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 9 ${target}/pcp_${n_val}_9.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 10 ${target}/pcp_${n_val}_10.off
	python n_m_cupolaic_prismatoids.py $args -s $n_val 11 ${target}/pcp_${n_val}_11.off
fi

