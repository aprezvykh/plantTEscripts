#!/bin/bash

gen=D_multi_fasta.fa

aa=$(cat $gen | grep "Fatima" | sed 's/[<>,]//g')
for f in $aa;
	do
		samtools faidx D_multi_fasta.fa $f >> D_RLG_Tura_Fatima_B_210J24-1.fa
	done 


aa=$(cat $gen | grep "Atau_Sabrina" | sed 's/[<>,]//g')
for f in $aa;
        do
                samtools faidx D_multi_fasta.fa $f >> D_RLG_Atau_Sabrina_AY368673-2.fa

        done


aa=$(cat $gen | grep "Tura_Sabrina" | sed 's/[<>,]//g')
for f in $aa;
        do
                samtools faidx D_multi_fasta.fa $f >> D_RLG_Tura_Sabrina_210J24-1.fa
        done


aa=$(cat $gen | grep "Tura_Angela" | sed 's/[<>,]//g')
for f in $aa;
        do
                samtools faidx D_multi_fasta.fa $f >> D_RLC_Tura_Angela_210J24-3.fa
        done


aa=$(cat $gen | grep "Atau_Angela" | sed 's/[<>,]//g')
for f in $aa;
        do
                samtools faidx D_multi_fasta.fa $f >> D_RLC_Atau_Angela_AF497474-3.fa
        done

