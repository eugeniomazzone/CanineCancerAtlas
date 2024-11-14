#! /bin/bash

GENOME='/genome/canFam3.fa.gz'
VCF='/vcf1/Filtred_Published1.vcf.gz'

NORMAL=$1
TUMOR=$2

gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '1_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr1 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '2_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr2 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '3_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr3 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '4_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr4 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '5_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr5 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '6_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr6 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '7_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr7 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '8_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr8 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '9_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr9 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '10_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr10 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '11_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr11 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '12_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr12 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '13_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr13 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '14_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr14 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '15_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr15 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.')  -max-mnp-distance 0 -O '16_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr16 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR -NORMAL $(echo $NORMAL | cut -f1 -d'.')  -max-mnp-distance 0 -O '17_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr17 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '18_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr18 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '19_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr19 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '20_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr20 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '21_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr21 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '22_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr22 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '23_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr23 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '24_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr24 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '25_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr25 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '26_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr26 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '27_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr27 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '28_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr28 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '29_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr29 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '30_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr30 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '31_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr31 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '32_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr32 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '33_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr33 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '34_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr34 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '35_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr35 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '36_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr36 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '37_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr37 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O '38_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chr38 & \
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL -I 'BQRS/'$TUMOR  -NORMAL $(echo $NORMAL | cut -f1 -d'.') -max-mnp-distance 0 -O 'UN_'$NORMAL'-'$TUMOR'.vcf.gz' --panel-of-NORMALs  'pon.vcf.gz' -L chrX 
wait
gatk FilterMutectCalls -R $GENOME -V '1_'$NORMAL'-'$TUMOR'.vcf.gz' -O '1fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '2_'$NORMAL'-'$TUMOR'.vcf.gz' -O '2fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '3_'$NORMAL'-'$TUMOR'.vcf.gz' -O '3fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '4_'$NORMAL'-'$TUMOR'.vcf.gz' -O '4fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '5_'$NORMAL'-'$TUMOR'.vcf.gz' -O '5fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '6_'$NORMAL'-'$TUMOR'.vcf.gz' -O '6fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '7_'$NORMAL'-'$TUMOR'.vcf.gz' -O '7fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '8_'$NORMAL'-'$TUMOR'.vcf.gz' -O '8fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '9_'$NORMAL'-'$TUMOR'.vcf.gz' -O '9fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '10_'$NORMAL'-'$TUMOR'.vcf.gz' -O '10fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '11_'$NORMAL'-'$TUMOR'.vcf.gz' -O '11fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '12_'$NORMAL'-'$TUMOR'.vcf.gz' -O '12fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '13_'$NORMAL'-'$TUMOR'.vcf.gz' -O '13fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '14_'$NORMAL'-'$TUMOR'.vcf.gz' -O '14fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '15_'$NORMAL'-'$TUMOR'.vcf.gz' -O '15fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '16_'$NORMAL'-'$TUMOR'.vcf.gz' -O '16fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '17_'$NORMAL'-'$TUMOR'.vcf.gz' -O '17fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '18_'$NORMAL'-'$TUMOR'.vcf.gz' -O '18fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '19_'$NORMAL'-'$TUMOR'.vcf.gz' -O '19fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '20_'$NORMAL'-'$TUMOR'.vcf.gz' -O '20fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '21_'$NORMAL'-'$TUMOR'.vcf.gz' -O '21fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '22_'$NORMAL'-'$TUMOR'.vcf.gz' -O '22fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '23_'$NORMAL'-'$TUMOR'.vcf.gz' -O '23fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '24_'$NORMAL'-'$TUMOR'.vcf.gz' -O '24fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '25_'$NORMAL'-'$TUMOR'.vcf.gz' -O '25fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '26_'$NORMAL'-'$TUMOR'.vcf.gz' -O '26fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '27_'$NORMAL'-'$TUMOR'.vcf.gz' -O '27fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '28_'$NORMAL'-'$TUMOR'.vcf.gz' -O '28fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '29_'$NORMAL'-'$TUMOR'.vcf.gz' -O '29fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '30_'$NORMAL'-'$TUMOR'.vcf.gz' -O '30fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '31_'$NORMAL'-'$TUMOR'.vcf.gz' -O '31fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '32_'$NORMAL'-'$TUMOR'.vcf.gz' -O '32fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '33_'$NORMAL'-'$TUMOR'.vcf.gz' -O '33fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '34_'$NORMAL'-'$TUMOR'.vcf.gz' -O '34fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '35_'$NORMAL'-'$TUMOR'.vcf.gz' -O '35fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '36_'$NORMAL'-'$TUMOR'.vcf.gz' -O '36fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '37_'$NORMAL'-'$TUMOR'.vcf.gz' -O '37fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V '38_'$NORMAL'-'$TUMOR'.vcf.gz' -O '38fi_'$NORMAL'-'$TUMOR'.vcf.gz'
gatk FilterMutectCalls -R $GENOME -V 'UN_'$NORMAL'-'$TUMOR'.vcf.gz' -O 'UNfi_'$NORMAL'-'$TUMOR'.vcf.gz'

