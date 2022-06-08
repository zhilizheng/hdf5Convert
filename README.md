# hdf5 convertor
```
# conert hdf5 to VCF format
# Only for a special format, don't try in your own data.
# h5ToVCF INPUT.hdf5  OUPUT.vcf.gz
h5ToVCF impute_runner_chr4.hdf5 test_chr4.vcf.gz
# convert to pgen
plink2 --vcf test_chr4.vcf.gz dosage=DS --make-pgen --out test_chr4
```
