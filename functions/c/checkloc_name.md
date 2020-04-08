# checkloc_name
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v2

__Overview__
This function checks for the presence of the variable ``loc_name``. ``loc_name`` is a SNP naming convention that concatenates the chromosome, base position and the IUPAC genotype code. The IUPAC codes are defaulted to ``R`` for ``G|A`` and ``C|T``, and  ``M`` for ``C|A`` and ``G|T`` SNPs. This prevents strand specific names. If ``loc_name`` is absent, it uses the variable ``chr``, ``bp`` and ``gt`` to create the variable.

The ``loc_name`` format is ``chr``:``bp``-``gt*``

|chr|snp|bp|a1|a2|gt|loc_name|
|--------|--------|--------|--------|--------|--------|--------|
|1|rs3131972|752721|A|G|R|``chr1:752721-R``
|1|rs11240777|798959|T|C|Y|``chr1:798959-R`` | - note the ``Y`` genotype is flagged as ``R`` in the ``loc_name``
|1|rs4970383|838555|A|C|M|``chr1:838555-M``


__Dependencies__
The function has the following dependencies are ;
``recodegenotype``

__Syntax__
```
checkloc_name
```

