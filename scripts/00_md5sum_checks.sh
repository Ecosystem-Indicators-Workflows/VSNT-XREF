# This works in mac with zsh or bash
## Create md5sums with
md5sum data/source-data/*/* > data/metadata/data-checksums.txt
md5sum data/processed_data/*.csv >> data/metadata/data-checksums.txt
md5sum data/xrefdata/*.csv >> data/metadata/data-checksums.txt
md5sum results/*/* >> data/metadata/data-checksums.txt

### check all md5sums with
## md5sum -c data/metadata/data-checksums.txt
