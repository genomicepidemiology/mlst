
#!/bin/env bash
#

PERLBREW='http://install.perlbrew.pl'
BLASTLINUX='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-x64-linux.tar.gz'
BLASTMAC='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-universal-macosx.tar.gz'
BLASTFOLDER=blast

# PerlBrew needs to be installed to manage isolated perl environemnts

# FIXME
if [ which `perlbrew` ]; then
    echo 'Perl brew installed'
else
    echo 'Installing Perl brew...'
    curl -L ${PERLBREW} | bash
    echo 'source ~/perl5/perlbrew/etc/bashrc' >> ~/.bash_profile
    perlbrew init
    perlbrew install perl-5.23.0
    perlbrew use perl-5.23.0
    perlbrew install-cpanm
fi

# Installing NCBI Blast tools
# # FIXME
if [ which `blastall` ]; then
    echo 'Blast tools installed'
else
    echo 'Installing Blast tools...'
    curl ${BLASTMAC} -o ${BLASTFOLDER}.tar.gz
    tar -zxvf ${BLASTFOLDER}.tar.gz ${BLASTFOLDER}
    export PATH=$PATH:blast-2.2.26/bin
fi
