
#!/bin/env bash
#

PERLBREW='http://install.perlbrew.pl'
BLASTLINUX='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-x64-linux.tar.gz'
BLASTMAC='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-universal-macosx.tar.gz'
BLASTFOLDER=blast

# PerlBrew needs to be installed to manage isolated perl environemnts if missing
command -v perlbrew >/dev/null 2>&1 || {
    echo 'Installing Perl brew...'
    curl -L ${PERLBREW} | bash
    echo 'source ~/perl5/perlbrew/etc/bashrc' >> ~/.bash_profile

    perlbrew init

    echo "Do you want to install a local perl? [Y]/[N]"
    read answer
    if  [ $answer == 'Y' ]; then
        echo 'Installing perl-5.23.0';
        perlbrew install perl-5.23.0
        perlbrew use perl-5.23.0
    else
        echo 'Local perl will be used';
    fi
}

echo "Do you want to install a cpanmin locally through perlbrew? [Y]/[N]"
read answer
if  [ $answer == 'Y' ]; then
    echo 'Installing perlbrew install-cpanm';
    perlbrew install-cpanm
else
  curl -L https://cpanmin.us | perl - --sudo App::cpanminus
fi

# Installing NCBI Blast tools if missing
command -v blastall >/dev/null 2>&1 || {
    echo 'Installing Blast tools...'
    curl ${BLASTMAC} -o ${BLASTFOLDER}.tar.gz
    tar -zxvf ${BLASTFOLDER}.tar.gz
    export PATH=$PATH:blast-2.2.26/bin
}
