
#!/bin/env bash
# PerlBrew needs to be installed to manage isolated perl environemnts

if [ perlbrew ]; then
  echo 'Perl brew installed'
else
  echo 'Installing Perl brew...'
  curl -L http://install.perlbrew.pl | bash
  echo 'source ~/perl5/perlbrew/etc/bashrc' >> ~/.bash_profile
  perlbrew init
  perlbrew install perl-5.23.0
  perlbrew use perl-5.23.0
  perlbrew install-cpanm
fi
