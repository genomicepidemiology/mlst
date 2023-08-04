  GNU nano 6.2                                                                                        download-virulence-db.sh
#!/usr/bin/env bash

echo "Downloading lastest version of the mlst database to current directory..."

mkdir mlst_db
cd mlst_db

wget https://bitbucket.org/genomicepidemiology/mlst_db/get/master.tar.gz
tar -xvf master.tar.gz --strip-components 1

echo "Installing the mlst database with KMA"
python INSTALL.py

echo "The mlst database has been downloaded and installed."

exit 0
