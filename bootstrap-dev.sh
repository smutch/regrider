echo "If you haven't already, please be sure and run\n\tspack compiler find && spack external find\n"

while true; do
    read -p "Do you wish to setup the spack env with all requirements?" yn
    case $yn in
        [Yy]* ) make install; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

spack env create -d ./spack-env spack.yaml
spack env activate ./spack-env
spack concretize
spack install
