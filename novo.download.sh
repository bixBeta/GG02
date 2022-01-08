#!/bin/sh

#SBATCH -J novoFTP
#SBATCH -o %x.out

usage(){

    echo "Novogene Download - @bixBeta"
    echo ""
    echo ""
    echo "Usage: bash" $0 "[-h arg] [-u arg] [-p args] [-m arg] "
    echo
    echo "---------------------------------------------------------------------------------------------------------------------------"
    echo ""
    echo "[-h] --> Display Help"
    echo "[-u] --> Comma sep list of Novogene ftp username and pw e.g. < username,password > "
    echo "[-p] --> Comma sep list of all projects in the ftp download e.g. < 1002D,1004A,1024S > "
    echo "[-m] --> check md5s < yes or no >"
    echo ""
    echo "---------------------------------------------------------------------------------------------------------------------------"
    echo ""


}


donwload(){
	echo ""
	echo "User = $USER, PW = $PW "
	echo ""
	wget -r -c ftp://${USER}:${PW}@usftp21.novogene.com:21/
}


mdcheck(){

	echo ""

	cd usftp21.novogene.com/
	md5sum raw_data/*/*.gz >> myMD5s.txt

	sort MD5.txt > sorted.MD5.txt
	sort myMD5s.txt > sorted.myMD5s.txt

	diff sorter.MD5.txt sorted.myMD5s.txt > diff.command.out

    cd ..

    mv usftp21.novogene.com $MASTER
}


# need username
# need password
# need md5 checks
#

while getopts "hu:p:m:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        echo
        exit 1

    ;;

    u )

        USERPW=$OPTARG

    ;;

    p )

        PROJ=$OPTARG

    ;;

    m)

        MD=$OPTARG

    ;;


    \?)
        echo
        echo
        echo
        usage

    ;;

    esac
done

# check if PROJECT parameter provided

if [[ ! -z "${PROJ+x}" ]]; then

	MASTER=`echo $PROJ | sed 's/,/-/g'`

fi

# split user and pw

if [[ ! -z "${USERPW+x}" ]]; then

	USER=`echo $USERPW | cut -d , -f1`
	PW=`echo $USERPW | cut -d , -f2 `
	donwload

fi


## check if md5sum True

if [[ ! -z "${MD+x}" ]]; then
    if [[ $MD = "yes" ]]; then
    mdcheck
  else
    echo "md5sum check is not True"
  fi
fi



if [[ -z $1 ]] || [[  $1 = "help"  ]] || [[  $1 = "--help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    echo
    exit 1

else
    echo >> ftp.log
    echo `date` >> ftp.log
    echo "Projects included in the download = " $MASTER >> ftp.log
    echo "md5s checked  = " $MD >> ftp.log



fi
