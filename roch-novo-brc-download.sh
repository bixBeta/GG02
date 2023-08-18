#!/bin/bash

#SBATCH -J brc_roch_novo_download
#SBATCH -o %x.out


# ------------------------------------------------------------------------------------------------------------------------
# Download Script for BRC Novogene and Rochester:
#	Downloads, Archives and Copy to TREx_shared/projects

# Archive Path: /home/RSCshare/RSC/Projects/dataArchive/

# Archive Dir Structure
# .
# ├── logfiles
# ├── raw
# │	└── < sequencing provider distribution date format YYMMDD > e.g. 230209
# └── reports



# 
# ------------------------------------------------------------------------------------------------------------------------

WD=`pwd`

usage(){

    echo "BRC // Rochester // Novogene Download - @bixBeta"
    echo ""
    echo "hu:p:b:d:l:"
    echo "Usage: bash" $0 "[-h arg] [-u arg] [-p args] [-b arg] [-d arg] [-l arg] "
    echo
    echo "---------------------------------------------------------------------------------------------------------------------------"
    echo ""
    echo "[-h] --> Display Help"
    echo "[-p] --> Comma sep list of all projects in the ftp download e.g. < 1002D,1004A,1024S > "
    echo "[-u] --> Comma sep list of Novogene ftp username and pw e.g. < username,password > "
    echo "[-b] --> BRC download < yes or no >"
    echo "[-d] --> Data distibution date from BRC, Roch or Novo < YYMMDD >"
    echo "[-l] --> Roch link"
    echo ""
    echo "---------------------------------------------------------------------------------------------------------------------------"
    echo ""


}


dirCheck(){

    # check if project dir exists
    
    cat $WD/projects.list | while read line
    do
            DIR="/home/RSCshare/RSC/Projects/dataArchive/500x/${line}*" 
            
            if [ -d /home/RSCshare/RSC/Projects/dataArchive/500x/$line*  ]; then

                echo "-----------------------------------------------------------------------------"
                echo "Project Directory Exists, adding reseqs"
                echo "-----------------------------------------------------------------------------"

                cd $DIR/raw
                mkdir -m 777 $DATE
                cd $DATE

                pwd >> $WD/archive.Path

            else

                echo "-----------------------------------------------------------------------------"
                echo "Project Directory Does not Exist!!! Creating dirs logfiles; raw/$DATE; reports"
                echo "-----------------------------------------------------------------------------"

                cd /home/RSCshare/RSC/Projects/dataArchive/500x/
                mkdir -m 777 -p $line/logfiles $line/raw $line/reports
                cd $line/raw
                mkdir -m 777 $DATE && cd $DATE

                pwd >> $WD/archive.Path
            fi

    done

}

donwload.novo(){
	echo "N O V O G E N E   D O W N L O A D"
	echo "User = $USER, PW = $PW "
	echo ""
	wget -r -c ftp://${USER}:${PW}@usftp21.novogene.com:21/

    cd usftp21.novogene.com

    md5sum 01.RawData/*/*.gz >> myMD5s.txt
    sort -k2 -n MD5.txt > sorted.MD5.txt
    sort -k2 -n myMD5s.txt > sorted.myMD5s.txt
    diff sorted.MD5.txt sorted.myMD5s.txt > diff.command.out
    
    cd 01.RawData
    mkdir -m 777 downloaded_fastqs
    mv */*.gz downloaded_fastqs

    cd downloaded_fastqs 
    cat $WD/projects.list | while read fq
        do  
            mkdir -m 777 $fq
            mv `find . -type f | grep $fq` $fq 
        done
    
    cat $WD/projects.list | while read fq
        do
            rsync -avR $fq `grep "$fq" $WD/archive.Path`
            rsync -arctuxzv --remove-source-files $fq /workdir/TREx_shared/projects/
        done

}


donwload.roch(){
    echo "R O C H E S T E R   D O W N L O A D"
    echo "Roch URL -- ${URL}"
    echo ""
    wget -r -c ${URL}
    md5sum *tar.gz >> TAR.md5

    tar -xzvf *tar.gz
    cd deliv*

    md5sum *.gz >> myMD5s.txt
    sort -k2 -n md5sums.log > sorted.MD5.txt
    sort -k2 -n myMD5s.txt > sorted.myMD5s.txt
    diff sorted.MD5.txt sorted.myMD5s.txt > diff.command.out

    cat $WD/projects.list | while read fq
        do  
            mkdir -m 777 $fq
            mv `find . -type f | grep $fq` $fq 
        done
    
    cat $WD/projects.list | while read fq
        do
            rsync -avR $fq `grep "$fq" $WD/archive.Path`
            rsync -arctuxzv --remove-source-files $fq /workdir/TREx_shared/projects/
        done

}


donwload.BRC(){

    echo "B R C   D O W N L O A D"
    `cat download.sh`

    cat $WD/projects.list | while read fq
        do  
            mkdir -m 777 $fq
            mv `find . -type f | grep $fq` $fq 
        done
    
    cat $WD/projects.list | while read fq
        do
            rsync -avR $fq `grep "$fq" $WD/archive.Path`
            rsync -arctuxzv --remove-source-files $fq /workdir/TREx_shared/projects/
        done


}




# PROJ=$1
# WD=`pwd`
# DATE="YYMMDD"
# dirCheck
# rm $WD/projects.list


#  `find ${DIR} -type d | wc -l` -gt 0


# ------------------------------------------------------------------------------------------------------------------------
# CASE ESAC FLOW

while getopts "hu:p:b:d:l:" opt; do
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

    b)

        BRC=$OPTARG

    ;;

    d) 
        DATE=$OPTARG

    ;;

    l)

        URL=$OPTARG

    ;;

    \?)
        echo
        echo
        echo
        usage

    ;;

    esac
done

# PARAMETER CHECKS 

if [[ ! -z "${PROJ+x}" ]]; then

    echo $PROJ | sed -e 's/,/\n/g' >> $WD/projects.list 
    dirCheck
fi

# split user and pw

if [[ ! -z "${USERPW+x}" ]]; then

    echo 
    echo

    USER=`echo $USERPW | cut -d , -f1`
    PW=`echo $USERPW | cut -d , -f2 `

    # dirCheck
    donwload.novo


fi


if [[ ! -z "${URL+x}" ]]; then

    echo 
    echo "Downloading data from Rochester"
    echo

    # dirCheck
    donwload.roch

fi

if [[ ! -z "${BRC+x}" ]]; then

    echo 
    echo "Downloading data from BRC"
    echo

    # dirCheck
    donwload.BRC

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
    echo >> download.log
    echo `date` >> download.log
    echo "Projects included in the download = " $PROJ >> download.log



fi