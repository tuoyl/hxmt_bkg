#!/bin/bash

# initiating HXMT environment
# description: in this script, three scripts generating HXMT backgroud
#              will be added to HXMT toolkit by creating three symlinks
#              from $HEADAS/refdata directory.
# Set environment paramter $REFPATH



echo "HEADAS DIRECTORY:" $HEADAS 

# check existence of bkgmap python scripts
FILE="./hebkgmap.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/hebkgmap ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/hebkgmap
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi

FILE="./mebkgmap.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/mebkgmap ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/mebkgmap
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi

FILE="./lebkgmap.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/lebkgmap ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/lebkgmap
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi


FILE="./hen2pi.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/hen2pi ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/hen2pi
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi

FILE="./hprint_detid.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/hprint_detid ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/hprint_detid
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi


FILE="./hspec_merge.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/hspec_merge ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/hspec_merge
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi


FILE="./legti.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/legti ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/legti
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi


FILE="./megti.py"
if [ -f $FILE ]; then
    if [[ -x "$FILE" ]]; then
        echo "symlink $FILE"
        if [ -f $HEADAS/bin/megti ]; then
            :
        else
            ln -s $FILE $HEADAS/bin/megti
        fi
    else
        echo "Permission Denied, File $FILE is not executable"
        exit 1
    fi
else
    echo "ERROR:$FILE does not exist:"
    exit 1
fi









