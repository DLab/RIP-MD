#! /bin/bash 

wget -q --tries=10 --timeout=20 --spider http://google.com
if [[ $? -eq 0 ]]; then
	echo "Looking for dependencies..."
	toInstallPip="0"
	hash pip || toInstallPip="1"
	if [[ $toInstallPip == "1" ]]; then
		echo "    Installing pip..."; 
		errorAPT="0"
		apt-get install python-pip python-dev build-essential || errorAPT="1"
		if [[ $errorAPT == "1" ]]; then
			echo "	Error while trying to install pip, are you in sudo mode?. Exiting"
			exit 1
		fi
	fi
	
	echo "Installing python-tk"
	apt-get install python-tk
	
	errorPIP="0"
	pip install --upgrade pip || errorPIP="1"
	if [[ $errorPIP == "1" ]]; then
		echo "	Error upgrading pip, are you in sudo mode?. Exiting"
		exit 1
	fi			
	pip install -r libraries || errorPIP="1"
	if [[ $errorPIP == "1" ]]; then
		echo "	Error installing pip libraries, are you in sudo mode?. Exiting"
		exit 1
	fi	
else
	echo "Can not connect to internet. Exiting."
	exit 1
fi
