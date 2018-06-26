#! /bin/bash 

################ Function to install RIP-MD standalone core #############
function standalone {
	echo "Installing standalone core..."

	echo "Checking for previously RIP-MD installation. Press enter to continue..."	
	read aux
	file="$HOME/.bashrc"
	while IFS= read -r var
	do
		if [[ "$(echo $var | grep 'RIP_MD' )" != "" ]]; then
			echo "It seems like RIP-MD is current installed on your system. Please uninstall it before to continue."
			exit 1
		fi				
	done < "$file"		
	
		
	installationFolder=""
	option=""
	while [[ $option == "" ]]
		do
			echo ""
			echo -n "Do you want to install the standalone core of RIP-MD in default folder ($HOME). Options [Y/N]: "
			read option
			if [[ $option == "Y" ]] || [[ $option == "y" ]]
				then
					installationFolder=$HOME
			elif [[ $option == "N" ]] || [[ $option == "n" ]]
				then
					installationFolder=""
					while [[ $installationFolder == "" ]]; do
						echo -n "Select a new installation folder: "
						read installationFolder
						if [ -d "$installationFolder" ]; then
							# Control will enter here if $DIRECTORY exists.
							if [ -d "$installationFolder/RIP-MD" ]; then
								echo "It seems like RIP-MD is current installed on your system. exiting..."
								exit 1
							fi						
						else
							echo "Destination folder does not exist."
							installationFolder=""
						fi
					done
			else
				echo "Incorrect option...."
				option=""
			fi
		done
	

	echo "Copying RIP-MD standalone core..."
	if [ -d "$installationFolder/RIP-MD" ]; then
		echo "It seems like RIP-MD is current installed on your system. exiting..."
		exit 1
	fi
	errorCP=""
	cp -r RIP-MD/ $installationFolder || errorCP="1"
	if [[ $errorCP == "1" ]]; then
		echo "Error while copying RIP-MD folder. Exiting"
		exit 1
	fi
	ripFinalFolder="$installationFolder/RIP-MD"

	if [ -d "$installationFolder/RIP-MD" ]; then
		chmod -R 777 $ripFinalFolder		
		echo "RIP-MD standalone core was installed succefully"
		echo "Creating environment variable of RIP-MD in .bashrc file"
		echo "export RIP_MD=$installationFolder/RIP-MD" >> $HOME/.bashrc	
	else
		echo "RIP-MD standalone core was not installed. Maybe you should run it in sudo mode"
	fi

}

################ VMD plugin interface #############

function vmd_installation {
	echo "Installing plugin interface. Press enter to continue..."	
	read aux
	file="$HOME/.bashrc"
	vmdpath="/usr/local/lib/vmd"
	vmdOption=""
	
	while [[ $vmdOption == "" ]]
	do
		echo -n "Is VMD installed in $vmdpath? [Y/N]: "
		read vmdOption
		if [[ $vmdOption != "Y" ]] && [[ $vmdOption != "y" ]] && [[ $vmdOption != "N" ]] && [[ $vmdOption != "n" ]]; then
			vmdOption=""
		else
			if [[ $vmdOption == "N" ]] || [[ $vmdOption == "n" ]]; then
				echo -n "Enter the VMD installation directory: "
				read vmdpath
				if [ -d "$vmdpath" ]; then
					# Control will enter here if $DIRECTORY exists.
					if [ -d "$vmdpath/RIP-MD" ]; then
							echo "It seems like RIP-MD is current installed on your system. exiting..."
							exit 1
					fi						
				else
					echo "Destination folder does not exist."
					vmdpath="/usr/local/lib/vmd"
					vmdOption=""
				fi
			fi
		fi
	done


	while IFS= read -r var
	do
		if [[ $var == 'vmd_install_extension ripmd ripmd_tk "Analysis/RIP-MD"' ]]; then
			echo "It seems like RIP-MD VMD plugin interface is current installed on your system. Please uninstall it before to continue."
			exit 1
		fi				
	done < "$vmdpath/.vmdrc"

	echo "Copying RIP-MD plugin interface core..."
	if [ -d "$vmdpath/RIP-MD" ]; then
		echo "It seems like RIP-MD is current installed on your system. exiting..."
		exit 1
	fi
	vmdRIPMDpath="$vmdpath/plugins/noarch/tcl/RIP-MD"
	
	errorRIPVMD=""
	cp -r RIP-MD_vmd/ $vmdRIPMDpath || errorRIPVMD="1"
	if [[ $errorRIPVMD == "1" ]]; then
		echo "Error while copying RIP-MD plugin interface. Exiting"
		exit 1
	fi
    
	if [ -d "$vmdRIPMDpath" ]; then
		chmod -R 777 $vmdRIPMDpath
		echo "RIP-MD VMD plugin interface was installed succefully"
		aux=""
		echo
		echo 'vmd_install_extension ripmd ripmd_tk "Analysis/RIP-MD"' >> "$vmdpath/.vmdrc" || echo "Error, please add the following line to $vmdpath/.vmdrc: vmd_install_extension ripmd ripmd_tk \"Analysis/RIP-MD\""

	else
		echo "RIP-MD VMD plugin interface was not installed. Maybe you should run it in sudo mode"
	fi
	
	
}

############# main ##########################
echo "####################################"
echo "##                                ##"
echo "##        RIP-MD INSTALLER        ##"
echo "##                                ##"
echo "####################################"
echo
echo "Options"
echo "[1] Install all"
echo "[2] Install Standalone  core"
echo "[3] Install Plugin interface"
echo
echo -n "Choose your option [1-3]: "
menu_option=""
while [[ $menu_option == "" ]]
do
	read menu_option
	if [[ $menu_option == "1" ]]; then
		standalone
		vmd_installation
		exit 0
	elif [[ $menu_option == "2" ]]; then
		standalone
		exit 0
	elif [[ $menu_option == "3" ]]; then
		vmd_installation
		exit 0
	else
		echo -n "Invalid option. Choose your option [1-6]: " 
		menu_option=""
	fi
done

echo "DONE"
