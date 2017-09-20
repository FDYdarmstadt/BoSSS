#!/bin/sh

OAuth=$2
tag="BoSSS-version$1"
url="https://api.github.com/repos/FDYdarmstadt/BoSSS/releases"

# set release
release_json=$(curl -H "Authorization: token $OAuth" -d '{"tag_name":"'$tag'","body":"This is an automatic release. All NUnit tests were successful."}' $url)

upload=$(echo $release_json | grep -Po '"upload_url": "\K.*?(?=")')
if [[ $? -eq 0 ]]; then
	echo Release published
else 
	echo Error publishing release!
	return
fi
	
upload_url=$(echo $upload | cut -d "{" -f1)
#echo $upload_url
release_id=$(echo $upload | grep -Po 'releases/\K.*?(?=/)')
#echo $release_id


# set assets (installer and handbook)
installer_path="InnoSetup/Output/BoSSS-setup-$1.exe"
installer_asset="$upload_url?name=BoSSS-setup-$1.exe"
download=$(curl -H "Authorization: token $OAuth" -H "Accept: application/vnd.github.v3+json" -H "Content-Type: application/octet-stream" --data-binary @$installer_path $installer_asset)
if [[ $? -eq 0 ]]; then
	echo $download | cut -d: -f2,3 | cut -d\" -f2
else
	echo Upload error!
fi

handbook_path="doc/handbook/BoSSShandbook.pdf"
handbook_asset="$upload_url?name=BoSSS-handbook.pdf"
download=$(curl -H "Authorization: token $OAuth" -H "Accept: application/vnd.github.v3+json" -H "Content-Type: application/pdf" --data-binary @$handbook_path $handbook_asset)
if [[ $? -eq 0 ]]; then
	echo $download | cut -d: -f2,3 | cut -d\" -f2
else
	echo Upload error!
fi