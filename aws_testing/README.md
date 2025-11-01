# Running RUFUS on AWS

Running RUFUS on AWS requires obtaining the RUFUS container and the RUFUS resource directory.

## Step 1: Download the RUFUS Container
```bash
curl "https://zenodo.org/records/13871423/files/rufus_stable" -o rufus.sif
```

## Step 2: Download the RUFUS Resource Directory
The RUFUS resource directory is approximately 150GB, and will need to be in the directory (referred to as DATA_DIR below) which gets mounted during the RUFUS container execution (see below).

### Option 1: Download with AWS-CLI
```bash
aws s3 sync s3://rufus.marth.lab/public_access_data/rufus_resources/ ${HOST_DATA_DIR}/rufus_resources --no-sign-request
```

### Option 2: Download with Rclone

```bash
rclone copy :s3:rufus.marth.lab/public_access_data/rufus_resources/ ${HOST_DATA_DIR}/rufus_resources/ --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P
```

### Option 3: Download with wget
```bash
# Download file manifest
curl "https://zenodo.org/records/13871423/files/rufus_resource_manifest.txt"

# wget all the files in the manifest
mkdir -p ${HOST_DATA_DIR}/rufus_resources
while IFS= read -r file; do
    wget -P ${HOST_DATA_DIR}/rufus_resources https://rufus.marth.lab.s3.us-east-1.amazonaws.com/"$file"
done < rufus_resource_manifest.txt

# Check file integrity with MD5 (recommended)
curl "https://zenodo.org/records/13871423/files/rufus_resource_checksums.md5"
md5sum -c rufus_resource_checksums.md5 
```

## Step 3: Run RUFUS
1. **Populate data directory** - place your sample file(s) and the RUFUS resources directory (`rufus_resources`) in the data directory that will be bound to your Docker container (`HOST_DATA_DIR`). Note, the bound data directory can contain soft links if it is undesireable to actually move things around.
2. **Fill out the rufus.env file** - the rufus.env file is located in the `rufus_resources` directory, and tells the launch script where to find required input files and run time arguments
3. **Run RUFUS** - RUFUS is meant to be run in a distributed fashion over 1MB regions of the genome. An example script can be downloaded and run:
    ```bash
    curl "https://zenodo.org/records/13871423/files/launch_parallel_rufus.sh"
    chmod u+x launch_parallel_rufus.sh
    ./launch_parallel_rufus.sh ${OPTIONAL_PATH_TO_ENV_FILE}
    ``` 



## Recommendations for multi-node EC2 running
1. Create attachable volume with 150GB (for RUFUS resource directory) + size of input sample file(s)
2. Attach volume to all EC2 instances
3. Update `launch_rufus.sh` script to only pull `# regions / # EC2 instances`, comment out post-processing command, then execute. 
4. Run post-processing command after all RUFUS runs complete.