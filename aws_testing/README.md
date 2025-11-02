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
aws s3 sync s3://rufus.marth.lab/public_access_data/rufus_resources/ ${DATA_DIR}/rufus_resources --no-sign-request
```

### Option 2: Download with Rclone

```bash
rclone copy :s3:rufus.marth.lab/public_access_data/rufus_resources/ ${DATA_DIR}/rufus_resources/ --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P
```

### Option 3: Download with wget
```bash
# Download file manifest
curl "https://zenodo.org/records/13871423/files/rufus_resource_manifest.txt"

# wget all the files in the manifest
mkdir -p ${DATA_DIR}/rufus_resources
while IFS= read -r file; do
    wget -P ${DATA_DIR}/rufus_resources https://rufus.marth.lab.s3.us-east-1.amazonaws.com/"$file"
done < rufus_resource_manifest.txt

# Check file integrity with MD5 (recommended)
curl "https://zenodo.org/records/13871423/files/rufus_resource_checksums.md5"
md5sum -c rufus_resource_checksums.md5 
```



1. Spin up EC2 instance(s) to run RUFUS. We recommend running ___ machines.
2. Create a data directory that will be mounted to the Docker container during exectution. This can be an external storage volume that mounts onto the EC2 instance. The data volume should be 150GB + whatever the size of the input file(s) going in to the RUFUS call.
3. Download the RUFUS resources directory into the data directory. This usually takes about 20 minutes. The resource directory can be easily downloaded with AWS-CLI, rclone, or wget. A convenience script can be downloaded and run to auto-detect the best way to download the resource directory with: