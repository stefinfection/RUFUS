# Running RUFUS on AWS

## Step 0: Prerequisites
Running RUFUS requires Docker. Installation instructions for getting Docker on your machine can be found [here](https://docs.docker.com/engine/install/).

## Step 1: Download the RUFUS Container and launch script
```bash
curl "https://zenodo.org/records/13871423/files/rufus_stable" -o rufus.sif
curl "https://zenodo.org/records/13871423/files/launch_rufus.sh" -o launch_rufus.sh
```

## Step 2: Create a data directory with all inputs
RUFUS requires all input files to be in a single directory that gets mounted during exection. One simple way to accomplish this is by soft linking your input file(s) within the directory:
```bash
mkdir ${HOST_DATA_DIR}
ln -s ${SUBJECT_FILE} ${HOST_DATA_DIR}
ln -s ${REFERENCE_FILE} ${HOST_DATA_DIR}
```

## Step 3: Download RUFUS resource directory
RUFUS requires a small resource directory containing an environment file, a regions file, and reference files that allow RUFUS to run quicker. This resource directory **must** be placed in the same data directory from step 2. There are a few ways to easily download these resources.

### Option A: Download with AWS-CLI
```bash
aws s3 sync s3://rufus.marth.lab/public_access_data/rufus_resources/ ${HOST_DATA_DIR}/rufus_resources --no-sign-request
```

### Option B: Download with Rclone
```bash
rclone copy :s3:rufus.marth.lab/public_access_data/rufus_resources/ ${HOST_DATA_DIR}/rufus_resources/ --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P
```

### Option C: Download with wget
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

## Step 4: Fill out the rufus.env file
The `rufus.env` file coordinates passing files into the RUFUS launch script. This environment file must remain `${HOST_DATA_DIR}/rufus_resources/rufus.env`. 

## Step 5: Launch RUFUS
```bash
chmod u+x launch_rufus.sh
./launch_rufus.sh
```

Upon completion of RUFUS, result files will be stored in `${HOST_DATA_DIR}/rufus_calls.vcf.gz`