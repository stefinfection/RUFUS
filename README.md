# Running RUFUS on AWS

## Step 0: Prerequisites
Running RUFUS requires Docker. Installation instructions for getting Docker on your machine can be found [here](https://docs.docker.com/engine/install/).

## Step 1: Fetch RUFUS image from Docker Hub and Setup Scripts from Zenodo
```bash
docker pull stefinfection/rufus:latest
mkdir rufus_resources
curl "https://zenodo.org/records/17510792/files/launch_rufus.sh?download=1" -o launch_rufus.sh
curl "https://zenodo.org/records/17510792/files/process_region_worker.sh?download=1" -o rufus_resources/process_region_worker.sh
curl "https://zenodo.org/records/17510792/files/grch38_1mb_regions.txt?download=1" -o rufus_resources/grch38_1mb_regions.txt
curl "https://zenodo.org/records/17510792/files/rufus.env?download=1" -o rufus_resources/rufus.env
```

## Step 2: Create a Data Directory With All Input Files
RUFUS requires all input files (including the above downloaded `rufus_resources`) to be in a single directory that gets mounted during exection. Since subject files can be very large, and copying them may not be ideal, one simple way to accomplish this is by soft linking. (WARNING: soft-linking only works at the top level of the mounted directory. If you want to use sub-directories within your HOST_DATA_DIR for your source files, be sure to hard link instead):
```bash
mkdir ${HOST_DATA_DIR}

# Move downloaded resrouces to data directory
mv rufus_resources ${HOST_DATA_DIR}

# Soft link any required large files to data directory
ln -s ${SUBJECT_FILE} ${HOST_DATA_DIR}
ln -s ${REFERENCE_FILE} ${HOST_DATA_DIR}
```

## Step 3: [*OPTIONAL*] Download Pre-Built RUFUS Hash Tables and Reference Indexes (~1TB)
RUFUS utilizes hash tables to identify unique kmers within a subject sample. Some of these hash tables have already been created and can be downloaded and stored locally, if desired. This is a good idea if your servers have firewalls or do not allow https traffic. **NOTE: this is approximately 1TB of data.** 

If these resources are not downloaded, RUFUS will automatically fetch them as needed during the run - each run only downloads ~500MB of this data at a time, and deletes after use as to not overwhelm a host directory. 

### Option A: Download with AWS-CLI
```bash
# Internal control hashes (for single sample mode)
aws s3 sync s3://rufus.marth.lab/public_access_data/control_hashes/ ${HOST_DATA_DIR}/rufus_resources/control_hashes --no-sign-request

# 1000G hashes
aws s3 sync s3://rufus.marth.lab/public_access_data/kg1_hashes/ ${HOST_DATA_DIR}/rufus_resources/kg1_hashes --no-sign-request
```

### Option B: Download with Rclone
```bash
# Internal control hashes (for single sample mode)
rclone copy :s3:rufus.marth.lab/public_access_data/control_hashes/ ${HOST_DATA_DIR}/rufus_resources/control_hashes --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P

# 1000G hashes
rclone copy :s3:rufus.marth.lab/public_access_data/control_hashes/ ${HOST_DATA_DIR}/rufus_resources/control_hashes --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P
```

## Step 4: Fill Out the RUFUS Environment File
The `rufus.env` file, which should now be found at `${HOST_DATA_DIR}/rufus_resources/rufus.env`, coordinates passing arguments into the RUFUS launch script. 

## Step 5: Launch RUFUS
```bash
chmod u+x launch_rufus.sh 
PATH_TO_RUFUS_ENV="${HOST_DATA_DIR}/rufus_resources/rufus.env"
./launch_rufus.sh $PATH_TO_RUFUS_ENV
```

Upon completion of RUFUS, the resulting vcf will be stored in `${HOST_DATA_DIR}/` and supplemental files in `${HOST_DATA_DIR}/rufus_supplementals`
