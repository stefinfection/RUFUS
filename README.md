# Running RUFUS on AWS

## Step 0: Prerequisites
Running RUFUS on AWS requires Docker. Installation instructions for getting Docker on your machine can be found [here](https://docs.docker.com/engine/install/).

## Step 1: Fetch RUFUS image from Docker Hub and Setup Files from Zenodo
```bash
docker pull stefinfection/rufus:latest
curl "https://s3.us-east-1.amazonaws.com/rufus.marth.lab/public_access_data/launch_resources/launch_rufus.sh" -o launch_rufus.sh
curl "https://s3.us-east-1.amazonaws.com/rufus.marth.lab/public_access_data/launch_resources/grch38_1mb_regions.txt" -o grch38_1mb_regions.txt
curl "https://s3.us-east-1.amazonaws.com/rufus.marth.lab/public_access_data/launch_resources/rufus.env" -o rufus.env
```
## Step 2: Fill Out the RUFUS Environment File
The `rufus.env` file, downloaded in Step 1, coordinates passing arguments into the RUFUS launch script. 

## Step 3: Launch RUFUS
```bash
chmod u+x launch_rufus.sh 
./launch_rufus.sh "${PATH_TO}/rufus.env"
```
RUFUS will automatically run in and write results to the current directory, unless '$WORKING_DIR' is set to otherwise in `rufus.env`.

## [*OPTIONAL*] Download Pre-Built RUFUS Hash Tables and Reference Indexes (~1TB)
RUFUS utilizes hash tables to identify unique kmers within a subject sample. Some of these hash tables have already been created and can be downloaded and stored locally, if desired. This is a good idea if your servers have firewalls or do not allow https traffic. **NOTE: this is approximately 1TB of data.** 

If these resources are not downloaded, RUFUS will automatically fetch them as needed during the run - each run only downloads ~500MB of this data at a time, and deletes after use as to not overwhelm a host directory. 

## Step 1: Download Hashes
#### Option A: Download with AWS-CLI
```bash
# Internal control hashes (for single sample mode)
aws s3 sync s3://rufus.marth.lab/public_access_data/control_hashes/ ${LOCAL_DESTINATION_DIR}/control_hashes --no-sign-request

# 1000G hashes
aws s3 sync s3://rufus.marth.lab/public_access_data/kg1_hashes/ ${LOCAL_DESTINATION_DIR}/kg1_hashes --no-sign-request
```

#### Option B: Download with Rclone
```bash
# Internal control hashes (for single sample mode)
rclone copy :s3:rufus.marth.lab/public_access_data/control_hashes/ ${LOCAL_DESTINATION_DIR}/control_hashes --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P

# 1000G hashes
rclone copy :s3:rufus.marth.lab/public_access_data/kg1_hashes/ ${LOCAL_DESTINATION_DIR}/kg1_hashes --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P
```

### Step 2: Update rufus.env Variables
Add or assign the following variables in `rufus.env`:
```bash
KG1_HASH_LOCAL_DIR=${PATH_TO_KG1_HASHES}
CONTROL_HASH_LOCAL_DIR=${PATH_TO_CONTROL_HASHES}

# Example if using download command from above
CONTROL_HASH_LOCAL_DIR=${LOCAL_DESTINATION_DIR}/control_hashes
KG1_HASH_LOCAL_DIR=${LOCAL_DESTINATION_DIR}/kg1_hashes
```