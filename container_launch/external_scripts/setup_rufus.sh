#!/bin/bash
# Sets up slurm or bash script for RUFUS run according to config_parser.yaml
# Checks user environment for Singularity or Docker as input in args.
# Spins up ephemeral container to build launch script.
# Stops container safely and prompts user how to run container.
# Run externally/directly on host.

# Args
container_type="${1:-}"             # "docker" or "singularity"
image="${2:-}"                      # image name or path
config="${3:-./rufus_config.yaml}"  # config path (host path)
workdir="$(pwd)"                    # host working dir to bind into containe

# Constants
DOCKER="docker"
SINGULARITY="singularity"
INSTANCE_NAME="rufus_setup_instance"
CONTAINER_NAME="rufus_setup"
BUILD_SCRIPT="/opt/RUFUS/container_launch/build_launch_script.sh"

# NOTE: these paths must match those in globals.env
RUN_DIR="/work/"
RUNTIME_TEMP_DIR="$RUN_DIR/rufus_temp/"
INTERNAL_CONFIG_PATH="${RUNTIME_TEMP_DIR}rufus_config.yaml"

# Check for required args
if [[ -z "$container_type" || -z "$image" ]]; then
  echo "Usage: $0 <container_type: docker|singularity> <image> [config]" >&2
  exit 2
fi

if [[ ! -f "$config" ]]; then
  echo "Config file not found: $config" >&2
  exit 3
fi

# Check that docker/singularity installed
if [[ "$container_type" == "$DOCKER" ]]; then
  command -v docker >/dev/null 2>&1 || { echo "docker not found"; exit 4; }
elif [[ "$container_type" == "$SINGULARITY" ]]; then
  command -v singularity >/dev/null 2>&1 || { echo "singularity not found"; exit 4; }
else
  echo "Unknown container_type: $container_type (expected: docker or singularity)" >&2
  exit 5
fi

# Set trap to stop instances if unexpected exit
cleanup_docker() {
  if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}\$"; then
    docker stop "${CONTAINER_NAME}" >/dev/null 2>&1 || true
  fi
}
cleanup_singularity() {
  if singularity instance list | grep -q "${INSTANCE_NAME}"; then
    singularity instance stop "${INSTANCE_NAME}" >/dev/null 2>&1 || true
  fi
}
trap 'cleanup_docker; cleanup_singularity' EXIT

# Ephemeral container spin up to parse yaml + pull out helper functions
if [ "$container_type" == "$DOCKER" ]; then
    launch_out="${workdir}/launch_rufus.sh"

    echo "[launcher] Starting ephemeral docker container to generate launch script..."  
    # Run detached container as same uid so files are written with correct ownership
    docker run -d --rm --name "${CONTAINER_NAME}" \
        -u "$(id -u):$(id -g)" \
        -v "${workdir}:${RUN_DIR}" \
        -v "${config}:${INTERNAL_CONFIG_PATH}" \
        --cap-add SYS_ADMIN \
        --entrypoint sleep \
        "$image" 3600 >/dev/null
        
    echo "[launcher] Generating bash script for RUFUS launch..."
    docker exec --user "$(id -u):$(id -g)" "${CONTAINER_NAME}" "$BUILD_SCRIPT" "$container_type" "$config"

    # TODO: is this the correct way to copy over file?
    # TODO: can I get access to env file here? If not, could return region mode from $BUILD_SCRIPT call above?
    if [ "$WINDOWED_MODE" ]; then
      echo "[launcher] Creating RUFUS temporary run files..."
      docker exec --user "$(id -u):$(id -g)" "${CONTAINER_NAME}" mkdir "${RUN_DIR}/rufus_temp"
      docker exec --user "$(id -u):$(id -g)" "${CONTAINER_NAME}" cp "$DOCKER_REG_WRAPPER" "${RUN_DIR}/rufus_temp"
    fi

    docker stop "${CONTAINER_NAME}" >/dev/null
    echo "[launcher] Launch script generated for RUFUS run. Review script in ${workdir} and then run with \"sh $launch_out\""
else 
    launch_out="${workdir}/launch_rufus.slurm"

    echo "[launcher] Starting singularity instance to run planner..."
    singularity instance start --bind "${workdir}:/work" --bind "${config}:${INTERNAL_CONFIG_PATH}" "$image" "${INSTANCE_NAME}"

    echo "[launcher] Generating launch script inside singularity instance..."
    singularity exec instance://"${INSTANCE_NAME}" "$BUILD_SCRIPT" "$container_type" "$config"
    singularity instance stop "${INSTANCE_NAME}"
    echo "[launcher] Slurm launch script generated for RUFUS run. Review script in ${workdir} and then run with \"sbatch $launch_out\""
fi