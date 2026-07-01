# Building and publishing RUFUS containers

RUFUS has **one** container definition — the root [`Dockerfile`](../Dockerfile). There is no
separate Singularity `.def`; SIFs are produced from the Docker image via `apptainer`. Builds,
tests, and publishing are automated by [`.github/workflows/build-publish.yml`](../.github/workflows/build-publish.yml).

## Submodule: modified jellyfish

RUFUS's custom jellyfish fork lives in its own repository
([stefinfection/modified-jellyfish](https://github.com/stefinfection/modified-jellyfish)) and is
pinned here as a git submodule at `src/modifiedJellyfish`. The container build compiles it from
that submodule. **Clone RUFUS with the submodule**, or the build (and any local `cmake`) will fail:
```bash
git clone --recursive https://github.com/stefinfection/RUFUS.git
# already cloned without --recursive:
git submodule update --init --recursive
```
To change the jellyfish source: edit/commit/tag the modified-jellyfish repo, then in RUFUS
`cd src/modifiedJellyfish && git checkout <new-tag> && cd ../.. && git add src/modifiedJellyfish`
and commit the updated pointer.

## Tag mapping (CI)

| Git event                | Docker Hub tag(s)                         | Stage |
|--------------------------|-------------------------------------------|-------|
| push to `docker` / `dev` | `stefinfection/rufus:dev`                 | dev   |
| push to `master`         | `stefinfection/rufus:stage`               | stage |
| push git tag `v*`        | `stefinfection/rufus:<version>` + `:latest` | prod |

Every build runs [`tests/smoke_test.sh`](../tests/smoke_test.sh) inside the freshly built image
and only pushes if it passes.

## Cutting a production release

1. Bump `RUFUS_VERSION` in [`resources/globals.txt`](../resources/globals.txt) and merge to `master`.
2. Tag and push:
   ```bash
   git tag v1.1.11
   git push origin v1.1.11
   ```
3. CI then automatically:
   - builds the image, smoke-tests it, and pushes `:<version>` + `:latest` to Docker Hub;
   - builds a SIF from that image with `apptainer`;
   - publishes it as a **new version of the existing Zenodo concept record** (shared concept DOI).

## Staging / using an image on HPC (CHPC)

Pull a published image straight into a SIF — no `sudo`, no manual `.def` build:
```bash
bash singularity/pull_staged_image.sh stage          # staging image
bash singularity/pull_staged_image.sh v1.1.11         # a specific prod release
```
This writes `rufus_<tag>.sif` into the zenodo_images dir and runs the smoke test against it.

## Required GitHub repo secrets

Set these on the `stefinfection/RUFUS` repository (Settings → Secrets and variables → Actions):

| Secret                     | Purpose                                            |
|----------------------------|----------------------------------------------------|
| `DOCKERHUB_USERNAME`       | Docker Hub login                                   |
| `DOCKERHUB_TOKEN`          | Docker Hub access token (push)                     |
| `ZENODO_TOKEN`             | Zenodo token (`deposit:write` + `deposit:actions`) |
| `ZENODO_CONCEPT_RECORD_ID` | Concept (all-versions) record id of the RUFUS Zenodo record |
