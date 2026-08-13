# Build de l'image
podman build -t python-with-git-lfs:1.0 -f -f Dockerfile-pipeline .

# Tag de l'image
podman tag python-with-git-lfs:1.0 docker.io/obs2co/python-with-git-lfs:1.0

# Push de l'image sur Docker Hub
podman login docker.io
podman  push docker.io/obs2co/python-with-git-lfs:1.0