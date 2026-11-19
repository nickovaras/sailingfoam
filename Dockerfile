# 1. Start from your desired base image (equivalent to docker pull)
FROM opencfd/openfoam-default:latest

USER root

# 2. Install your dependencies (equivalent to your manual apt-get steps)
RUN apt-get update && apt-get install -y paraview && rm -rf /var/lib/apt/lists/*

ENV QT_X11_NO_MITSHM=1

# Set your working directory inside the container
WORKDIR /home/openfoam/sailingfoam

# (Optional) Copy your local project files into the container
# COPY . /app

# 3. Keep the container alive or run your application 
# (Replaces 'sleep infinity' if you just want it to stay running, or use your app command)
CMD ["sleep", "infinity"]
