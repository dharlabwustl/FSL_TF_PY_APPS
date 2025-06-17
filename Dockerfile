# Use NVIDIA TensorFlow base image with CUDA
FROM nvcr.io/nvidia/tensorflow:23.06-tf2-py3

# Set working directory
WORKDIR /workspace/deepreg

# Copy DeepReg code into container
COPY DeepReg/ /workspace/deepreg/DeepReg

# Install system dependencies
# Add NeuroDebian repository and install FSL
# Install basic system tools and dependencies


RUN apt-get update && apt-get install -y \
    curl \
    gnupg \
    unzip \
    git \
    vim \
    zip \
    tree \
    libgl1 \
    libglib2.0-0 \
    dcm2niix \
    bc \
    python3-pyqt5 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
#######################

# Prepare folders for scripts and data
RUN mkdir -p /callfromgithub
RUN chmod 755 /callfromgithub
RUN chmod 755 /workspace/deepreg/DeepReg
COPY downloadcodefromgithub.sh /callfromgithub/
RUN chmod +x /callfromgithub/downloadcodefromgithub.sh
RUN mkdir -p /workspace/deepreg/DeepReg/demos/classical_mr_prostate_nonrigid/dataset
RUN chmod 755 /workspace/deepreg/DeepReg/demos/classical_mr_prostate_nonrigid/dataset

# Install Python dependencies
RUN pip3 install \
  nibabel \
  numpy \
  xmltodict \
  pandas \
  requests \
  pydicom \
  python-gdcm \
  glob2 \
  scipy \
  pypng \
  PyGithub \
  SimpleITK \
  h5py \
  webcolors \
  antspyx \
  SQLAlchemy \
  mysql-connector-python==8.0.27 \
  opencv-python \
  scikit-image
# Install LaTeX (with common packages)
#RUN apt-get update && apt-get install -y \
#    texlive-latex-base \
#    texlive-latex-extra \
#    texlive-fonts-recommended \
#    texlive-fonts-extra \
#    texlive-science \
#    texlive-luatex \
#    texlive-xetex \
#    latexmk \
#    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Environment variables for REDCap / database connections (customize as needed)
ENV REDCAP_API='36F3BA05DE0507BEBDFB94CC5DA13F93'
ENV GOOGLE_MYSQL_DB_IP='34.58.59.235'
ENV GOOGLE_MYSQL_DB_PASS='dharlabwustl1!'



################################
# Install FSL from official FMRIB script
RUN curl -sSL https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/fslinstaller.py -o fslinstaller.py && \
    python3 fslinstaller.py -V 6.0.7 -d /usr/local/fsl -q && \
    rm fslinstaller.py

# Set FSL environment
#ENV FSLDIR=/usr/local/fsl
#ENV PATH=${FSLDIR}/bin:${PATH}
#ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV DEBIAN_FRONTEND=noninteractive

#ENV FSLDIR=/usr/share/fsl/6.0
#ENV PATH=${FSLDIR}/bin:${PATH}
#ENV FSLOUTPUTTYPE=NIFTI_GZ

# Set FSL environment variables
#ENV FSLDIR=/usr/local/fsl
#ENV PATH=${FSLDIR}/bin:${PATH}
#ENV FSLOUTPUTTYPE=NIFTI_GZ
#ENV DEBIAN_FRONTEND=noninteractive

# Upgrade pip
RUN pip install --upgrade pip

# Default command
CMD ["/bin/bash"]
