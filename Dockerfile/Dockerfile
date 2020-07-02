FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu16.04
MAINTAINER CSuperlei<csuperlei@163.com>

ADD sources.list /etc/apt/sources.list

RUN apt upgrade -y && apt update -y && \
    mkdir -p /root/.cpan && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

RUN apt install -y wget curl net-tools iputils-ping locales  \
    zip unzip bzip2 apt-utils zlib1g zlib1g-dev \
    tmux screen \
    gcc g++ \
    vim git htop xclip cmake sudo tree jq \
    build-essential gfortran automake bash-completion \
    libapparmor1 libedit2 libc6 \
    psmisc rrdtool libzmq3-dev libtool apt-transport-https \
    && locale-gen en_US.UTF-8 && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# bash && ctags && cscope && gtags
RUN apt install cscope libncurses5-dev -y && \
    cd /tmp && \
    curl https://ftp.gnu.org/gnu/bash/bash-5.0.tar.gz -o bash-5.0.tar.gz && \
    tar xzf bash-5.0.tar.gz && cd bash-5.0 && ./configure && make && make install && \
    cd /tmp && \
#   git clone --depth 1 https://github.com/universal-ctags/ctags.git && \
#   cd ctags && ./autogen.sh && ./configure && make && make install && \
#   cd /tmp && \
#   curl http://ftp.vim.org/ftp/gnu/global/global-6.6.3.tar.gz -o global.tar.gz && \
#   tar xzf global.tar.gz && cd global-6.6.3 && ./configure --with-sqlite3 && make && make install && \
#   cd /tmp && \
#   curl https://ftp.gnu.org/pub/gnu/libiconv/libiconv-1.16.tar.gz -o libiconv.tar.gz && \
#   tar xzf libiconv.tar.gz && cd libiconv-1.16 && ./configure && make && make install && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# R language
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    apt upgrade -y && apt update -y && \
    apt install -y r-base-dev r-base r-base-core && \
    apt install openjdk-8-jdk xvfb libswt-gtk-4-java -y && \
    R CMD javareconf && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# rstudio
RUN cd /tmp && \
    curl https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.5033-amd64.deb -o rstudio.deb && \
    gdebi -n rstudio.deb && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# language
RUN apt update && \
    apt install -y language-pack-zh-hans && locale-gen zh_CN.UTF-8 && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# conda
ADD .condarc /root
ENV PATH=/opt/miniconda3/bin:$PATH
RUN cd /tmp && \
    rm -f /bin/bash && ln -s /usr/local/bin/bash /bin/bash && \
    curl https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda3.sh && \
    bash miniconda3.sh -b -p /opt/miniconda3 && \
    conda update -n base -c defaults conda pip && \
    conda clean -a -y && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# conda install
RUN conda install -n base -c conda-forge time libxml2 libxslt libssh2 krb5 ripgrep lazygit yarn nodejs=12.16 jupyterlab=2.0.1 && \
    /opt/miniconda3/bin/pip install --no-cache-dir -i http://pypi.douban.com/simple --trusted-host pypi.douban.com pynvim neovim-remote flake8 pygments ranger-fm msgpack-python jedi==0.15.2 && \
    /opt/miniconda3/bin/pip install --no-cache-dir -i http://pypi.douban.com/simple --trusted-host pypi.douban.com python-language-server && \
    conda clean -a -y && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# machine learning
RUN /opt/miniconda3/bin/pip install --no-cache-dir -i http://pypi.douban.com/simple --trusted-host pypi.douban.com pandas scikit-learn numpy matplotlib scipy seaborn ggplot plotly xgboost && \
    conda clean -a -y && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

## vim
RUN cd /tmp && \
    git clone --depth 1 https://github.com/vim/vim.git && \
    cd vim && \
    export LDFLAGS='-L/opt/miniconda3/lib -Wl,-rpath,/opt/miniconda3/lib' && \
    ./configure --with-features=huge \
      --enable-multibyte \
      --enable-python3interp=yes \
      --with-python3-config-dir=/opt/miniconda3/lib/python3.7/config-3.7m-x86_x64-linux-gnu \
      --prefix=/usr/local && \
    make -j16 && make install && \
    rm -rf /tmp/* /var/tmp/* /root/.cpan/*

# conda env deeplearning
#RUN conda create -n lei_env python=3.6 && \
#    source activate lei_env && \
#    conda install pytorch torchvision -c pytorch && \
#    conda install tensorflow-gpu keras && \
#    conda clean -a -y && \
#    conda deactivate && \
#    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*
